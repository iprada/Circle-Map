import numpy as np
from numba import njit, prange
from typing import List, Tuple, Optional
import pandas as pd
import pysam
from multiprocessing import Pool, cpu_count
from functools import partial
from config import BASE_MAP, REALIGN_FLANK, PROB_CUTOFF, MIN_SPLIT_SUPPORT


def compute_bg_freqs(ref_seq_str: str) -> np.ndarray:
    """
    Compute nucleotide background frequencies from a reference sequence.

    Follows Circle-Map's approach: count A/C/G/T occurrences and normalize.
    Returns [freq_A, freq_C, freq_G, freq_T] matching BASE_MAP encoding
    (A=0, C=1, G=2, T=3). N bases are ignored.

    When fewer than 10 non-N bases are present, falls back to uniform 0.25
    to avoid degenerate frequency vectors.
    """
    seq = ref_seq_str.upper()
    n_a = seq.count('A')
    n_c = seq.count('C')
    n_g = seq.count('G')
    n_t = seq.count('T')
    total = n_a + n_c + n_g + n_t
    if total < 10:
        return np.array([0.25, 0.25, 0.25, 0.25], dtype=np.float64)
    return np.array([n_a, n_c, n_g, n_t], dtype=np.float64) / total


@njit(parallel=True, cache=True)
def _batch_realign_kernel(
    read_seqs: np.ndarray,
    read_probs: np.ndarray,
    ref_seq: np.ndarray,
    clip_is_head: np.ndarray,
    bg_freqs: np.ndarray,
) -> np.ndarray:
    """
    Vectorized semi-global realignment scoring kernel with directional constraint
    and GC-aware null model.

    For each offset in the reference window, the log-likelihood ratio is:

        LLR = log2(P(read | aligned at offset) / P(read | null))

    where:
      - P(read | aligned) uses base-quality-dependent match/mismatch probabilities
      - P(read | null) uses the local nucleotide composition of the reference
        interval (freq_A, freq_C, freq_G, freq_T) instead of uniform 0.25

    This follows Circle-Map's PSSM approach and prevents inflated scores in
    GC-rich or AT-rich regions where uniform 0.25 overestimates evidence.

    Head clips (soft-clipped 5'): aligned portion extends RIGHT of junction
    Tail clips (soft-clipped 3'): aligned portion extends LEFT of junction
    """
    MIN_OVERLAP = 5
    NEG_INF = -1e9
    EPS = 1e-10

    n_reads = read_seqs.shape[0]
    read_len = read_seqs.shape[1]
    ref_len = ref_seq.shape[0]
    scores = np.full(n_reads, NEG_INF, dtype=np.float64)

    if ref_len == 0 or read_len == 0:
        return scores

    mid = ref_len // 2  # approximate junction position

    for i in prange(n_reads):
        best_for_read = NEG_INF
        is_head = clip_is_head[i]

        if is_head:
            offset_start = mid
            offset_end = ref_len - 1
        else:
            offset_start = 0
            offset_end = mid

        for offset in range(offset_start, min(offset_end + 1, ref_len)):
            window_score = 0.0
            overlaps = 0
            k_max = min(read_len, ref_len - offset)
            for k in range(k_max):
                r_base = read_seqs[i, k]
                ref_base = ref_seq[k + offset]
                if r_base == 4 or ref_base == 4:
                    continue
                p_corr = read_probs[i, k]
                # GC-aware null: use base frequency from reference window
                null_freq = max(bg_freqs[ref_base], EPS)
                if r_base == ref_base:
                    window_score += np.log2(max(p_corr, EPS) / null_freq)
                else:
                    window_score += np.log2(max((1.0 - p_corr) / 3.0, EPS) / null_freq)
                overlaps += 1
            if overlaps >= MIN_OVERLAP and window_score > best_for_read:
                best_for_read = window_score
        scores[i] = best_for_read

    return scores


def encode_batch(
    sequences: List[str], qualities: List[str], max_len: int
) -> Tuple[np.ndarray, np.ndarray]:
    """Convert lists of strings to padded numpy matrices."""
    n = len(sequences)
    seq_mat = np.full((n, max_len), 4, dtype=np.int8)
    prob_mat = np.full((n, max_len), 0.1, dtype=np.float32)

    for i, (seq, qual) in enumerate(zip(sequences, qualities)):
        l = min(len(seq), max_len)
        for j in range(l):
            seq_mat[i, j] = BASE_MAP.get(seq[j].upper(), 4)
            q = ord(qual[j]) - 33
            if q < 0:
                q = 0
            elif q > 40:
                q = 40
            prob_mat[i, j] = 1.0 - (10.0 ** (-q / 10.0))

    return seq_mat, prob_mat


def _per_read_probabilities(log2_scores: np.ndarray) -> np.ndarray:
    llr = np.asarray(log2_scores, dtype=np.float64)
    llr = np.where(np.isfinite(llr), llr, -700.0)
    return 1.0 / (1.0 + np.power(2.0, -llr))


def combine_read_probabilities(log2_scores: np.ndarray) -> float:
    if log2_scores.size == 0:
        return 0.0
    p_per_read = _per_read_probabilities(log2_scores)
    p_clipped = np.clip(p_per_read, 0.0, 1.0 - 1e-15)
    log_complement = np.log(1.0 - p_clipped).sum()
    return float(1.0 - np.exp(log_complement))


def _process_junction_batch(args):
    """
    Process a BATCH of junctions in a single worker call.
    Opens the FastaFile ONCE per batch instead of once per junction.

    Args:
        args: (ref_fasta, [(chrom, pos, seq_list, qual_list, clip_type_list), ...])

    Returns:
        List of dicts with junction results.
    """
    ref_fasta, junctions = args
    results = []
    try:
        genome = pysam.FastaFile(ref_fasta)
    except Exception:
        return []

    for chrom, pos, seqs, quals, clip_types in junctions:
        try:
            start = max(0, pos - REALIGN_FLANK)
            end = pos + REALIGN_FLANK
            ref_seq_str = genome.fetch(chrom, start, end)
        except Exception:
            continue

        ref_enc = np.array(
            [BASE_MAP.get(b, 4) for b in ref_seq_str.upper()], dtype=np.int8
        )
        if ref_enc.shape[0] == 0:
            continue

        bg_freqs = compute_bg_freqs(ref_seq_str)
        max_len = max(len(s) for s in seqs)
        r_seq, r_prob = encode_batch(seqs, quals, max_len)

        clip_is_head = np.array(
            [1 if ct == 'head' else 0 for ct in clip_types], dtype=np.int8
        )

        scores = _batch_realign_kernel(r_seq, r_prob, ref_enc, clip_is_head, bg_freqs)

        finite_scores = scores[np.isfinite(scores) & (scores > -1e8)]
        prob = combine_read_probabilities(scores)

        results.append({
            'chr': chrom,
            'pos': pos,
            'support_reads': len(seqs),
            'mean_score': float(finite_scores.mean()) if finite_scores.size else 0.0,
            'posterior_prob': prob,
            'is_valid': prob > PROB_CUTOFF,
        })

    genome.close()
    return results


class ProbabilisticRealigner:
    def __init__(self, genome_fasta: str, n_threads: int = 0):
        self.genome_fasta = genome_fasta
        self.n_threads = n_threads if n_threads > 0 else cpu_count()
        self.genome = pysam.FastaFile(genome_fasta)

    def realign_junctions(self, candidates: pd.DataFrame) -> pd.DataFrame:
        if candidates.empty:
            return pd.DataFrame(
                columns=['chr', 'pos', 'support_reads', 'mean_score',
                         'posterior_prob', 'is_valid']
            )

        grouped = candidates.groupby(['chr', 'pos'])
        total_raw = len(grouped)
        filtered_groups = [(n, g) for n, g in grouped if len(g) >= MIN_SPLIT_SUPPORT]
        total_kept = len(filtered_groups)
        print(f"   Processing {total_kept} junctions "
              f"(filtered {total_raw - total_kept} with <{MIN_SPLIT_SUPPORT} support)",
              flush=True)

        if total_kept == 0:
            return pd.DataFrame(
                columns=['chr', 'pos', 'support_reads', 'mean_score',
                         'posterior_prob', 'is_valid']
            )

        # --- Batch junctions to reduce pickle + FastaFile overhead ---
        # Each batch opens pysam.FastaFile ONCE instead of once per junction.
        BATCH_SIZE = 2000
        batches = []
        batch_items = []
        for (chrom, pos), group in filtered_groups:
            batch_items.append((
                chrom, pos,
                group['seq'].tolist(),
                group['qual'].tolist(),
                group['type'].tolist(),
            ))
            if len(batch_items) >= BATCH_SIZE:
                batches.append((self.genome_fasta, batch_items))
                batch_items = []
        if batch_items:
            batches.append((self.genome_fasta, batch_items))

        n_batches = len(batches)
        n_workers = min(self.n_threads, n_batches)
        print(f"   {n_batches} batches of ~{BATCH_SIZE}, {n_workers} workers", flush=True)

        if n_workers <= 1:
            results = []
            for i, ba in enumerate(batches):
                for r in _process_junction_batch(ba):
                    results.append(r)
                if (i + 1) % 10 == 0:
                    print(f"   Batch {i + 1}/{n_batches} done ({len(results)} junctions)", flush=True)
        else:
            with Pool(n_workers) as pool:
                results = []
                for i, batch_out in enumerate(pool.imap_unordered(
                    _process_junction_batch, batches, chunksize=1
                )):
                    results.extend(batch_out)
                    if (i + 1) % 10 == 0:
                        print(f"   Batch {i + 1}/{n_batches} done ({len(results)} junctions)", flush=True)

        print(f"   Realignment complete: {len(results)} junctions processed", flush=True)
        return pd.DataFrame(results)