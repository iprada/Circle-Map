"""
Read Extractor Module for ecDNA Detection

This module implements the first stage of the CircleDetect pipeline: extracting 
candidate circular DNA reads from paired-end sequencing data. It uses bwa-mem2 
for alignment, samblaster for extracting split and discordant reads, and pysam 
for efficient BAM parsing.

Key Features:
- Paired-end read alignment using bwa-mem2
- Split read extraction (soft-clipped reads indicating junction breakpoints)
- Discordant read pair extraction (long-range connectivity evidence)
- Quality filtering based on MAPQ and clip length thresholds
- Preservation of sequence and quality scores for downstream realignment

Author: CircleDetect Team
Based on: Circle-Map ReadExtractor methodology
"""

import subprocess
import pysam
import pandas as pd
import numpy as np
import os
from typing import Tuple, List, Dict, Optional
from config import (
    BWA_MEM2_PATH, SAMBLASTER_PATH, SAMTOOLS_PATH,
    MIN_MAPQ, MIN_MAPQ_DISCORDANT, MIN_CLIP_LEN, MAX_INSERT_SIZE, OUTPUT_DIR
)


class ReadExtractor:
    """
    Extract candidate ecDNA reads from paired-end sequencing data.

    Two modes:
    1. **FASTQ mode** (default): Aligns reads with bwa-mem2, then extracts split
       and discordant reads via samblaster.
    2. **Pre-aligned BAM mode**: Parses existing split-reads BAM and discordant
       BAM files produced by a previous samblaster run, plus a coordinate-sorted,
       indexed full-alignment BAM for depth/continuity analysis.

    Attributes:
        ref (str): Path to reference FASTA file
        fastq_r1 (str): Path to Read 1 FASTQ file (None in BAM mode)
        fastq_r2 (str): Path to Read 2 FASTQ file (None in BAM mode)
        bam_path (str): Path to coordinate-sorted, indexed BAM (full alignment)

    Example:
        >>> # FASTQ mode
        >>> extractor = ReadExtractor("ref.fasta", "reads_1.fastq", "reads_2.fastq")
        >>> split_df, discord_df = extractor.align_and_extract()

        >>> # Pre-aligned BAM mode
        >>> extractor = ReadExtractor.from_bams(
        ...     "ref.fasta", "split.bam", "discordant.bam", "aligned.bam"
        ... )
        >>> split_df, discord_df = extractor.parse_existing_bams()
    """

    def __init__(
        self,
        ref_fasta: str,
        fastq_r1: Optional[str] = None,
        fastq_r2: Optional[str] = None,
        bam_path: Optional[str] = None,
        align_threads: int = 8,
        two_pass: bool = False,
    ) -> None:
        """
        Initialize the ReadExtractor.

        Args:
            ref_fasta: Path to reference genome FASTA (must be indexed).
            fastq_r1: Path to Read 1 FASTQ (None in BAM mode).
            fastq_r2: Path to Read 2 FASTQ (None in BAM mode).
            bam_path: Path to coordinate-sorted, indexed BAM (None in FASTQ mode).
        """
        self.ref = ref_fasta
        self.fastq_r1 = fastq_r1
        self.fastq_r2 = fastq_r2
        self.align_threads = align_threads
        self.two_pass = two_pass
        self.bam_path = bam_path if bam_path else os.path.join(OUTPUT_DIR, "aligned_sorted.bam")
        self._split_bam_path: Optional[str] = None
        self._discord_bam_path: Optional[str] = None

    @classmethod
    def from_bams(
        cls,
        ref_fasta: str,
        split_bam: str,
        discordant_bam: str,
        aligned_bam: str,
    ) -> "ReadExtractor":
        """
        Create a ReadExtractor from pre-existing BAM files (skip alignment).

        Args:
            ref_fasta: Path to reference genome FASTA.
            split_bam: Path to split-reads BAM (from samblaster --splitterFile).
            discordant_bam: Path to discordant-reads BAM (from samblaster --discordantFile).
            aligned_bam: Path to coordinate-sorted, INDEXED full-alignment BAM.

        Returns:
            ReadExtractor configured for BAM-mode parsing.
        """
        extractor = cls(ref_fasta, bam_path=aligned_bam)
        extractor._split_bam_path = split_bam
        extractor._discord_bam_path = discordant_bam
        return extractor

    def parse_existing_bams(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Parse pre-existing split-reads and discordant BAMs without running alignment.

        Requires that the extractor was created via ``from_bams()``.

        Returns:
            Tuple of (split_df, discord_df) — same format as ``align_and_extract()``.
        """
        if self._split_bam_path is None or self._discord_bam_path is None:
            raise ValueError(
                "parse_existing_bams() requires a ReadExtractor created via from_bams(). "
                "Use align_and_extract() for FASTQ mode."
            )
        print(f">> Parsing split reads: {self._split_bam_path}")
        split_df = self._parse_split_reads(self._split_bam_path)

        print(f">> Parsing discordant reads: {self._discord_bam_path}")
        discord_df = self._parse_discordant_reads(self._discord_bam_path)

        return split_df, discord_df
        
    def align_and_extract(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Execute the complete read alignment and extraction pipeline.

        Two modes (controlled by ``two_pass`` flag):

        **Single-pass** (default, ``two_pass=False``):
        bwa-mem2 ``-5`` (split-read aware, keeps PE pairing)
        → samblaster extracts both split and discordant reads
        → coordinate-sort for confidence-scoring BAM.

        Good compromise: ``-5`` gives junction sensitivity without breaking
        samblaster's PE pairing detection.

        **Two-pass** (``two_pass=True``):
        Pass 1 — bwa-mem2 ``-5SP`` (max split-read sensitivity)
        → samblaster ``--splitterFile`` only.
        Pass 2 — bwa-mem2 standard mode (max PE accuracy)
        → samblaster ``--discordantFile`` + coordinate-sort for BAM.

        2× alignment time but maximises both split-read and discordant-pair
        detection sensitivity.
        """
        split_out = os.path.join(OUTPUT_DIR, "split_reads.bam")
        discord_out = os.path.join(OUTPUT_DIR, "discordant_reads.bam")
        aligned_sam = os.path.join(OUTPUT_DIR, "aligned.sam")
        coord_bam = self.bam_path
        T = str(self.align_threads)

        if self.two_pass:
            # ── Pass 1: -5SP mode → split reads only ──
            print(">> [Pass 1/2] bwa-mem2 -5SP → samblaster (split reads)...",
                  flush=True)
            cmd_pass1 = (
                f"set -o pipefail; "
                f"{BWA_MEM2_PATH} mem -t {T} -5SP -Y -k 19 "
                f"'{self.ref}' '{self.fastq_r1}' '{self.fastq_r2}' | "
                f"{SAMBLASTER_PATH} "
                f"--splitterFile '{split_out}' "
                f"-o /dev/null"
            )
            subprocess.run(cmd_pass1, shell=True, check=True, executable='/bin/bash')

            # ── Pass 2: standard mode → discordant + full BAM ──
            print(">> [Pass 2/2] bwa-mem2 (standard) → samblaster (discordant) "
                  "+ aligned.sam...", flush=True)
            cmd_pass2 = (
                f"set -o pipefail; "
                f"{BWA_MEM2_PATH} mem -t {T} -Y -k 19 "
                f"'{self.ref}' '{self.fastq_r1}' '{self.fastq_r2}' | "
                f"{SAMBLASTER_PATH} "
                f"--discordantFile '{discord_out}' "
                f"> '{aligned_sam}'"
            )
            subprocess.run(cmd_pass2, shell=True, check=True, executable='/bin/bash')

        else:
            # ── Single pass: -5 mode (compromise) ──
            print(">> Aligning: bwa-mem2 -5 + samblaster (split + discordant)...",
                  flush=True)
            cmd = (
                f"set -o pipefail; "
                f"{BWA_MEM2_PATH} mem -t {T} -5 -Y -k 19 "
                f"'{self.ref}' '{self.fastq_r1}' '{self.fastq_r2}' | "
                f"{SAMBLASTER_PATH} "
                f"--splitterFile '{split_out}' "
                f"--discordantFile '{discord_out}' "
                f"> '{aligned_sam}'"
            )
            subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')

        # ── Coordinate-sort the full alignment for confidence scoring ──
        print(">> Sorting aligned.sam by coordinate...", flush=True)
        subprocess.run(
            [SAMTOOLS_PATH, "sort", "-@", T, "-o", coord_bam, aligned_sam],
            check=True,
        )
        subprocess.run([SAMTOOLS_PATH, "index", coord_bam], check=True)

        # Clean up intermediate SAM
        if os.path.exists(aligned_sam):
            os.remove(aligned_sam)

        print(">> Parsing Split Reads...", flush=True)
        split_df = self._parse_split_reads(split_out)

        print(">> Parsing Discordant Reads...", flush=True)
        discord_df = self._parse_discordant_reads(discord_out)

        return split_df, discord_df

    def _parse_split_reads(self, bam_path: str) -> pd.DataFrame:
        """
        Stream through split-reads BAM and collect junction candidates.

        Uses a single pass with flat list storage (not dict-of-lists) to balance
        memory and I/O. Progress is reported every 5M reads.
        """
        # Flat lists, one element per clip that passes filters.
        seqs: List[str] = []
        quals: List[str] = []
        chrs: List[str] = []
        poss: List[int] = []
        ctypes: List[str] = []
        read_names: List[str] = []
        strands: List[str] = []
        mate_chrs: List[Optional[str]] = []
        mate_poss: List[int] = []
        mate_unmapped_list: List[bool] = []
        proper_list: List[bool] = []

        total_reads = 0
        kept_clips = 0
        report_interval = 5_000_000

        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                total_reads += 1
                if total_reads % report_interval == 0:
                    print(f"   Parsed {total_reads:,} reads, {kept_clips:,} clips "
                          f"kept ({len(seqs):,} total)", flush=True)

                if read.is_unmapped or read.mapping_quality < MIN_MAPQ:
                    continue

                clips = self._get_soft_clips(read)
                if not clips:
                    continue

                mate_unmapped = read.mate_is_unmapped if read.is_paired else True
                strand = '-' if read.is_reverse else '+'
                mate_chr = read.next_reference_name if not mate_unmapped else None
                mate_pos = read.next_reference_start if not mate_unmapped else -1
                proper = read.is_paired and read.is_proper_pair
                chrom = read.reference_name
                qname = read.query_name

                for clip_pos, clip_seq, clip_qual, clip_type in clips:
                    if len(clip_seq) < MIN_CLIP_LEN:
                        continue
                    kept_clips += 1
                    seqs.append(clip_seq)
                    quals.append(clip_qual)
                    chrs.append(chrom)
                    poss.append(clip_pos)
                    ctypes.append(clip_type)
                    read_names.append(qname)
                    strands.append(strand)
                    mate_chrs.append(mate_chr)
                    mate_poss.append(mate_pos)
                    mate_unmapped_list.append(mate_unmapped)
                    proper_list.append(proper)

        if total_reads == 0:
            return pd.DataFrame()

        print(f"   Parsing complete: {total_reads:,} reads scanned, "
              f"{kept_clips:,} clips kept.", flush=True)

        df = pd.DataFrame({
            'read_name': read_names,
            'chr': chrs,
            'pos': poss,
            'strand': strands,
            'seq': seqs,
            'qual': quals,
            'type': ctypes,
            'mate_chr': mate_chrs,
            'mate_pos': mate_poss,
            'is_proper_pair': proper_list,
            'mate_unmapped': mate_unmapped_list,
        })
        df['is_same_chr'] = df['chr'] == df['mate_chr']

        mem_gb = df.memory_usage(deep=True).sum() / 1e9
        print(f"   DataFrame: {len(df):,} rows, {mem_gb:.1f} GB", flush=True)
        return df

    @staticmethod
    def _revcomp(seq: str) -> str:
        """Reverse-complement a nucleotide string (N stays N)."""
        comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'N': 'N', 'n': 'n'}
        return ''.join(comp.get(b, 'N') for b in reversed(seq))

    def _parse_discordant_reads(self, bam_path: str) -> pd.DataFrame:
        """
        Parse BAM file to extract discordant read pairs supporting long-range connections.
        
        Discordant reads are paired-end reads where:
        - The mates map to different chromosomes (inter-chromosomal)
        - The mates map too far apart on the same chromosome (> MAX_INSERT_SIZE)
        - The orientation is unexpected (not properly paired)
        
        These provide supporting evidence for connections between distant genomic regions,
        which is crucial for detecting multi-segment ecDNA circles.
        
        Args:
            bam_path (str): Path to BAM file containing discordant reads
            
        Returns:
            pd.DataFrame: DataFrame with columns:
                - 'read_name': Read identifier
                - 'chr1', 'pos1': First mate's position
                - 'chr2', 'pos2': Second mate's position
                - 'strand1', 'strand2': Strand orientations
        """
        records: List[Dict[str, any]] = []
        seen_pairs: set = set()
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped or read.mate_is_unmapped:
                    continue
                if read.mapping_quality < MIN_MAPQ_DISCORDANT:
                    continue

                chr1, pos1 = read.reference_name, read.reference_start
                chr2, pos2 = read.next_reference_name, read.next_reference_start

                # Deduplicate: both R1 and R2 appear in the discordant BAM.
                a = (chr1, pos1)
                b = (chr2, pos2)
                key = (min(a, b), max(a, b), read.query_name)
                if key in seen_pairs:
                    continue
                seen_pairs.add(key)

                if chr1 == chr2:
                    dist = abs(pos1 - pos2)
                    if dist <= MAX_INSERT_SIZE:
                        continue

                records.append({
                    'read_name': read.query_name,
                    'chr1': chr1, 'pos1': pos1,
                    'chr2': chr2, 'pos2': pos2,
                    'strand1': '-' if read.is_reverse else '+',
                    'strand2': '-' if read.mate_is_reverse else '+',
                })
        return pd.DataFrame(records)

    def _get_soft_clips(self, read: pysam.AlignedSegment) -> List[Tuple[int, str, str, str]]:
        """
        Extract soft-clipped sequences and their quality scores from a read.
        
        Parses the CIGAR string to identify soft-clipped (S) operations and extracts
        the corresponding sequence and quality substrings.
        
        Args:
            read (pysam.AlignedSegment): Aligned read from pysam
            
        Returns:
            List[Tuple[int, str, str, str]]: List of tuples containing:
                - clip_pos: Genomic position at junction
                - clip_seq: Clipped nucleotide sequence
                - clip_qual: Quality scores as Phred string
                - clip_type: 'head' (5') or 'tail' (3') clip
                
        Note:
            - Head clips: occur at start of read (s_pos == 0)
            - Tail clips: occur after aligned portion
            - Coordinates approximate the junction position
        """
        clips: List[Tuple[int, str, str, str]] = []
        cigartypes = dict(zip(range(8), "MIDNSHP=X"))
        seq = read.query_sequence
        qual = read.query_qualities
        
        if seq is None or qual is None:
            return []
            
        pos = read.reference_start
        s_pos = 0  # Position in read sequence
        
        for op, length in read.cigartuples:
            op_char = cigartypes[op]
            if op_char == 'S':
                clip_seq = seq[s_pos:s_pos+length]
                # Convert quality integers to Phred string format
                clip_qual = "".join([chr(q + 33) for q in qual[s_pos:s_pos+length]])
                
                # Determine coordinate based on clip position:
                # Head clip: junction approximately at reference start
                # Tail clip: junction at end of aligned block
                if s_pos == 0:  # Head clip (5' end)
                    coord = pos
                    clip_type = 'head'
                else:  # Tail clip (3' end)
                    coord = pos
                    clip_type = 'tail'
                
                clips.append((coord, clip_seq, clip_qual, clip_type))
            
            # Update positions based on CIGAR operation.
            # Reference-consuming ops: M D N = X (advance reference coordinate)
            # Read-consuming ops:       M I S = X (advance read coordinate)
            if op_char in "MDN=X":
                pos += length
            if op_char in "MIS=X":
                s_pos += length

        return clips


def estimate_insert_size(
    bam_path: str,
    sample_size: int = 100_000,
    mapq_cutoff: int = 60,
    max_insert: int = 5000,
) -> Tuple[float, float]:
    """
    Estimate the paired-end insert size distribution (mean, std) from a
    coordinate-sorted BAM by sampling concordant proper pairs.

    Following Circle-Map's approach: collect `sample_size` R1 reads from
    same-chromosome proper pairs with MAPQ >= `mapq_cutoff`, record their
    template lengths, then robustly estimate mean and std after removing
    outliers beyond 3 standard deviations.

    Args:
        bam_path: Path to a coordinate-sorted, indexed BAM.
        sample_size: Number of concordant pairs to sample.
        mapq_cutoff: Minimum mapping quality for inclusion.
        max_insert: Hard cap to discard chimeric / outlier templates.

    Returns:
        (mean_insert, std_insert). Falls back to (350.0, 50.0) if too few
        pairs are found.
    """
    sizes: List[int] = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.mate_is_unmapped:
                continue
            if not read.is_paired or not read.is_proper_pair:
                continue
            if not read.is_read1:
                continue
            if read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < mapq_cutoff:
                continue
            if read.reference_name != read.next_reference_name:
                continue
            insert = abs(read.template_length)
            if 0 < insert <= max_insert:
                sizes.append(insert)
            if len(sizes) >= sample_size:
                break

    if len(sizes) < 100:
        return (350.0, 50.0)

    arr = np.asarray(sizes, dtype=np.float64)
    mean = float(arr.mean())
    std = float(arr.std())
    # One round of outlier trimming at 3 sigma for robustness.
    keep = np.abs(arr - mean) <= 3.0 * std
    if keep.sum() >= 100:
        arr = arr[keep]
        mean = float(arr.mean())
        std = float(arr.std())

    return (mean, std)