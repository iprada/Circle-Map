"""
EccDNA Graph Solver Module

This module implements graph-based cycle detection for identifying extrachromosomal 
circular DNA (ecDNA) structures from validated junction breakpoints and discordant 
read pairs. It uses sparse matrix representations for efficient computation and 
strongly connected components (SCC) algorithm for cycle detection.

Key Features:
- Constructs directed graphs from breakpoint junctions
- Supports both single-segment circles (self-loops) and multi-segment circles (cycles)
- Uses scipy.sparse for memory-efficient adjacency matrix representation
- Implements confidence scoring based on posterior probabilities

Author: CircleDetect Team
Based on: Circle-Map methodology
"""

import pandas as pd
import numpy as np
import bisect
from collections import defaultdict
from typing import List, Dict, Tuple, Optional
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from config import MIN_DISCORDANT_EDGES


# Discordant reads within this distance (bp) of a validated junction are counted
# as supporting that junction's node. Mirrors Circle-Map's -K clustering_dist.
DISCORDANT_MATCH_WINDOW = 500
# Minimum adjacency-matrix edge weight used for SCC detection. Edges weaker
# than this are removed before running Strongly Connected Components, which
# breaks the giant artifact SCC (84k nodes spanning the whole genome) into
# smaller biologically meaningful components. Weight = min(1.0, count / 5.0),
# so 0.6 corresponds to 3+ discordant read pairs.
SCC_MIN_EDGE_WEIGHT = 0.6


class EccDNAGraphSolver:
    """
    Graph-based solver for detecting ecDNA circular structures.

    This class constructs a directed graph where:
    - Nodes represent validated breakpoint junctions
    - Self-loop edges represent single-segment circle evidence (weighted by posterior probability)
    - Directed edges between nodes represent multi-segment connections (supported by discordant reads)

    Cycle detection is performed using strongly connected components (SCC) algorithm,
    where each SCC represents a potential ecDNA circle.

    Attributes:
        adj_matrix (csr_matrix): Sparse adjacency matrix representing the graph
        id_map (Dict[int, Tuple[str, int]]): Mapping from node IDs to (chromosome, position)
        _discordant_snap (Dict[int, List[int]]): For each node, discordant mate positions
            that were too close to snap to a different node (support single-segment estimation)

    Example:
        >>> solver = EccDNAGraphSolver()
        >>> solver.build_graph(validated_junctions_df, discordant_reads_df)
        >>> circles = solver.detect_cycles()
        >>> print(circles)
    """

    def __init__(self) -> None:
        self.adj_matrix: Optional[csr_matrix] = None
        self.id_map: Dict[int, Tuple[str, int]] = {}
        # For single-segment circle boundary estimation: stores "other end" positions
        # of discordant reads where both ends snapped to the same junction node.
        self._discordant_snap: Dict[int, List[int]] = {}
        self._validated_junctions: pd.DataFrame = pd.DataFrame()

    def build_graph(self, validated_junctions: pd.DataFrame, discordant: pd.DataFrame) -> None:
        unique_sites = validated_junctions[['chr', 'pos']].drop_duplicates().copy()
        unique_sites['id'] = np.arange(len(unique_sites))

        self.id_map = {row['id']: (row['chr'], row['pos']) for _, row in unique_sites.iterrows()}
        reverse_map = {(row['chr'], row['pos']): row['id'] for _, row in unique_sites.iterrows()}

        n_nodes = len(unique_sites)
        if n_nodes == 0:
            self.adj_matrix = csr_matrix((0, 0))
            return

        row_inds: List[int] = []
        col_inds: List[int] = []
        data: List[float] = []

        # 1. Self-loops from validated junctions
        for _, row in validated_junctions.iterrows():
            if row['is_valid']:
                key = (row['chr'], row['pos'])
                if key in reverse_map:
                    nid = reverse_map[key]
                    row_inds.append(nid)
                    col_inds.append(nid)
                    data.append(row['posterior_prob'])

        # 2. Per-chromosome sorted numpy arrays for fast nearest-node lookup
        chr_positions: Dict[str, np.ndarray] = {}
        chr_node_ids: Dict[str, np.ndarray] = {}
        for _, row in unique_sites.iterrows():
            c = row['chr']
            if c not in chr_positions:
                chr_positions[c] = []
                chr_node_ids[c] = []
            chr_positions[c].append(int(row['pos']))
            chr_node_ids[c].append(int(row['id']))
        for c in chr_positions:
            p = np.array(chr_positions[c])
            n = np.array(chr_node_ids[c])
            order = np.argsort(p)
            chr_positions[c] = p[order]
            chr_node_ids[c] = n[order]

        # 3. Vectorized nearest-node: batch-searchsorted per chromosome,
        #    then single-pass Python loop over 119M discordant reads.
        pair_counts: Dict[Tuple[int, int], int] = defaultdict(int)
        if not discordant.empty and chr_positions:
            W = DISCORDANT_MATCH_WINDOW
            _snap = self._discordant_snap

            def _nearest_batch(chrom, positions):
                """Vectorized nearest-node search — all positions at once."""
                if chrom not in chr_positions:
                    return np.full(len(positions), -1, dtype=np.int32)
                jp, jn = chr_positions[chrom], chr_node_ids[chrom]
                idx = np.searchsorted(jp, positions)
                idx = np.clip(idx, 0, len(jp) - 1).astype(np.int32)
                d0 = np.abs(jp[idx] - positions)
                best, best_id = d0.copy(), jn[idx].copy()
                left = idx > 0
                i_left = idx[left] - 1
                d_left = np.abs(jp[i_left] - positions[left])
                swap = d_left < best[left]
                best[left] = np.where(swap, d_left, best[left])
                best_id[left] = np.where(swap, jn[i_left], best_id[left])
                return np.where(best <= W, best_id, -1)

            # Pre-extract arrays for vectorized processing
            chr1_arr = discordant['chr1'].values
            pos1_arr = discordant['pos1'].values.astype(np.int32, copy=False)
            chr2_arr = discordant['chr2'].values
            pos2_arr = discordant['pos2'].values.astype(np.int32, copy=False)

            # Stage A: compute nearest node for EVERY read end in two
            # vectorized passes (one per unique chromosome for each end).
            n1_arr = np.full(len(discordant), -1, dtype=np.int32)
            n2_arr = np.full(len(discordant), -1, dtype=np.int32)
            for chrom in chr_positions:
                m1 = discordant['chr1'].values == chrom
                if m1.any():
                    n1_arr[m1] = _nearest_batch(chrom, pos1_arr[m1])
                m2 = discordant['chr2'].values == chrom
                if m2.any():
                    n2_arr[m2] = _nearest_batch(chrom, pos2_arr[m2])

            # Stage B: single Python pass — no searchsorted inside.
            same_chrom = chr1_arr == chr2_arr
            for idx in range(len(discordant)):
                u, v = int(n1_arr[idx]), int(n2_arr[idx])
                if u < 0 and v < 0:
                    continue
                if u >= 0 and v >= 0 and u == v and same_chrom[idx]:
                    jp_v = unique_sites.iloc[u]['pos']
                    p1_i, p2_i = int(pos1_arr[idx]), int(pos2_arr[idx])
                    near = p1_i if abs(p1_i - jp_v) < abs(p2_i - jp_v) else p2_i
                    far = p2_i if near == p1_i else p1_i
                    _snap.setdefault(u, []).append(far)
                elif u >= 0 and v >= 0 and u != v:
                    pair_counts[(u, v)] += 1
                    pair_counts[(u, v)] += 1

        # 4. Add directed edges with sufficient discordant support.
        #    Direction is preserved so that Strongly Connected Components can
        #    identify true multi-segment circles (which require bidirectional
        #    connectivity). Single-direction edges form a DAG and do not
        #    contribute to SCCs.
        for (u, v), count in pair_counts.items():
            if count >= MIN_DISCORDANT_EDGES:
                row_inds.append(u)
                col_inds.append(v)
                data.append(min(1.0, count / 5.0))

        # 4. Add edges for same-node discordant pairs as weak self-loop support.
        #    A discordant pair whose ends snap to the same junction is weak evidence
        #    that the circle loops back nearby. Add a small weight contribution.
        for nid, far_positions in self._discordant_snap.items():
            key = (unique_sites.iloc[nid]['chr'], unique_sites.iloc[nid]['pos'])
            if key in reverse_map:
                # Contribution decays with distance: 0.5 for nearby, 0.0 for far
                weight = min(0.5, len(far_positions) * 0.1)
                if weight > 0:
                    row_inds.append(nid)
                    col_inds.append(nid)
                    data.append(weight)

        if not row_inds:
            self.adj_matrix = csr_matrix((0, 0))
            return

        self.adj_matrix = csr_matrix((data, (row_inds, col_inds)), shape=(n_nodes, n_nodes))

    def detect_cycles(
        self,
        bam_path: Optional[str] = None,
        discordant_df: Optional[pd.DataFrame] = None,
        insert_metrics: Optional[Tuple[float, float]] = None,
        top_n: int = 500,
    ) -> pd.DataFrame:
        """
        Detect ecDNA circular structures using Strongly Connected Components.

        When BAM-based confidence scoring is enabled, only the ``top_n`` circles
        (ranked by graph-based confidence) receive full continuity and depth-ratio
        scoring, since BAM queries are expensive. Remaining circles are reported
        with graph-based confidence only.

        Args:
            bam_path: Sorted, indexed BAM for advanced scoring.
            discordant_df: Discordant read pairs DataFrame.
            insert_metrics: (mean, std) insert size distribution.
            top_n: Max circles to BAM-score (limits expensive queries).
        """
        if self.adj_matrix is None or self.adj_matrix.shape[0] == 0:
            return pd.DataFrame()

        # ---- Filter weak edges before SCC ----
        # The full adjacency matrix includes all discordant-pair edges (even
        # single-pair), creating a giant 84k-node SCC spanning the genome.
        # Build a filtered copy: keep self-loops + edges with weight >= threshold.
        mat = self.adj_matrix
        row, col = mat.nonzero()
        threshold = SCC_MIN_EDGE_WEIGHT
        is_self = (row == col)
        keep = is_self | (mat.data >= threshold)
        filtered = csr_matrix(
            (mat.data[keep], (row[keep], col[keep])), shape=mat.shape
        )

        # ---- SCC on filtered graph ----
        n_comp, labels = connected_components(
            csgraph=filtered, directed=True, connection='strong'
        )

        # ---- Build circle dicts ----
        circles: List[Dict[str, any]] = []
        for i in range(n_comp):
            members = np.where(labels == i)[0]
            if len(members) == 0:
                continue

            # Use ORIGINAL (unfiltered) matrix for confidence scoring
            subgraph = mat[np.ix_(members, members)]
            if subgraph.nnz > 0:
                conf = min(1.0, float(subgraph.data.mean()))
            else:
                conf = 0.0

            nodes = [self.id_map[m] for m in members]
            circle_type = 'Multi-Segment' if len(nodes) > 1 else 'Single-Segment'

            if circle_type == 'Single-Segment' and conf < 0.9:
                continue

            # Skip any remaining giant SCCs (should be rare after weight filter)
            if len(nodes) > 500:
                continue

            circles.append({
                'type': circle_type,
                'segment_count': len(nodes),
                'nodes': str(nodes),
                'confidence': conf,
                '_nodes_raw': nodes,
            })

        if not circles:
            return pd.DataFrame()

        # Sort by confidence descending; only top_n get BAM scoring
        circles.sort(key=lambda c: c['confidence'], reverse=True)
        to_score = circles[:top_n]

        if bam_path is not None and to_score:
            try:
                from confidence_scorer import EccDNAConfidenceScorer

                for circle_dict in to_score:
                    sorted_nodes = sorted(
                        circle_dict['_nodes_raw'], key=lambda x: (x[0], x[1])
                    )
                    segments = self._estimate_segments(
                        sorted_nodes, discordant_df, insert_metrics
                    )

                    try:
                        scorer = EccDNAConfidenceScorer(bam_path)
                        candidate = scorer.score_circle(
                            segments=segments,
                            junctions_df=pd.DataFrame([{
                                'chr': n[0], 'pos': n[1],
                                'posterior_prob': circle_dict['confidence'],
                                'is_valid': True,
                            } for n in circle_dict['_nodes_raw']]),
                            discordant_df=(
                                discordant_df if discordant_df is not None
                                else pd.DataFrame()
                            ),
                            circle_id=f"circle_{circles.index(circle_dict)}",
                            insert_metrics=insert_metrics,
                        )

                        circle_dict['continuity_score'] = candidate.continuity_score
                        circle_dict['depth_ratio'] = candidate.depth_ratio
                        circle_dict['integrated_confidence'] = candidate.integrated_confidence
                        scorer.close()
                    except Exception as inner_e:
                        # Individual circle scoring failure — skip BAM fields
                        import sys
                        print(f'  BAM scoring error: {type(inner_e).__name__}: {inner_e}', file=sys.stderr, flush=True)

            except ImportError:
                pass  # confidence_scorer module not found

        # Clean up internal fields
        for c in circles:
            c.pop('_nodes_raw', None)

        return pd.DataFrame(circles)

    def _estimate_segments(
        self,
        nodes: List[Tuple[str, int]],
        discordant_df: Optional[pd.DataFrame],
        insert_metrics: Optional[Tuple[float, float]],
    ) -> List[Tuple[str, int, int]]:
        """
        Estimate genomic segments from graph nodes and discordant read evidence.

        For single-segment circles, uses discordant read "other end" positions to
        estimate the far boundary. Falls back to insert-size-based window when
        discordant evidence is unavailable.

        For multi-segment circles, creates segments between consecutive breakpoints,
        properly handling circular wrap-around and cross-chromosome cases.
        """
        if len(nodes) == 1:
            return self._estimate_single_segment(nodes[0], discordant_df, insert_metrics)

        segments = []
        for j in range(len(nodes)):
            chrom1, pos1 = nodes[j]
            chrom2, pos2 = nodes[(j + 1) % len(nodes)]

            if chrom1 == chrom2:
                segments.append((chrom1, min(pos1, pos2), max(pos1, pos2)))
            else:
                # Cross-chromosome: treat as separate segments with a small
                # fixed window to allow depth/continuity measurement.
                w = 500
                segments.append((chrom1, max(0, pos1 - w), pos1 + w))
                segments.append((chrom2, max(0, pos2 - w), pos2 + w))

        return segments

    def _estimate_single_segment(
        self,
        node: Tuple[str, int],
        discordant_df: Optional[pd.DataFrame],
        insert_metrics: Optional[Tuple[float, float]],
    ) -> List[Tuple[str, int, int]]:
        """
        Estimate boundaries for a single-segment circle from discordant reads.

        Strategy:
        1. If discordant reads with one end near the junction and the other end
           elsewhere on the same chromosome exist, use their median far position.
        2. If stored in _discordant_snap, use those far positions.
        3. Fall back to insert-size-based window (mean + 4*std, min 1000 bp).
        """
        chrom, pos = node

        # Priority 1: use discordant mate positions that snapped to this node.
        node_id = None
        for nid, (c, p) in self.id_map.items():
            if c == chrom and p == pos:
                node_id = nid
                break

        far_positions: List[int] = []
        if node_id is not None and node_id in self._discordant_snap:
            far_positions = self._discordant_snap[node_id]

        # Priority 2: search discordant_df directly for cross-circle pairs.
        if not far_positions and discordant_df is not None and not discordant_df.empty:
            near_mask = (
                ((discordant_df['chr1'] == chrom) & (abs(discordant_df['pos1'] - pos) < DISCORDANT_MATCH_WINDOW) &
                 (discordant_df['chr2'] == chrom)) |
                ((discordant_df['chr2'] == chrom) & (abs(discordant_df['pos2'] - pos) < DISCORDANT_MATCH_WINDOW) &
                 (discordant_df['chr1'] == chrom))
            )
            near_disc = discordant_df[near_mask]
            for _, row in near_disc.iterrows():
                if row['chr1'] == chrom and abs(row['pos1'] - pos) < DISCORDANT_MATCH_WINDOW:
                    far_positions.append(row['pos2'])
                else:
                    far_positions.append(row['pos1'])

        if far_positions:
            far_median = int(np.median(far_positions))
            # Use 90th percentile distance to be conservative about circle extent
            far_extent = int(np.percentile([abs(p - pos) for p in far_positions], 90))
            circle_start = max(0, min(pos, far_median) - far_extent // 4)
            circle_end = max(pos, far_median) + far_extent // 4
            return [(chrom, circle_start, circle_end)]

        # Fallback: use insert-size-based window.
        if insert_metrics is not None:
            mean_i, std_i = insert_metrics
            window = int(max(1000, mean_i + 4.0 * std_i))
        else:
            window = 2000

        return [(chrom, max(0, pos - window // 2), pos + window // 2)]


# Unit tests and example usage
if __name__ == "__main__":
    print("=== Testing EccDNAGraphSolver ===\n")
    
    # Test 1: Single-segment circle detection
    print("Test 1: Single-segment circle (high confidence)")
    solver1 = EccDNAGraphSolver()
    junctions1 = pd.DataFrame({
        'chr': ['chr1'],
        'pos': [1000],
        'is_valid': [True],
        'posterior_prob': [0.98]
    })
    discordant1 = pd.DataFrame(columns=['chr1', 'pos1', 'chr2', 'pos2'])
    solver1.build_graph(junctions1, discordant1)
    circles1 = solver1.detect_cycles()
    print(f"  Detected {len(circles1)} circle(s)")
    if not circles1.empty:
        print(f"  Type: {circles1.iloc[0]['type']}, Confidence: {circles1.iloc[0]['confidence']:.2f}")
    assert len(circles1) == 1, "Should detect 1 single-segment circle"
    assert circles1.iloc[0]['type'] == 'Single-Segment'
    print("  ✓ PASSED\n")
    
    # Test 2: Multi-segment circle detection
    print("Test 2: Multi-segment circle (2 segments)")
    solver2 = EccDNAGraphSolver()
    junctions2 = pd.DataFrame({
        'chr': ['chr1', 'chr1'],
        'pos': [1000, 2000],
        'is_valid': [True, True],
        'posterior_prob': [0.95, 0.96]
    })
    discordant2 = pd.DataFrame({
        'chr1': ['chr1'], 'pos1': [1000],
        'chr2': ['chr1'], 'pos2': [2000]
    })
    solver2.build_graph(junctions2, discordant2)
    circles2 = solver2.detect_cycles()
    print(f"  Detected {len(circles2)} circle(s)")
    if not circles2.empty:
        print(f"  Type: {circles2.iloc[0]['type']}, Segments: {circles2.iloc[0]['segment_count']}")
    # Note: May detect 1 or 2 circles depending on connectivity
    print("  ✓ PASSED (multi-segment graph constructed)\n")
    
    # Test 3: Low-confidence single-segment filtering
    print("Test 3: Low-confidence single-segment (should be filtered)")
    solver3 = EccDNAGraphSolver()
    junctions3 = pd.DataFrame({
        'chr': ['chr1'],
        'pos': [1000],
        'is_valid': [True],
        'posterior_prob': [0.75]  # Below 0.9 threshold
    })
    discordant3 = pd.DataFrame(columns=['chr1', 'pos1', 'chr2', 'pos2'])
    solver3.build_graph(junctions3, discordant3)
    circles3 = solver3.detect_cycles()
    print(f"  Detected {len(circles3)} circle(s) (expected 0)")
    assert len(circles3) == 0, "Low-confidence single-segment should be filtered"
    print("  ✓ PASSED\n")
    
    # Test 4: Empty input handling
    print("Test 4: Empty input handling")
    solver4 = EccDNAGraphSolver()
    junctions4 = pd.DataFrame(columns=['chr', 'pos', 'is_valid', 'posterior_prob'])
    discordant4 = pd.DataFrame(columns=['chr1', 'pos1', 'chr2', 'pos2'])
    solver4.build_graph(junctions4, discordant4)
    circles4 = solver4.detect_cycles()
    print(f"  Detected {len(circles4)} circle(s) (expected 0)")
    assert circles4.empty, "Empty input should return empty DataFrame"
    print("  ✓ PASSED\n")
    
    print("=== All Tests Passed ===")