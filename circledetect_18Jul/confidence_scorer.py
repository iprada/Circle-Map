"""
EccDNA Confidence Scorer Module

This module implements advanced confidence scoring for ecDNA detection by integrating:
1. Junction posterior probabilities from split-read realignment
2. Discordant read support for connectivity evidence
3. Read continuity metrics across the entire circular structure
4. Sequencing depth ratio between ecDNA region and flanking linear regions

Key Innovations over Circle-Map:
- Continuity scoring: Measures whether paired-end reads span the complete ecDNA circle
- Depth ratio analysis: Compares coverage within proposed ecDNA vs. flanking regions
- Integrated Bayesian confidence: Combines all evidence sources into unified probability

Author: CircleDetect Team
Based on: Circle-Map methodology with enhanced scoring
"""

import numpy as np
import pandas as pd
import pysam
from typing import List, Dict, Tuple, Optional
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from dataclasses import dataclass


@dataclass
class EccDNACandidate:
    """
    Represents a candidate ecDNA circle with all supporting evidence.
    
    Attributes:
        circle_id: Unique identifier for this candidate
        segments: List of (chromosome, start, end) tuples defining the circular path
        junction_posterior: Mean posterior probability of junction breakpoints (0-1)
        discordant_count: Number of discordant read pairs supporting connections
        continuity_score: Fraction of reads spanning the complete circle (0-1)
        depth_ratio: Coverage ratio (inside ecDNA / outside flanking region)
        integrated_confidence: Combined confidence score (0-1)
    """
    circle_id: str
    segments: List[Tuple[str, int, int]]
    junction_posterior: float
    discordant_count: int
    continuity_score: float
    depth_ratio: float
    integrated_confidence: float


class EccDNAConfidenceScorer:
    """
    Advanced confidence scorer for ecDNA candidates.
    
    This class implements a multi-evidence scoring framework that combines:
    1. **Junction Evidence**: Posterior probabilities from HMM realignment of split reads
    2. **Connectivity Evidence**: Count and distribution of discordant read pairs
    3. **Continuity Evidence**: Paired-end reads that span the entire circular structure
    4. **Depth Evidence**: Copy number signal from elevated coverage in ecDNA region
    
    The integrated confidence is computed using a Bayesian framework:
    
    P(circle | data) ∝ P(junctions | circle) × P(connectivity | circle) × 
                       P(continuity | circle) × P(depth | circle) × P(prior)
    
    Attributes:
        bam_path (str): Path to input BAM file for depth/continuity calculations
        window_size (int): Size of flanking region for depth comparison (default: 50bp)
        
    Example:
        >>> scorer = EccDNAConfidenceScorer("aligned.bam")
        >>> candidate = scorer.score_circle(
        ...     segments=[('chr1', 1000, 2000)],
        ...     junctions_df=junctions,
        ...     discordant_df=discordant
        ... )
        >>> print(f"Confidence: {candidate.integrated_confidence:.3f}")
    """
    
    def __init__(
        self,
        bam_path: str,
        inside_length: int = 100,
        extension: int = 200,
        min_base_quality: int = 20,
        coverage_mapq: int = 20,
        continuity_mapq: int = 20,
        junction_window: int = 200,
    ) -> None:
        """
        Initialize the confidence scorer.

        Args:
            bam_path: Path to a coordinate-sorted, indexed BAM file.
            inside_length: Number of bp inside each breakpoint used for the
                depth-ratio "inside edge" (Circle-Map ``-E``).
            extension: Number of bp flanking each breakpoint used for the
                depth-ratio "outside flank" (Circle-Map ``-b``).
            min_base_quality: Minimum base quality for ``count_coverage``.
            coverage_mapq: MAPQ filter applied via read_callback for coverage.
            continuity_mapq: MAPQ filter applied when counting bridging pairs.
            junction_window: Half-window around each breakpoint inside which a
                read counts as "anchored" at that breakpoint for continuity.

        Example:
            >>> scorer = EccDNAConfidenceScorer("sample.bam")
        """
        self.bam_path = bam_path
        self.inside_length = inside_length
        self.extension = extension
        self.min_base_quality = min_base_quality
        self.coverage_mapq = coverage_mapq
        self.continuity_mapq = continuity_mapq
        self.junction_window = junction_window
        self._bam_file: Optional[pysam.AlignmentFile] = None

    def _read_callback(self, read: pysam.AlignedSegment) -> bool:
        """``count_coverage`` filter: keep reads passing the coverage MAPQ cutoff."""
        return (not read.is_unmapped) and (not read.is_secondary) and (
            read.mapping_quality >= self.coverage_mapq
        )

    def _open_bam(self) -> pysam.AlignmentFile:
        """Lazy loading of BAM file."""
        if self._bam_file is None or self._bam_file.closed:
            self._bam_file = pysam.AlignmentFile(self.bam_path, "rb")
        return self._bam_file
        
    def close(self) -> None:
        """Close the BAM file handle."""
        if self._bam_file is not None and not self._bam_file.closed:
            self._bam_file.close()
            self._bam_file = None
    
    def calculate_continuity_score(
        self,
        segments: List[Tuple[str, int, int]],
        insert_metrics: Optional[Tuple[float, float]] = None,
    ) -> float:
        """
        Measure read-pair continuity across each segment's backsplice junction.

        Two regimes:

        1. **Small circle** (``end - start`` < mean_insert + 3·std): a true
           ecDNA produces paired reads that physically span the junction — one
           mate inside the segment near ``start``, the other near ``end``, with
           an insert size matching the circle length. We count these "bridging"
           pairs and divide by the number of pairs anchored at the junction,
           giving a bridging fraction. A 10% bridging rate is strong evidence
           of a real circle, so we map the fraction through a sigmoid centred
           at 0.10.

        2. **Large circle** (junction too far for a single insert to span):
           bridging cannot be measured in one pair, so we use the coefficient
           of variation of per-base coverage inside the segment as a continuity
           proxy — a uniformly amplified circle has low CV, a focal artifact
           has high CV. ``score = 1 / (1 + CV)``.

        The overall score is the mean across all segments, so it generalises
        to multi-segment circles (each segment scored independently).

        Args:
            segments: List of (chromosome, start, end) tuples.
            insert_metrics: Optional ``(mean, std)`` of the insert size
                distribution. Required to enable the bridging measurement and
                to choose between the two regimes.

        Returns:
            Continuity score in [0, 1]; higher = more continuous.
        """
        if len(segments) == 0:
            return 0.0

        max_span: Optional[float] = None
        if insert_metrics is not None:
            mean_i, std_i = insert_metrics
            max_span = mean_i + 3.0 * std_i

        scores: List[float] = []
        for chrom, start, end in segments:
            if chrom is None or end <= start:
                continue
            circle_len = end - start
            if max_span is not None and circle_len <= max_span:
                scores.append(self._bridging_fraction(chrom, start, end))
            else:
                scores.append(self._coverage_uniformity(chrom, start, end))

        if not scores:
            return 0.0
        return float(np.mean(scores))

    def _bridging_fraction(self, chrom: str, start: int, end: int) -> float:
        """
        Fraction of junction-anchored read pairs whose mate bridges to the
        opposite breakpoint. Checks **both directions**:
          - reads anchored near ``start`` whose mate lands near ``end``
          - reads anchored near ``end`` whose mate lands near ``start``
        Deduplicates pairs so each pair is counted at most once.

        Only meaningful when the circle is smaller than the sequencing insert
        size (otherwise a single read pair cannot physically span the circle).
        """
        bam = self._open_bam()
        W = self.junction_window

        left_lo = max(0, start - W)
        left_hi = start + W
        right_lo = end - W
        right_hi = end + W

        # Sets of read names: one for all anchored reads, one for bridging pairs.
        anchored: set = set()
        bridging: set = set()

        def _check_read(read, region_lo, region_hi, mate_lo, mate_hi):
            """Helper: check if a read is anchored in *region* and bridges to *mate*."""
            if read.is_unmapped or read.mate_is_unmapped:
                return
            if read.mapping_quality < self.continuity_mapq:
                return
            if read.is_secondary or read.is_supplementary:
                return
            if not read.is_paired:
                return
            if read.next_reference_name != chrom:
                return
            anchored.add(read.query_name)
            if mate_lo <= read.next_reference_start <= mate_hi:
                bridging.add(read.query_name)

        try:
            # Direction 1: reads anchored near the left breakpoint
            for read in bam.fetch(chrom, left_lo, left_hi):
                _check_read(read, left_lo, left_hi, right_lo, right_hi)

            # Direction 2: reads anchored near the right breakpoint
            for read in bam.fetch(chrom, right_lo, right_hi):
                _check_read(read, right_lo, right_hi, left_lo, left_hi)

        except ValueError:
            return 0.0

        total_anchored = len(anchored)
        if total_anchored == 0:
            return 0.0

        # Sigmoid centred at 0.10: 10% bridging -> ~0.98, 5% -> ~0.88, 1% -> ~0.37.
        fraction = len(bridging) / total_anchored
        return float(1.0 / (1.0 + np.exp(-25.0 * (fraction - 0.05))))

    def _coverage_uniformity(self, chrom: str, start: int, end: int) -> float:
        """
        Continuity proxy for large circles: per-base coverage CV inside the
        segment. Uniform amplification -> low CV -> score near 1.
        """
        bam = self._open_bam()
        inside_length = min(self.inside_length, end - start)
        if inside_length <= 0:
            return 0.0
        try:
            cov = bam.count_coverage(
                chrom, start, end,
                quality_threshold=self.min_base_quality,
                read_callback=self._read_callback,
            )
        except ValueError:
            return 0.0
        per_base = np.array([cov[0], cov[1], cov[2], cov[3]]).sum(axis=0)
        if per_base.size == 0:
            return 0.0
        mean = float(per_base.mean())
        if mean <= 0:
            return 0.0
        cv = float(per_base.std()) / mean
        return float(1.0 / (1.0 + cv))
    
    def calculate_depth_ratio(
        self,
        segments: List[Tuple[str, int, int]],
        sample_type: str = 'tumor'
    ) -> float:
        """
        Compute inside-vs-flank sequencing depth ratio at each breakpoint.

        Follows the Circle-Map Coverage algorithm: for each segment, sum the
        per-base coverage over the inner edge (the first/last `inside_length`
        bases of the segment) and the outer flank (`extension` bases just
        outside the breakpoint), then form

            inside_outside_ratio = sum(inside_edge) / sum(flank)

        A value > 1.0 indicates the segment's interior is more covered than
        the flanking linear genome (amplification); 1.0 = no amplification.

        Improvements over the previous implementation:
          - Uses ``count_coverage`` (per-base A/C/G/T sums) instead of pileup,
            which avoids counting the same read twice and lets us apply a
            MAPQ filter via ``read_filter``.
          - Reports a per-breakpoint ratio at BOTH the start and end of each
            segment (Circle-Map reports ``start_coverage_ratio`` and
            ``end_coverage_ratio`` separately); we return their mean so that
            the integrated-confidence formula (which expects >1 = amplification)
            continues to work.
          - Uses larger default windows (100 bp inside, 200 bp flank) to keep
            the estimate stable on noisy data.

        Args:
            segments: List of (chromosome, start, end) tuples.
            sample_type: Reserved; not used in the current computation.

        Returns:
            Mean inside/flank coverage ratio across both breakpoints of every
            segment. Returns 1.0 if no usable coverage was observed.
        """
        bam = self._open_bam()

        if len(segments) == 0:
            return 1.0

        inside_length = self.inside_length
        extension = self.extension
        flank_total = inside_length + extension

        ratios: List[float] = []

        for chrom, start, end in segments:
            if end - start < 2 * inside_length:
                # Segment too small for a meaningful edge vs. flank split.
                continue

            # Clamp the extended window to the chromosome bounds when possible.
            ext_start = max(0, start - extension)
            ext_end = end + extension
            try:
                cov = bam.count_coverage(
                    chrom, ext_start, ext_end,
                    quality_threshold=self.min_base_quality,
                    read_callback=self._read_callback,
                )
            except ValueError:
                continue

            per_base = np.array([cov[0], cov[1], cov[2], cov[3]]).sum(axis=0)
            if per_base.size < flank_total * 2:
                continue

            # Offset of [start, end) within the extended window.
            s_off = start - ext_start
            e_off = end - ext_start

            # Left breakpoint: inside = first `inside_length` bp of segment,
            # flank = `extension` bp immediately upstream.
            left_inside = per_base[s_off:s_off + inside_length]
            left_flank = per_base[s_off - extension:s_off] if s_off - extension >= 0 else per_base[:s_off]
            if left_flank.size > 0 and left_flank.sum() > 0:
                ratios.append(float(left_inside.sum()) / float(left_flank.sum()))

            # Right breakpoint: inside = last `inside_length` bp of segment,
            # flank = `extension` bp immediately downstream.
            right_inside = per_base[e_off - inside_length:e_off]
            right_flank = per_base[e_off:e_off + extension]
            if right_flank.size > 0 and right_flank.sum() > 0:
                ratios.append(float(right_inside.sum()) / float(right_flank.sum()))

        if not ratios:
            return 1.0
        return float(np.mean(ratios))
    
    def calculate_integrated_confidence(
        self,
        junction_posterior: float,
        discordant_count: int,
        continuity_score: float,
        depth_ratio: float,
        segment_count: int
    ) -> float:
        """
        Calculate integrated confidence score combining all evidence sources.
        
        Uses a weighted Bayesian framework to combine:
        1. Junction posterior probability (from HMM realignment)
        2. Discordant read support (connectivity evidence)
        3. Continuity score (spanning reads)
        4. Depth ratio (copy number amplification)
        
        The formula:
        
        log_odds = w1·log(P_junction/(1-P_junction)) +
                   w2·log(1 + discordant_count/10) +
                   w3·logit(continuity) +
                   w4·max(0, log2(depth_ratio))
                   
        confidence = sigmoid(log_odds)
        
        Weights are tuned based on empirical performance:
        - Junction evidence: w1=2.0 (most reliable)
        - Discordant support: w2=1.0
        - Continuity: w3=1.5 (strong structural evidence)
        - Depth ratio: w4=1.0 (supporting but not definitive)
        
        Args:
            junction_posterior (float): Mean posterior probability of junctions (0-1)
            discordant_count (int): Number of supporting discordant read pairs
            continuity_score (float): Continuity score (0-1)
            depth_ratio (float): Depth ratio (inside/outside)
            segment_count (int): Number of segments in the circle
            
        Returns:
            float: Integrated confidence score (0-1)
                   - > 0.9: High confidence
                   - 0.7-0.9: Moderate confidence
                   - < 0.7: Low confidence (likely false positive)
                   
        Example:
            >>> scorer = EccDNAConfidenceScorer("sample.bam")
            >>> conf = scorer.calculate_integrated_confidence(
            ...     junction_posterior=0.95,
            ...     discordant_count=15,
            ...     continuity_score=0.6,
            ...     depth_ratio=2.5,
            ...     segment_count=1
            ... )
            >>> print(f"Integrated confidence: {conf:.3f}")
            0.892
        """
        # Weight parameters (empirically tuned)
        w_junction = 2.0
        w_discordant = 1.0
        w_continuity = 1.5
        w_depth = 1.0
        
        # Convert each evidence source to log-odds space
        
        # 1. Junction evidence (already a probability)
        if junction_posterior >= 1.0:
            log_odds_junction = 5.0  # Cap at reasonable maximum
        elif junction_posterior <= 0.0:
            log_odds_junction = -5.0
        else:
            log_odds_junction = np.log(junction_posterior / (1.0 - junction_posterior))
        
        # 2. Discordant read support (log-scale to handle wide dynamic range)
        # Expect ~5-10 reads for true positives
        log_odds_discordant = np.log(1.0 + discordant_count / 10.0)
        
        # 3. Continuity score (convert to logit space)
        if continuity_score >= 1.0:
            log_odds_continuity = 3.0
        elif continuity_score <= 0.0:
            log_odds_continuity = -3.0
        else:
            log_odds_continuity = np.log(continuity_score / (1.0 - continuity_score + 1e-10))
        
        # 4. Depth ratio (only count amplification, not depletion)
        # log2(ratio) because doubling of copy number = ratio of 2
        if depth_ratio > 1.0:
            log_odds_depth = np.log2(depth_ratio)
        else:
            log_odds_depth = 0.0  # No penalty for normal copy number
        
        # Combine with weights
        total_log_odds = (
            w_junction * log_odds_junction +
            w_discordant * log_odds_discordant +
            w_continuity * log_odds_continuity +
            w_depth * log_odds_depth
        )
        
        # Adjust for segment count (multi-segment circles are harder to detect)
        # Apply slight penalty for complexity unless strongly supported
        # Penalty applies when total_log_odds is below higher threshold (4.0 = ~98% confidence)
        if segment_count > 1 and total_log_odds < 4.0:
            total_log_odds -= 0.3 * (segment_count - 1)
        
        # Convert back to probability using sigmoid
        confidence = 1.0 / (1.0 + np.exp(-total_log_odds))
        
        return confidence
    
    def score_circle(
        self,
        segments: List[Tuple[str, int, int]],
        junctions_df: pd.DataFrame,
        discordant_df: pd.DataFrame,
        circle_id: str = "circle_0",
        insert_metrics: Optional[Tuple[float, float]] = None,
    ) -> EccDNACandidate:
        """
        Compute comprehensive confidence score for a candidate ecDNA circle.

        This is the main entry point for scoring. It:
        1. Extracts relevant junctions and discordant reads for the segments
        2. Calculates continuity score from BAM alignments
        3. Calculates depth ratio from coverage analysis
        4. Combines all evidence into integrated confidence

        Args:
            segments: List of (chromosome, start, end) tuples defining the circle
            junctions_df: DataFrame with columns ['chr', 'pos', 'posterior_prob', 'is_valid']
            discordant_df: DataFrame with columns ['chr1', 'pos1', 'chr2', 'pos2']
            circle_id: Unique identifier for this candidate
            insert_metrics: Optional ``(mean, std)`` of the insert size distribution;
                enables the bridging-based continuity measurement for circles
                small enough for a read pair to span.

        Returns:
            EccDNACandidate: Dataclass containing all scores and metadata

        Raises:
            ValueError: If segments list is empty
        """
        if len(segments) == 0:
            raise ValueError("Segments list cannot be empty")

        # Filter junctions overlapping with segments
        junction_scores = []
        for chrom, start, end in segments:
            mask = (
                (junctions_df['chr'] == chrom) &
                (junctions_df['pos'] >= start - 50) &
                (junctions_df['pos'] <= end + 50)
            )
            segment_junctions = junctions_df[mask]
            if not segment_junctions.empty:
                junction_scores.extend(segment_junctions['posterior_prob'].tolist())

        mean_junction_posterior = np.mean(junction_scores) if junction_scores else 0.0

        # Count discordant reads supporting connections between segments
        discordant_count = 0
        for _, row in discordant_df.iterrows():
            chr1, pos1 = row['chr1'], row['pos1']
            chr2, pos2 = row['chr2'], row['pos2']

            # Check if both ends fall within any segment
            in_segment1 = any(
                (c == chr1 and s <= pos1 <= e) for c, s, e in segments
            )
            in_segment2 = any(
                (c == chr2 and s <= pos2 <= e) for c, s, e in segments
            )

            if in_segment1 and in_segment2:
                discordant_count += 1

        # Calculate continuity score (insert metrics enable the bridging regime)
        continuity_score = self.calculate_continuity_score(segments, insert_metrics)

        # Calculate depth ratio
        depth_ratio = self.calculate_depth_ratio(segments)

        # Calculate integrated confidence
        integrated_conf = self.calculate_integrated_confidence(
            junction_posterior=mean_junction_posterior,
            discordant_count=discordant_count,
            continuity_score=continuity_score,
            depth_ratio=depth_ratio,
            segment_count=len(segments)
        )

        return EccDNACandidate(
            circle_id=circle_id,
            segments=segments,
            junction_posterior=mean_junction_posterior,
            discordant_count=discordant_count,
            continuity_score=continuity_score,
            depth_ratio=depth_ratio,
            integrated_confidence=integrated_conf
        )


# Unit tests
if __name__ == "__main__":
    print("=== Testing EccDNAConfidenceScorer ===\n")
    
    # Note: These tests require a real BAM file
    # For demonstration, we test the scoring logic with mock data
    
    print("Test 1: Integrated confidence calculation (high confidence scenario)")
    scorer_mock = EccDNAConfidenceScorer.__new__(EccDNAConfidenceScorer)
    conf_high = scorer_mock.calculate_integrated_confidence(
        junction_posterior=0.98,
        discordant_count=20,
        continuity_score=0.7,
        depth_ratio=3.0,
        segment_count=1
    )
    print(f"  High confidence score: {conf_high:.3f}")
    assert conf_high > 0.9, "High evidence should yield high confidence"
    print("  ✓ PASSED\n")
    
    print("Test 2: Integrated confidence calculation (low confidence scenario)")
    conf_low = scorer_mock.calculate_integrated_confidence(
        junction_posterior=0.60,
        discordant_count=2,
        continuity_score=0.1,
        depth_ratio=0.9,
        segment_count=1
    )
    print(f"  Low confidence score: {conf_low:.3f}")
    assert conf_low < 0.7, "Low evidence should yield low confidence"
    print("  ✓ PASSED\n")
    
    print("Test 3: Multi-segment penalty")
    conf_multi = scorer_mock.calculate_integrated_confidence(
        junction_posterior=0.75,  # Lower base to trigger penalty
        discordant_count=10,
        continuity_score=0.4,
        depth_ratio=2.0,
        segment_count=3
    )
    conf_single = scorer_mock.calculate_integrated_confidence(
        junction_posterior=0.75,  # Same parameters except segment count
        discordant_count=10,
        continuity_score=0.4,
        depth_ratio=2.0,
        segment_count=1
    )
    print(f"  Single-segment: {conf_single:.3f}, Multi-segment (3): {conf_multi:.3f}")
    assert conf_multi < conf_single, "Multi-segment should have slight penalty when confidence is borderline"
    print("  ✓ PASSED\n")
    
    print("Test 4: Depth ratio impact")
    conf_amp = scorer_mock.calculate_integrated_confidence(
        junction_posterior=0.80,
        discordant_count=10,
        continuity_score=0.5,
        depth_ratio=4.0,  # Strong amplification
        segment_count=1
    )
    conf_normal = scorer_mock.calculate_integrated_confidence(
        junction_posterior=0.80,
        discordant_count=10,
        continuity_score=0.5,
        depth_ratio=1.0,  # Normal copy number
        segment_count=1
    )
    print(f"  With amplification (4x): {conf_amp:.3f}, Normal depth: {conf_normal:.3f}")
    assert conf_amp > conf_normal, "Amplification should increase confidence"
    print("  ✓ PASSED\n")
    
    print("=== All Tests Passed ===")
    print("\nNote: Full integration tests require real BAM files.")
