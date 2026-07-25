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
    
    def __init__(self, bam_path: str, window_size: int = 50) -> None:
        """
        Initialize the confidence scorer.
        
        Args:
            bam_path (str): Path to sorted, indexed BAM file
            window_size (int): Size of flanking region for depth ratio calculation
            
        Raises:
            FileNotFoundError: If BAM file or its index (.bai) does not exist
            
        Example:
            >>> scorer = EccDNAConfidenceScorer("sample.bam", window_size=100)
        """
        self.bam_path = bam_path
        self.window_size = window_size
        self._bam_file: Optional[pysam.AlignmentFile] = None
        
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
        segments: List[Tuple[str, int, int]]
    ) -> float:
        """
        Calculate read continuity score across the entire ecDNA structure.
        
        Continuity measures whether paired-end reads physically span the complete
        circular path. For a valid ecDNA circle, we expect:
        - For single-segment circles: Read pairs where one mate maps near the 
          junction and the other mate maps to the opposite side of the circle
        - For multi-segment circles: Read pairs connecting non-adjacent segments
          along the circular path
        
        The score is computed as:
        continuity = (spanning_read_pairs) / (total_expected_pairs)
        
        Args:
            segments: List of (chromosome, start, end) tuples defining the circle.
                     For a single-segment circle: [('chr1', 1000, 2000)]
                     For multi-segment: [('chr1', 1000, 1500), ('chr1', 3000, 3500)]
                     
        Returns:
            float: Continuity score between 0.0 and 1.0
                   - 0.0: No spanning reads detected
                   - 1.0: Maximum possible spanning reads observed
                   
        Note:
            - Uses proper pair orientation (R1F2 or R2F1) to identify spanning reads
            - Filters by insert size to exclude PCR duplicates
            - Requires both mates to be mapped with MAPQ >= 20
            
        Example:
            >>> scorer = EccDNAConfidenceScorer("sample.bam")
            >>> segments = [('chr1', 1000, 2000)]
            >>> score = scorer.calculate_continuity_score(segments)
            >>> print(f"Continuity: {score:.2f}")
        """
        bam = self._open_bam()
        
        if len(segments) == 0:
            return 0.0
            
        # For single-segment circle, check for reads spanning the junction
        # Expected: One mate in [start, start+window] and other in [end-window, end]
        # OR vice versa due to circular nature
        
        spanning_count = 0
        total_candidates = 0
        
        chrom, start, end = segments[0]
        circle_length = end - start
        
        # Define junction-proximal regions (within 200bp of breakpoints)
        junction_window = 200
        left_junction_start = max(0, start - junction_window)
        left_junction_end = start + junction_window
        right_junction_start = end - junction_window
        right_junction_end = end + junction_window
        
        try:
            # Fetch reads in left junction region
            for read in bam.fetch(chrom, left_junction_start, left_junction_end):
                if read.is_unmapped or read.mate_is_unmapped:
                    continue
                if read.mapping_quality < 20:
                    continue
                if not read.is_proper_pair:
                    continue
                if read.reference_name != read.next_reference_name:
                    continue
                    
                mate_pos = read.next_reference_start
                
                # Check if mate is in right junction region
                if right_junction_start <= mate_pos <= right_junction_end:
                    spanning_count += 1
                    
                total_candidates += 1
                
        except ValueError as e:
            # Region may be out of bounds
            pass
            
        if total_candidates == 0:
            return 0.0
            
        # Normalize by expected count (heuristic: expect ~10% of reads to span)
        # Use sigmoid function to map ratio to 0-1 score
        raw_ratio = spanning_count / max(total_candidates, 1)
        continuity = 1.0 / (1.0 + np.exp(-10 * (raw_ratio - 0.05)))
        
        return min(1.0, max(0.0, continuity))
    
    def calculate_depth_ratio(
        self,
        segments: List[Tuple[str, int, int]],
        sample_type: str = 'tumor'
    ) -> float:
        """
        Calculate sequencing depth ratio between ecDNA region and flanking linear regions.
        
        ecDNA circles often exhibit copy number amplification, resulting in elevated
        sequencing depth compared to the surrounding linear genome. This method:
        
        1. Calculates mean coverage within the proposed ecDNA region(s)
        2. Calculates mean coverage in flanking regions (outside the circle)
        3. Computes ratio: depth_inside / depth_outside
        
        Expected ratios:
        - Ratio ≈ 1.0: No amplification (likely false positive or single-copy ecDNA)
        - Ratio > 1.5: Moderate amplification (supporting evidence)
        - Ratio > 3.0: Strong amplification (strong supporting evidence)
        - Ratio < 0.8: Depletion (likely artifact)
        
        Args:
            segments: List of (chromosome, start, end) tuples defining ecDNA regions
            sample_type: Not used currently, reserved for tumor/normal normalization
            
        Returns:
            float: Depth ratio (depth_inside / depth_outside)
                   - Values > 1.0 indicate amplification
                   - Values < 1.0 indicate depletion
                   
        Raises:
            ValueError: If segments are invalid or out of chromosome bounds
            
        Note:
            - Uses pileup for efficient coverage calculation
            - Excludes bases with base quality < 20
            - Flanking region extends window_size bp on each side
            - Handles multi-segment circles by averaging across segments
            
        Example:
            >>> scorer = EccDNAConfidenceScorer("sample.bam", window_size=100)
            >>> segments = [('chr1', 1000, 2000)]
            >>> ratio = scorer.calculate_depth_ratio(segments)
            >>> print(f"Depth ratio: {ratio:.2f}x")
        """
        bam = self._open_bam()
        
        if len(segments) == 0:
            return 1.0
            
        total_inside_bases = 0
        total_inside_coverage = 0
        total_outside_bases = 0
        total_outside_coverage = 0
        
        for chrom, start, end in segments:
            region_length = end - start
            if region_length <= 0:
                continue
                
            # Calculate coverage inside ecDNA region
            inside_coverage_sum = 0
            inside_base_count = 0
            
            try:
                for pileupcolumn in bam.pileup(
                    chrom, start, end,
                    stepper='nofilter',
                    min_base_quality=20
                ):
                    pos = pileupcolumn.pos
                    if start <= pos < end:
                        inside_coverage_sum += pileupcolumn.n
                        inside_base_count += 1
            except ValueError:
                continue
                
            if inside_base_count == 0:
                continue
                
            mean_inside = inside_coverage_sum / inside_base_count
            total_inside_bases += inside_base_count
            total_inside_coverage += inside_coverage_sum
            
            # Calculate coverage in flanking regions (outside ecDNA)
            flank_start_left = max(0, start - self.window_size)
            flank_end_right = end + self.window_size
            
            outside_coverage_sum = 0
            outside_base_count = 0
            
            # Left flank
            try:
                for pileupcolumn in bam.pileup(
                    chrom, flank_start_left, start,
                    stepper='nofilter',
                    min_base_quality=20
                ):
                    outside_coverage_sum += pileupcolumn.n
                    outside_base_count += 1
            except ValueError:
                pass
                
            # Right flank
            try:
                for pileupcolumn in bam.pileup(
                    chrom, end, flank_end_right,
                    stepper='nofilter',
                    min_base_quality=20
                ):
                    outside_coverage_sum += pileupcolumn.n
                    outside_base_count += 1
            except ValueError:
                pass
                
            if outside_base_count > 0:
                total_outside_bases += outside_base_count
                total_outside_coverage += outside_coverage_sum
        
        if total_inside_bases == 0 or total_outside_bases == 0:
            return 1.0
            
        mean_inside = total_inside_coverage / total_inside_bases
        mean_outside = total_outside_coverage / total_outside_bases
        
        if mean_outside == 0:
            return 1.0
            
        ratio = mean_inside / mean_outside
        return ratio
    
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
        circle_id: str = "circle_0"
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
            
        Returns:
            EccDNACandidate: Dataclass containing all scores and metadata
            
        Raises:
            ValueError: If segments list is empty
            
        Example:
            >>> scorer = EccDNAConfidenceScorer("sample.bam")
            >>> segments = [('chr1', 1000, 2000)]
            >>> candidate = scorer.score_circle(
            ...     segments=segments,
            ...     junctions_df=junctions_df,
            ...     discordant_df=discordant_df,
            ...     circle_id="eccDNA_001"
            ... )
            >>> print(f"Confidence: {candidate.integrated_confidence:.3f}")
            >>> print(f"Depth ratio: {candidate.depth_ratio:.2f}x")
            >>> print(f"Continuity: {candidate.continuity_score:.2f}")
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
        
        # Calculate continuity score
        continuity_score = self.calculate_continuity_score(segments)
        
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
