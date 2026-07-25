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
import os
from typing import Tuple, List, Dict, Optional
from config import (
    BWA_MEM2_PATH, SAMBLASTER_PATH, SAMTOOLS_PATH,
    MIN_MAPQ, MIN_CLIP_LEN, MAX_INSERT_SIZE, OUTPUT_DIR
)


class ReadExtractor:
    """
    Extract candidate ecDNA reads from paired-end FASTQ files.
    
    This class orchestrates the alignment and extraction pipeline:
    1. Aligns reads to reference genome using bwa-mem2
    2. Sorts alignments by query name (required for samblaster)
    3. Extracts split reads (soft-clipped) and discordant read pairs
    4. Parses BAM files to create DataFrames for downstream processing
    
    Attributes:
        ref (str): Path to reference FASTA file
        fastq_r1 (str): Path to Read 1 FASTQ file
        fastq_r2 (str): Path to Read 2 FASTQ file
        bam_path (str): Path to output sorted BAM file
        
    Example:
        >>> extractor = ReadExtractor("ref.fasta", "reads_1.fastq", "reads_2.fastq")
        >>> split_df, discord_df = extractor.align_and_extract()
        >>> print(f"Found {len(split_df)} split reads, {len(discord_df)} discordant pairs")
    """
    
    def __init__(self, ref_fasta: str, fastq_r1: str, fastq_r2: str) -> None:
        """
        Initialize the ReadExtractor with input files.
        
        Args:
            ref_fasta (str): Path to reference genome FASTA file (must be indexed)
            fastq_r1 (str): Path to Read 1 FASTQ file (paired-end)
            fastq_r2 (str): Path to Read 2 FASTQ file (paired-end)
            
        Raises:
            FileNotFoundError: If any input file does not exist
            
        Example:
            >>> extractor = ReadExtractor("/path/to/ref.fasta", "r1.fastq", "r2.fastq")
        """
        self.ref = ref_fasta
        self.fastq_r1 = fastq_r1
        self.fastq_r2 = fastq_r2
        self.bam_path = os.path.join(OUTPUT_DIR, "aligned_sorted.bam")
        
    def align_and_extract(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Execute the complete read alignment and extraction pipeline.
        
        Pipeline Steps:
        1. **Alignment**: Map paired-end reads to reference using bwa-mem2
           - Uses `-Y` flag for soft-clipping (important for junction detection)
           - Uses `-k 19` minimum seed length for sensitivity
        2. **Sorting**: Name-sort BAM for samblaster processing
        3. **Extraction**: Use samblaster to split into:
           - Split reads (soft-clipped, potential junction evidence)
           - Discordant reads (abnormal insert size or orientation)
        4. **Parsing**: Convert BAM files to pandas DataFrames
        
        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]: 
                - split_df: DataFrame of split reads with columns:
                  ['read_name', 'chr', 'pos', 'strand', 'seq', 'qual', 'type', 
                   'mate_chr', 'mate_pos', 'is_proper_pair', 'is_same_chr']
                - discord_df: DataFrame of discordant pairs with columns:
                  ['read_name', 'chr1', 'pos1', 'chr2', 'pos2', 'strand1', 'strand2']
                  
        Raises:
            subprocess.CalledProcessError: If alignment or extraction fails
            
        Example:
            >>> extractor = ReadExtractor("ref.fa", "r1.fq", "r2.fq")
            >>> split_df, discord_df = extractor.align_and_extract()
            >>> print(f"Split reads: {len(split_df)}, Discordant: {len(discord_df)}")
        """
        print(">> Step 1: Aligning Paired-End reads with bwa-mem2...")
        # BWA-MEM2 accepts two FASTQ files for PE data
        # -M: Mark shorter split hits as secondary (for Picard compatibility)
        # -Y: Use soft-clipping for supplementary alignments
        # -k 19: Minimum seed length (lower = more sensitive)
        bwa_cmd = [
            BWA_MEM2_PATH, "mem", "-t", "8", "-M", "-Y", "-k", "19",
            self.ref, self.fastq_r1, self.fastq_r2
        ]
        
        print(">> Step 2: Sorting and Extracting with samblaster...")
        split_out = os.path.join(OUTPUT_DIR, "split_reads.bam")
        discord_out = os.path.join(OUTPUT_DIR, "discordant_reads.bam")
        
        # Pipe: bwa -> samtools sort (name) -> samblaster
        # Name sorting is required for samblaster to properly identify read pairs
        full_cmd = (
            f"{' '.join(bwa_cmd)} | "
            f"{SAMTOOLS_PATH} sort -n -@ 8 - | "
            f"{SAMBLASTER_PATH} "
            f"--splitFile {split_out} "
            f"--discordantFile {discord_out} "
            f"-o /dev/null"
        )
        
        # Execute pipeline
        subprocess.run(full_cmd, shell=True, check=True, executable='/bin/bash')
        
        print(">> Step 3: Parsing Split Reads (Junction Candidates)...")
        split_df = self._parse_split_reads(split_out)
        
        print(">> Step 4: Parsing Discordant Reads (Supporters)...")
        discord_df = self._parse_discordant_reads(discord_out)
        
        return split_df, discord_df

    def _parse_split_reads(self, bam_path: str) -> pd.DataFrame:
        """
        Parse BAM file to extract soft-clipped reads indicating potential junction breakpoints.
        
        Soft-clipped reads occur when part of a read aligns to the reference while another
        part (the clipped portion) does not. In ecDNA detection, these often indicate
        circular junction points where the read spans the backsplice junction.
        
        Args:
            bam_path (str): Path to sorted BAM file containing split reads
            
        Returns:
            pd.DataFrame: DataFrame with columns:
                - 'read_name': Read identifier
                - 'chr': Chromosome of aligned portion
                - 'pos': Genomic position at clip junction
                - 'strand': Strand orientation ('+' or '-')
                - 'seq': Clipped sequence (for realignment)
                - 'qual': Quality scores as Phred string (for realignment)
                - 'type': Clip type ('head' or 'tail')
                - 'mate_chr': Mate's chromosome
                - 'mate_pos': Mate's alignment position
                - 'is_proper_pair': Whether read is in proper pair
                - 'is_same_chr': Whether mate is on same chromosome
                
        Note:
            - Filters by MIN_MAPQ and MIN_CLIP_LEN thresholds from config
            - Preserves sequence and quality for downstream probabilistic realignment
        """
        records: List[Dict[str, any]] = []
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                # Filter by mapping quality
                if read.is_unmapped or read.mapping_quality < MIN_MAPQ:
                    continue
                
                # Extract soft-clipped regions
                clips = self._get_soft_clips(read)
                if not clips:
                    continue
                
                for clip_pos, clip_seq, clip_qual, clip_type in clips:
                    # Filter by minimum clip length
                    if len(clip_seq) < MIN_CLIP_LEN:
                        continue
                        
                    records.append({
                        'read_name': read.query_name,
                        'chr': read.reference_name,
                        'pos': clip_pos,
                        'strand': '-' if read.is_reverse else '+',
                        'seq': clip_seq,  # Preserve for realignment
                        'qual': clip_qual,  # Preserve quality scores
                        'type': clip_type,
                        'mate_chr': read.next_reference_name,
                        'mate_pos': read.next_reference_start,
                        'is_proper_pair': not read.is_paired or read.is_proper_pair
                    })
        
        df = pd.DataFrame(records)
        if df.empty:
            return df
            
        # Add same-chromosome flag for single-segment circle hypothesis
        df['is_same_chr'] = df['chr'] == df['mate_chr']
        return df

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
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                # Skip if either mate is unmapped
                if read.is_unmapped or read.mate_is_unmapped:
                    continue
                    
                if read.reference_name == read.next_reference_name:
                    # Same chromosome: check insert size
                    dist = abs(read.reference_start - read.next_reference_start)
                    if dist > MAX_INSERT_SIZE:
                        records.append({
                            'read_name': read.query_name,
                            'chr1': read.reference_name,
                            'pos1': read.reference_start,
                            'chr2': read.next_reference_name,
                            'pos2': read.next_reference_start,
                            'strand1': '-' if read.is_reverse else '+',
                            'strand2': '-' if read.mate_is_reverse else '+'
                        })
                else:
                    # Inter-chromosomal: always counts as discordant
                    records.append({
                        'read_name': read.query_name,
                        'chr1': read.reference_name,
                        'pos1': read.reference_start,
                        'chr2': read.next_reference_name,
                        'pos2': read.next_reference_start,
                        'strand1': '-' if read.is_reverse else '+',
                        'strand2': '-' if read.mate_is_reverse else '+'
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
            
            # Update positions based on CIGAR operation
            if op_char in "M=X":
                pos += length  # Advance reference position
            if op_char in "MIS=X":
                s_pos += length  # Advance read position
                
        return clips
                
        return clips