import argparse
import sys
from extractor import ReadExtractor
from realigner import ProbabilisticRealigner
from graph_solver import EccDNAGraphSolver
from config import OUTPUT_DIR
import pandas as pd

def run_pipeline(ref: str, fastq_r1: str, fastq_r2: str, bam_path: str = None):
    """
    Execute the complete CircleDetect pipeline for ecDNA detection.
    
    Pipeline stages:
    1. Read alignment and extraction (split reads + discordant pairs)
    2. Probabilistic realignment of junction candidates (HMM-based)
    3. Graph construction from validated junctions
    4. Cycle detection with advanced confidence scoring
    5. Output generation with comprehensive metrics
    
    Args:
        ref (str): Path to reference FASTA file
        fastq_r1 (str): Path to Read 1 FASTQ file
        fastq_r2 (str): Path to Read 2 FASTQ file
        bam_path (str, optional): Path to sorted BAM for advanced scoring
        
    Returns:
        pd.DataFrame: Detected ecDNA circles with confidence scores
    """
    print(f"=== Starting EccDNA Detection Pipeline (Paired-End) ===")
    print(f"Reference: {ref}")
    print(f"Read 1: {fastq_r1}")
    print(f"Read 2: {fastq_r2}")
    
    # 1. Extraction (Handles PE alignment internally)
    extractor = ReadExtractor(ref, fastq_r1, fastq_r2)
    split_df, discord_df = extractor.align_and_extract()
    
    if split_df.empty:
        print("No split reads found. Exiting.")
        return pd.DataFrame()

    # 2. Probabilistic Realignment (Matrix-Based HMM)
    print(">> Running Matrix-Based Probabilistic Realignment (with Gap Penalties)...")
    realigner = ProbabilisticRealigner(ref)
    validated_df = realigner.realign_junctions(split_df)
    
    valid_count = validated_df['is_valid'].sum()
    print(f"Validated {valid_count} high-confidence junctions.")
    
    if valid_count == 0:
        print("No high-confidence junctions found.")
        return pd.DataFrame()

    # 3. Graph Construction & Cycle Detection
    print(">> Building Breakpoint Graph and Detecting Cycles...")
    solver = EccDNAGraphSolver()
    solver.build_graph(validated_df, discord_df)
    
    # Use BAM path for advanced scoring if available
    # If not provided, use the aligned BAM from extraction step
    if bam_path is None:
        bam_path = extractor.bam_path
    
    circles = solver.detect_cycles(bam_path=bam_path)
    
    # 4. Output
    output_file = f"{OUTPUT_DIR}/eccdna_calls.csv"
    circles.to_csv(output_file, index=False)
    print(f"=== Detection Complete. Results saved to {output_file} ===")
    if not circles.empty:
        print("\nDetected ecDNA circles:")
        print(circles.to_string())
        
        # Print summary statistics
        print(f"\nSummary:")
        print(f"  Total circles detected: {len(circles)}")
        print(f"  Single-segment: {(circles['type'] == 'Single-Segment').sum()}")
        print(f"  Multi-segment: {(circles['type'] == 'Multi-Segment').sum()}")
        
        if 'integrated_confidence' in circles.columns and circles['integrated_confidence'].notna().any():
            high_conf = (circles['integrated_confidence'] >= 0.9).sum()
            print(f"  High confidence (≥0.9): {high_conf}")
            
        if 'depth_ratio' in circles.columns and circles['depth_ratio'].notna().any():
            mean_ratio = circles['depth_ratio'].mean()
            print(f"  Mean depth ratio: {mean_ratio:.2f}x")
    else:
        print("No circular DNA structures detected.")
    
    return circles

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="High-Performance Multi-Segment ecDNA Detector")
    parser.add_argument("-r", "--reference", required=True, help="Path to reference FASTA")
    parser.add_argument("-1", "--fastq1", required=True, help="Path to Read 1 (FASTQ)")
    parser.add_argument("-2", "--fastq2", required=True, help="Path to Read 2 (FASTQ)")
    parser.add_argument("-b", "--bam", required=False, help="Path to sorted BAM (optional, for advanced scoring)")
    
    args = parser.parse_args()
    run_pipeline(args.reference, args.fastq1, args.fastq2, args.bam)