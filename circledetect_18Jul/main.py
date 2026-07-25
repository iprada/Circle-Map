import argparse
import sys
from extractor import ReadExtractor, estimate_insert_size
from realigner import ProbabilisticRealigner
from graph_solver import EccDNAGraphSolver
from config import OUTPUT_DIR
import pandas as pd


def run_pipeline(
    ref: str,
    fastq_r1: str = None,
    fastq_r2: str = None,
    bam_path: str = None,
    split_bam: str = None,
    discord_bam: str = None,
    threads: int = 16,
    two_pass: bool = False,
):
    """
    Execute the CircleDetect pipeline for ecDNA detection.

    Two input modes:
    - **FASTQ mode**: Provide ``fastq_r1`` and ``fastq_r2``. Runs bwa-mem2 + samblaster.
    - **Pre-aligned BAM mode**: Provide ``split_bam``, ``discord_bam``, and ``bam_path``.
      Skips alignment, parses existing samblaster output.

    Pipeline stages:
    1. Read extraction (split reads + discordant pairs)
    2. Probabilistic realignment of junction candidates (matrix-based)
    3. Insert size distribution estimation from concordant pairs
    4. Graph construction and cycle detection
    5. Continuity scoring + depth ratio analysis
    6. Integrated confidence calculation

    Args:
        ref: Path to reference FASTA.
        fastq_r1: Read 1 FASTQ (FASTQ mode).
        fastq_r2: Read 2 FASTQ (FASTQ mode).
        bam_path: Coordinate-sorted, indexed full-alignment BAM.
        split_bam: Split-reads BAM (pre-aligned mode).
        discord_bam: Discordant-reads BAM (pre-aligned mode).

    Returns:
        pd.DataFrame: Detected ecDNA circles with confidence scores.
    """
    print(f"=== Starting EccDNA Detection Pipeline ===")
    print(f"Reference: {ref}")

    # 1. Extraction
    if fastq_r1 and fastq_r2:
        print(">> Mode: FASTQ alignment with bwa-mem2 + samblaster")
        if two_pass:
            print("   Two-pass: -5SP for split reads, standard for discordant pairs")
        extractor = ReadExtractor(ref, fastq_r1, fastq_r2, align_threads=threads, two_pass=two_pass)
        split_df, discord_df = extractor.align_and_extract()
        if bam_path is None:
            bam_path = extractor.bam_path
    elif split_bam and discord_bam and bam_path:
        print(">> Mode: pre-aligned BAM files")
        extractor = ReadExtractor.from_bams(ref, split_bam, discord_bam, bam_path)
        split_df, discord_df = extractor.parse_existing_bams()
    else:
        raise ValueError(
            "Provide either (fastq_r1 + fastq_r2) for FASTQ mode, "
            "or (split_bam + discord_bam + bam_path) for pre-aligned BAM mode."
        )

    if split_df.empty:
        print("No split reads found. Exiting.")
        return pd.DataFrame()

    # 2. Probabilistic Realignment
    print(">> Running Matrix-Based Probabilistic Realignment...")
    realigner = ProbabilisticRealigner(ref, n_threads=threads)
    validated_df = realigner.realign_junctions(split_df)

    valid_count = validated_df['is_valid'].sum()
    print(f"Validated {valid_count} high-confidence junctions.")

    if valid_count == 0:
        print("No high-confidence junctions found.")
        return pd.DataFrame()

    # 3. Estimate insert size distribution
    print(">> Estimating insert size distribution...")
    insert_metrics = estimate_insert_size(bam_path)
    print(f"   Mean insert: {insert_metrics[0]:.1f} bp, SD: {insert_metrics[1]:.1f} bp")

    # 4. Graph Construction & Cycle Detection
    print(">> Building Breakpoint Graph and Detecting Cycles...")
    solver = EccDNAGraphSolver()
    solver.build_graph(validated_df, discord_df)

    circles = solver.detect_cycles(
        bam_path=bam_path,
        discordant_df=discord_df,
        insert_metrics=insert_metrics,
    )

    # 5. Output
    output_file = f"{OUTPUT_DIR}/eccdna_calls.csv"
    circles.to_csv(output_file, index=False)
    print(f"=== Detection Complete. Results saved to {output_file} ===")

    if not circles.empty:
        print("\nDetected ecDNA circles:")
        print(circles.to_string())

        print(f"\nSummary:")
        print(f"  Total circles detected: {len(circles)}")
        print(f"  Single-segment: {(circles['type'] == 'Single-Segment').sum()}")
        print(f"  Multi-segment: {(circles['type'] == 'Multi-Segment').sum()}")

        if 'integrated_confidence' in circles.columns and circles['integrated_confidence'].notna().any():
            high_conf = (circles['integrated_confidence'] >= 0.9).sum()
            print(f"  High confidence (>=0.9): {high_conf}")

        if 'depth_ratio' in circles.columns and circles['depth_ratio'].notna().any():
            mean_ratio = circles['depth_ratio'].mean()
            print(f"  Mean depth ratio: {mean_ratio:.2f}x")

        if 'continuity_score' in circles.columns and circles['continuity_score'].notna().any():
            mean_cont = circles['continuity_score'].mean()
            print(f"  Mean continuity score: {mean_cont:.3f}")
    else:
        print("No circular DNA structures detected.")

    return circles


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="ecDNA Detector — Fast, matrix-based junction detection"
    )

    # Reference is always required
    parser.add_argument("-r", "--reference", required=True, help="Reference FASTA (indexed)")

    # FASTQ mode (mutually exclusive with BAM mode)
    parser.add_argument("-1", "--fastq1", help="Read 1 FASTQ (paired-end)")
    parser.add_argument("-2", "--fastq2", help="Read 2 FASTQ (paired-end)")

    # Pre-aligned BAM mode
    parser.add_argument("-sb", "--split-bam", help="Split-reads BAM from samblaster")
    parser.add_argument("-db", "--discordant-bam", help="Discordant-reads BAM from samblaster")
    parser.add_argument("-b", "--bam", help="Coordinate-sorted, INDEXED full-alignment BAM")
    parser.add_argument("-t", "--threads", type=int, default=16, help="Worker threads for realignment")
    parser.add_argument("--two-pass", action="store_true",
                        help="Two bwa-mem2 passes: -5SP for split reads, standard for discordant pairs (2x alignment time)")

    args = parser.parse_args()
    run_pipeline(
        ref=args.reference,
        fastq_r1=args.fastq1,
        fastq_r2=args.fastq2,
        bam_path=args.bam,
        split_bam=args.split_bam,
        discord_bam=args.discordant_bam,
        threads=args.threads,
        two_pass=args.two_pass,
    )