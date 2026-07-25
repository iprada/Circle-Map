import sys
sys.path.insert(0, '/home/renhb/test/circledetect/circledetect_V2/Circle-Map-master/circledetect')
from graph_solver import EccDNAGraphSolver
from extractor import estimate_insert_size
import pandas as pd
import traceback

OUT = '/home/renhb/test/circledetect/circledetect_V2/Circle-Map-master/circledetect/eccdna_results'
FULL_BAM = '/home/renhb/test/out.sorted.bam'

validated_df = pd.read_pickle(f'{OUT}/validated_df.pkl')
discord_df = pd.read_pickle(f'{OUT}/discord_df.pkl')
print(f'Loaded {validated_df.is_valid.sum()} valid junctions, {len(discord_df)} discordant', flush=True)

valid_only = validated_df[validated_df['is_valid']].copy()
print(f'Building graph from {len(valid_only)} valid junctions...', flush=True)

insert_metrics = estimate_insert_size(FULL_BAM)
print(f'Insert size: mean={insert_metrics[0]:.1f} std={insert_metrics[1]:.1f}', flush=True)

solver = EccDNAGraphSolver()
solver.build_graph(valid_only, discord_df)

print('Detecting cycles (top 500 BAM-scored)...', flush=True)
try:
    circles = solver.detect_cycles(
        bam_path=FULL_BAM,
        discordant_df=discord_df,
        insert_metrics=insert_metrics,
        top_n=500,
    )
except Exception as e:
    print(f'ERROR in detect_cycles: {e}', flush=True)
    traceback.print_exc()
    sys.exit(1)

print(f'Total circles: {len(circles)}', flush=True)
if circles.empty:
    print('No circles detected', flush=True)
    sys.exit(0)

singles = (circles['type'] == 'Single-Segment').sum()
multis = (circles['type'] == 'Multi-Segment').sum()
print(f'Single: {singles}, Multi: {multis}', flush=True)

# Save immediately
circles.to_csv(f'{OUT}/eccdna_calls.csv', index=False)
print(f'Saved to {OUT}/eccdna_calls.csv', flush=True)

# Check BAM-scored columns
has_bam = 'integrated_confidence' in circles.columns
print(f'Has BAM-scored columns: {has_bam}', flush=True)
if has_bam:
    bam_scored = circles['integrated_confidence'].notna().sum()
    print(f'BAM-scored: {bam_scored}', flush=True)
    with_bam = circles[circles['integrated_confidence'].notna()]
    if not with_bam.empty:
        print(f'\nTop 20 BAM-scored circles:', flush=True)
        sorted_bam = with_bam.sort_values('integrated_confidence', ascending=False)
        for _, r in sorted_bam.head(20).iterrows():
            cont = r.get('continuity_score', 'N/A')
            depth = r.get('depth_ratio', 'N/A')
            print(f'  int={r["integrated_confidence"]:.4f} cont={cont} depth={depth} graph={r["confidence"]:.4f} segs={r["segment_count"]} {r["type"]}', flush=True)

print('\nDone.', flush=True)