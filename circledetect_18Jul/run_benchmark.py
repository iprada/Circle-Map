import sys, ast
sys.path.insert(0, '.')
from graph_solver import EccDNAGraphSolver
from extractor import estimate_insert_size
import pandas as pd, numpy as np
from collections import Counter
from scipy.sparse.csgraph import connected_components

OUT = 'eccdna_results/'
validated_df = pd.read_pickle(f'{OUT}/validated_df.pkl')
discord_df = pd.read_pickle(f'{OUT}/discord_df.pkl')

solver = EccDNAGraphSolver()
solver.build_graph(validated_df[validated_df['is_valid']], discord_df)

n_comp, labels = connected_components(solver.adj_matrix, directed=True, connection='strong')
sizes = Counter(labels)
single = sum(1 for c in range(n_comp) if sizes[c] == 1)
multi = sum(1 for c in range(n_comp) if sizes[c] > 1)
print(f'Graph: {solver.adj_matrix.shape}, nnz={solver.adj_matrix.nnz:,}')
print(f'SCCs: {n_comp} ({single} single, {multi} multi)')
for sz in sorted(set(sizes.values())):
    if sz > 1:
        cnt = sum(1 for c in range(n_comp) if sizes[c]==sz)
        print(f'  size={sz}: {cnt}')

circles = solver.detect_cycles(bam_path=None, top_n=0)
print(f'\nCircleDetect: {len(circles)} circles', flush=True)
if not circles.empty:
    singles = (circles['type']=='Single-Segment').sum()
    multis = (circles['type']=='Multi-Segment').sum()
    print(f'  Single: {singles}, Multi: {multis}', flush=True)
    circles.to_csv(f'{OUT}/eccdna_calls_v2.csv', index=False)

cm = pd.read_csv('/home/renhb/test/circle.bed', sep='\t', header=None,
                 names=['chr','start','end','discordants','sc','score','mean','std','start_ratio','end_ratio','continuity'])
print(f'\nCircle-Map: {len(cm)} circles', flush=True)

cm_near = 0
for _, r in cm.iterrows():
    nearby = validated_df[(validated_df['chr']==r['chr']) & (validated_df['pos']>=r['start']-500) & (validated_df['pos']<=r['end']+500) & validated_df['is_valid']]
    if len(nearby) > 0: cm_near += 1
print(f'CM circles with nearby junction: {cm_near}/{len(cm)} ({cm_near/len(cm)*100:.1f}%)', flush=True)

if 'Multi-Segment' in circles['type'].values:
    print(f'\nMulti-segment circles:', flush=True)
    for _, r in circles[circles['type']=='Multi-Segment'].iterrows():
        print(f'  segs={r["segment_count"]} conf={r["confidence"]:.4f}', flush=True)

print(f'\n--- Circle-Map column stats ---', flush=True)
for col in cm.columns:
    vals = cm[col].dropna()
    if vals.dtype in ['float64','int64']:
        print(f'  {col}: mean={vals.mean():.4f} median={vals.median():.4f}', flush=True)

print(f'\n--- Discrepancy Analysis ---', flush=True)
print(f'Circle-Map uses coverage ratio filter (-r 0.6):', flush=True)
print(f'  start_ratio>=0.6: {(cm["start_ratio"]>=0.6).sum()} ({(cm["start_ratio"]>=0.6).mean()*100:.1f}%)', flush=True)
print(f'Circle-Map split read count filter:', flush=True)
print(f'  sc>=2: {(cm["sc"]>=2).sum()}/{len(cm)}, mean sc={cm["sc"].mean():.1f}', flush=True)
print(f'  sc>=5: {(cm["sc"]>=5).sum()}/{len(cm)}', flush=True)
print(f'\nCircleDetect uses matrix posterior prob (GC-aware null model):', flush=True)
print(f'  Mean posterior: {validated_df[validated_df["is_valid"]]["posterior_prob"].mean():.4f}', flush=True)
print(f'  Junctions with >0.99 prob: {(validated_df[validated_df["is_valid"]]["posterior_prob"]>0.99).sum()}', flush=True)

print('\nDone.', flush=True)