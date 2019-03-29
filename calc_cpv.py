from Coverage import coverage
import pybedtools as bt
import os
import pysam as ps

os.chdir("/home/iprada/faststorage/projects/circle_map_benchmark")
circ_explorer = bt.BedTool("2sv_circexplorer.bed")


#parse circexplorer
parsed = []
for interval in circ_explorer:
    parsed.append([interval.chrom,interval.start,interval.end,0,interval[4],0])


print(bt.BedTool(parsed).sort())
circlexplorer_cov = coverage("coord_SRR6315430.bam",bt.BedTool(parsed),200,0,100,"/home/iprada/faststorage/projects/circle_map_benchmark")

output = circlexplorer_cov.compute_coverage(circlexplorer_cov.get_wg_coverage())

unparsed_pd = output.to_dataframe(
    names=['chrom', 'start', 'end', 'discordants', 'soft-clipped', 'score', 'mean', 'std', 'start_ratio', 'end_ratio',
           'continuity'])
unparsed_pd.to_csv("CIRCexplorer_w_coverage.bed", header=None, index=None, sep='\t', mode='w')

print(output)