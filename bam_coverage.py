import pysam as ps
import os
import pybedtools as bt
import numpy as np
import time

begin = time.time()

os.chdir("/isdata/kroghgrp/xsh723/projects/circle_map/test_data/circle_map_coverage")

bam = ps.AlignmentFile('sorted_B06.bam','rb')
bed = bt.BedTool('coverage.bed')
print(len(bed))
reference_contigs = bam.header['SQ']


dict = {}


for reference in reference_contigs:
    dict[reference['SN']] = reference['LN']

print(dict)


output = []
for interval in bed:
    if interval.start < 0:
        start = 0
    else:
        start = interval.start

    if interval.end > dict[interval.chrom]:
        end = interval.end
    else:
        end = interval.end


    cov = bam.count_coverage(contig=interval.chrom,start=start,end=end,quality_threshold=0)
    summarized_cov = np.array([cov[0],cov[1],cov[2],cov[3]]).sum(axis=0)
    extended_cov = summarized_cov
    interval_cov = summarized_cov[100:-100]
    mean = np.mean(summarized_cov)
    sd = np.std(summarized_cov)
    coverage_ratio = np.sum(interval_cov)/np.sum(extended_cov)
    zero_frac = np.count_nonzero( summarized_cov == 0)/len(summarized_cov)

    interval.append(str(mean))
    interval.append(str(sd))
    interval.append(str(coverage_ratio))
    interval.append(str(zero_frac))
    print(interval)
    output.append(interval)

print(bt.BedTool(output))

end = time.time()

print((end-begin)/60)