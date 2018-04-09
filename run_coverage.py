from test_coverage import coverage
import os
import pybedtools as bt
import time

begin = time.time()

os.chdir("/isdata/kroghgrp/xsh723/projects/13_resequencing_4_aged/whole_merge/circle_bams/B06/")


bam = "sorted_B06.bam"
bed= bt.BedTool('B06_conservative.bed')

object = coverage(bam,bed,200,0,300)


coverage_dict,header_dict = object.get_wg_coverage()

output = object.compute_coverage(coverage_dict,header_dict)
print(output)

output.saveas("data_with_coverage.bed")
end = time.time()

print((end-begin)/60)