import pysam as ps
import pybedtools as bt
import os
import time


os.chdir("/isdata/kroghgrp/xsh723/projects/13_resequencing_4_aged/whole_merge/circle_bams/BA6")

bed = bt.BedTool("BA6_conservative.bed")

bam = ps.AlignmentFile("sorted_BA6.bam", "rb")


begin = time.time()
for pileupcolumn in bam.fetch("chrXII",450728,1072005):
    a = 0

    #print(interval)
    #print(pileupcolumn)
    break

end = time.time()

print((end-begin)/60)

