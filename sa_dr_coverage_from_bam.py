#!/data/xsh723/anaconda/bin/python3.6

import os
import pysam as ps
import pybedtools as bt
import time
from utils import *

begin = time.time()

os.chdir("/isdata/kroghgrp/xsh723/circle_map/test_data/real_data_t01/")
print(os.listdir())

eccdna_bam = ps.AlignmentFile("circle_qname_sorted_t08.bam", "rb")

circle_sv_reads = ps.AlignmentFile("sa_dr_bam.bam", "wb", template=eccdna_bam)


soft_clipped = 0
sa = 0
hard_clipped = 0
dr = 0
read1 = ''
for read in eccdna_bam:
    if read.is_read1:
        read1 = read
    else:
        if read.is_read2 and read.qname == read1.qname:
            # both reads in memory
            read2 = read

            # both reads need to be mapped
            if read1.is_unmapped == False and read2.is_unmapped == False:

                if read2.is_reverse and read1.is_reverse == False:

                    # read2 leftmost mapping position smaller than read1 leftmost mapping position
                    if read2.reference_start < read1.reference_start:

                        # aligned to the same chromosome
                        if read1.reference_id == read2.reference_id:

                            if is_soft_clipped(read1) == True or is_soft_clipped(read2) == True:

                                continue

                            elif is_hard_clipped(read1) == True or is_hard_clipped(read2) == True:

                                continue

                            else:

                                circle_sv_reads.write(read1)
                                circle_sv_reads.write(read2)
                                dr +=2

                        else:
                            continue





end = time.time()


wall_time = (end-begin)/60
print("elapsed time ",wall_time)

print("Soft-clipped reads ",soft_clipped)
print("Supplementary alignments: ",sa)
print("Hard-clipped: ", hard_clipped)
print("Discordant reads ",dr)