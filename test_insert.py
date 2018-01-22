import numpy as np

import pysam as ps
import os

os.chdir("/isdata/kroghgrp/xsh723/circle_map/test_data/test_insert_size/")

raw_bam = ps.AlignmentFile('paired_end_sim_aln.bam', "rb")

counter = 0

insert_size = []
read1 = ''
for read in raw_bam:

    if read.is_read1:
        read1 = read
    else:
        if read.is_read2 and read.qname == read1.qname:
            read2 = read
            #both reads in memory
            if read1.mapq >= 60 and read2.mapq >= 60:
                if read2.is_reverse == True and read1.is_reverse == False:
                    if (read2.tlen + read1.infer_query_length())< 0   and (read1.tlen - read2.infer_query_length()) > 0:
                        insert_size.append(read1.tlen)
                        counter +=1

    if counter >= 100000:
        break
    else:
        continue







#insert_data = np.loadtxt("inserts.txt")
#insert_data = np.absolute(insert_data)
mean = np.mean(insert_size)
std = np.std(insert_size)
print(mean,std)



