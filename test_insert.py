import numpy as np

import pysam as ps
import os
from utils import *

os.chdir("/isdata/kroghgrp/xsh723/circle_map/test_data/real_data_t01/")

raw_bam = ps.AlignmentFile('qname_t08.bam', "rb")

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
                if read1.is_proper_pair:
                    if is_hard_clipped(read1) == False and is_hard_clipped(read2) == False:
                        if is_soft_clipped(read1) == False and is_soft_clipped(read2) == False:
                            if read1.is_reverse == False and read2.is_reverse == True:
                                if read1.tlen > 0:
                                    insert_size.append(read1.tlen)
                                    counter +=1

    if counter >= 1000000:
        break
    else:
        continue







#insert_data = np.loadtxt("inserts.txt")
#insert_data = np.absolute(insert_data)
print(insert_size)
mean = np.mean(insert_size)
std = np.std(insert_size)
print(mean,std)



