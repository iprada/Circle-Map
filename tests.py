import pysam as ps
import os
import matplotlib.pyplot as plt

os.chdir("/isdata/kroghgrp/xsh723/rasmus/B02_05_samples/")

mapping_qualities = []
bamfile02 = ps.AlignmentFile("BWAB03.bam")

i = 0

for read in bamfile02:
    while i < 10000000:
        print(i)
        mapping_qualities.append(read.mapq)
        i +=1
    else:
        plt.hist(mapping_qualities)
        plt.show()

        exit()