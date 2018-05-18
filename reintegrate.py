import pysam as ps
import pybedtools as bt
import os

os.chdir('/isdata/kroghgrp/xsh723/projects/scbrain')
bed = bt.BedTool('deletions_eccdna_neuron_095.bed')
bam = ps.AlignmentFile('sorted_neuron_sra.bam','rb')

for interval in bed:

    print(interval[0],interval[1],interval[2])
