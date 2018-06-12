import pysam as ps
import pybedtools as bt
import os
from utils import rightmost_from_sa

os.chdir('/isdata/kroghgrp/xsh723/projects/scbrain')
bed = bt.BedTool('deletions_eccdna_neuron_095.bed')
bam = ps.AlignmentFile('sorted_neuron_sra.bam','rb')


supplementaries = []

for interval in bed:
    supplementaries = []
    print(interval)
    try:

        for read in bam.fetch(interval[0],int(interval[1]),int(interval[2])):
            if read.has_tag('SA'):
                sa = read.get_tag('SA')
                chrom=sa.split(',')[0]
                start=sa.split(',')[1]
                end=rightmost_from_sa(sa.split(',')[1],sa.split(',')[2])

                if int(start)<int(end):
                    a=0
                    #supplementaries.append([chrom,start,end,'sa'])
                else:
                    a=0
                    #supplementaries.append([chrom,end,start,'sa'])

            elif read.reference_id != read.next_reference_id:

                supplementaries.append([bam.get_reference_name(read.next_reference_id),
                                        read.next_reference_start-300,read.next_reference_start+300,'dr'])



        bed= bt.BedTool(supplementaries).sort().merge(c=[1],o='count')
        for sa_interval in bed:
            if int(sa_interval[3]) >=10 and interval.chrom != sa_interval.chrom:
                print(sa_interval)

    except BaseException as e:
        print(e)





