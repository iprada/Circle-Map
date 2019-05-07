import pysam as ps
import pybedtools as bt
import os
import time
from utils import merge_coverage_bed,rightmost_from_sa


class repeat:
    """Class for indentifying repeat derived eccDNA by looking of the reads with two alignments"""


    def __init__(self,bam,directory,mismatch,fraction,read_number):
        self.bam = bam
        self.dir = directory
        self.mismatch = mismatch
        self.fraction = fraction
        self.number = read_number

    def find_circles(self):

        begin = time.time()

        os.chdir("%s" % self.dir)

        bam = ps.AlignmentFile("%s" % self.bam,'rb')



        print("Iterating trough the bam file")

        output = []
        for read in bam:

            try:
                if read.has_tag('XA'):
                    tag = read.get_tag('XA').split(';')[:-1]

                    read_edit_distance = read.get_tag('NM')

                    if read_edit_distance <= self.mismatch and len(tag) ==1:

                        read_chrom = bam.get_reference_name(read.reference_id)
                        chrom = tag[0].split(',')[0]


                        if chrom == read_chrom:


                            aln = int(tag[0].split(',')[1][1:])

                            if aln < read.reference_start:

                                interval = [chrom,aln,read.reference_start+ read.infer_read_length(),1]

                                output.append(interval)

                            else:

                                interval = [chrom,read.reference_start,rightmost_from_sa(aln,tag[0].split(',')[2]),1]
                                output.append(interval)


            except BaseException as e:
                print(e)



        bed = merge_coverage_bed(output,self.fraction,self.number)

        #add dots to read metrics stats

        with_dot = []
        for interval in bed:
            interval.append(".")
            with_dot.append(interval)

        bed= bt.BedTool(with_dot)

        return(bed)
