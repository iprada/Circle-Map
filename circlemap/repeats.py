#MIT License
#
#Copyright (c) 2019 Iñigo Prada Luengo
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

import pysam as ps
import pybedtools as bt
import os
import time
from circlemap.utils import merge_coverage_bed,rightmost_from_sa


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
