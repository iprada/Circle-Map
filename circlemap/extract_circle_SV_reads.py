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
import os
from src.utils import *
import time
import sys
import warnings


class readExtractor:
    """Class for managing the read extracting part of circle map"""
    def __init__(self,sorted_bam,output_bam,working_dir,mapq_cutoff,extract_discordant,extract_soft_clipped,extract_hard_clipped,
                 verbose,parser
                 ):
        #input-output
        self.sorted_bam = sorted_bam
        self.output_bam = output_bam
        #working place
        self.working_dir = working_dir

        #read options
        self.no_discordants = extract_discordant
        self.no_soft_clipped = extract_soft_clipped
        self.no_hard_clipped = extract_hard_clipped

        #mapq cutoff

        self.mapq_cutoff = mapq_cutoff

        #verbose level
        self.verbose = int(verbose)
        #parser options
        self.parser = parser

    def extract_sv_circleReads(self):

        """Function that extracts Structural Variant reads that indicate circular DNA,
        The programme with extract soft-clipped reads and R2F1 (<- ->) oriented reads"""

        os.chdir(self.working_dir)

        #input
        raw_bam = ps.AlignmentFile(self.working_dir + "/" + self.sorted_bam, "rb")

        #HD the tag for the header line. SO indicates sorting order of the alignements
        if 'HD' in raw_bam.header:

            if raw_bam.header['HD']['SO'] != 'queryname':
                sys.stderr.write(
                    "The input bam header says that bam is not sorted by queryname. It is sorted by %s\n\n" % (raw_bam.header['HD']['SO']))
                sys.stderr.write(
                    "Sort your bam file queryname with the following command:\n\n\tsamtools sort -n -o output.bam input.bam")

                time.sleep(0.01)

                self.parser.print_help()
                sys.exit(1)
        else:

            if self.verbose >=2:
                warnings.warn("WARNING:Circle-Map does not know if the input bam is queryname sorted\n Please check that, the output would be unexpected otherwise")
                print("As sanity check, sort your bam file queryname with the following command:\n\n\tsamtools sort -n -o output.bam input.bam")




        circle_sv_reads = ps.AlignmentFile(self.working_dir + "/" + self.output_bam , "wb", template=raw_bam)


        #modify the tag to unsorted
        if 'HD' in raw_bam.header == True:
            circle_sv_reads.header['HD']['SO'] = 'unsorted'

        if self.verbose >=3:
            print("Extracting circular structural variants")

        #timing
        begin = time.time()


        #cache read1. operate in read2. this speed-ups the search
        read1 = ''

        #counter for processed reads
        processed_reads = 0


        for read in raw_bam:

            if self.verbose >=3:
                processed_reads +=1


                if (processed_reads/1000000).is_integer() == True:
                    partial_timer = time.time()
                    partial_time = (partial_timer - begin)/60
                    print("Processed %s reads in %s mins" % (processed_reads,round(partial_time,3)))

            if read.is_read1:
                read1 = read
            else:
                if read.is_read2 and read.qname == read1.qname:
                    # both reads in memory
                    read2 = read


                    #both reads need to be mapped
                    if read1.is_unmapped == False and read2.is_unmapped == False:


                        if read2.is_reverse and read1.is_reverse == False:

                            # read2 leftmost mapping position smaller than read1 leftmost mapping position
                            if read2.reference_start < read1.reference_start:

                                #aligned to the same chromosome
                                if read1.reference_id == read2.reference_id:

                                    if read1.mapq >= self.mapq_cutoff and read2.mapq >= self.mapq_cutoff:

                                        #is discordant extraction turn off?

                                        if self.no_discordants == False:
                                            #add mate mapping quality info
                                            read1.tags += [('MQ',read2.mapq)]
                                            read2.tags += [('MQ', read1.mapq)]

                                            circle_sv_reads.write(read1)
                                            circle_sv_reads.write(read2)
                                        else:
                                            pass
                                    else:
                                        #extract soft-clipped if the mapq is high enough
                                        if read1.mapq >= self.mapq_cutoff:
                                            if is_soft_clipped(read1) == True:
                                                if self.no_soft_clipped == False:

                                                    if read1.mapq >= self.mapq_cutoff:
                                                        read1.tags += [('MQ', read2.mapq)]
                                                        circle_sv_reads.write(read1)

                                                else:

                                                    pass

                                            else:

                                                if is_hard_clipped(read1) == True:

                                                    if self.no_hard_clipped == False:

                                                        # gets its on mapq since mate is unmapped
                                                        if read1.mapq >= self.mapq_cutoff:
                                                            read1.tags += [('MQ', read1.mapq)]
                                                            circle_sv_reads.write(read1)

                                                    else:

                                                        pass
                                        if read2.mapq >= self.mapq_cutoff:
                                            if is_soft_clipped(read2) == True:

                                                if self.no_soft_clipped == False:

                                                    # gets its on mapq since mate is unmapped
                                                    if read2.mapq >= self.mapq_cutoff:
                                                        read2.tags += [('MQ', read2.mapq)]
                                                        circle_sv_reads.write(read2)

                                                else:
                                                    pass
                                            if is_hard_clipped(read2) == True:

                                                if self.no_hard_clipped == False:

                                                    # gets its on mapq since mate is unmapped
                                                    if read2.mapq >= self.mapq_cutoff:
                                                        read2.tags += [('MQ', read2.mapq)]
                                                        circle_sv_reads.write(read1)

                                                else:

                                                    pass


                                else:
                                    if is_soft_clipped(read1) == True:

                                        if self.no_soft_clipped == False:

                                            if read1.mapq >= self.mapq_cutoff:
                                                read1.tags += [('MQ', read2.mapq)]
                                                circle_sv_reads.write(read1)

                                        else:
                                            pass

                                    else:

                                        if is_hard_clipped(read1) == True:

                                            if self.no_hard_clipped == False:

                                                if read1.mapq >= self.mapq_cutoff:
                                                    read1.tags += [('MQ', read2.mapq)]
                                                    circle_sv_reads.write(read1)
                                            else:
                                                pass

                                    if is_soft_clipped(read2) == True:

                                        if self.no_soft_clipped == False:

                                            if read2.mapq >= self.mapq_cutoff:
                                                read2.tags += [('MQ', read1.mapq)]
                                                circle_sv_reads.write(read2)

                                        else:
                                            pass

                                    else:

                                        if is_hard_clipped(read2) == True:

                                            if self.no_hard_clipped == False:

                                                if read2.mapq >= self.mapq_cutoff:
                                                    read2.tags += [('MQ', read1.mapq)]
                                                    circle_sv_reads.write(read2)

                                            else:

                                                pass

                            else:
                                #if the leftmost mapping condition is not met check if they are soft-clipped
                                if is_soft_clipped(read1) == True:


                                    if self.no_soft_clipped == False:

                                        if read1.mapq >= self.mapq_cutoff:

                                            read1.tags += [('MQ', read2.mapq)]
                                            circle_sv_reads.write(read1)

                                    else:
                                        pass

                                else:

                                    if is_hard_clipped(read1) == True:


                                        if self.no_hard_clipped == False:

                                            if read1.mapq >= self.mapq_cutoff:

                                                read1.tags += [('MQ', read2.mapq)]
                                                circle_sv_reads.write(read1)
                                        else:
                                            pass


                                if is_soft_clipped(read2) == True:

                                    if self.no_soft_clipped == False:

                                        if read2.mapq >= self.mapq_cutoff:

                                            read2.tags += [('MQ', read1.mapq)]
                                            circle_sv_reads.write(read2)

                                    else:
                                        pass

                                else:


                                    if is_hard_clipped(read2) == True:


                                        if self.no_hard_clipped == False:

                                            if read2.mapq >= self.mapq_cutoff:

                                                read2.tags += [('MQ', read1.mapq)]
                                                circle_sv_reads.write(read2)

                                        else:

                                            pass


                        else:
                            #check soft-clipped if R2F1 orientation is not True
                            if is_soft_clipped(read1) == True:


                                if self.no_soft_clipped == False:

                                    if read1.mapq  >= self.mapq_cutoff:

                                        read1.tags += [('MQ', read2.mapq)]
                                        circle_sv_reads.write(read1)

                                else:

                                    pass

                            else:

                                if is_hard_clipped(read1) == True:



                                    if self.no_hard_clipped == False:

                                        if read1.mapq >= self.mapq_cutoff:
                                            read1.tags += [('MQ', read2.mapq)]
                                            circle_sv_reads.write(read1)

                                    else:

                                        pass


                            if is_soft_clipped(read2) == True:


                                if self.no_soft_clipped == False:

                                    if read2.mapq >= self.mapq_cutoff:

                                        read2.tags += [('MQ', read1.mapq)]
                                        circle_sv_reads.write(read2)


                                else:
                                    pass
                    else:
                        if read1.is_unmapped == False:
                            if is_soft_clipped(read1) == True:
                                if self.no_soft_clipped == False:

                                    if read1.mapq >= self.mapq_cutoff:
                                        read1.tags += [('MQ', read2.mapq)]
                                        circle_sv_reads.write(read1)

                                else:

                                    pass

                            else:

                                if is_hard_clipped(read1) == True:

                                    if self.no_hard_clipped == False:

                                        #gets its on mapq since mate is unmapped
                                        if read1.mapq >= self.mapq_cutoff:
                                            read1.tags += [('MQ', read1.mapq)]
                                            circle_sv_reads.write(read1)

                                    else:

                                        pass
                        if read2.is_unmapped == False:
                            if is_soft_clipped(read2) == True:

                                if self.no_soft_clipped == False:

                                    #gets its on mapq since mate is unmapped
                                    if read2.mapq >= self.mapq_cutoff:
                                        read2.tags += [('MQ', read2.mapq)]
                                        circle_sv_reads.write(read2)

                                else:
                                    pass
                            if is_hard_clipped(read2) == True:

                                if self.no_hard_clipped == False:

                                    # gets its on mapq since mate is unmapped
                                    if read2.mapq >= self.mapq_cutoff:
                                        read2.tags += [('MQ', read2.mapq)]
                                        circle_sv_reads.write(read1)

                                else:

                                    pass



                else:
                    # reads are not queryname sorted and cannot be processed in paired mode
                    warnings.warn("Unpaired reads found. Is your bam file queryname sorted?")


        end = time.time()

        circle_sv_reads.close()



        if self.verbose >=3:


            print("finished extracting reads. Elapsed time:", (end - begin) / 60, "mins")

            print("Thanks for using Circle-Map")