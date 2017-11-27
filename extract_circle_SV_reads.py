#!/data/xsh723/anaconda/bin/python3.6
from __future__ import division
#Author Inigo Prada Luengo
#email: inigo.luengo@bio.ku.dk



import pysam as ps
import os
import time
import sys
import warnings
from utils import is_soft_clipped,is_hard_clipped

class readExtractor:
    """Class for managing the read extracting part of circle map"""
    def __init__(self,sorted_bam,output_bam,working_dir,extract_discordant,extract_soft_clipped,extract_hard_clipped,
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
        #verbose level
        self.verbose = int(verbose)
        #parser options
        self.parser = parser

    def extract_sv_circleReads(self):

        """Function that extracts Structural Variant reads that indicate circular DNA,
        The programme with extract soft-clipped reads and R2F1 (<- ->) oriented reads"""

        os.chdir(self.working_dir)

        #open input bam and create output bam
        raw_bam = ps.AlignmentFile(self.working_dir + self.sorted_bam, "rb")

        #check that the tag is present
        if 'HD' in raw_bam.header:

            #check if bam is qname sorted
            if raw_bam.header['HD']['SO'] != 'queryname':

                #give some errors if not
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




        circle_sv_reads = ps.AlignmentFile(self.working_dir + self.output_bam, "wb", template=raw_bam)
        #change the SO tag of the header to unsorted

        #modify the tag to unsorted
        if 'HD' in raw_bam.header == True:
            circle_sv_reads.header['HD']['SO'] = 'unsorted'

        if self.verbose >=3:
            print("Extracting circular structural variants")
        #timing
        begin = time.time()


        #cache read1. this speed-ups the search
        read1 = ''

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
                    # both reads in memory. Now operate
                    read2 = read
                    #check that both are mapped

                    #both reads need to be mapped
                    if read1.is_unmapped == False and read2.is_unmapped == False:

                        # Check read 2 aligned to the reverse strand, read1 aligned to the forward strand
                        if read2.is_reverse and read1.is_reverse == False:

                            # read2 leftmost mapping position smaller than read1 leftmost mapping position
                            if read2.reference_start < read1.reference_start:

                                #aligned to the same chromosome
                                if read1.reference_id == read2.reference_id:

                                    #check that extraction is not turn off
                                    if self.no_discordants == False:
                                        circle_sv_reads.write(read1)
                                        circle_sv_reads.write(read2)
                                    else:
                                        continue
                            else:
                                #if the leftmost mapping condition is not met check if they are soft-clipped
                                if is_soft_clipped(read1) == True:


                                    # check that extraction is not turn off
                                    if self.no_soft_clipped == False:

                                        circle_sv_reads.write(read1)
                                        circle_sv_reads.write(read2)
                                    else:
                                        continue

                                else:
                                    #check hard-clipped
                                    if is_hard_clipped(read1) == True:

                                        # check that extraction is not turn off
                                        if self.no_hard_clipped == False:

                                            circle_sv_reads.write(read1)
                                        else:
                                            continue


                                if is_soft_clipped(read2) == True:

                                    # check that extraction is not turn off
                                    if self.no_soft_clipped == False:

                                        circle_sv_reads.write(read2)
                                        circle_sv_reads.write(read1)
                                    else:
                                        continue

                                else:

                                    #check hard-clipped
                                    if is_hard_clipped(read2) == True:

                                        # check that extraction is not turn off
                                        if self.no_hard_clipped == False:

                                            circle_sv_reads.write(read2)

                                        else:

                                            continue


                        else:
                            #check soft-clipped if R2F1 orientation is not met
                            if is_soft_clipped(read1) == True:

                                # check that extraction is not turn off
                                if self.no_soft_clipped == False:

                                    circle_sv_reads.write(read1)
                                    circle_sv_reads.write(read2)

                                else:

                                    continue

                            else:
                                #check hard-clipped
                                if is_hard_clipped(read1) == True:


                                    # check that extraction is not turn off
                                    if self.no_hard_clipped == False:
                                        circle_sv_reads.write(read1)

                                    else:

                                        continue


                            if is_soft_clipped(read2) == True:

                                # check that extraction is not turn off
                                if self.no_soft_clipped == False:

                                    circle_sv_reads.write(read2)
                                    circle_sv_reads.write(read1)

                                else:

                                    continue

        end = time.time()

        circle_sv_reads.close()



        if self.verbose >=3:


            print("finished extracting reads. Elapsed time:", (end - begin) / 60, "mins")

            print("Thanks for using Circle-Map")








