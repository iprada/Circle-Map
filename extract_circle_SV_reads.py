#!/data/xsh723/anaconda/bin/python3.6
from __future__ import division
#Author Inigo Prada Luengo
#email: inigo.luengo@bio.ku.dk
import pysam as ps
import os
import time
from utils import is_soft_clipped

class readExtractor:
    def __init__(self,sorted_bam,output_bam,working_dir):
        self.sorted_bam = sorted_bam
        self.output_bam = output_bam
        self.working_dir = working_dir

    def extract_sv_circleReads(self):
        """Function that extracts Structural Variant reads that indicate circular DNA,
        The programme with extract soft-clipped reads and R2F1 (<- ->) oriented reads"""
        os.chdir(self.working_dir)

        #open input bam and create output bam
        raw_bam = ps.AlignmentFile(self.working_dir + self.sorted_bam, "rb")
        circle_sv_reads = ps.AlignmentFile(self.working_dir + self.output_bam, "wb", template=raw_bam)

        print("Extracting circular structural variants")
        #timing
        begin = time.time()
        #cache read1. this speed-ups the search
        read1 = ''
        for read in raw_bam:
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
                                    circle_sv_reads.write(read1)
                                    circle_sv_reads.write(read2)
                            else:
                                #if the leftmost mapping condition is not met check if they are soft-clipped
                                if is_soft_clipped(read1) == True:
                                    circle_sv_reads.write(read1)
                                    circle_sv_reads.write(read2)
                                if is_soft_clipped(read2) == True:
                                    circle_sv_reads.write(read2)
                                    circle_sv_reads.write(read1)

                        else:
                            #check soft-clipped if R2F1 orientation is not met
                            if is_soft_clipped(read1) == True:
                                circle_sv_reads.write(read1)
                                circle_sv_reads.write(read2)
                            if is_soft_clipped(read2) == True:
                                circle_sv_reads.write(read2)
                                circle_sv_reads.write(read1)
                    else:
                        #check if they are soft-clipped if both reads are not mapped
                        if is_soft_clipped(read1) == True:
                            circle_sv_reads.write(read1)
                            circle_sv_reads.write(read2)
                        if is_soft_clipped(read2) == True:
                            circle_sv_reads.write(read2)
                            circle_sv_reads.write(read1)



        end = time.time()

        print("finished extracting reads. Elapsed time:",(end-begin)/60,"mins")




        #cache read1
#read1=''
#for read in without_plasmids:
    #check that is read 1 and cache it
#    if read.is_read1:
#        read1=read


#    else:
        #now, assert that reads are pairs
#        if read.is_read2 and read1.query_name == read.query_name:
#            read2 = read
#            #assert that it is paired
#            if read1.is_paired and read1.is_read1:
#                #assert that mate is R and read1 is F
#                if read1.mate_is_reverse and read1.is_reverse == False:
#                    # check that R2 reference start is < than F1
#                    if read1.reference_start > read1.next_reference_start:
#                        #only save then if they are mapped to the same chromosome
#                        if read1.reference_id == read2.reference_id:
#                            circle_reads.write(read1)
#                            circle_reads.write(read2)
#            else:
#                # if they are not paired, check it r1 and r2 are soft clipped
#                if is_soft_clipped(read1) == True:
#                    circle_reads.write(read1)
#                if is_soft_clipped(read):
#                    circle_reads.write(read)



#circle_reads.close()
#print("finished extracting reads")











