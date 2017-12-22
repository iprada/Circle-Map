#!/data/xsh723/anaconda/bin/python3.6
from __future__ import division


import subprocess as sb
import os
import pysam as ps
import pybedtools as bt
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import time
import sys
from utils import *

class realignment:
    """Class for managing the realignment and eccDNA indetification of circle-map"""

    def __init__(self, input_bam,genome_fasta,directory,mapq_cutoff,ncores):
        self.input_bam = input_bam
        self.genome_fa = genome_fasta
        self.directory = directory
        self.mapq_cutoff = mapq_cutoff
        self.cores = ncores



    def realignment(self):
        """Function that will iterate trough the bam file containing reads indicating eccDNA structural variants and
        will output a bed file containing the soft-clipped reads, the discordant and the coverage within the interval"""

        begin = time.time()

        os.chdir(self.directory)

        eccdna_bam = ps.AlignmentFile("%s" % self.input_bam, "rb")


        circ_peaks,sorted_bam = bam_circ_sv_peaks(eccdna_bam,self.input_bam,self.cores)





        #Find a mate interval for every interval
        for interval in circ_peaks:

            #this holds the intervals where the discordant and SA mates have mapped
            candidate_mates = []
            for read in sorted_bam.fetch(interval.chrom, interval.start, interval.end,multiple_iterators=True):

                if read.mapq >= self.mapq_cutoff:

                    # create mate interval based on the soft-clipped SA alignments
                    if is_soft_clipped(read) == True and read.has_tag('SA'):

                        read_chr = sorted_bam.get_reference_name(read.reference_id)
                        suplementary = read.get_tag('SA')


                        # [chr, left_most start, "strand,CIGAR,mapq, edit_distance]
                        supl_info = [x.strip() for x in suplementary.split(',')]




                        if read_chr == supl_info[0] and int(supl_info[4]) >= self.mapq_cutoff:

                            #split read with the same orientation
                            if (read.is_reverse == True and supl_info[2] == '-') or (read.is_reverse == False and supl_info[2] == '+'):


                                #SA is downstream, the interval is start, start+read length

                                #Future Inigo, this part of the code is over complicated. you can create a function of this
                                if read.reference_start > int(supl_info[1]):

                                    ref_alignment_length = genome_alignment_from_cigar(supl_info[3])


                                    #ref_alignment_length * 2 is done for extending the realignment region

                                    mate_interval = [interval.chrom,int(supl_info[1])-(ref_alignment_length*2),(int(supl_info[1]) + (ref_alignment_length*2)),"SA"]

                                    candidate_mates.append(mate_interval)


                                #SA is upstream, the interval is end - read length, end
                                elif read.reference_start < int(supl_info[1]):

                                    ref_alignment_length = genome_alignment_from_cigar(supl_info[3])


                                    # ref_alignment_length * 2 is done for extending the realignment region

                                    mate_interval = [interval.chrom,(int(supl_info[1])-(ref_alignment_length*2)),int(supl_info[1])+(ref_alignment_length*2),"SA"]

                                    candidate_mates.append(mate_interval)








                    # check discordant reads (R2F1 orientation)
                    elif read.is_unmapped  == False and read.mate_is_unmapped == False:

                        #check R2F1 orientation,when the R2 read
                        if read.is_reverse == True and read.mate_is_reverse == False:

                            #R2F1 order
                            if read.reference_start < read.next_reference_start:

                                if read.reference_id == read.next_reference_id:

                                    continue

                                    #create mate interval

                                    #print("R2F1 hit, for the R2")
                                    #print(read)

                                    #sys.exit(1)

                                    #mate_interval = [interval.chrom, read.next_reference_start - read_length,read.next_reference_start + read_length]
                                    #candidate_mates.append(bt.create_interval_from_list(mate_interval))


                        #R2F1 when iterating trough F1 read
                        elif read.is_reverse == False and read.mate_is_reverse == False:
                            if read.next_reference_start < read.reference_start:
                                if read.reference_id == read.next_reference_id:
                                    continue
                                    #a = 0
                                    #print("R2F1 hit, for the F1")
                                    #print(read)
                                    #sys.exit(1)
                                    # create mate interval


                    else:
                        #hard clips, soft clipped without SA
                        a = 0
                else:
                    #low mapping quality, do nothing
                    a=0


            print(candidate_mates)







        end = time.time()

        eccdna_bam.close()

        print((end-begin)/60)
