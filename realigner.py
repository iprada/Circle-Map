#!/data/xsh723/anaconda/bin/python3.6
from __future__ import division


import subprocess as sb
import warnings
import os
import pysam as ps
import pybedtools as bt
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import time
import sys
import numpy as np
import pandas as pd
from utils import *

class realignment:
    """Class for managing the realignment and eccDNA indetification of circle-map"""

    def __init__(self, input_bam,qname_bam,genome_fasta,directory,mapq_cutoff,insert_size_mapq,std_extension,
                 insert_size_sample_size,n_hits,prob_cutoff,ncores):
        self.input_bam = input_bam
        self.qname_bam = qname_bam
        self.genome_fa = ps.FastaFile(genome_fasta)
        self.directory = directory
        self.mapq_cutoff = mapq_cutoff
        self.insert_size_mapq = insert_size_mapq
        self.std_extenstion = std_extension
        self.insert_sample_size = insert_size_sample_size
        self.cores = ncores
        self.n_hits = n_hits
        self.prob_cutoff = prob_cutoff
        self.phreds_to_probs = np.vectorize(phred_to_prob)





    def realignment(self):
        """Function that will iterate trough the bam file containing reads indicating eccDNA structural variants and
        will output a bed file containing the soft-clipped reads, the discordant and the coverage within the interval"""

        begin = time.time()

        os.chdir(self.directory)

        eccdna_bam = ps.AlignmentFile("%s" % self.input_bam, "rb")





        circ_peaks,sorted_bam = bam_circ_sv_peaks(eccdna_bam,self.input_bam,self.cores)

        print("Computing insert size and standard deviation from %s F1R2 reads with a mapping quality of %s" %
              (self.insert_sample_size,self.insert_size_mapq))

        insert_metrics = insert_size_dist(self.insert_sample_size,self.insert_size_mapq,self.qname_bam)



        print("The computed insert size is %f with a standard deviation of %s" % (insert_metrics[0],insert_metrics[1]))

        #define realignment extension interval

        extension = insert_metrics[0] + self.std_extenstion*insert_metrics[1]

        #Find a mate interval for every interval
        iteration = 0




        for interval in circ_peaks:




            candidate_mates = get_mate_intervals(sorted_bam,interval,self.mapq_cutoff)

            if len(candidate_mates) > 0:
                iteration +=1
                print(iteration)


                # sort merge and extend
                realignment_interval_extended = get_realignment_intervals(candidate_mates,extension)


                #add counters
                realignment_intervals = realignment_intervals_with_counter(realignment_interval_extended)


                for mate_interval in realignment_intervals:

                    print(mate_interval)


                    #sample realignment intervals
                    plus_coding_interval = self.genome_fa.fetch(mate_interval.chrom,mate_interval.start,mate_interval.end).upper()
                    minus_coding_interval = str(Seq(plus_coding_interval).complement())


                    #note that I am getting the reads of the interval. Not the reads of the mates
                    for read in sorted_bam.fetch(interval.chrom,interval.start,interval.end):


                        if is_soft_clipped(read):

                            if read.mapq >= self.mapq_cutoff:

                                # no need to realignment
                                if read.has_tag('SA'):


                                    #check realignment from SA tag

                                    if circle_from_SA(read,self.mapq_cutoff,mate_interval) == True:

                                        mate_interval[3] = int(mate_interval[3]) + 1

                                else:
                                    #realignment

                                    realignment_bases = get_longest_soft_clipped_bases(read)

                                    #realign(read,self.n_hits,plus_coding_interval,minus_coding_interval)





                            else:
                                a = 0
                                #Nothing! no more loops
                        else:

                            # check discordancy if the read is discordant
                            a = 0


        end = time.time()

        eccdna_bam.close()

        print((end-begin)/60)
