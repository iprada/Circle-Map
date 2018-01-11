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
                 insert_size_sample_size,ncores):
        self.input_bam = input_bam
        self.qname_bam = qname_bam
        self.genome_fa = ps.FastaFile(genome_fasta)
        self.directory = directory
        self.mapq_cutoff = mapq_cutoff
        self.insert_size_mapq = insert_size_mapq
        self.std_extenstion = std_extension
        self.insert_sample_size = insert_size_sample_size
        self.cores = ncores
        self.phred_to_prob = np.vectorize(phred_to_prob)





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

            #clean temp files
            #bt.helpers.cleanup()

            #
            candidate_mates = get_mate_intervals(sorted_bam,interval,self.mapq_cutoff)

            if len(candidate_mates) > 0:
                iteration +=1
                print(iteration)

                # sort by column 4
                sorted_candidate_mates = sorted(candidate_mates, key=lambda x: x[3])
                candidate_mates = bt.BedTool(sorted_candidate_mates)
                # group by column 4. and get the minimum start point and the maximum end point
                grouped = candidate_mates.groupby(g=4, c=[1, 2, 3, 5], o=['distinct', 'min', 'max','distinct'])

                # reformat to fit pybedtools requirements and get pandas object. Usefull for calculations
                grouped_pandas = pd.read_table(grouped.fn, names=['read_type', 'chrom', 'start', 'stop','orientation'])
                grouped = grouped_pandas[['chrom', 'start', 'stop', 'read_type','orientation']]
                grouped = bt.BedTool.from_dataframe(grouped)

                realignment_intervals = get_realignment_interval(grouped,grouped_pandas,extension,sorted_bam)

                # add counter for discordanst and soft-clipped to realign
                realignment_intervals = realignment_intervals_with_counter(realignment_intervals)

                for mate_interval in realignment_intervals:

                    #sample realignment intervals
                    plus_coding_interval = self.genome_fa.fetch(mate_interval.chrom,mate_interval.start,mate_interval.end).upper()
                    minus_coding_interval = Seq(plus_coding_interval).reverse_complement()


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





                                    #bases = np.char.array(list(read.query_sequence))

                                    base_qualities = np.vstack((read.query_qualities,read.query_qualities,read.query_qualities
                                                                ,read.query_qualities))



                                    #print(phred_to_prob(base_qualities))

                                    #print(test)



                                    #exit()



                            else:
                                a = 0
                                #Nothing! no more loops
                        else:

                            # check discordancy if the read is discordant
                            a = 0


        end = time.time()

        eccdna_bam.close()

        print((end-begin)/60)
