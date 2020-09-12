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

from __future__ import division


import os
import sys
from Bio.Seq import Seq
import time
from circlemap.utils import *
import pandas as pd
import traceback
import multiprocessing as mp
import warnings
import datetime



class bam2bam:
    """Class for managing the realignment and eccDNA indetification of circle-map"""

    queue = mp.Manager().Queue()

    def __init__(self, input_bam,output,qname_bam,genome_fasta,directory,mapq_cutoff,insert_size_mapq,std_extension,
                 insert_size_sample_size,gap_open,gap_ext,n_hits,prob_cutoff,min_soft_clipped_length,
                 interval_p_cut,ncores,locker,verbose,pid,edit_distance_frac,
                 remap_splits,only_discordants,score,insert_size,manager):
        #I/O
        self.edit_distance_frac = edit_distance_frac
        self.ecc_dna_str = input_bam
        self.output = output
        self.qname_bam = qname_bam
        self.directory = directory
        self.genome_fa = genome_fasta

        #realignment parameters

        # probabilistic realignment options
        self.n_hits = n_hits
        self.prob_cutoff = prob_cutoff
        self.min_sc_length = min_soft_clipped_length
        self.mapq_cutoff = mapq_cutoff
        self.interval_p = interval_p_cut
        self.remap = remap_splits
        self.only_discordants = only_discordants
        self.score = score
        self.insert = insert_size

        # affine gap scoring options
        self.gap_open = gap_open
        self.gap_ext = gap_ext


        #insert size stimation parameters
        self.insert_size_mapq = insert_size_mapq
        self.std_extenstion = std_extension
        self.insert_sample_size = insert_size_sample_size

        #regular options
        self.cores = ncores
        self.verbose = verbose
        self.lock = locker



        #for instances running on the same directoiry

        self.pid = pid

        #parallel enviroment
        self.read_list = manager.list()
        self.read_count = manager.Value('i', 0)
        self.write_round = manager.Value('i', 0)







    def listener_writer(self,bam):

        f = open('test.sam',"w")

        header = bam.header

        while True:

            # Read from the queue and do nothing
            read = self.queue.get()


            if read == "DONE":
               f.close()
               print("breaking")
               bam.close()
               break
            else:

                pysam_read = ps.AlignedSegment.fromstring(read,bam.header)
                f.write(read + "\n")
                bam.write(pysam_read)

    def kill(self):
        print("KILLING")
        self.queue.put("DONE")

    def beta_version_warning(self):
        """Warn the user that this is experimental"""
        print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S: You are using a beta version feature"))
        warnings.warn("The bam2bam feature on Circle-Map is experimental. The development of this feature is active, but"
                      " have in mind that it might produce unintended results. Check https://github.com/iprada/Circle-Map"
                      " for the development status.")




    def realign(self,peaks):
        """Function that will iterate trough the bam file containing reads indicating eccDNA structural variants and
        will output a bed file containing the soft-clipped reads, the discordant and the coverage within the interval"""

        #open files for every process
        try:
            peaks_pd = pd.DataFrame.from_records(peaks,columns=['chrom', 'start', 'end'])
            genome_fa = ps.FastaFile(self.genome_fa)
            ecc_dna = ps.AlignmentFile(self.ecc_dna_str,"rb")

            begin = time.time()








            # compute insert size distribution

            insert_metrics = self.insert


            #define realignment extension interval
            extension = insert_metrics[0] + self.std_extenstion*insert_metrics[1]


            iteration = 0



            for index,interval in peaks_pd.iterrows():




                try:


                    #find out the prior distribution (mate alignment positions).
                    candidate_mates = get_mate_intervals(ecc_dna,interval,self.mapq_cutoff,self.verbose,self.only_discordants)






                    if len(candidate_mates) > 0 or candidate_mates != None:


                        realignment_interval_extended = get_realignment_intervals(candidate_mates,extension,self.interval_p,
                                                                                  self.verbose)
                        

                        if realignment_interval_extended is None:
                            continue



                        iteration_results = []
                        for index,mate_interval in realignment_interval_extended.iterrows():

                            iteration += 1



                            #sample realignment intervals
                            #fasta file fetch is 1 based that why I do +1

                            plus_coding_interval = genome_fa.fetch(str(mate_interval['chrom']),int(int(mate_interval['start'])+1),int(int(mate_interval['end'])+1)).upper()
                            interval_length = len(plus_coding_interval)
                            minus_coding_interval = str(Seq(plus_coding_interval).complement())

                            # precompute the denominators of the error model. They will be constants for every interval
                            plus_base_freqs = background_freqs(plus_coding_interval)

                            minus_base_freqs = {'T':plus_base_freqs['A'],'A':plus_base_freqs['T'],
                                                'C':plus_base_freqs['G'],'G':plus_base_freqs['C']}

                            minus_base_freqs = np.array([plus_base_freqs['T'],plus_base_freqs['A'],plus_base_freqs['G'],plus_base_freqs['C']])
                            plus_base_freqs = np.array([plus_base_freqs['A'],plus_base_freqs['T'],plus_base_freqs['C'],plus_base_freqs['G']])


                            #note that I am getting the reads of the interval. Not the reads of the mates

                            for read in ecc_dna.fetch(interval['chrom'],int(interval['start']),int(interval['end']),multiple_iterators=True):


                                if is_soft_clipped(read):

                                    if read.mapq >= self.mapq_cutoff:

                                        # no need to realignment
                                        if read.has_tag('SA') and self.remap != True:

                                            # check realignment from SA tag
                                            support = circle_from_SA(read, self.mapq_cutoff, mate_interval)

                                            if support is None:
                                                pass

                                            else:

                                                if support['support'] == True:
                                                    self.queue.put(read.to_string())

                                                else:
                                                    # uninformative read
                                                    pass



                                        else:
                                            #sc length
                                            sc_len = len(get_longest_soft_clipped_bases(read)['seq'])


                                            if non_colinearity(int(read.cigar[0][0]),int(read.cigar[-1][0]),int(read.pos),
                                                               int(mate_interval.start),int(mate_interval.end)) == True:
                                                if sc_len >= self.min_sc_length:
                                                    edits_allowed = adaptative_myers_k(sc_len, self.edit_distance_frac)
                                                #realignment

                                                    realignment_dict = realign(read,self.n_hits,plus_coding_interval,minus_coding_interval,
                                                                               plus_base_freqs,minus_base_freqs,self.gap_open,self.gap_ext,self.verbose,edits_allowed)


                                                    if realignment_dict == None:

                                                        pass

                                                    else:
                                                        #calc edit distance allowed
                                                        prob = realignment_probability(realignment_dict,interval_length)
                                                        if prob >= self.prob_cutoff and realignment_dict['alignments'][1][3] <= edits_allowed:

                                                            # here I have to retrieve the nucleotide mapping positions. Which should be the
                                                            # the left sampling pysam coordinate - edlib coordinates

                                                            read_end = rightmost_from_read(read)

                                                            #aln start on the reference
                                                            soft_clip_start = int(mate_interval['start'])+ int(realignment_dict['alignments'][1][0][0])

                                                            soft_clip_end = int(mate_interval['start']) + int(realignment_dict['alignments'][1][0][1])

                                                            score = sc_len*prob


                                                            # I store the read name to the output, so that a read counts as 1 no matter it is SC in 2 pieces
                                                            # Soft-clipped aligned upstream. Primary aligned downstream
                                                            if read.reference_start < int(mate_interval['start']) + int(
                                                                    realignment_dict['alignments'][1][0][0]):
                                                                    # construct tag
                                                                    sa_tag = realignment_read_to_SA_string(realignment_dict,
                                                                                                  prob, interval['chrom'],
                                                                                                  soft_clip_start)


                                                                    #read.tags += [('SA', sa_tag)]

                                                                    self.queue.put(read.to_string())



                                                            # soft-clipped aligned downstream primary alignment is upstream
                                                            elif read.reference_start + int(mate_interval['start']) + int(
                                                                    realignment_dict['alignments'][1][0][0]):

                                                                sa_tag = realignment_read_to_SA_string(realignment_dict,
                                                                                                       prob, interval[
                                                                                                           'chrom'],
                                                                                                       soft_clip_start)

                                                                read.tags += [('SA', sa_tag)]

                                                                self.queue.put(read.to_string())

                                                            else:
                                                                # uninformative read
                                                                pass



                                                        else:
                                                            pass
                                            else:

                                                pass

                                    else:
                                        pass
                                else:
                                    pass


                except BaseException as e:
                    traceback.print_exc(file=sys.stdout)
                    warnings.warn(
                        "Failed on interval %s due to the error %s" % (
                            str(interval), str(e)))
                    return([1,1])




            ecc_dna.close()
            genome_fa.close()

        except:
            print("Failed on cluster:")
            print(traceback.print_exc(file=sys.stdout))
            return([1,1])

        genome_fa.close()
        ecc_dna.close()


        return([0,0])