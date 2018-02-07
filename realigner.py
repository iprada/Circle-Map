#!/data/xsh723/anaconda/bin/python3.6
from __future__ import division


import os
from Bio.Seq import Seq
import time
from utils import *


class realignment:
    """Class for managing the realignment and eccDNA indetification of circle-map"""

    def __init__(self, input_bam,qname_bam,genome_fasta,directory,mapq_cutoff,insert_size_mapq,std_extension,
                 insert_size_sample_size,gap_open,gap_ext,n_hits,prob_cutoff,ncores,min_soft_clipped_length,overlap_frac,
                 interval_p_cut):
        #I/O
        self.input_bam = input_bam
        self.qname_bam = qname_bam
        self.directory = directory
        self.genome_fa = ps.FastaFile(genome_fasta)

        #realignment parameters

        # probabilistic realignment options
        self.n_hits = n_hits
        self.prob_cutoff = prob_cutoff
        self.min_sc_length = min_soft_clipped_length
        self.mapq_cutoff = mapq_cutoff
        self.interval_p = interval_p_cut

        # affine gap scoring options
        self.gap_open = gap_open
        self.gap_ext = gap_ext


        #insert size stimation parameters
        self.insert_size_mapq = insert_size_mapq
        self.std_extenstion = std_extension
        self.insert_sample_size = insert_size_sample_size


        #output options

        self.overlap_fraction = overlap_frac

        #regular options
        self.cores = ncores
        self.phreds_to_probs = np.vectorize(phred_to_prob)






    def realignment(self):
        """Function that will iterate trough the bam file containing reads indicating eccDNA structural variants and
        will output a bed file containing the soft-clipped reads, the discordant and the coverage within the interval"""

        partial_time = 0
        begin = time.time()

        os.chdir(self.directory)

        eccdna_bam = ps.AlignmentFile("%s" % self.input_bam, "rb")




        #compute genome coverage of the eccDNA SV reads and sort bam
        circ_peaks,sorted_bam = bam_circ_sv_peaks(eccdna_bam,self.input_bam,self.cores)


        # compute insert size distribution
        print("Computing insert size and standard deviation from %s F1R2 reads with a mapping quality of %s" %
              (self.insert_sample_size,self.insert_size_mapq))

        insert_metrics = insert_size_dist(self.insert_sample_size,self.insert_size_mapq,self.qname_bam)



        print("The computed insert size is %f with a standard deviation of %s" % (insert_metrics[0],insert_metrics[1]))


        #define realignment extension interval
        extension = insert_metrics[0] + self.std_extenstion*insert_metrics[1]


        iteration = 0


        results = []

        for interval in circ_peaks:

            print(interval)

            interval_sc = []
            interval_sa = []
            interval_dr = []


            #find out the prior distribution (mate alignment positions).
            candidate_mates = get_mate_intervals(sorted_bam,interval,self.mapq_cutoff)
            print("priors",candidate_mates)



            #check that the output is not empty
            if len(candidate_mates) > 0:

                iteration +=1
                print(iteration)


                # sort merge and extend
                realignment_interval_extended = get_realignment_intervals(candidate_mates,extension,self.interval_p)

                print("extended intervals")
                print(realignment_interval_extended)

                exit()

                if realignment_interval_extended == None:
                    continue




                for mate_interval in realignment_interval_extended:



                    #sample realignment intervals
                    #fasta file fetch is 1 based that why I do +1

                    plus_coding_interval = self.genome_fa.fetch(mate_interval.chrom,mate_interval.start+1,mate_interval.end+1).upper()
                    interval_length = len(plus_coding_interval)
                    minus_coding_interval = str(Seq(plus_coding_interval).complement())

                    # precompute the denominators of the error model. They will be constants for every interval
                    plus_base_freqs = background_freqs(plus_coding_interval)
                    minus_base_freqs = {'T':plus_base_freqs['A'],'A':plus_base_freqs['T'],
                                        'C':plus_base_freqs['G'],'G':plus_base_freqs['C']}


                    #note that I am getting the reads of the interval. Not the reads of the mates
                    for read in sorted_bam.fetch(interval.chrom,interval.start,interval.end):


                        if is_soft_clipped(read):

                            if read.mapq >= self.mapq_cutoff:

                                # no need to realignment
                                if read.has_tag('SA'):


                                    #check realignment from SA tag
                                    support = circle_from_SA(read, self.mapq_cutoff, mate_interval)



                                    if support['support'] == True:

                                        #compute mapping positions

                                        read_end = rightmost_from_read(read)

                                        supplementary_end = rightmost_from_sa(support['leftmost'],support['cigar'])



                                        # I store the read name to the output, so that a read counts as 1 no matter it is SC in 2 pieces
                                        if read.reference_start < support['leftmost']:

                                            interval_sa.append([interval.chrom,read.reference_start,(supplementary_end-1),read.qname,'SA'])

                                        elif read.reference_start > support['leftmost']:

                                            interval_sa.append(
                                                [interval.chrom, (support['leftmost']-1), read_end, read.qname, 'SA'])

                                        else:
                                            #uninformative read
                                            continue



                                else:
                                    #sc length
                                    sc_len = len(get_longest_soft_clipped_bases(read)['seq'])


                                    if sc_len >= self.min_sc_length:
                                    #realignment

                                        realignment_dict = realign(read,self.n_hits,plus_coding_interval,minus_coding_interval,
                                                                   plus_base_freqs,minus_base_freqs,self.gap_open,self.gap_ext)


                                        if realignment_dict == None:

                                            continue

                                        else:


                                            if realignment_probability(realignment_dict,interval_length) >= (1 - self.prob_cutoff):

                                                # here I have to retrieve the nucleotide mapping positions. Which should be the
                                                # the left sampling pysam coordinate - edlib coordinates

                                                read_end = rightmost_from_read(read)


                                                soft_clip_start = mate_interval.start + int(realignment_dict['alignments'][1][0][0])

                                                soft_clip_end = mate_interval.start + int(realignment_dict['alignments'][1][0][1])




                                                # I store the read name to the output, so that a read counts as 1 no matter it is SC in 2 pieces
                                                if read.reference_start < mate_interval.start + int(
                                                        realignment_dict['alignments'][1][0][0]):

                                                    interval_sc.append([interval.chrom, read.reference_start, soft_clip_end+1, read.qname, 'SC'])

                                                elif read.reference_start + mate_interval.start + int(
                                                        realignment_dict['alignments'][1][0][0]):

                                                    interval_sc.append([interval.chrom, soft_clip_start, read_end, read.qname, 'SC'])

                                                else:
                                                    # uninformative read
                                                    continue



                                            else:
                                                continue

                            else:
                                continue
                        else:
                            #discordant reads
                            #R2F1 oriented when iterating trough R2
                            if read.is_reverse == True and read.mate_is_reverse == False:
                                if read.is_read2:
                                    if read.reference_start < read.next_reference_start:
                                        # discordant read
                                        interval_dr.append([interval.chrom,read.reference_start,read.next_reference_start + read.infer_query_length(),0,read.qname])




                            #R2F1 when iterating trough F1
                            elif read.is_reverse == False and read.mate_is_reverse ==  True:
                                if read.is_read2 == False:
                                    if read.next_reference_start < read.reference_start:
                                        interval_dr.append([interval.chrom, read.next_reference_start,read.reference_start + read.infer_query_length(),0,read.qname])

            time_sc_sa = time.time()


            #OPTIMIZING
            sc_sa_dict = intersect_sa_sc(interval_sa,interval_sc,self.overlap_fraction)

            end_sc_sa = time.time()

            partial_time += (time_sc_sa - end_sc_sa) / 60
            print("timing_sc_sa",partial_time)

            if len(interval_dr) > 0:

                if sc_sa_dict != None:

                    if len(sc_sa_dict['data']) > 0:

                        interval_output = add_discordants(sc_sa_dict,interval_dr)

                        for interval in interval_output:
                            results.append(interval)
                            if len(interval) !=5:
                                exit()
                    else:
                        continue

            else:
                if sc_sa_dict != None:

                    if len(sc_sa_dict['data']) > 0:
                        for interval in sc_sa_dict['data']:
                            output_test = [interval.chrom,interval.start,interval.end,interval[3],0]
                            results.append(output_test)
                            if len(output_test) !=5:
                                exit()












            end = time.time()
            print((end - begin) / 60)

        eccdna_bam.close()

        #the end

        unparsed_bed = bt.BedTool(results)

        grouped_bed = unparsed_bed.sort().groupby(g=[1,2,3],c=[4,5],o=['sum','sum'])

        grouped_bed.saveas("circle_map_results.bed")

        print((end-begin)/60)


