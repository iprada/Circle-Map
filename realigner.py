#!/data/xsh723/anaconda/bin/python3.6
from __future__ import division


import os
from Bio.Seq import Seq
import time
from utils import *

class realignment:
    """Class for managing the realignment and eccDNA indetification of circle-map"""

    def __init__(self, input_bam,qname_bam,genome_fasta,directory,mapq_cutoff,insert_size_mapq,std_extension,
                 insert_size_sample_size,gap_open,gap_ext,n_hits,prob_cutoff,ncores):
        self.input_bam = input_bam
        self.qname_bam = qname_bam
        self.genome_fa = ps.FastaFile(genome_fasta)
        self.directory = directory
        self.mapq_cutoff = mapq_cutoff
        self.insert_size_mapq = insert_size_mapq
        self.std_extenstion = std_extension
        self.insert_sample_size = insert_size_sample_size
        self.cores = ncores
        self.gap_open = gap_open
        self.gap_ext = gap_ext
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


        #list that will store the output
        results = []

        for interval in circ_peaks:

            interval_sc_sa = []
            interval_dr = []


            candidate_mates = get_mate_intervals(sorted_bam,interval,self.mapq_cutoff)

            if len(candidate_mates) > 0:
                iteration +=1
                print(iteration)


                # sort merge and extend
                realignment_interval_extended = get_realignment_intervals(candidate_mates,extension)

                if realignment_interval_extended == None:

                    continue


                for mate_interval in realignment_interval_extended:



                    #sample realignment intervals
                    plus_coding_interval = self.genome_fa.fetch(mate_interval.chrom,mate_interval.start,mate_interval.end).upper()
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

                                        read_end = read.reference_start  + aligned_bases(read)

                                        supplementary_end = support['leftmost'] + aligned_bases_from_sa(support['cigar'])

                                        # I store the read name to the output, so that a read counts as 1 no matter it is SC in 2 pieces
                                        if read.reference_start < support['leftmost'] and read_end < supplementary_end:

                                            interval_sc_sa.append([interval.chrom,read.reference_start,supplementary_end,read.qname,'SA'])

                                        elif read.reference_start > support['leftmost'] and read_end > supplementary_end:

                                            interval_sc_sa.append(
                                                [interval.chrom, (support['leftmost']-1), read_end, read.qname, 'SA'])

                                        else:
                                            #uninformative read
                                            continue



                                else:
                                    #realignment

                                    realignment_dict = realign(read,self.n_hits,plus_coding_interval,minus_coding_interval,
                                                               plus_base_freqs,minus_base_freqs,self.gap_open,self.gap_ext)


                                    if realignment_dict == None:

                                        continue

                                    else:

                                        if realignment_probability(realignment_dict,interval_length) >= (1 - self.prob_cutoff):

                                            # here I have to retrieve the nucleotide mapping positions. Which should be the
                                            # the left sampling pysam coordinate - edlib coordinates

                                            read_end = read.reference_start + aligned_bases(read)


                                            soft_clip_start = mate_interval.start + int(realignment_dict['alignments'][1][0][0])

                                            soft_clip_end = mate_interval.start + int(realignment_dict['alignments'][1][0][1])


                                            # I store the read name to the output, so that a read counts as 1 no matter it is SC in 2 pieces
                                            if read.reference_start < mate_interval.start + int(
                                                    realignment_dict['alignments'][1][0][0]) and read_end < soft_clip_end:

                                                interval_sc_sa.append([interval.chrom, read.reference_start, soft_clip_end, read.qname, 'SC'])

                                            elif read.reference_start > soft_clip_start and read_end > soft_clip_end:

                                                interval_sc_sa.append([interval.chrom, soft_clip_start, read_end, read.qname, 'SC'])

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

            print(bt.BedTool(interval_sc_sa))
            print(bt.BedTool(interval_dr))






        end = time.time()

        eccdna_bam.close()

        print((end-begin)/60)
