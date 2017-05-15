import pysam as ps
import os
import numpy as np
import pybedtools as bt
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


#TO-DO
# remove loop that checks for N reads in realignment step




class alignment:
    def __init__(self, bam, working_dir,genome_fa_dir,genome_name):
        self.bam = bam
        self.circ_bam = "circ_calls.bam"
        self.circ_bedgraph = "circ_dna_boundaries.bedGraph"
        self.working_dir = working_dir
        self.all_bam = "sorted_paired_end_sim_aln.bam"
        self.number_of_cores = 3
        self.circ_boundaries = bt.BedTool("circ_support_calls.bed")
        self.coverage = bt.BedTool("read_coverage_merged.bed")
        self.genome_dir = genome_fa_dir
        self.genome = genome_name

    def remove_concordant_pairs(self):
        """Function that removes concordant pairs"""
        os.system("samtools view -hb -F 2 %s > %s" % (self.bam, self.circ_bam))

    def query_name_sorted_circs(self):
        os.system("samtools sort -n -o query_name_sorted_circs.bam circ_calls.bam")
        os.system("samtools index circ_calls.bam ")
        sorted_circs = ps.AlignmentFile("query_name_sorted_circs.bam", "rb")
        return (sorted_circs)

    def check_mapping_range(self, interval, read):
        if read.reference_start >= (interval.start - 50):
            if read.reference_start <= (interval.end + 50):
                return (True)

        else:

            return (False)

    def is_soft_clipped(self, read):
        for cigar in read.cigar:
            if cigar[0] == 4:
                return (True)
            else:
                return (False)

    def generate_bed_from_circ_calls(self):
        # bedtools genomecov -bg -ibam circ_calls.bam | mergeBed > circ_support_calls.bed
        return (None)

    def split_to_cores(self,len):
        list = np.arange(len)
        cores_split = np.array_split(list,self.number_of_cores)



        return(cores_split)

    def filter_reads(self):
        """Function that aims to filters the reads based on mapq and soft-clipped/discordant reads"""




        # file of the sorted circles
        os.system("samtools index circ_calls.bam")
        circs = ps.AlignmentFile("circ_calls.bam", "rb")

        # bam file of all reads
        all_bam = ps.AlignmentFile(self.all_bam, "rb")

        # mapq_filter
        mapq_filtered = []
        each_one = 0
        for interval in self.coverage:
            if interval.end - interval.start < 1000:
                interval_mapq = []
                for read in all_bam.fetch(interval.chrom, interval.start, interval.end):
                    interval_mapq.append(read.mapq)

                if np.mean(interval_mapq) < 1:
                    continue
                else:
                    mapq_filtered.append(interval)
                    each_one += 1

            else:
                each_one += 1
                mapq_filtered.append(interval)
                continue

            print(each_one)

        mapq_filtered = bt.BedTool(mapq_filtered)
        mapq_filtered.saveas("mapq_filter.bed")

        filtered_intervals = []
        # soft_clip_filtering
        second_filters = 0
        for interval in mapq_filtered:
            suports = circs.count(interval.chrom, interval.start, interval.end)
            if (interval.end - interval.start) < 300:
                if suports > 4:
                    filtered_intervals.append(interval)
                    second_filters += 1
                    print(second_filters)
            else:
                if suports > 2:
                    filtered_intervals.append(interval)
                    second_filters += 1
                    print(second_filters)

        filtered_intervals = bt.BedTool(filtered_intervals)

        circ_boundaries = bt.BedTool("circ_support_calls.bed")

        length_filtered_intervals = len(filtered_intervals)



        # save the first file
        filtered_intervals.saveas("filtered_cutoff.bed")

        return(filtered_intervals)


    def realign_with_one_boundary(self,interval,circ_intervals,circ_bam,fastafile):
        split_read = 0
        discordant = 0

        for circ_interval in circ_intervals:


            # count reads in the interval
            count_reads = circ_bam.count(circ_interval.chrom, circ_interval.start, circ_interval.end)

            for read in circ_bam.fetch(circ_interval.chrom, circ_interval.start, circ_interval.end):

                # loop trough read
                if read.mapq > 10:
                    # filter by mapq

                    if alignment.is_soft_clipped(self, read) == True:
                        try:
                            # if there is SA (split alignment) we do not need to realign
                            suplementary = read.get_tag('SA')
                            # split the list to get the info
                            # this list will have the following information [chr,left_most start,"strand,CIGAR,mapq, edit_distance]
                            supl_info = [x.strip() for x in suplementary.split(',')]

                            # check the chromosome
                            if supl_info[0] == interval.chrom:
                                # aligned in the left boundary, supplementary on right
                                if read.reference_start == interval.start and interval.end - read.query_length <= supl_info[
                                    1] <= interval.end:
                                    split_read += 1
                                # aligned on the right boundary, supplementary alignment in the left
                                elif read.reference_end == interval.end and interval.start <= supl_info[
                                    1] <= interval.start + read.query_length:
                                    split_read += 1



                        except:
                            # realignment

                            if read.cigar[0][0] == 4:
                                # get the nucleotides of the beginning
                                nucleotides = read.cigar[0][1]
                                soft_clip_seq = read.seq[0:nucleotides]

                                if 'N' in soft_clip_seq:
                                    continue
                                else:

                                    if read.is_reverse == True:
                                        # get support interval seq
                                        seq = fastafile.fetch(circ_interval.chrom, circ_interval.start, circ_interval.end)
                                        # reverse it
                                        seq = Seq(seq, generic_dna)
                                        seq = seq.reverse_complement()

                                        # smith watermann
                                        pairwise_al = pairwise2.align.localms(seq.upper(), soft_clip_seq, 1, -1, -.5, -.1)


                                        #if the score is equal to read length

                                        if int(pairwise_al[0][2]) == len(soft_clip_seq):
                                            split_read += 1


                                    else:
                                        #read mapped to + strand
                                        seq = fastafile.fetch(circ_interval.chrom, circ_interval.start, circ_interval.end)
                                        #smith watermann
                                        pairwise_al = pairwise2.align.localms(seq.upper(), soft_clip_seq, 1, -1, -.5, -.1)

                                        if int(pairwise_al[0][2]) == len(soft_clip_seq):
                                            split_read += 1


                            elif read.cigar[-1][0] == 4:
                                # get nucleotides in the end
                                nucleotides = read.cigar[-1][1]
                                soft_clip_seq = read.seq[nucleotides:]

                                if 'N' in soft_clip_seq:
                                    continue

                                else:

                                    if read.is_reverse == True:
                                        # get support interval seq
                                        seq = fastafile.fetch(circ_interval.chrom, circ_interval.start,
                                                              circ_interval.end)
                                        # reverse it
                                        seq = Seq(seq, generic_dna)
                                        seq = seq.reverse_complement()
                                        #smith watermann
                                        pairwise_al = pairwise2.align.localms(seq.upper(), soft_clip_seq, 1, -1, -.5, -.1)
                                        if int(pairwise_al[0][2]) == len(soft_clip_seq):
                                            split_read += 1


                                    else:
                                        #get seq
                                        seq = fastafile.fetch(circ_interval.chrom, circ_interval.start,
                                                              circ_interval.end)

                                        #smith watermann
                                        pairwise_al = pairwise2.align.localms(seq.upper(), soft_clip_seq, 1, -1, -.5, -.1)

                                        if int(pairwise_al[0][2]) == len(soft_clip_seq):
                                            split_read += 1


                    else:
                        # check the mapping of the discordants
                        discordant += 1

                    return(split_read,discordant,circ_interval)


    def realign_with_two_boundaries(self,interval,circ_intervals,circ_bam,fastafile):
        """Function that realigns the reads of a pair of intervals"""
        split_read = 0
        discordant = 0
        print(circ_intervals[0])

        for read in circ_bam.fetch(circ_intervals[0].chrom, circ_intervals[0].start, circ_intervals[0].end):
            # loop trough read
            if read.mapq > 10:
            # filter by mapq

                if alignment.is_soft_clipped(self, read) == True:
                    try:
                        # if there is SA (split alignment) we do not need to realign
                        suplementary = read.get_tag('SA')
                        # split the list to get the info
                        # this list will have the following information [chr,left_most start,"strand,CIGAR,mapq, edit_distance]
                        supl_info = [x.strip() for x in suplementary.split(',')]

                        # check the chromosome
                        if supl_info[0] == circ_intervals[1].chrom:
                            # aligned in the left boundary, supplementary on right
                            if  circ_intervals[1].start<= supl_info[1] <= circ_intervals[1].end:
                                split_read += 1

                    except:
                        #realignment
                        if read.cigar[0][0] == 4:
                            # get the nucleotides of the beginning
                            nucleotides = read.cigar[0][1]
                            soft_clip_seq = read.seq[0:nucleotides]

                            if 'N' in soft_clip_seq:
                                continue
                            else:

                                if read.is_reverse == True:
                                    # get support interval seq
                                    seq = fastafile.fetch(circ_intervals[1].chrom, circ_intervals[1].start, circ_intervals[1].end)
                                    # reverse it
                                    seq = Seq(seq, generic_dna)
                                    seq = seq.reverse_complement()
                                    pairwise_al = pairwise2.align.localms(seq.upper(), soft_clip_seq, 1, -1, -.5, -.1)

                                    if int(pairwise_al[0][2]) == len(soft_clip_seq):
                                        split_read += 1


                                else:
                                    seq = fastafile.fetch(circ_intervals[1].chrom, circ_intervals[1].start, circ_intervals[1].end)
                                    pairwise_al = pairwise2.align.localms(seq.upper(), soft_clip_seq, 1, -1, -.5, -.1)

                                    if int(pairwise_al[0][2]) == len(soft_clip_seq):
                                        split_read += 1


                        elif read.cigar[-1][0] == 4:
                            # get nucleotides in the end
                            nucleotides = read.cigar[-1][1]
                            soft_clip_seq = read.seq[nucleotides:]

                            if 'N' in soft_clip_seq:
                                continue

                            else:

                                if read.is_reverse == True:
                                    # get support interval seq
                                    seq = fastafile.fetch(circ_intervals[1].chrom, circ_intervals[1].start,
                                                          circ_intervals[1].end)
                                    # reverse it
                                    seq = Seq(seq, generic_dna)
                                    seq = seq.reverse_complement()
                                    pairwise_al = pairwise2.align.localms(seq.upper(), soft_clip_seq, 1, -1, -.5, -.1)
                                    if int(pairwise_al[0][2]) == len(soft_clip_seq):
                                        split_read += 1


                                else:
                                    seq = fastafile.fetch(circ_intervals[1].chrom, circ_intervals[1].start,
                                                          circ_intervals[1].end)

                                    pairwise_al = pairwise2.align.localms(seq.upper(), soft_clip_seq, 1, -1, -.5, -.1)

                                    if int(pairwise_al[0][2]) == len(soft_clip_seq):
                                        split_read += 1


                else:
                    discordant +=1

            split_read = 0
            discordant = 0


        for read in circ_bam.fetch(circ_intervals[1].chrom, circ_intervals[1].start, circ_intervals[1].end):
            # loop trough read
            if read.mapq > 10:
                # filter by mapq

                if alignment.is_soft_clipped(self, read) == True:
                    try:
                        # if there is SA (split alignment) we do not need to realign
                        suplementary = read.get_tag('SA')
                        # split the list to get the info
                        # this list will have the following information [chr,left_most start,"strand,CIGAR,mapq, edit_distance]
                        supl_info = [x.strip() for x in suplementary.split(',')]

                        # check the chromosome
                        if supl_info[0] == circ_intervals[0].chrom:
                            # aligned in the left boundary, supplementary on right
                            if circ_intervals[0].start <= supl_info[0] <= circ_intervals[0].end:
                                split_read += 1

                    except:
                        # realignment
                        if read.cigar[0][0] == 4:
                            # get the nucleotides of the beginning
                            nucleotides = read.cigar[0][1]
                            soft_clip_seq = read.seq[0:nucleotides]

                            if 'N' in soft_clip_seq:
                                continue
                            else:

                                if read.is_reverse == True:
                                    # get support interval seq
                                    seq = fastafile.fetch(circ_intervals[0].chrom, circ_intervals[0].start,
                                                          circ_intervals[0].end)
                                    # reverse it
                                    seq = Seq(seq, generic_dna)
                                    seq = seq.reverse_complement()
                                    pairwise_al = pairwise2.align.localms(seq.upper(), soft_clip_seq, 1, -1, -.5,
                                                                          -.1)

                                    if int(pairwise_al[0][2]) == len(soft_clip_seq):
                                        split_read += 1


                                else:

                                    seq = fastafile.fetch(circ_intervals[0].chrom, circ_intervals[0].start,
                                                          circ_intervals[0].end)

                                    pairwise_al = pairwise2.align.localms(seq.upper(), soft_clip_seq, 1, -1, -.5,
                                                                          -.1)


                                    if int(pairwise_al[0][2]) == len(soft_clip_seq):
                                        split_read += 1


                        elif read.cigar[-1][0] == 4:
                            # get nucleotides in the end
                            nucleotides = read.cigar[-1][1]
                            soft_clip_seq = read.seq[nucleotides:]

                            if 'N' in soft_clip_seq:
                                continue

                            else:

                                if read.is_reverse == True:
                                    # get support interval seq
                                    seq = fastafile.fetch(circ_intervals[0].chrom, circ_intervals[0].start,
                                                          circ_intervals[0].end)
                                    # reverse it
                                    seq = Seq(seq, generic_dna)
                                    seq = seq.reverse_complement()
                                    pairwise_al = pairwise2.align.localms(seq.upper(), soft_clip_seq, 1, -1, -.5,
                                                                          -.1)
                                    if int(pairwise_al[0][2]) == len(soft_clip_seq):
                                        split_read += 1


                                else:
                                    seq = fastafile.fetch(circ_intervals[0].chrom, circ_intervals[0].start,
                                                          circ_intervals[0].end)

                                    pairwise_al = pairwise2.align.localms(seq.upper(), soft_clip_seq, 1, -1, -.5,
                                                                          -.1)

                                    if int(pairwise_al[0][2]) == len(soft_clip_seq):
                                        split_read += 1


                else:
                    discordant += 1

        return(split_read,discordant,interval)


    def realignment(self,circ_boundaries,list):

        circles = 0


        circ_bam = ps.AlignmentFile(self.circ_bam, "rb")

        filtered_intervals = bt.BedTool("filtered_cutoff.bed")

        os.chdir(self.genome_dir)

        fastafile = ps.FastaFile("%s" % self.genome)

        os.chdir(self.working_dir)

        for i in range(list[0],list[-1]):
            # each interval to realign/analyze
            interval = filtered_intervals[i]
            # number of boundaries in the interval
            boundaries = circ_boundaries.count_hits(interval)

            # get the intervals of the boundaries
            each_overlapping_interval = circ_boundaries.all_hits(interval)

            intervals = bt.BedTool(each_overlapping_interval)

            results = []

            f = open('test_status.txt', 'a')

            # analize the intervals that have only one overlapping boundary
            if boundaries == 1:
                #realign_one_interval(self,interval,circ_intervals,circ_bam,fastafile):
                realigned_interval = alignment.realign_with_one_boundary(self,interval,intervals,circ_bam,fastafile)

                try:

                    if realigned_interval[0] + realigned_interval[1] >= 4:

                        circles +=1
                        results.append(realigned_interval[2])
                        each_line = ["a", "\n"]
                        join_lines = ' '.join(map(str, each_line))
                        f.write(join_lines)
                except:
                    continue

            elif boundaries ==2:
                realigned_boundaries = alignment.realign_with_two_boundaries(self,interval,intervals,circ_bam,fastafile)

                try:
                    if realigned_boundaries[0] + realigned_boundaries[1] >= 4:
                        circles +=1
                        results.append(realigned_boundaries[2])
                        each_line = ["a", "\n"]
                        join_lines = ' '.join(map(str, each_line))
                        f.write(join_lines)


                except:
                    continue


            f.close()








        return(None)


