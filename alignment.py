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
    def __init__(self, bam,circ_bam, working_dir,genome_fa_dir,genome_name,number_of_cores):
        self.init_bam = circ_bam
        self.circ_bam = "sorted_" + circ_bam
        self.working_dir = working_dir
        self.sort_circ_reads()
        self.all_bam = bam
        self.index_bams()
        self.number_of_cores = number_of_cores
        self.generate_bed_from_bams()
        self.circ_boundaries = bt.BedTool("circ_supports.bed")
        self.coverage = bt.BedTool("read_coverage_merged.bed")
        self.genome_dir = genome_fa_dir
        self.genome = genome_name

    def sort_circ_reads(self):
        os.system("samtools sort %s > %s" % (self.init_bam,self.circ_bam))

    def index_bams(self):
        os.system("samtools index %s ; samtools index %s" % (self.circ_bam, self.all_bam))


    def is_soft_clipped(self, read):
        for cigar in read.cigar:
            if cigar[0] == 4:
                return (True)
            else:
                return (False)

    def generate_bed_from_bams(self):
        os.system("bedtools genomecov -bg -ibam %s | mergeBed | sortBed > %s" % (self.circ_bam,"circ_supports.bed"))
        os.system("bedtools genomecov -bg -ibam %s | mergeBed | sortBed > %s" % (self.all_bam, "read_coverage_merged.bed"))
        return (None)

    def split_to_cores(self,len):
        list = np.arange(len)
        cores_split = np.array_split(list,self.number_of_cores)



        return(cores_split)

    def filter_reads(self):
        """Function that aims to filters the reads based on mapq and soft-clipped/discordant reads"""



        # file of the sorted circles
        circs = ps.AlignmentFile(self.circ_bam, "rb")

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

        # save the first file


        return(mapq_filtered)


    def realign_with_one_boundary(self,interval,circ_intervals,circ_bam,fastafile):
        split_read = 0
        discordant = 0

        for circ_interval in circ_intervals:


            # count reads in the interval
            count_reads = circ_bam.count(circ_interval.chrom, circ_interval.start, circ_interval.end)

            for read in circ_bam.fetch(circ_interval.chrom, circ_interval.start, circ_interval.end):

                # loop trough read
                if read.mapq > 0:
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
            if read.mapq > 0:
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
            if read.mapq > 0:
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




    def realignment(self,circ_boundaries,list,temp_bed):

        f = open("test_status.txt","a")

        circ_bam = ps.AlignmentFile(self.circ_bam, "rb")

        filtered_intervals = bt.BedTool("mapq_filter.bed")



        os.chdir(self.genome_dir)

        fastafile = ps.FastaFile("%s" % self.genome)

        results = []

        os.chdir(self.working_dir)

        for i in range(list[0],list[-1]):
            # each interval to realign/analyze
            interval = filtered_intervals[i]
            # number of boundaries in the interval
            boundaries = circ_boundaries.count_hits(interval)

            # get the intervals of the boundaries
            each_overlapping_interval = circ_boundaries.all_hits(interval)

            intervals = bt.BedTool(each_overlapping_interval)






            # analize the intervals that have only one overlapping boundary
            if boundaries == 1:

                #realign_one_interval(self,interval,circ_intervals,circ_bam,fastafile):
                realigned_interval = alignment.realign_with_one_boundary(self,interval,intervals,circ_bam,fastafile)

                try:

                    if realigned_interval[0] + realigned_interval[1] >= 4:


                        results.append(realigned_interval[2])
                        each_line = ["a", "\n"]
                        join_lines = ' '.join(map(str, each_line))
                        f.write(join_lines)
                        print("one boundary correct")
                except:
                    print("error in one boundary")
                    continue

            elif boundaries ==2:

                realigned_boundaries = alignment.realign_with_two_boundaries(self,interval,intervals,circ_bam,fastafile)

                try:
                    if realigned_boundaries[0] + realigned_boundaries[1] >= 4:
                        results.append(realigned_boundaries[2])
                        each_line = ["a", "\n"]
                        join_lines = ' '.join(map(str, each_line))
                        f.write(join_lines)
                        print("two boundaries correct")


                except:
                    print("error in two boundaries")
                    continue

            else:


                while len(intervals) > 1:




                    for each_interval in intervals:

                        #create scoring list

                        interval_list = []
                        scoring_list = []
                        for each_interval in intervals:
                            interval_list.append(each_interval)
                            scoring_list.append(0)






                        for read in circ_bam.fetch(each_interval.chrom, each_interval.start, each_interval.end):

                            next_read = read.next_reference_id

                            if next_read == -1:
                                continue

                            else:
                                mate_map_chr = circ_bam.get_reference_name(next_read)
                                mate_map_pos = read.next_reference_start

                                for element in range(0,len(interval_list)):
                                    if mate_map_chr == interval_list[element].chrom:
                                        if interval_list[element].start < mate_map_pos < interval_list[element].end:
                                            scoring_list[element] +=1






                    #print(scoring_list)
                    index_number = scoring_list.index(max(scoring_list))
                    index = [0,index_number]
                    pair = [ interval_list[i] for i in index]
                    interval_list = [i for j, i in enumerate(interval_list) if j not in index]
                    intervals = bt.BedTool(interval_list)
                    pair = bt.BedTool(pair)
                    realigned_boundaries = alignment.realign_with_two_boundaries(self,interval,pair,circ_bam,fastafile)

                    try:
                        if realigned_boundaries[0] + realigned_boundaries[1] >= 4:
                            results.append(realigned_boundaries[2])
                            each_line = ["a", "\n"]
                            join_lines = ' '.join(map(str, each_line))
                            f.write(join_lines)
                            print("multiple boundaries correct")


                    except:
                       continue


                else:
                    realigned_interval = alignment.realign_with_one_boundary(self, interval, intervals, circ_bam,
                                                                             fastafile)

                    try:

                        if realigned_interval[0] + realigned_interval[1] >= 4:

                            results.append(realigned_interval[2])
                            each_line = ["a", "\n"]
                            join_lines = ' '.join(map(str, each_line))
                            f.write(join_lines)
                            print("1 boundary correct")
                    except:
                        print("one boundary not correct")
                        continue



        f.close()

        results = bt.BedTool(results)

        results.saveas(str(temp_bed))




        return(None)


