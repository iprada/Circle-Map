import pysam as ps
import os
import pandas as pd
import numpy as np
import pybedtools as bt
import multiprocessing as mp


class alignment:
    def __init__(self, bam, working_dir):
        self.bam = bam
        self.circ_bam = "circ_calls.bam"
        self.circ_bedgraph = "circ_dna_boundaries.bedGraph"
        self.working_dir = working_dir
        self.all_bam = "sorted_paired_end_sim_aln.bam"
        self.number_of_cores = 3
        self.circ_boundaries = bt.BedTool("read_coverage_merged.bed")

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

    def realignment(self,circ_boundaries,list):

        filtered_intervals = bt.BedTool("filtered_cutoff.bed")


        for i in range(list[0],list[-1]):
            interval = filtered_intervals[i]
            boundaries = circ_boundaries.count_hits(interval)
            each_overlapping_interval = circ_boundaries.all_hits(interval)
            i += 1
            print(i)



        return("a")



    def call_circles(self):
        """Function that aims to call circles based on the boundaries identified
        using the discordant reads and soft clipped reads"""




        # file of the sorted circles
        os.system("samtools index circ_calls.bam")
        circs = ps.AlignmentFile("circ_calls.bam", "rb")

        # bam file of all reads
        all_bam = ps.AlignmentFile(self.all_bam, "rb")

        # mapq_filter
        mapq_filtered = []
        each_one = 0
        for interval in self.circ_boundaries:
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
                if suports > 0:
                    filtered_intervals.append(interval)
                    second_filters += 1
                    print(second_filters)
            else:
                if suports > 0:
                    filtered_intervals.append(interval)
                    second_filters += 1
                    print(second_filters)

        filtered_intervals = bt.BedTool(filtered_intervals)

        circ_boundaries = bt.BedTool("circ_support_calls.bed")

        length_filtered_intervals = len(filtered_intervals)



        # save the first file
        filtered_intervals.saveas("filtered_cutoff.bed")

        return(filtered_intervals)