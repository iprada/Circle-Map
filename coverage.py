import os
import pysam as ps
import pybedtools as bt
import numpy as np


class coverage:
    """Class for managing the coverage metrics of circle-map"""

    def __init__(self,sorted_bam,eccdna_bed,extension,mapq,inside_length):

        self.bam = ps.AlignmentFile(sorted_bam, "rb")
        self.bed = eccdna_bed

        #length of out
        self.ext = extension
        self.mapq = mapq

        #length of the region for the ratio
        self.ilen = inside_length

    def get_wg_coverage(self):
        """Function that takes as input a sorted bam and a merged bam of the circles in the whole genome and returns a numpy
        array for every interval with the coverage"""

        print("Computing whole genome coverage")

        reference_contigs = self.bam.header['SQ']

        header_dict = {}
        for reference in reference_contigs:
            header_dict[reference['SN']] = reference['LN']

        merged_bed = self.bed.sort().merge()

        print("%s merged intervals for coverage computation " % len(merged_bed))


        coverage_dict = {}
        for interval in merged_bed:

            print(interval)

            if interval.start - self.ext < 0:
                start = 0

            else:
                start = interval.start - self.ext

            if header_dict[interval.chrom] < (interval.end + self.ext):
                end = interval.end + self.ext
            else:
                end = interval.end

            cov = self.bam.count_coverage(contig=interval.chrom, start=start, end=end, quality_threshold=self.mapq)
            summarized_cov = np.array([cov[0], cov[1], cov[2], cov[3]]).sum(axis=0)

            # save memory, convert to uint32.
            summ_cov = np.uint32(summarized_cov)

            coverage_dict[bt.Interval(interval.chrom, start, end)] = summ_cov

        return (coverage_dict,header_dict)


    def compute_coverage(self,cov_dict,header_dict):

        """Function that takes as input a sorted bam file and a bed file a bed file with summarized statistics of the cove-
        rage within the intervals"""

        print("Computing the coverage of the identified eccDNA")

        output = []
        for key,value in cov_dict.items():

            print(key)


            overlaps = bt.BedTool(self.bed.all_hits(key))


            for interval in overlaps:

                # compute array slicing indices
                start = interval.start -key.start
                end = interval.end - key.start


                if start - self.ext < 0:
                    ext_start = 0
                else:
                    ext_start = start - self.ext

                if header_dict[interval.chrom] < (end+ self.ext):
                    ext_end = header_dict[interval.chrom]
                else:
                    ext_end = end + self.ext

                # slice extended array and coverage array
                ext_array = value[ext_start:ext_end]
                region_array = value[start:end]


                mean = np.mean(region_array)
                sd = np.std(region_array)


                #compute ratios
                start_coverage_ratio = np.sum(region_array[0:self.ilen]) / np.sum(ext_array[0:(self.ilen+self.ext)])
                end_coverage_ratio = np.sum(region_array[-self.ilen:]) / np.sum(ext_array[-(self.ilen + self.ext):])


                zero_frac = np.count_nonzero(region_array == 0) / len(region_array)

                interval.append(str(mean))
                interval.append(str(sd))

                interval.append(str(start_coverage_ratio))
                interval.append(str(end_coverage_ratio))

                interval.append(str(zero_frac))

                output.append(interval)

        return(bt.BedTool(output))

