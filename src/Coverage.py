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
import os
import pysam as ps
import pybedtools as bt
import numpy as np


class coverage:
    """Class for managing the coverage metrics of circle-map"""

    def __init__(self,sorted_bam,eccdna_bed,extension,mapq,inside_length,directory):

        self.bam = ps.AlignmentFile(directory + "/" + sorted_bam, "rb")
        self.bed = eccdna_bed

        #length of out
        self.ext = extension
        self.mapq = mapq

        #length of the region for the ratio
        self.ilen = inside_length

        def print_parameters(self):
            print("Running coverage computations \n")



    def get_wg_coverage(self):
        """Generator that takes as input a sorted bam and a merged bam of the circles in the whole genome and returns a numpy
        array for every interval with the coverage"""



        reference_contigs = self.bam.header['SQ']

        header_dict = {}
        for reference in reference_contigs:
            header_dict[reference['SN']] = reference['LN']

        merged_bed = self.bed.sort().merge()

        for interval in merged_bed:

            coverage_dict = {}
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

            print("Computing coverage on interval %s:%s-%s" % (interval.chrom,interval.start,interval.end))
            coverage_dict[bt.Interval(interval.chrom, start, end)] = summ_cov

            yield(coverage_dict,header_dict)


    def compute_coverage(self,cov_generator):


        """Function that takes as input generator returning  coverage numpy arrays and  file with summarized statistics
        of the coverage within the intervals"""

        print("Computing the coverage of the identified eccDNA")
        print("Merging intervals for coverage computation")

        output = []
        for cov_dict,header_dict in cov_generator:
            for key,value in cov_dict.items():


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






                    try:

                        mean = np.mean(region_array)
                        sd = np.std(region_array)

                        interval.append(str(mean))
                        interval.append(str(sd))

                    except:

                        interval.append('NA')
                        interval.append('NA')


                    # compute ratios

                    try:

                        start_coverage_ratio = np.sum(region_array[0:self.ilen]) / np.sum(
                            ext_array[0:(self.ilen + self.ext)])
                        end_coverage_ratio = np.sum(region_array[-self.ilen:]) / np.sum(ext_array[-(self.ilen + self.ext):])

                        interval.append(str(start_coverage_ratio))
                        interval.append(str(end_coverage_ratio))

                    except:

                        interval.append('NA')
                        interval.append('NA')

                    try:


                        zero_frac = np.count_nonzero(region_array == 0) / len(region_array)
                        interval.append(str(zero_frac))

                    except:

                        interval.append('NA')
                    output.append(interval)


        return(bt.BedTool(output))
