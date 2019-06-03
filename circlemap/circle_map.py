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


import argparse
import sys
from functools import partial
import os
import time
import pandas as pd
from circlemap.extract_circle_SV_reads import readExtractor
from circlemap.realigner import realignment
from circlemap.repeats import repeat
from circlemap.utils import merge_final_output, filter_by_ratio, start_realign, start_simulate, mutate, insert_size_dist
from circlemap.Coverage import coverage
import multiprocessing as mp
import pybedtools as bt
from circlemap.simulations import sim_ecc_reads
import subprocess as sp
import glob
from tqdm import *


class circle_map:

    def __getpid__(self):

        pid = os.getpid()
        return (pid)

    def __init__(self):
        self.parser = argparse.ArgumentParser(
            description='Circle-Map',
            usage='''Circle-Map <subprogram> [options]

version=1.0
contact= https://github.com/iprada/Circle-Map/issues

The Circle-Map suite

Commands:

   ReadExtractor   Extracts circular DNA read candidates
   Realign         Realign circular DNA read candidates
   Repeats         Identify circular DNA from repetitive regions
   Simulate        Simulate circular DNA

''')
        subparsers = self.parser.add_subparsers()

        self.readextractor = subparsers.add_parser(
            name="ReadExtractor",
            description='Extracts circular DNA read candidates',
            prog="Circle-Map ReadExtractor",
            usage='''Circle-Map ReadExtractor [options]'''

        )

        self.realigner = subparsers.add_parser(
            name="Realign",
            description='Realign circular DNA read candidates',
            prog="Circle-Map Realign",
            usage='''Circle-Map Realign [options]'''

        )

        self.repeats = subparsers.add_parser(
            name="Repeats",
            description='Identify circular DNA from repetitive regions',
            prog="Circle-Map Repeats",
            usage='''Circle-Map Repeats [options]'''

        )
        self.simulate = subparsers.add_parser(
            name="Simulate",
            description='Simulate eccDNA NGS datastes',
            prog="Circle-Map Reepeats",
            usage='''Circle-Map Simulate [options]'''

        )

        if len(sys.argv) <= 1:
            self.parser.print_help()
            time.sleep(0.01)
            sys.stderr.write("\nNo argument given to Circle-Map"
                             "\nExiting\n")
            sys.exit(1)

        else:
            if sys.argv[1] == "ReadExtractor":


                self.subprogram = self.args_readextractor()
                self.args = self.subprogram.parse_args(sys.argv[2:])

                object = readExtractor(self.args.i, self.args.output, self.args.directory, self.args.quality,
                                       self.args.nodiscordant,
                                       self.args.nohardclipped, self.args.nosoftclipped, self.args.verbose,
                                       self.subprogram)
                object.extract_sv_circleReads()

            elif sys.argv[1] == "Realign":
                self.subprogram = self.args_realigner()
                self.args = self.subprogram.parse_args(sys.argv[2:])

                os.chdir(self.args.directory)
                # get clusters
                splitted, sorted_bam, begin = start_realign(self.args.i, self.args.output, self.args.threads,
                                                            self.args.verbose, self.__getpid__(),
                                                            self.args.clustering_dist)




                sorted_bam.close()
                #get global insert size prior
                metrics = insert_size_dist(self.args.sample_size, self.args.insert_mapq, self.args.qbam)


                # pool based parallel of religment
                m = mp.Manager()

                lock = m.Lock()

                object = realignment(self.args.i, self.args.qbam, self.args.sbam, self.args.fasta,
                                     self.args.directory,
                                     self.args.mapq,
                                     self.args.insert_mapq, self.args.std, self.args.sample_size,
                                     self.args.gap_open,
                                     self.args.gap_ext, self.args.nhits, self.args.cut_off, self.args.min_sc,
                                     self.args.merge_fraction, self.args.interval_probability, self.args.output,
                                     self.args.threads, self.args.allele_frequency, lock, self.args.split,
                                     self.args.ratio, self.args.verbose, self.__getpid__(),
                                     self.args.edit_distance_fraction, self.args.remap_splits,
                                     self.args.only_discordants, self.args.split,
                                     self.args.split_quality, metrics,self.args.number_of_discordants)


                pool = mp.Pool(processes=self.args.threads)

                #progress bar
                with tqdm(total=len(splitted)) as pbar:
                    for i,exits in tqdm(enumerate(pool.imap_unordered(object.realign, splitted))):
                        pbar.update()

                pbar.close()
                pool.close()
                pool.join()
                output = merge_final_output(self.args.sbam, self.args.output, begin, self.args.split,
                                            self.args.directory,
                                            self.args.merge_fraction, self.__getpid__())

                # compute coverage statistics
                if self.args.no_coverage == False:

                    coverage_object = coverage(self.args.sbam, output,
                                               self.args.bases, self.args.cmapq, self.args.extension,
                                               self.args.directory)

                    #Generator function for the coverage calculations
                    output = coverage_object.compute_coverage(coverage_object.get_wg_coverage())
                    filtered_output = filter_by_ratio(output, self.args.ratio)
                    filtered_output.to_csv(r'%s' % self.args.output, header=None, index=None, sep='\t', mode='w')

                else:
                    output.saveas("%s" % self.args.output)

            elif sys.argv[1] == "Repeats":

                self.subprogram = self.args_repeats()
                self.args = self.subprogram.parse_args(sys.argv[2:])


                object = repeat(self.args.i, self.args.directory, self.args.mismatch, self.args.fraction,
                                self.args.read_number)

                bed = object.find_circles()

                coverage_object = coverage(self.args.sbam, bed,
                                           self.args.bases, self.args.cmapq, self.args.extension,
                                           self.args.directory)

                output = coverage_object.compute_coverage(coverage_object.get_wg_coverage())

                filtered_output = filter_by_ratio(output, self.args.ratio)
                filtered_output.saveas("%s" % self.args.output)


            elif sys.argv[1] == "Simulate":

                self.subprogram = self.args_simulate()
                self.args = self.subprogram.parse_args(sys.argv[2:])

                sim_pid = start_simulate(self.__getpid__())

                lock = mp.Lock()

                paired_end_fastq_1 = open("%s_1.fastq" % self.args.base_name, "w")
                paired_end_fastq_2 = open("%s_2.fastq" % self.args.base_name, "w")
                paired_end_fastq_1.close()
                paired_end_fastq_2.close()

                # mutate reference genome
                if self.args.variants == True:
                    mutate(self.args.g, sim_pid, self.args.Indels, self.args.substitution, self.args.java_memory)


                manager = mp.Manager()
                # Shared memory object
                circle_list = manager.list()
                skipped_circles = mp.Value('i', 0)
                correct_circles = mp.Value('i', 0)
                jobs = []
                # init the processes

                for i in range(self.args.processes):
                    p = mp.Process(target=sim_ecc_reads,
                                   args=(self.args.g, self.args.read_length, self.args.directory,
                                         int(round(self.args.read_number / self.args.processes)),
                                         self.args.skip_region, self.args.base_name,
                                         self.args.mean_insert_size, self.args.error,
                                         self.args.mean_coverage, lock, i, circle_list,
                                         "%s_1.fastq" % self.args.base_name, "%s_2.fastq" % self.args.base_name,
                                         skipped_circles,
                                         correct_circles, self.args.insRate, self.args.insRate2, self.args.delRate,
                                         self.args.delRate2, sim_pid,))
                    jobs.append(p)
                    p.start()
                # kill the process
                for p in jobs:
                    p.join()
                print("Skipped %s circles, that overlapped the provided regions to exclude" % skipped_circles.value)
                print("Simulated %s circles across %s parallel processes" % (
                correct_circles.value, self.args.processes))
                print("Writting to disk bed file containing the simulated circle coordinates")

                bt.BedTool(list(circle_list)).saveas(self.args.output)

            else:
                self.parser.print_help()
                time.sleep(0.01)
                sys.stderr.write("\nWrong argument given to Circle-Map"
                                 "\nExiting\n")
                sys.exit(1)

    def args_readextractor(self):

        parser = self.readextractor

        parser._action_groups.pop()
        required = parser.add_argument_group('required arguments')
        optional = parser.add_argument_group('optional arguments')
        # prefixing the argument with -- means it's optional
        # input and output

        required.add_argument('-i', metavar='', help="Input: query name sorted bam file")

        if "-i" in sys.argv:
            optional.add_argument('-o', '--output', metavar='',
                                  help="Ouput: Reads indicating circular DNA structural variants",
                                  default="circle_%s" % sys.argv[sys.argv.index("-i") + 1])

            optional.add_argument('-dir', '--directory', metavar='',
                                  help="Working directory, default is the working directory",
                                  default=os.getcwd())

            # mapping quality cutoff

            optional.add_argument('-q', '--quality', type=int, metavar='',
                                  help="bwa-mem mapping quality cutoff. Default value 10",
                                  default=10)

            # read extraction options
            # extract discordant reads
            optional.add_argument('-nd', '--nodiscordant', help="Turn off discordant (R2F1 oriented) read extraction",
                                  action='store_true')

            # soft-clipped argument
            optional.add_argument('-nsc', '--nosoftclipped', help="Turn off soft-clipped read extraction",
                                  action='store_true')
            # extract hard-clippped reads
            optional.add_argument('-nhc', '--nohardclipped', help="Turn off hard-clipped read extraction",
                                  action='store_true')

            # verbose level

            optional.add_argument('-v', '--verbose', type=int, metavar='',
                                  help='Verbose level, 1=error,2=warning, 3=message',
                                  choices=[1, 2, 3], default=3)

        else:
            optional.add_argument('-o', '--output', metavar='',
                                  help="Ouput: Reads indicating circular DNA structural variants")

            optional.add_argument('-dir', '--directory', metavar='',
                                  help="Working directory, default is the working directory",
                                  default=os.getcwd())

            # mapping quality cutoff
            optional.add_argument('-q', '--quality', type=int, metavar='',
                                  help="bwa-mem mapping quality cutoff. Default value 10",
                                  default=10)

            # read extraction options
            # extract discordant reads
            optional.add_argument('-nd', '--nodiscordant', help="Turn off discordant (R2F1 oriented) read extraction",
                                  action='store_true')

            # soft-clipped argument
            optional.add_argument('-nsc', '--nosoftclipped', help="Turn off soft-clipped read extraction",
                                  action='store_true')
            # extract hard-clippped reads
            optional.add_argument('-nhc', '--nohardclipped', help="Turn off hard-clipped read extraction",
                                  action='store_true')

            # verbose level

            optional.add_argument('-v', '--verbose', type=int, metavar='',
                                  help='Verbose level, 1=error,2=warning, 3=message. Default=3',
                                  choices=[1, 2, 3], default=3)

            parser.print_help()

            time.sleep(0.01)
            sys.stderr.write(
                "\nNo input or output input given to readExtractor, be sure that you are providing the flags'-i' and '-o'"
                "\nExiting\n")
            sys.exit(1)

        # parse the commands

        if len(sys.argv[2:]) == 0:
            parser.print_help()
            time.sleep(0.01)
            sys.stderr.write("\nNo arguments given to read extractor. Exiting\n")
            sys.exit(1)

        return (parser)

    def args_realigner(self):
        parser = self.realigner

        # declare the different groups for the parser
        parser._action_groups.pop()
        io_options = parser.add_argument_group('Input/Output options')
        alignment_options = parser.add_argument_group('Alignment options')
        out_decision = parser.add_argument_group('eccDNA output options')
        i_size_estimate = parser.add_argument_group('Insert size estimation options')
        interval = parser.add_argument_group('Interval processing options')
        coverage_metrics = parser.add_argument_group('Coverage metrics options')
        running = parser.add_argument_group('Running options')

        io_options.add_argument('-i', metavar='',
                                help="Input: bam file containing the reads extracted by ReadExtractor")
        io_options.add_argument('-qbam', metavar='', help="Input: query name sorted bam file")
        io_options.add_argument('-sbam', metavar='', help="Input: coordinate sorted bam file")
        io_options.add_argument('-fasta', metavar='', help="Input: Reference genome fasta file")

        if "-i" and "-qbam" and "-fasta" in sys.argv:
            # output

            io_options.add_argument('-o', '--output', metavar='', help="Output filename",
                                    default="circle_%s.bed" % sys.argv[sys.argv.index("-i") + 1])

            # alignment
            alignment_options.add_argument('-n', '--nhits', type=int, metavar='',
                                           help="Number of realignment attempts. Default: 10",
                                           default=10)

            alignment_options.add_argument('-p', '--cut_off', type=float, metavar='',
                                           help="Probability cut-off for considering a soft-clipped as realigned: Default: 0.99",
                                           default=0.99)

            alignment_options.add_argument('-m', '--min_sc', type=float, metavar='',
                                           help="Minimum soft-clipped length to attempt the realignment. Default: 6",
                                           default=6)

            alignment_options.add_argument('-g', '--gap_open', type=int, metavar='',
                                           help="Gap open penalty in the position specific scoring matrix. Default: 5",
                                           default=5)

            alignment_options.add_argument('-e', '--gap_ext', type=int, metavar='',
                                           help="Gap extension penalty in the position specific scoring matrix. Default: 1",
                                           default=1)

            alignment_options.add_argument('-q', '--mapq', type=int, metavar='',
                                           help="Minimum mapping quality allowed in the supplementary alignments. Default: 20",
                                           default=20)

            alignment_options.add_argument('-d', '--edit_distance-fraction', type=float, metavar='',
                                           help="Maximum edit distance fraction allowed in the first realignment. Default (0.05)",
                                           default=0.05)

            alignment_options.add_argument('-Q', '--split_quality', type=float, metavar='',
                                           help="Minium split score to output an interval. Default (0.0)",
                                           default=0.0)

            alignment_options.add_argument('-R', '--remap_splits', help="Remap probabilistacally the split reads",
                                           action='store_true')

            # insert size

            i_size_estimate.add_argument('-iq', '--insert_mapq', type=int, metavar='',
                                         help="Mapq cutoff for stimating the insert size distribution. Default 60",
                                         default=60)

            i_size_estimate.add_argument('-sd', '--std', type=int, metavar='',
                                         help="Standard deviations of the insert size to extend the intervals. Default 5",
                                         default=4)

            i_size_estimate.add_argument('-s', '--sample_size', type=int, metavar='',
                                         help="Number of concordant reads (R2F1) to use for estimating the insert size distribution. Default 100000",
                                         default=100000)

            # Interval options

            interval.add_argument('-f', '--merge_fraction', type=float, metavar='',
                                  help="Merge intervals reciprocally overlapping by a fraction. Default 0.95",
                                  default=0.95)

            interval.add_argument('-P', '--interval_probability', type=float, metavar='',
                                  help="Skip edges of the graph with a probability below the threshold. Default: 0.01",
                                  default=0.01)
            interval.add_argument('-K', '--clustering_dist', type=int, metavar='',
                                  help="Cluster reads that are K nucleotides appart in the same node. Default: 500",
                                  default=500)

            interval.add_argument('-D', '--only_discordants', help="Use only discordant reads to build the graph",
                                  action='store_false')
            interval.add_argument('-F', '--allele_frequency', type=float, metavar='',
                                  help="Minimum allele frequency required to report the circle interval. Default (0.1)",
                                  default=0.1)
            # When to call a circle

            out_decision.add_argument('-S', '--split', type=int, metavar='',
                                      help="Number of required split reads to output a eccDNA. Default: 0",
                                      default=0)

            out_decision.add_argument('-O', '--number_of_discordants', type=int, metavar='',
                                      help="Number of required discordant reads for intervals with only discordants. Default: 3",
                                      default=3)
            out_decision.add_argument('-r', '--ratio', type=float, metavar='',
                                      help="Minimum in/out required coverage ratio. Default: 0.0",
                                      default=0.0)

            # coverage metrics

            coverage_metrics.add_argument('-N', '--no_coverage', help="Don't compute coverage statistics",
                                          action='store_true')

            coverage_metrics.add_argument('-b', '--bases', type=int, metavar='',
                                          help="Number of bases to extend for computing the coverage ratio. Default: 200",
                                          default=200)

            coverage_metrics.add_argument('-cq', '--cmapq', type=int, metavar='',
                                          help="Minimum mapping quality treshold for coverage computation. Default: 0",
                                          default=0)

            coverage_metrics.add_argument('-E', '--extension', type=int, metavar='',
                                          help="Number of bases inside the eccDNA breakpoint coordinates to compute the ratio. Default: 100",
                                          default=100)

            # run options

            running.add_argument('-t', '--threads', type=int, metavar='',
                                 help="Number of threads to use.Default 1",
                                 default=1)

            running.add_argument('-dir', '--directory', metavar='',
                                 help="Working directory, default is the working directory",
                                 default=os.getcwd())

            running.add_argument('-v', '--verbose', type=int, metavar='',
                                 help='Verbose level, 1=error,2=warning, 3=message',
                                 choices=[1, 2, 3], default=3)



        else:

            # output

            io_options.add_argument('-o', metavar='', help="Output filename")

            alignment_options.add_argument('-n', '--nhits', type=int, metavar='',
                                           help="Number of realignment attempts. Default: 10",
                                           default=10)

            alignment_options.add_argument('-p', '--cut_off', type=float, metavar='',
                                           help="Probability cut-off for considering a soft-clipped as realigned: Default: 0.99",
                                           default=0.99)

            alignment_options.add_argument('-m', '--min_sc', type=float, metavar='',
                                           help="Minimum soft-clipped length to attempt the realignment. Default: 6",
                                           default=6)

            alignment_options.add_argument('-g', '--gap_open', type=int, metavar='',
                                           help="Gap open penalty in the position specific scoring matrix. Default: 5",
                                           default=5)

            alignment_options.add_argument('-e', '--gap_ext', type=int, metavar='',
                                           help="Gap extension penalty in the position specific scoring matrix. Default: 1",
                                           default=1)

            alignment_options.add_argument('-q', '--mapq', type=int, metavar='',
                                           help="Minimum mapping quality allowed in the supplementary alignments. Default: 20",
                                           default=20)

            alignment_options.add_argument('-d', '--edit_distance-fraction', type=float, metavar='',
                                           help="Maximum edit distance fraction allowed in the first realignment. Default (0.05)",
                                           default=0.05)

            alignment_options.add_argument('-Q', '--split_quality', type=float, metavar='',
                                           help="Minium split score to output an interval. Default (0.0)",
                                           default=0.0)
            alignment_options.add_argument('-R', '--remap_splits', help="Remap probabilistacally bwa-mem split reads",
                                           action='store_true')

            # insert size

            i_size_estimate.add_argument('-iq', '--insert_mapq', type=int, metavar='',
                                         help="Mapq cutoff for stimating the insert size distribution. Default 60",
                                         default=60)

            i_size_estimate.add_argument('-sd', '--std', type=int, metavar='',
                                         help="Standard deviations of the insert size to extend the intervals. Default 5",
                                         default=5)

            i_size_estimate.add_argument('-s', '--sample_size', type=int, metavar='',
                                         help="Number of concordant reads (R2F1) to use for estimating the insert size distribution. Default 100000",
                                         default=100000)

            # Interval options

            interval.add_argument('-f', '--merge_fraction', type=float, metavar='',
                                  help="Merge intervals reciprocally overlapping by a fraction. Default 0.95",
                                  default=0.95)

            interval.add_argument('-P', '--interval_probability', type=float, metavar='',
                                  help="Skip edges of the graph with a probability below the threshold. Default: 0.01",
                                  default=0.01)
            interval.add_argument('-K', '--clustering_dist', type=int, metavar='',
                                  help="Cluster reads that are K nucleotides appart in the same node. Default: 500",
                                  default=500)
            interval.add_argument('-D', '--only_discordants', help="Use only discordant reads to build the graph",
                                  action='store_true')
            interval.add_argument('-F', '--allele_frequency', type=float, metavar='',
                                  help="Minimum allele frequency required to report the circle interval. Default (0.1)",
                                  default=0.1)

            # When to call a circle

            out_decision.add_argument('-S', '--split', type=int, metavar='',
                                      help="Number of required split reads to output a eccDNA. Default: 0",
                                      default=0)
            out_decision.add_argument('-O', '--number_of_discordants', type=int, metavar='',
                                      help="Number of required discordant reads for intervals with only discordants. Default: 3",
                                      default=3)

            out_decision.add_argument('-r', '--ratio', type=float, metavar='',
                                      help="Minimum in/out required coverage ratio. Default: 0.0",
                                      default=0.0)

            # coverage metrics

            coverage_metrics.add_argument('-N', '--no_coverage', help="Don't compute coverage statistics",
                                          action='store_true')

            coverage_metrics.add_argument('-b', '--bases', type=int, metavar='',
                                          help="Number of bases to extend for computing the coverage ratio. Default: 200",
                                          default=200)

            coverage_metrics.add_argument('-cq', '--cmapq', type=int, metavar='',
                                          help="Minimum mapping quality treshold for coverage computation. Default: 0",
                                          default=0)

            coverage_metrics.add_argument('-E', '--extension', type=int, metavar='',
                                          help="Number of bases inside the eccDNA breakpoint coordinates to compute the ratio. Default: 100",
                                          default=100)

            # Running options

            running.add_argument('-t', '--threads', type=int, metavar='',
                                 help="Number of threads to use.Default 1",
                                 default=1)

            running.add_argument('-dir', '--directory', metavar='',
                                 help="Working directory, default is the working directory",
                                 default=os.getcwd())

            running.add_argument('-v', '--verbose', type=int, metavar='',
                                 help='Verbose level, 1=error,2=warning, 3=message',
                                 choices=[1, 2, 3], default=3)

            # find out which arguments are missing

            parser.print_help()

            time.sleep(0.01)
            sys.stderr.write("\nInput does not match. Check that you provide the -i, -qbam and -fasta options"
                             "\nExiting\n")
            sys.exit(1)

        if len(sys.argv[2:]) == 0:
            parser.print_help()
            time.sleep(0.01)
            sys.stderr.write("\nNo arguments given to Realign. Exiting\n")
            sys.exit(1)

        return (parser)

    def args_repeats(self):

        parser = self.repeats

        parser._action_groups.pop()
        required = parser.add_argument_group('required arguments')
        optional = parser.add_argument_group('optional arguments')
        # prefixing the argument with -- means it's optional
        # input and output

        required.add_argument('-i', metavar='', help="Input: coordinate name sorted bam file")

        if "-i" in sys.argv:

            optional.add_argument('-o', '--output', metavar='',
                                  help="Ouput: Reads indicating circular DNA structural variants from repeat regions",
                                  default="circle_repeats_%s" % sys.argv[sys.argv.index("-i") + 1])

            optional.add_argument('-dir', '--directory', metavar='',
                                  help="Working directory, default is the working directory",
                                  default=os.getcwd())

            # coverage metrics
            optional.add_argument('-m', '--mismatch', metavar='',
                                  help="Number of mismatches allowed on the reads",
                                  default=2)

            optional.add_argument('-b', '--bases', type=int, metavar='',
                                  help="Number of bases to extend for computing the coverage ratio. Default: 200",
                                  default=200)

            optional.add_argument('-cq', '--cmapq', type=int, metavar='',
                                  help="Minimum mapping quality treshold for coverage computation. Default: 0",
                                  default=0)

            optional.add_argument('-E', '--extension', type=int, metavar='',
                                  help="Number of bases inside the eccDNA coordinates to compute the ratio. Default: 100",
                                  default=100)

            optional.add_argument('-r', '--ratio', type=float, metavar='',
                                  help="Minimum in/out required ratio. Default: 0.6",
                                  default=0.6)

            optional.add_argument('-f', '--fraction', type=float, metavar='',
                                  help="Required fraction to merge the intervals of the double mapped reads. Default 0.8",
                                  default=0.8)

            optional.add_argument('-n', '--read_number', metavar='',
                                  help="Minimum number of reads required to output",
                                  default=20)



        else:

            optional.add_argument('-o', '--output', metavar='',
                                  help="Ouput: Reads indicating circular DNA structural variants",
                                  )

            optional.add_argument('-dir', '--directory', metavar='',
                                  help="Working directory, default is the working directory",
                                  default=os.getcwd())

            # coverage metrics

            optional.add_argument('-m', '--mismatch', metavar='',
                                  help="Number of mismatches allowed on the reads",
                                  default=2)

            optional.add_argument('-b', '--bases', type=int, metavar='',
                                  help="Number of bases to extend for computing the coverage ratio. Default: 200",
                                  default=200)

            optional.add_argument('-cq', '--cmapq', type=int, metavar='',
                                  help="Minimum mapping quality treshold for coverage computation. Default: 0",
                                  default=0.6)

            optional.add_argument('-E', '--extension', type=int, metavar='',
                                  help="Number of bases inside the eccDNA coordinates to compute the ratio. Default: 100",
                                  default=100)

            optional.add_argument('-r', '--ratio', type=float, metavar='',
                                  help="Minimum in/out required ratio. Default: 0.6",
                                  default=0.6)

            optional.add_argument('-f', '--fraction', type=float, metavar='',
                                  help="Required fraction to merge the intervals of the double mapped reads. Default 0.8",
                                  default=0.8)

            optional.add_argument('-n', '--read_number', metavar='',
                                  help="Minimum number of reads required to output",
                                  default=20)

            parser.print_help()

            time.sleep(0.01)
            sys.stderr.write("\nNo input input given to Repeats, be sure that you are providing the flag '-i'"
                             "\nExiting\n")
            sys.exit(1)

        # parse the commands

        if len(sys.argv[2:]) == 0:
            parser.print_help()
            time.sleep(0.01)
            sys.stderr.write("\nNo arguments given to Repeats. Exiting\n")
            sys.exit(1)

        return (parser)

    def args_simulate(self):

        parser = self.simulate

        parser._action_groups.pop()
        required = parser.add_argument_group('required arguments')
        optional = parser.add_argument_group('optional arguments')
        # prefixing the argument with -- means it's optional
        # input and output

        if "-g" and "-N" in sys.argv:
            required.add_argument('-g', metavar='',
                                  help="Genome fasta file (Needs to be indexed with samtools faidx)")
            required.add_argument('-N', '--read-number', type=int, metavar='',
                                  help="Number of reads to simulate")
            optional.add_argument('-o', '--output', default='simulated.bed',
                                  help="Output file name")
            optional.add_argument('-dir', '--directory', metavar='',
                                  help="Working directory, default is the working directory",
                                  default=os.getcwd())
            optional.add_argument('-b', '--base-name', metavar='', default='simulated',
                                  help="Fastq output basename")
            optional.add_argument('-s', '--skip-region', metavar='', default=None,
                                  help="Regions of the genome to skip the simulation. The input needs to be in bed format")
            optional.add_argument('-r', '--read-length', metavar='', type=int, default=150,
                                  help="Read length to simulate")
            optional.add_argument('-m', '--mean-insert-size', metavar='', type=int, default=300,
                                  help="Mean of the insert size distribution")
            optional.add_argument('-c', '--mean-coverage', metavar='', type=int, default=30,
                                  help="Mean sequencing coverage within the eccDNA coordinates")
            optional.add_argument('-p', '--processes', metavar='', type=int, default=1,
                                  help="Mean sequencing coverage within the eccDNA coordinates")

            optional.add_argument('-v', '--variants', action='store_true',
                                  help="If set to true, introduce mutations in the reference genome prior to simulating"
                                       "reads.")
            optional.add_argument('-S', '--substitution', metavar='', type=float, default=0.0001,
                                  help="Fraction of base substitutions to introduce on the genome. Default: 0.0001")

            optional.add_argument('-I', '--Indels', metavar='', type=float, default=0.001,
                                  help="Fraction of indels to introduce on the genome. Default: 0.001")
            optional.add_argument('-J', '--java_memory', metavar='', type=str, default="-Xmx16g",
                                  help="Java memory allocation, required for mutating the genome. Default: -Xmx16g")

            optional.add_argument('-e', '--error', action='store_true',
                                  help="Introduce sequencing errors ( Uses ART on the background)")

            optional.add_argument('-i', '--instrument', metavar='', type=str, default="HS25",
                                  help="Illumina sequecing instrument to simulate reads from (Default HiSeq 2500)")

            optional.add_argument('-ir', '--insRate', metavar='', type=float, default=0.00009,
                                  help="the first-read insertion rate (default: 0.00009)")
            optional.add_argument('-ir2', '--insRate2', metavar='', type=float, default=0.00015,
                                  help="the second-read insertion rate (default: 0.00015)")
            optional.add_argument('-dr', '--delRate', metavar='', type=float, default=0.00011,
                                  help="the first-read deletion rate (default:  0.00011)")
            optional.add_argument('-dr2', '--delRate2', metavar='', type=float, default=0.00023,
                                  help="the second-read deletion rate (default: 0.00023)")
        else:
            required.add_argument('-g', metavar='',
                                  help="Genome fasta file (Needs to be indexed with samtools faidx)")
            required.add_argument('-N', '--read-number', type=int, metavar='',
                                  help="Number of reads to simulate")
            optional.add_argument('-o', '--output', default='simulated.bed',
                                  help="Output file name")
            optional.add_argument('-dir', '--directory', metavar='',
                                  help="Working directory, default is the working directory",
                                  default=os.getcwd())
            optional.add_argument('-b', '--base-name', metavar='', default='simulated',
                                  help="Fastq output basename")
            optional.add_argument('-s', '--skip-region', metavar='', default=None,
                                  help="Regions of the genome to skip the simulation. The input needs to be in bed format")
            optional.add_argument('-r', '--read-length', metavar='', type=int, default=150,
                                  help="Read length to simulate")
            optional.add_argument('-m', '--mean-insert', metavar='', type=int, default=300,
                                  help="Mean of the insert size distribution")

            optional.add_argument('-c', '--mean-coverage', metavar='', type=int, default=30,
                                  help="Mean sequencing coverage within the eccDNA coordinates")

            optional.add_argument('-p', '--processes', metavar='', type=int, default=1,
                                  help="Number of parallel processes to use")

            optional.add_argument('-v', '--variants', action='store_true',
                                  help="If set to true, introduce mutations in the reference genome prior to simulating"
                                       "reads.")
            optional.add_argument('-S', '--substitution', metavar='', type=float, default=0.0001,
                                  help="Fraction of base substitutions to introduce on the genome. Default: 0.0001")

            optional.add_argument('-I', '--Indels', metavar='', type=float, default=0.001,
                                  help="Fraction of indels to introduce on the genome. Default: 0.001")
            optional.add_argument('-J', '--java_memory', metavar='', type=str, default="-Xmx16g",
                                  help="Java memory allocation, required for mutating the genome. Default: -Xmx16g")

            optional.add_argument('-e', '--error', action='store_true',
                                  help="Introduce sequencing errors ( Uses ART on the background)")

            optional.add_argument('-i', '--instrument', metavar='', type=str, default="HS25",
                                  help="Illumina sequecing instrument to simulate reads from (Default HiSeq 2500)")
            optional.add_argument('-ir', '--insRate', metavar='', type=float, default=0.00009,
                                  help="the first-read insertion rate (default: 0.00009). Requires -e")
            optional.add_argument('-ir2', '--insRate2', metavar='', type=float, default=0.00015,
                                  help="the second-read insertion rate (default: 0.00015). Requires -e")
            optional.add_argument('-dr', '--delRate', metavar='', type=float, default=0.00011,
                                  help="the first-read deletion rate (default:  0.00011). Requires -e")
            optional.add_argument('-dr2', '--delRate2', metavar='', type=float, default=0.00023,
                                  help="the second-read deletion rate (default: 0.00023). Requires -e")

            parser.print_help()

            time.sleep(0.01)
            sys.stderr.write(
                "\nNo input input given to Simulate, be sure that you are providing the flags '-g' and '-N'"
                "\nExiting\n")

            sys.exit(1)


        if len(sys.argv[2:]) == 0:
            parser.print_help()
            time.sleep(0.01)
            sys.stderr.write("\nNo arguments given to Simulate. Exiting\n")




        return (parser)

def main():
    run = circle_map()
    pid = run.__getpid__()
    # clean
    os.system("rm -rf temp_files_%s" % pid)

if __name__ == '__main__':
    main()

