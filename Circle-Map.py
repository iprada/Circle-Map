#!/home/iprada/bin/miniconda3/bin/python3.6
from __future__ import division
from utils import start_realign
#Author Inigo Prada Luengo
#email: inigo.luengo@bio.ku.dk


import argparse
import sys
import os
import time
import pandas as pd
from extract_circle_SV_reads import readExtractor
from realigner import realignment
from repeats import repeat
from utils import merge_final_output,filter_by_ratio
from coverage import coverage
import multiprocessing as mp
import pybedtools as bt






class circle_map:

    def __init__(self):
        self.parser = argparse.ArgumentParser(
            description='Circle-Map',
            usage='''CircleMap <subprogram> [options]

Author= Inigo Prada-Luengo
version=1.0
contact= https://github.com/iprada/Circle-Map/issues
The Regenberg laboratory

The Circle-Map suite

Commands:

   ReadExtractor     Extracts the reads indicating extrachromosomal circular DNA structural variants
   Realign         Realign the soft-clipped reads and identify circular DNA
   Repeats         Identify eccDNA coming from repeat regions
''')
        subparsers = self.parser.add_subparsers()

        self.readextractor = subparsers.add_parser(
            name="ReadExtractor",
            description='Extracts the reads indicating extrachromosomal circular DNA structural variants',
            prog="CircleMap ReadExtractor",
            usage='''CircleMap ReadExtractor [options]

                     Author= Inigo Prada-Luengo
                     version=1.0
                     contact= https://github.com/iprada/Circle-Map/issues
                     The Regenberg laboratory
                     '''

                                            )

        self.realigner = subparsers.add_parser(
            name="Realign",
            description='Realigns the soft-clipped reads and indentifies circular DNA',
            prog="CircleMap Realign",
            usage='''CircleMap Realign [options]

                             Author= Inigo Prada-Luengo
                             version=1.0
                             contact= https://github.com/iprada/Circle-Map/issues
                             The Regenberg laboratory
                             '''

        )

        self.repeats = subparsers.add_parser(
            name="Repeats",
            description='Identify cluster of reads with two alignments',
            prog="CircleMap Reepeats",
            usage='''CircleMap Repeats [options]

                             Author= Inigo Prada-Luengo
                             version=1.0
                             contact= https://github.com/iprada/Circle-Map/issues
                             The Regenberg laboratory
                             '''

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



                object = readExtractor(self.args.i,self.args.output,self.args.directory,self.args.quality,self.args.nodiscordant,
                                       self.args.nohardclipped,self.args.nosoftclipped,self.args.verbose,self.subprogram)
                object.extract_sv_circleReads()

            elif sys.argv[1] == "Realign":
                self.subprogram = self.args_realigner()
                self.args = self.subprogram.parse_args(sys.argv[2:])



                splitted,sorted_bam,begin = start_realign(self.args.i,self.args.output,self.args.threads,self.args.verbose)





                if __name__ == '__main__':

                    lock = mp.Lock()

                    processes = []

                    for core in range(0,self.args.threads):

                        object = realignment(sorted_bam, self.args.qbam,self.args.sbam, self.args.fasta, self.args.directory,
                                             self.args.mapq,
                                             self.args.insert_mapq, self.args.std, self.args.sample_size,
                                             self.args.gap_open,
                                             self.args.gap_ext, self.args.nhits, self.args.cut_off, self.args.min_sc,
                                             self.args.merge_fraction, self.args.interval_probability, self.args.output,
                                             self.args.threads, splitted[core],lock,self.args.split,self.args.ratio,self.args.verbose)

                        processes.append(object)

                    jobs = []
                    # init the processes

                    processes[0].print_parameters()
                    for i in range(0,len(processes)):

                        p = mp.Process(target=processes[i].realign)
                        jobs.append(p)
                        p.start()
                    # kill the process
                    for p in jobs:
                        p.join()

                    output = merge_final_output("%s" % self.args.output, begin,self.args.split,self.args.directory,self.args.merge_fraction)


                    # compute coverage statistics

                    if self.args.no_coverage == False:

                        coverage_object = coverage(self.args.sbam,output,self.args.bases,self.args.cmapq,self.args.extension,self.args.directory)
                        coverage_dict, header_dict = coverage_object.get_wg_coverage()
                        output = coverage_object.compute_coverage(coverage_dict, header_dict)
                        filtered_output = filter_by_ratio(output,self.args.ratio)
                        filtered_output.saveas("%s" % self.args.output)

                    else:
                        output.saveas("%s" % self.args.output)

            elif sys.argv[1] == "Repeats":



                self.subprogram = self.args_repeats()
                self.args = self.subprogram.parse_args(sys.argv[2:])

                if self.args.directory[-1:] != "/":
                    self.args.dir = self.args.directory + "/"

                object = repeat(self.args.i,self.args.directory,self.args.mismatch,self.args.fraction,self.args.read_number)
                bed = object.find_circles()
                coverage_object = coverage(self.args.i, bed, self.args.bases, self.args.cmapq,
                                           self.args.extension, self.args.directory)

                coverage_dict, header_dict = coverage_object.get_wg_coverage()
                output = coverage_object.compute_coverage(coverage_dict, header_dict)
                filtered_output = filter_by_ratio(output, self.args.ratio)
                filtered_output.saveas("%s" % self.args.output)







    def args_readextractor(self):

        parser = self.readextractor


        parser._action_groups.pop()
        required = parser.add_argument_group('required arguments')
        optional = parser.add_argument_group('optional arguments')
        # prefixing the argument with -- means it's optional
        # input and output


        required.add_argument('-i', metavar='', help="Input: query name sorted bam file")
        required.add_argument('-o', '--output', metavar='',
                              help="Ouput: Reads indicating circular DNA structural variants",
                              default="circle_%s" % sys.argv[sys.argv.index("-i") + 1])


        if "-i" and "-o" in sys.argv:


            optional.add_argument('-dir', '--directory',metavar='',
                                  help="Working directory, default is the working directory",
                                  default=os.getcwd())

            #mapping quality cutoff

            optional.add_argument('-q', '--quality',type=int, metavar='',
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

            #verbose level


            optional.add_argument('-v', '--verbose',type=int, metavar='', help='Verbose level, 1=error,2=warning, 3=message',
                                  choices=[1, 2, 3],default=3)

        else:

            optional.add_argument('-dir', '--directory', metavar='',help="Working directory, default is the working directory",
                                  default=os.getcwd())


            #mapping quality cutoff
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


            #verbose level

            optional.add_argument('-v', '--verbose',type=int, metavar='', help='Verbose level, 1=error,2=warning, 3=message. Default=3',
                                  choices=[1, 2, 3],default=3)


            parser.print_help()

            time.sleep(0.01)
            sys.stderr.write("\nNo input or output input given to readExtractor, be sure that you are providing the flags'-i' and '-o'"
                             "\nExiting\n")
            sys.exit(1)



        #parse the commands



        if len(sys.argv[2:]) == 0:
            parser.print_help()
            time.sleep(0.01)
            sys.stderr.write("\nNo arguments given to read extractor. Exiting\n")
            sys.exit(1)

        return(parser)

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


        io_options.add_argument('-i', metavar='', help="Input: bam file containing the reads extracted by ReadExtractor")
        io_options.add_argument('-qbam', metavar='',help="Input: query name sorted bam file")
        io_options.add_argument('-sbam', metavar='', help="Input: coordinate sorted bam file")
        io_options.add_argument('-fasta', metavar='', help="Input: file contaning the genome fasta")




        if "-i" and "-qbam" and "-fasta" in sys.argv:

            #output

            io_options.add_argument('-c', metavar='', help="Genome Coverage file",default=None)

            io_options.add_argument('-o','--output', metavar='', help="Output filename",
                                    default="circle_%s.bed" % sys.argv[sys.argv.index("-i") + 1])

            #alignment
            alignment_options.add_argument('-n', '--nhits', type=int, metavar='',
                                  help="Number of alignment to compute the probability of the realignment. Default: 10",
                                  default=10)

            alignment_options.add_argument('-p', '--cut_off', type=float, metavar='',
                                           help="Probability cut-off for considering a soft-clipped as realigned: Default: 0.99",
                                           default=0.99)

            alignment_options.add_argument('-m', '--min_sc', type=float, metavar='',
                                           help="Minimum soft-clipped length to attempt the realignment. Default: 14",
                                           default=14)

            alignment_options.add_argument('-g', '--gap_open', type=int, metavar='',
                                           help="Gap open penalty in the position specific scoring matrix. Default: 17",
                                           default=17)

            alignment_options.add_argument('-e', '--gap_ext', type=int, metavar='',
                                           help="Gap extension penalty in the position specific scoring matrix. Default: 5",
                                           default=5)

            alignment_options.add_argument('-q', '--mapq', type=int, metavar='',
                                           help="Minimum mapping quality allowed in the supplementary alignments. Default: 20",
                                           default=20)


            #insert size

            i_size_estimate.add_argument('-iq', '--insert_mapq', type=int, metavar='',
                                           help="Mapq cutoff for stimating the insert size distribution. Default 60",
                                           default=60)

            i_size_estimate.add_argument('-sd', '--std', type=int, metavar='',
                                         help="Number of standard deviations of the insert size to extend the intervals. Default 4",
                                         default=4)


            i_size_estimate.add_argument('-s', '--sample_size', type=int, metavar='',
                                         help="Number of concordant reads (R2F1) to use for estimating the insert size distribution. Default 100000",
                                         default=100000)

            #Interval options

            interval.add_argument('-f', '--merge_fraction', type=float, metavar='',
                                         help="Fraction to merge the SC and SA called intervals. Default 0.95",
                                         default=0.95)

            interval.add_argument('-P', '--interval_probability', type=float, metavar='',
                                  help="Skip intervals where the probability of been the mate is less than the cutoff. Default: 0.01",
                                  default=0.01)

            #When to call a circle

            out_decision.add_argument('-S', '--split', type=int, metavar='',
                                         help="Number of required split reads to output a eccDNA. Default: 2",
                                         default=2)

            out_decision.add_argument('-r', '--ratio', type=float, metavar='',
                                      help="Minimum in/out required ratio. Default: 0.0",
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
                                          help="Number of bases inside the eccDNA coordinates to compute the ratio. Default: 100",
                                          default=100)


            #run options

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

            #output

            io_options.add_argument('-c', metavar='', help="Genome Coverage file", default=None)

            io_options.add_argument('-o', metavar='', help="Output filename")


            alignment_options.add_argument('-n', '--nhits', type=int, metavar='',
                                           help="Number of alignment to compute the probability of the realignment. Default: 10",
                                           default=10)

            alignment_options.add_argument('-p', '--cut_off', type=float, metavar='',
                                           help="Probability cut-off for considering a soft-clipped as realigned: Default: 0.99",
                                           default=0.99)

            alignment_options.add_argument('-m', '--min_sc', type=float, metavar='',
                                           help="Minimum soft-clipped length to attempt the realignment. Default: 14",
                                           default=14)

            alignment_options.add_argument('-g', '--gap_open', type=int, metavar='',
                                           help="Gap open penalty in the position specific scoring matrix. Default: 17",
                                           default=17)

            alignment_options.add_argument('-e', '--gap_ext', type=int, metavar='',
                                           help="Gap extension penalty in the position specific scoring matrix. Default: 5",
                                           default=5)

            alignment_options.add_argument('-q', '--mapq', type=int, metavar='',
                                           help="Minimum mapping quality allowed in the supplementary alignments. Default: 20",
                                           default=20)

            # insert size

            i_size_estimate.add_argument('-iq', '--insert_mapq', type=int, metavar='',
                                         help="Mapq cutoff for stimating the insert size distribution. Default 60",
                                         default=60)

            i_size_estimate.add_argument('-sd', '--std', type=int, metavar='',
                                         help="Number of standard deviations of the insert size to extend the intervals. Default 4",
                                         default=4)

            i_size_estimate.add_argument('-s', '--sample_size', type=int, metavar='',
                                         help="Number of concordant reads (R2F1) to use for estimating the insert size distribution. Default 100000",
                                         default=100000)

            # Interval options


            interval.add_argument('-f', '--merge_fraction', type=float, metavar='',
                                        help="Fraction to merge the SC and SA called intervals. Default 0.95",
                                        default=0.95)

            interval.add_argument('-P', '--interval_probability', type=float, metavar='',
                                  help="Skip intervals where the probability of been the mate is less than the cutoff. Default: 0.01",
                                  default=0.99)

            # When to call a circle

            out_decision.add_argument('-S', '--split', type=int, metavar='',
                                      help="Number of required split reads to output a eccDNA. Default: 2",
                                      default=2)

            out_decision.add_argument('-r', '--ratio', type=float, metavar='',
                                      help="Minimum in/out required ratio. Default: 0.0",
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
                                          help="Number of bases inside the eccDNA coordinates to compute the ratio. Default: 100",
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


            #find out which arguments are missing

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


        return(parser)


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

            optional.add_argument('-dir', '--directory',metavar='',
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



        #parse the commands



        if len(sys.argv[2:]) == 0:
            parser.print_help()
            time.sleep(0.01)
            sys.stderr.write("\nNo arguments given to Repeats. Exiting\n")
            sys.exit(1)

        return(parser)






if __name__ == '__main__':
    circle_map()
    #clean
    os.system("rm -rf temp_files")
