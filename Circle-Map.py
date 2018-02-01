#!/data/xsh723/anaconda/bin/python3.6
from __future__ import division
#Author Inigo Prada Luengo
#email: inigo.luengo@bio.ku.dk


import argparse
import sys
import os
import time
from extract_circle_SV_reads import readExtractor


class circle_map(object):
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

                if self.args.directory[-1:] != "/":
                    self.args.dir = self.args.directory + "/"

                object = readExtractor(self.args.i,self.args.output,self.args.dir,self.args.quality,self.args.nodiscordant,
                                       self.args.nohardclipped,self.args.nosoftclipped,self.args.verbose,self.subprogram)
                object.extract_sv_circleReads()

            elif sys.argv[1] == "Realign":
                self.subprogram = self.args_realigner()
                self.args = self.subprogram.parse_args(sys.argv[2:])

                if self.args.directory[-1:] != "/":
                    self.args.dir = self.args.directory + "/"

                object = (self.args.i,self.args.qbam,self.args.fasta,self.args.directory,self.args.q,
                          self.args.iq,self.args.sd,self.args.s,self.args.g,
                          self.args.e,self.args.n,self.args.p,self.args.t,self.args.m,
                          self.args.f)










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
            optional.add_argument('-o', '--output', metavar='',
                                  help="Ouput: Reads indicating circular DNA structural variants")

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
            sys.stderr.write("\nNo input input given to readExtractor, be sure that you are providing the flag '-i'"
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
        i_size_estimate = parser.add_argument_group('Insert size estimation options')
        output_merging = parser.add_argument_group('Output merging options')
        parallel = parser.add_argument_group('Parallele options')

        io_options.add_argument('-i', metavar='', help="Input: bam file containing the reads extracted by ReadExtractor")
        io_options.add_argument('-qbam', metavar='',help="Input: query name sorted bam file")
        io_options.add_argument('-fasta', metavar='', help="Input: file contaning the genome fasta")

        if "-i" and "qbam" and "fasta" in sys.argv:
            #alignment
            alignment_options.add_argument('-n', '--nhits', type=int, metavar='',
                                  help="Number of alignment to compute the probability of the realignment. Default: 10",
                                  default=10)

            alignment_options.add_argument('-p', '--cut_off', type=float, metavar='',
                                           help="Probability cut-off for considering a soft-clipped as realigned: Default: 0.99",
                                           default=0.99)

            alignment_options.add_argument('-m', '--min_sc', type=float, metavar='',
                                           help="Minimum soft-clipped length to attempt the realignment. Default: 8",
                                           default=8)

            alignment_options.add_argument('-g', '--gap_open', type=int, metavar='',
                                           help="Gap open penalty in the position specific scoring matrix. Default: 17",
                                           default=17)

            alignment_options.add_argument('-e', '--gap_ext', type=int, metavar='',
                                           help="Gap extension penalty in the position specific scoring matrix. Default: 1",
                                           default=1)

            alignment_options.add_argument('-q', '--mapq', type=int, metavar='',
                                           help="Minimum mapping quality allowed in the supplementary alignments. Default: 10",
                                           default=10)


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

            #overlap merging

            output_merging.add_argument('-f', '--merge_fraction', type=float, metavar='',
                                         help="Fraction to merge the SC and SA called intervals. Default 0.99",
                                         default=0.99)

            parallel.add_argument('-t', '--threads', type=int, metavar='',
                                        help="Number of threads to use.Default 1",
                                        default=1)



        else:

            alignment_options.add_argument('-n', '--nhits', type=int, metavar='',
                                           help="Number of alignment to compute the probability of the realignment. Default: 10",
                                           default=10)

            alignment_options.add_argument('-p', '--cut_off', type=float, metavar='',
                                           help="Probability cut-off for considering a soft-clipped as realigned: Default: 0.99",
                                           default=0.99)

            alignment_options.add_argument('-m', '--min_sc', type=float, metavar='',
                                           help="Minimum soft-clipped length to attempt the realignment. Default: 8",
                                           default=8)

            alignment_options.add_argument('-g', '--gap_open', type=int, metavar='',
                                           help="Gap open penalty in the position specific scoring matrix. Default: 17",
                                           default=17)

            alignment_options.add_argument('-e', '--gap_ext', type=int, metavar='',
                                           help="Gap extension penalty in the position specific scoring matrix. Default: 1",
                                           default=1)

            alignment_options.add_argument('-q', '--mapq', type=int, metavar='',
                                           help="Minimum mapping quality allowed in the supplementary alignments. Default: 10",
                                           default=10)

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

            # overlap merging

            output_merging.add_argument('-f', '--merge_fraction', type=float, metavar='',
                                        help="Fraction to merge the SC and SA called intervals. Default 0.99",
                                        default=0.99)

            parallel.add_argument('-t', '--threads', type=int, metavar='',
                                  help="Number of threads to use.Default 1",
                                  default=1)


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






if __name__ == '__main__':
    circle_map()