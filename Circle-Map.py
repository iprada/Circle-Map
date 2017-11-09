#!/data/xsh723/anaconda/bin/python3.6


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

                object = readExtractor(self.args.i,self.args.output,self.args.dir,self.args.nodiscordant,
                                       self.args.nohardclipped,self.args.nosoftclipped,self.args.verbose,self.subprogram)
                object.extract_sv_circleReads()

            elif sys.argv[1] == "realign":
                a = 0








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

    def fetch(self):
        parser = argparse.ArgumentParser(
            description='Download objects and refs from another repository')
        # NOT prefixing the argument with -- means it's not optional
        parser.add_argument('repository')
        args = parser.parse_args(sys.argv[2:])
        print('Running git fetch, repository=%s' % args.repository)


if __name__ == '__main__':
    circle_map()