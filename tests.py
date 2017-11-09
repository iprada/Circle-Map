#!/data/xsh723/anaconda/bin/python3.6


import argparse
import sys
import os


class circle_map(object):


    def __init__(self):
        parser = argparse.ArgumentParser(
            description='Circle-Map',
            usage='''CircleMap <subprogram> [options]
            
version=1.0
contact= https://github.com/iprada/Circle-Map/issues

The Circle-Map suite

Commands:

   ReadExtractor     Extracts the reads indicating extrachromosomal circular DNA structural variants
   Realigner         Realign the soft-clipped reads and identify circular DNA
''')
        parser.add_argument('command', help='Subcommand to run')
        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print('Unrecognized command')
            parser.print_help()
            exit(1)
        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()

    def ReadExtractor(self):
        from extract_circle_SV_reads import readExtractor
        parser = argparse.ArgumentParser(
            description='Extracts the reads indicating extrachromosomal circular DNA structural variants',
            prog="CircleMap ReadExtractor",
            usage = '''CircleMap ReadExtractor [options]
            
        Author= Inigo Prada-Luengo
        version=1.0
        contact= https://github.com/iprada/Circle-Map/issues
        The Regenberg laboratory
        ''')

        parser._action_groups.pop()
        required = parser.add_argument_group('required arguments')
        optional = parser.add_argument_group('optional arguments')
        # prefixing the argument with -- means it's optional
        #input and output

        internal_parser = argparse.ArgumentParser()
        internal_parser.add_argument('-i',metavar='', help="Input: query name sorted bam file")
        internal_args = internal_parser.parse_args(sys.argv[2:])




        required.add_argument('-i',metavar='', help="Input: query name sorted bam file")

        optional.add_argument('-o',metavar='', help="Ouput: Reads indicating circular DNA structural variants",
                              default="circle_%s" % internal_args.i)

        optional.add_argument('-dir', help="Working directory, default is the working directory",
                              default=os.getcwd())

        #read extraction options
        #extract discordant reads
        optional.add_argument('-nd', help="Turn off discordant (R2F1 oriented) read extraction",metavar='',default=True)

        #soft-clipped argument
        optional.add_argument('-sc', '--nosoftclipped', help="Turn off soft-clipped read extraction",metavar='', default=True)
        #extract hard-clippped reads
        optional.add_argument('-hc', '--nohardclipped', help="Turn off hard-clipped read extraction",metavar='', default=True)

        #
        # now that we're inside a subcommand, ignore the first
        args = parser.parse_args(sys.argv[2:])


        if len(sys.argv[2:]) == 0:
            parser.print_help()
            sys.stderr.write("\nNo arguments given to read extractor. Exiting\n")
            sys.exit(1)

        if args.dir[-1:] != "/":
            args.dir = args.dir + "/"

        readExtractor_object = readExtractor(args.i,args.o,args.dir)
        readExtractor_object.extract_sv_circleReads()




    def fetch(self):
        parser = argparse.ArgumentParser(
            description='Download objects and refs from another repository')
        # NOT prefixing the argument with -- means it's not optional
        parser.add_argument('repository')
        args = parser.parse_args(sys.argv[2:])
        print('Running git fetch, repository=%s' % args.repository)




if __name__ == '__main__':
    circle_map()

