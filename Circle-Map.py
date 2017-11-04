#!/data/xsh723/anaconda/bin/python3.6

#Author Inigo Prada Luengo
#email: inigo.luengo@bio.ku.dk


__author__ = [
    'Inigo Prada-Luengo (inigo.luengo@bio.ku.dk)',
]

__uni__ = [
    'University of Copenhagen (Department of Biology)',
]

import argparse
import os

def main():
    parser = argparse.ArgumentParser()
    #check this




def readextractor():
    from extract_circle_SV_reads import readExtractor
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help="Query name sorted bam file")
    parser.add_argument('-o', '--output', help="Output bam file with the reads that indicate circular DNA",default=parser.output_name)
    parser.add_argument('-dir', '--directory', help="Working directory, default is the working directory",default=os.getcwd())
    args = parser.parse_args()
    object = readExtractor(args.input,args.output,args.directory)
    object.extract_sv_circleReads()


def map_circles():
    a = 0

def annotate():
    a = 0

# this is done to know if this file is been run directly by python or imported
if __name__ == '__main__':
    main()