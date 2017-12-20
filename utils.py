#!/data/xsh723/anaconda/bin/python3.6
#Author: Inigo Prada Luengo
#email: inigo.luengo@bio.ku.dk

import pysam as ps
import pybedtools as bt
import warnings
from itertools import groupby
import numpy as np

def is_soft_clipped(read):
    """Function that checks the CIGAR string of the sam file and returns true if the read is soft-clipped"""

    # cigar 4 equals to S in pysam sam representation
    match = 0
    for cigar in read.cigar:
        if cigar[0] == 4:
            match +=1
        else:
            continue

    if match > 0:
        return(True)

    else:
        return(False)

def is_hard_clipped(read):
    """Function that checks the CIGAR string of the sam file and returns true if the read is hard-clipped"""

    # cigar 5 equals to H in pysam sam representation
    match = 0
    for cigar in read.cigar:
        if cigar[0] == 5:
            match += 1
        else:
            continue

    if match > 0:
        return (True)

    else:
        return (False)

def aligned_bases(read):
    """Function that counts the number of aligned bases from the CIGAR string and returns and integer"""

    aligned = 0

    for cigar in read.cigar:
        if cigar[0] == 0:
            aligned += cigar[1]
        else:
            continue

    return(aligned)

def aligned_bases_from_sa(sa_cigar):
    """Function that gets as input the SA tag and reports the number of bases that where matched to the genome"""

    cigar = [''.join(g) for _, g in groupby(sa_cigar, str.isalpha)]


    match_index = cigar.index('M')

    aligned = 0
    if type(match_index) == int:
        aligned += int(cigar[match_index -1])

    else:
        assert type(match_index) == list

        for index in match_index:
            aligned += int(cigar[index - 1])

    return(aligned)





def bam_circ_sv_peaks(bam,input_bam_name,cores):
    """Function that takes as input a bam file and returns a merged bed file of the genome covered by the bam, it will create
    and index too"""

    # check bam header for sorting state



    #check the header of the bam file for the sorting state, and sort if necessary

    if 'HD' in bam.header:
        if bam.header['HD']['SO'] == 'queryname':
            print("Bam is sorted by queryname, sorting bam by position")

            bam.close()
            ps.sort("-@","%s" % cores,"-o", "coordinate_%s" % input_bam_name, "%s" % input_bam_name)

            sorted_bam = ps.AlignmentFile("coordinate_%s" %input_bam_name)

            ps.index("coordinate_%s" % input_bam_name)

            bam = bt.BedTool("coordinate_%s" % input_bam_name)



        elif bam.header['HD']['SO'] == 'unsorted':

            bam.close()

            print("Bam is unsorted, sorting bam by position")

            ps.sort("-@","%s" % cores,"-o", "coordinate_%s" % input_bam_name, "%s" % input_bam_name)


            sorted_bam = ps.AlignmentFile("coordinate_%s" % input_bam_name)

            ps.index("coordinate_%s" % input_bam_name)

            bam = bt.BedTool("coordinate_%s" % input_bam_name)

        elif bam.header['HD']['SO'] == 'coordinate':

            bam.close()
            sorted_bam = ps.AlignmentFile("%s" % input_bam_name)

            ps.index("%s" % input_bam_name)

            bam = bt.BedTool("%s" % input_bam_name)

        else:

            warnings.warn(
                "WARNING: the bam file does not have an SO tag.\nCircle-Map cannot check if the bam file is sorted by coordinate.\n If the bam file is not sorted by coordinate the program will file")
            print(
                "As sanity check, sort your bam file coordinate with the following command:\n\n\tsamtools sort -o output.bam input.bam")


    else:
        warnings.warn(
            "WARNING: the bam file does not have an HD tag.\nCircle-Map cannot check if the bam file is sorted by coordinate.\n If the bam file is not sorted by coordinate the program will file")
        print(
            "As sanity check, sort your bam file coordinate with the following command:\n\n\tsamtools sort -o output.bam input.bam")


    #from bam to BedGraph

    peak_coverage = bam.genome_coverage(bg=True)

    #sort (sanity) and merge, BedGraph to bed
    sorted_peak_coverage = peak_coverage.sort()
    merged_peak_coverage = sorted_peak_coverage.merge()

    return(merged_peak_coverage,sorted_bam)













