#!/data/xsh723/anaconda/bin/python3.6
#Author: Inigo Prada Luengo
#email: inigo.luengo@bio.ku.dk

import pysam as ps
import pybedtools as bt
import warnings
import numpy as np
from itertools import groupby
import matplotlib.pyplot as plt
#remove this afterwards



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

    """Function that gets as input the SA tag CIGAR and reports the number of bases that where matched to the genome"""

    cigar = [''.join(g) for _, g in groupby(sa_cigar, str.isalpha)]


    match_index = cigar.index('M')

    aligned = 0
    #if only one hit
    if type(match_index) == int:
        aligned += int(cigar[match_index -1])

    #when there are more than 1 hits
    else:
        assert type(match_index) == list

        for index in match_index:
            aligned += int(cigar[index - 1])

    return(aligned)


def genome_alignment_from_cigar(sa_cigar):

    """Function that gets as input the SA tag CIGAR and returns the length of the alignment interval in the genome it
    will look at the number of matches and deletions in the CIGAR, as they are the elements that will explain the genome
    alignment
    """
    aligned = 0

    cigar = [''.join(g) for _, g in groupby(sa_cigar, str.isalpha)]

    #do it for the matches
    match_index = cigar.index('M')


    # if only one hit
    if type(match_index) == int:
        aligned += int(cigar[match_index - 1])

    # when there are more than 1 hits
    else:
        assert type(match_index) == list

        for index in match_index:
            aligned += int(cigar[index - 1])


    if 'D' in cigar == True:

        deletion_index = cigar.index('D')

        # if only one hit
        if type(deletion_index) == int:
            aligned += int(cigar[deletion_index - 1])

        # when there are more than 1 hits
        else:
            assert type(deletion_index) == list

            for index in deletion_index:
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

    ####################################################################################################################
    ####################################################################################################################
    ####################################################################################################################
    ################################### CHANGE THIS PART OF THE CODE ONCE DEVELOPMENT IS FINISHED ######################
    ####################################################################################################################
    ####################################################################################################################
    ####################################################################################################################

    #peak_coverage = bam.genome_coverage(bg=True)

    #sort (sanity) and merge, BedGraph to bed
    #sorted_peak_coverage = peak_coverage.sort()
    #merged_peak_coverage = sorted_peak_coverage.merge()
    #merged_peak_coverage.saveas("peak_coverage.bed")

    merged_peak_coverage = bt.BedTool("peak_coverage.bed")

    return(merged_peak_coverage,sorted_bam)



def insert_size_dist(sample_size,mapq_cutoff,qname_bam):
    """Function that takes as input a queryname sorted bam and computes the mean insert a size and
    the standard deviation from. This number is computed from the F1R2 read with a user defined sample size,
     using a user defined mapping quality cutoff in both reads"""

    #TO-DO Warning if the bam file is not query name sorted

    whole_bam = ps.AlignmentFile(qname_bam, "rb")

    counter = 0
    insert_length = []
    read1 = ''

    # this is similar to the code of read extractor. I save the first read in memory and then I operate
    # in both reads together
    for read in whole_bam:

        if read.is_read1:
            read1 = read
        else:
            if read.is_read2 and read.qname == read1.qname:
                read2 = read
                # both reads in memory
                if read1.mapq >= mapq_cutoff and read2.mapq >= mapq_cutoff:

                    if read1.is_reverse == False and read2.is_reverse == True:

                        if (read2.tlen + read1.infer_query_length()) < 0 and (read1.tlen - read2.infer_query_length()) > 0:



                            insert_length.append(read1.tlen)
                            counter += 1

        if counter >= sample_size:
            break
        else:
            continue

    mean = np.mean(insert_length)
    std = np.std(insert_length)

    return(mean, std)

def get_realignment_interval(grouped,grouped_pd,interval_extension):
    """Function that takes as input the insert metricsa grouped realignment interval and a pandas grouped one and will
    return the interval to perform the probabilistic realignment"""

    #interval definition


    read_types = grouped_pd.read_type.unique()



    if np.any(read_types == 'SC') == False:
        #print(grouped)
        grouped = grouped.sort()

        #complete realignment interval
        grouped = grouped.merge()

    elif (np.any(read_types == 'SC') == True) and (np.any(read_types == 'DR')== True or np.any(read_types == 'SA')== True):

        #remove the 'SC'
        no_sc_grouped = []
        for interval in grouped:
            if interval[3] != 'SC':
                no_sc_grouped.append(interval)

        grouped = bt.BedTool(no_sc_grouped)
        grouped = grouped.sort()
        grouped = grouped.merge()


    #orientation extension


    #only one extension

    extension_orientation = grouped_pd.orientation.unique()


    extended_grouped_L = []


    # check if the interval should be left extended
    if np.any(extension_orientation == 'L') == True:

        for interval in grouped:
            start = interval.start - interval_extension

            # in case that start is smaller than chromosome length
            if start < 0:
                extended_grouped_L.append([interval.chrom,str(0),interval.end])

            else:
                extended_grouped_L.append([interval.chrom, int(round(start)), interval.end])

        grouped = bt.BedTool(extended_grouped_L)



    extended_grouped = []

    #Check if the interval should be left extended
    if np.any(extension_orientation == 'R') == True:

        for interval in grouped:
            end = interval.end + interval_extension

            # in case that start is smaller than chromosome length
            #check chromosome length
            if end < 250000000:
                extended_grouped.append([interval.chrom,interval.start,int(round(end))])

            else:
                extended_grouped.append([interval.chrom, interval.start,int(round(end))])

        grouped = bt.BedTool(extended_grouped_L)

    return(grouped)













