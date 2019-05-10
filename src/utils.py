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

import pysam as ps
import pybedtools as bt
import warnings
import numpy as np
import pandas as pd
import itertools as it
import edlib
import os
import subprocess as sp
import glob
import time
import sys
from scipy import stats as st
import random
import re
from numba import jit
import math





def is_soft_clipped(read):

    """Function that checks the CIGAR string of the sam file and returns true if the read is soft-clipped"""

    # cigar 4 equals to S in pysam sam representation
    match = 0
    for cigar in read.cigar:
        if cigar[0] == 4:
            match +=1
        else:
            pass

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
            pass

    if match > 0:
        return (True)

    else:
        return (False)

def rightmost_from_read(read):
    """Function that takes as input a read a returns its rightmost mapping position"""

    rightmost = 0

    #matches, deletions and ref skip consume reference
    for cigar in read.cigar:

        if cigar[0] == 0:
            rightmost += cigar[1]

        elif cigar[0] == 2:
            rightmost += cigar[1]

        elif cigar[0] == 3:
            rightmost += cigar[1]


    return(read.reference_start + rightmost)

def rightmost_from_sa(leftmost,sa_cigar):
    """Function that takes as input the leftmost position of a supplementary alignment and returns it rightmost mapping
    position"""


    #the SA alignment is 1 based
    rightmost = int(leftmost)-1

    cigar = [''.join(g) for _, g in it.groupby(sa_cigar, str.isalpha)]
    # matches, deletions and ref skip consume reference
    match_index = [x for x in range(len(cigar)) if cigar[x] == 'M']
    deletion_index = [x for x in range(len(cigar)) if cigar[x] == 'D']
    ambiguous_index = [x for x in range(len(cigar)) if cigar[x] == 'N']


    for index in match_index:
        rightmost += int(cigar[index-1])

    for index in deletion_index:
        rightmost += int(cigar[index-1])


    for index in ambiguous_index:
        rightmost += int(cigar[index-1])

    assert rightmost >= (int(leftmost)-1)

    return(rightmost)





def aligned_bases(read):

    """Function that counts the number of aligned bases from the CIGAR string and returns and integer"""

    aligned = 0

    for cigar in read.cigar:
        if cigar[0] == 0:
            aligned += cigar[1]
        else:
            pass
    assert aligned >= 0
    return(aligned)

def aligned_bases_from_sa(sa_cigar):

    """Function that gets as input the SA tag CIGAR and reports the number of bases that where matched to the genome"""

    cigar = [''.join(g) for _, g in it.groupby(sa_cigar, str.isalpha)]


    match_index =  [x for x in range(len(cigar)) if cigar[x]=='M']

    aligned = 0
    #if only one hit
    if type(match_index) == int:
        aligned += int(cigar[match_index -1])

    #when there are more than 1 hits
    else:
        assert type(match_index) == list

        for index in match_index:
            aligned += int(cigar[index - 1])

    assert aligned >=0
    return(aligned)


def genome_alignment_from_cigar(sa_cigar):

    """Function that gets as input the SA tag CIGAR and returns the length of the alignment interval in the genome it
    will look at the number of matches and deletions in the CIGAR, as they are the elements that will explain the genome
    alignment
    """
    aligned = 0

    cigar = [''.join(g) for _, g in it.groupby(sa_cigar, str.isalpha)]

    #do it for the matches
    match_index = [x for x in range(len(cigar)) if cigar[x]=='M']


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

    assert aligned >=0
    return(aligned)





def bam_circ_sv_peaks(bam,input_bam_name,cores,verbose,pid,clusters):
    """Function that takes as input a bam file and returns a merged bed file of the genome covered by the bam, it will create
    and index too"""

    # check bam header for sorting state



    #check the header of the bam file for the sorting state, and sort if necessary

    if 'HD' in bam.header:
        if bam.header['HD']['SO'] == 'queryname':
            print("Bam is sorted by queryname, sorting bam by position")

            bam.close()
            ps.sort("-@","%s" % cores,"-o", "coordinate_%s" % input_bam_name, "%s" % input_bam_name)

            sorted_bam = ps.AlignmentFile("coordinate_%s" % input_bam_name)

            ps.index("coordinate_%s" % input_bam_name)





        elif bam.header['HD']['SO'] == 'unsorted':

            bam.close()

            print("Bam is unsorted, sorting bam by position")

            ps.sort("-@","%s" % cores,"-o", "coordinate_%s" % input_bam_name, "%s" % input_bam_name)


            sorted_bam = ps.AlignmentFile("coordinate_%s" % input_bam_name)

            ps.index("coordinate_%s" % input_bam_name)



        elif bam.header['HD']['SO'] == 'coordinate':

            bam.close()
            sorted_bam = ps.AlignmentFile("%s" % input_bam_name)

            ps.index("%s" % input_bam_name)



        else:
            if verbose < 2:

                warnings.warn(
                    "WARNING: the bam file does not have an SO tag.\nCircle-Map cannot check if the bam file is sorted by coordinate.\n If the bam file is not sorted by coordinate the program will file")
                print(
                    "As sanity check, sort your bam file coordinate with the following command:\n\n\tsamtools sort -o output.bam input.bam")


    else:

        if verbose < 2:
            warnings.warn(
                "WARNING: the bam file does not have an HD tag.\nCircle-Map cannot check if the bam file is sorted by coordinate.\n If the bam file is not sorted by coordinate the program will file")
            print(
                "As sanity check, sort your bam file coordinate with the following command:\n\n\tsamtools sort -o output.bam input.bam")


    #from bam to BedGraph

    sp.call("bedtools genomecov -bg -ibam %s | sort -T temp_files_%s -k 1,1 -k2,2n | mergeBed -d %s -c 4 -o mean | sort -r -n -k 4,4 > temp_files_%s/peaks.bed" %
            (input_bam_name,pid,clusters,pid),shell=True)

    split_peaks = []
    for interval in bt.BedTool("temp_files_%s/peaks.bed" % pid):

        if int(interval[2])-int(interval[1])>500:
            w_start = int(interval[1])
            while w_start < int(interval[2]):
                splitted = [interval.chrom,str(w_start),str(w_start+300)]
                w_start+=300
                split_peaks.append(splitted)
        else:
            split_peaks.append([interval.chrom,str(interval.start),str(interval.end)])


    return(sorted_bam,split_peaks)


def get_mate_intervals(sorted_bam,interval,mapq_cutoff,verbose,only_discordants):

    """Function that takes as input a sorted bam, an interval and the mapq cutoff and returns the mate alignment positions
        (the realignment prior) intervals"""

    try:



        candidate_mates = []
        for read in sorted_bam.fetch(interval['chrom'], int(interval['start']), int(interval['end']),multiple_iterators=True):

            if read.mapq >= mapq_cutoff:

                # create mate interval based on the soft-clipped SA alignments
                if is_soft_clipped(read) == True and read.has_tag('SA'):
                    if only_discordants != True:

                        read_chr = sorted_bam.get_reference_name(read.reference_id)
                        suplementary = read.get_tag('SA')

                        # [chr, left_most start, "strand,CIGAR,mapq, edit_distance]
                        supl_info = [x.strip() for x in suplementary.split(',')]

                        if read_chr == supl_info[0] and int(supl_info[4]) >= mapq_cutoff:

                            # split read with the same orientation
                            if (read.is_reverse == True and supl_info[2] == '-') or (
                                    read.is_reverse == False and supl_info[2] == '+'):

                                # SA is downstream, the interval is start, start+read length

                                if read.reference_start > int(supl_info[1]):

                                    ref_alignment_length = genome_alignment_from_cigar(supl_info[3])

                                    # ref_alignment_length * 2 is done for extending the realignment region
                                    # "SA" means that the realignment prior has been generated by a supplementary alignment
                                    # L means that the SA is aligned to to a rightmost part.

                                    mate_interval = [interval['chrom'], int(supl_info[1]) - (ref_alignment_length),
                                                     (int(supl_info[1]) + (ref_alignment_length)), "SA", "L",str(
                                            1-phred_to_prob(np.array(int(supl_info[4]),dtype=np.float64)))]

                                    candidate_mates.append(mate_interval)


                                # SA is upstream, the interval is end - read length, end
                                elif read.reference_start < int(supl_info[1]):

                                    ref_alignment_length = genome_alignment_from_cigar(supl_info[3])

                                    # ref_alignment_length * 2 is done for extending the realignment region, "SA" means that the realignment prior has been generated
                                    # by a supplementary alignment. R means that the SA is aligned to to a rightmost part.

                                    mate_interval = [interval['chrom'], (int(supl_info[1]) - (ref_alignment_length)),
                                                     int(supl_info[1]) + (ref_alignment_length), "SA", "R",str(1-phred_to_prob(np.array(int(supl_info[4]),dtype=np.float64)))]

                                    candidate_mates.append(mate_interval)
                    else:
                        pass




                # check discordant reads (R2F1 orientation)
                elif read.is_unmapped == False and read.mate_is_unmapped == False:

                    # check R2F1 orientation,when the R2 read
                    if read.is_reverse == True and read.mate_is_reverse == False:

                        # R2F1 order
                        if read.reference_start < read.next_reference_start:

                            if read.reference_id == read.next_reference_id:
                                # create mate interval
                                read_length = read.infer_query_length()

                                # DR means that the realignment prior has been generated by the discordants. R means
                                # that the mate has been aligned to a rightmost part





                                mate_interval = [interval['chrom'], read.next_reference_start,
                                                 (read.next_reference_start + read_length), "DR",
                                                 "R",str(1-phred_to_prob(np.array(read.get_tag('MQ'),dtype=np.float64)))]
                                candidate_mates.append(mate_interval)


                    # R2F1 when iterating trough F1 read
                    elif read.is_reverse == False and read.mate_is_reverse == True:

                        if read.next_reference_start < read.reference_start:

                            if read.reference_id == read.next_reference_id:
                                # create mate interval
                                read_length = read.infer_query_length()

                                # L means that the mate is aligned to a leftmost part

                                mate_interval = [interval['chrom'], read.next_reference_start,
                                                 (read.next_reference_start + read_length),
                                                 "DR", "L",str(1-phred_to_prob(np.array(read.get_tag('MQ'),dtype=np.float64)))]
                                candidate_mates.append(mate_interval)
                    else:

                        if only_discordants != True:
                            # soft clipped without and SA and hard clipped reads (secondary)


                            if is_soft_clipped(read) == True and read.has_tag('SA') == False:
                                # mate interval is whole chromosome

                                if 'SQ' in sorted_bam.header:

                                    for reference in sorted_bam.header['SQ']:

                                        if reference['SN'] == sorted_bam.get_reference_name(read.reference_id):
                                            # LR is added just not to crash the program

                                            mate_interval = [interval['chrom'], 1, reference['LN'], "SC", "LR",0]

                                            candidate_mates.append(mate_interval)


                                else:

                                    if verbose < 2:

                                        warnings.warn(
                                            "WARNING: the bam file does not have a SQ tag. Circle-Map cannot check the reference length for realigning\n"
                                            "soft clipped reads without a SA tag, hence, skipping. Please, check if your bam file is truncated")

                            elif is_hard_clipped(read):

                                # all hard clipped reads have SA tag with bwa, but just as sanity

                                if read.has_tag('SA'):

                                    read_chr = sorted_bam.get_reference_name(read.reference_id)

                                    suplementary = read.get_tag('SA')

                                    # [chr, left_most start, "strand,CIGAR,mapq, edit_distance]
                                    supl_info = [x.strip() for x in suplementary.split(',')]

                                    if read_chr == supl_info[0] and int(supl_info[4]) >= mapq_cutoff:

                                        # SA alignment with the same orientation
                                        if (read.is_reverse == True and supl_info[2] == '-') or (
                                                read.is_reverse == False and supl_info[2] == '+'):

                                            # SA is downstream, the interval is start, start+read length

                                            if read.reference_start > int(supl_info[1]):

                                                ref_alignment_length = genome_alignment_from_cigar(supl_info[3])

                                                # ref_alignment_length * 2 is done for extending the realignment region
                                                # "SA" means that the realignment prior has been generated by a supplementary alignment
                                                # L means that the SA is in a downstream region

                                                mate_interval = [interval['chrom'], int(supl_info[1]) - (ref_alignment_length * 2),
                                                                 (int(supl_info[1]) + (ref_alignment_length * 2)), "SA", "L",str(1-phred_to_prob(int(supl_info[4])))]

                                                candidate_mates.append(mate_interval)


                                            # SA is upstream, the interval is end - read length, end
                                            elif read.reference_start < int(supl_info[1]):

                                                ref_alignment_length = genome_alignment_from_cigar(supl_info[3])

                                                # ref_alignment_length * 2 is done for extending the realignment region, "SA" means that the realignment prior has been generated
                                                # by a supplementary alignment
                                                # R means that the SA is in a upstream region

                                                mate_interval = [interval['chrom'],
                                                                 (int(supl_info[1]) - (ref_alignment_length * 2)),
                                                                 int(supl_info[1]) + (ref_alignment_length * 2), "SA", "R",str(1-phred_to_prob(int(supl_info[4])))]



                                                candidate_mates.append(mate_interval)
                        else:
                            pass



            else:
                # low mapping quality reads, do nothing
                pass

        #this function should return the candidate mates (realignment prior, discordant intervals/split read intervals and soft-clipped reads)
        return(candidate_mates)

    except BaseException as e:

        if verbose < 2:

            warnings.warn(
                "WARNING: Could not get mate interval priors for the interval %s due to the following error %s \n Skipping interval" % (str(interval),str(e)))











def insert_size_dist(sample_size,mapq_cutoff,qname_bam):
    """Function that takes as input a queryname sorted bam and computes the mean insert a size and
    the standard deviation from. This number is computed from the F1R2 read with a user defined sample size,
     using a user defined mapping quality cutoff in both reads"""


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
                    if read1.is_proper_pair:
                        if is_hard_clipped(read1) == False and is_hard_clipped(read2) == False:
                            if is_soft_clipped(read1) == False and is_soft_clipped(read2) == False:
                                if read1.is_reverse == False and read2.is_reverse == True:
                                    if read1.tlen > 0:
                                        insert_length.append(read1.tlen)
                                        counter += 1

        if counter >= sample_size:
            break
        else:
            pass


    mean = np.mean(insert_length)
    std = np.std(insert_length)
    return(mean, std)

def get_realignment_intervals(bed_prior,interval_extension,interval_p_cutoff,verbose):


    """Function that takes as input a bed file with the read type information and will remove the soft-clipped if there
            are more informative priors (DR,SA). If there are only soft-clipped reads, they will be saved to a bed file to attemp
            lonely soft-clipped read rescue"""

    try:

        labels = ['chrom', 'start', 'end', 'read_type', 'orientation','probability']
        candidate_mates_dataframe = pd.DataFrame.from_records(bed_prior, columns=labels)

        read_types = candidate_mates_dataframe.read_type.unique()
        orientation = candidate_mates_dataframe.orientation.unique()





        #this contains the sumatory over all probabilities
        sum = 0


        if np.any(read_types == 'SC') == False:


            # nothing. Sort and merge


            candidate_mates_dataframe = candidate_mates_dataframe.sort_values(by=['chrom', 'start','end'],ascending=[True,True,True])
            candidate_mates_dataframe['probability'] = candidate_mates_dataframe.probability.astype(float)

            candidate_mates = candidate_mates_dataframe.groupby((candidate_mates_dataframe.end.shift()-candidate_mates_dataframe.start).lt(0).cumsum()).agg({'chrom':'first','start':'first','end':'last','probability':'sum'})

            sum = np.sum(float(x[3]) for index, x in candidate_mates.iterrows())


        elif np.any(read_types == 'SC') == True and (np.any(read_types == 'DR') == True or np.any(read_types == 'SA') == True):
            #remove lines with sc

            candidate_mates_no_sc = candidate_mates_dataframe.drop(candidate_mates_dataframe[candidate_mates_dataframe.read_type == 'SC'].index)
            candidate_mates_dataframe = candidate_mates_no_sc.sort_values(by=['chrom', 'start', 'end'],ascending=[True, True, True])
            candidate_mates_dataframe['probability'] = candidate_mates_dataframe.probability.astype(float)


            candidate_mates  = candidate_mates_dataframe.groupby((candidate_mates_dataframe.end.shift()-candidate_mates_dataframe.start).lt(0).cumsum()).agg({'chrom':'first','start':'first','end':'last','probability':'sum'})


            sum = np.sum(float(x[3]) for index,x in candidate_mates.iterrows())



        else:
            #only soft clipped

            return(None)

        extended = []


        #if argmax is turn on interval_p is 0
        if interval_p_cutoff == 0:
            #argmax(probability)

            candidate_mates =  candidate_mates.loc[candidate_mates['probability'] == candidate_mates['probability'].max()]

            for item,row in candidate_mates.iterrows():
                if ('LR' in orientation) or ('L' and 'R' in orientation):


                    start = row['start'] - interval_extension

                    end = row['end'] + interval_extension

                    if start < 0:
                        extended.append([row['chrom'], str(0), int(round(end))])

                    else:
                        extended.append([row['chrom'], int(round(start)), int(round(end))])

                elif 'L' in orientation:

                    start = row['start'] - interval_extension

                    if start < 0:
                        extended.append([row['chrom'], str(0), row['end']])

                    else:
                        extended.append([row['chrom'], int(round(start)), row['end']])

                elif 'R' in orientation:

                    end = row['end'] + interval_extension

                    extended.append([row['chrom'], row['start'], int(round(end))])

                return (pd.DataFrame.from_records(extended, columns=['chrom', 'start', 'end']))

        else:

            for index,interval in candidate_mates.iterrows():

                #small pseudocount to denominator to avoid div by zero

                if interval['probability']/(sum+0.00000001) >= interval_p_cutoff:

                    if ('LR' in orientation) or ('L' and 'R' in orientation):


                        start = interval['start'] - interval_extension

                        end = interval['end'] + interval_extension

                        if start < 0:
                            extended.append([interval['chrom'], str(0), int(round(end))])

                        else:
                            extended.append([interval['chrom'], int(round(start)), int(round(end))])

                    elif 'L' in orientation:

                        start = interval['start'] - interval_extension

                        if start < 0:
                            extended.append([interval['chrom'], str(0), interval['end']])

                        else:
                            extended.append([interval['chrom'], int(round(start)), interval['end']])

                    elif 'R' in orientation:

                        end = interval['end'] + interval_extension

                        extended.append([interval['chrom'], interval['start'], int(round(end))])


            return(pd.DataFrame.from_records(extended,columns=['chrom','start','end']))


    except BaseException as e:

        if verbose < 2:

            warnings.warn(
                "WARNING: Could not compute the probability for the mate interval priors %s due to the following error %s \n Skipping intervals" % (
                str(bed_prior), str(e)))




def circle_from_SA(read,mapq_cutoff,mate_interval):

    """Function that takes as input a read (soft-clipped) with a Suplementary alignment the mapping quality cut-off
    and the mate intervals and checks if it fits the conditions to call a circle. Will return True if the supplementary
    alignment matches the interval"""

    suplementary = read.get_tag('SA')

    #this list will have the following information [chr,left_most start,"strand,CIGAR,mapq, edit_distance]

    supl_info = [x.strip() for x in suplementary.split(',')]

    #mapq filter
    if int(supl_info[4]) > mapq_cutoff:
        #chromosome filter
        if supl_info[0] == mate_interval['chrom']:
            #aligned to the mate interval
            if int(mate_interval['start']) < int(supl_info[1]) < int(mate_interval['end']):

                #orientation
                if read.is_reverse == True and supl_info[2] == '-':
                    return{'support' : True, 'leftmost': int(supl_info[1]), 'cigar' : supl_info[3]}

                elif read.is_reverse == False and supl_info[2] == '+':

                    return{'support' : True, 'leftmost' : int(supl_info[1]), 'cigar' : supl_info[3]}

            else:

                return{'support' : False}

        else:

            return{'support' : False}

    else:
        return{'support' : False}

def number_encoding(seq):
    """Function that takes as input a DNA sequence, and encodes the sequence to numbers, so that it can be accelerated
    with numba"""
    encoded = []
    for i in seq:
        if i == "A":
            encoded.append(1)
        elif i == "T":
            encoded.append(2)
        elif i == "C":
            encoded.append(3)
        elif i == "G":
            encoded.append(4)
    return(np.array(encoded))


def check_alphabet(sequence):
    """Function that takes as input a sequence and it will check that there is at least a letter matching the alphabet
     in the sequence, returning true."""

    code = "ATCG"

    for base in sequence:
        if base in code:
            return(True)
    return(False)

def check_compatibility(seq1,seq2):
    """Function that takes as input two DNA sequence and checks whether their alphabets have at least one element
    in common. This due to an old bug in edlib"""

    for base in seq1:

        for base2 in seq2:

            if base == base2:

                return(True)

    return(False)

@jit(nopython=True)
def phred_to_prob(values):
    """Function that takes as input a numpy array with phred base quality scores and returns an array with base probabi-
    lity scores"""
    return(10**((values*-1)/10))

def get_longest_soft_clipped_bases(read):
    """Function that takes as input the cigar string and returns a dictionary containing the longest soft-clipped part of
     the read, the quality values and  the read mapping quality"""

    read_cigar = read.cigar


    #get index of the soft-clipped in the cigar
    match_index = [x for x in range(len(read_cigar)) if read_cigar[x][0] == 4]


    # soft-clipped in only one side
    if len(match_index) == 1:

        #return first n soft-clipped
        if match_index == [0]:
            return{'seq': read.seq[0:read_cigar[0][1]],'qual': read.query_qualities[0:read_cigar[0][1]],'mapq':read.mapq}

        #return last n nucleotides
        elif match_index[0] == (len(read_cigar)-1):

            return {'seq':read.seq[-read_cigar[match_index[0]][1]:],
                    'qual':read.query_qualities[-read_cigar[match_index[0]][1]:],'mapq':read.mapq}






    #soft-clipped in both sides of the read
    else:

        #make sure that is soft-clipped on both sides

        try:

            assert read_cigar[0][0] == 4 and read_cigar[-1][0] == 4

            # longest soft-clipped are first n nucleotides
            if read_cigar[0][1] >= read_cigar[-1][1]:



                return {'seq': read.seq[0:read_cigar[0][1]],'qual': read.query_qualities[0:read_cigar[0][1]],
                        'mapq':read.mapq}

            else:

                return{'seq':read.seq[-read_cigar[-1][1]:],'qual': read.query_qualities[-read_cigar[-1][1]:],
                       'mapq': read.mapq}

        except AssertionError as e:

            print(e)

def background_freqs(seq):
    """Function that takes as input the sequence of the nucletide frequencies in the realignment interval"""

    return{nucleotide: seq.count(nucleotide)/len(seq) for nucleotide in 'ATCG'}





def realign(read,n_hits,plus_strand,minus_strand,plus_base_freqs,minus_base_freqs,gap_open,gap_extend,verbose,max_edit):


    """Function that takes as input a read, the number of hits to find and the plus and minus strand and will return
    the number of hits, the sequencing qualities for that read and the g+c content of the realignment interval"""


    #get soft-clipped read
    soft_clipped_read = get_longest_soft_clipped_bases(read)

    #encoding of DNA and operations A,T,C,G,=,X,DI. THis is done for Numba
    nuc_and_ops =  np.array([1,2,3,4,5,6,7])
    encoded_nucs = number_encoding(soft_clipped_read['seq'])

    hits = 0

    min_score = len(soft_clipped_read['seq'])


    top_hits = {}


    if read.is_reverse:

        while hits < n_hits and min_score >= -10:

            alignment = edlib.align(soft_clipped_read['seq'], minus_strand, mode='HW', task='path')
            if hits ==0:
                if alignment['editDistance'] > max_edit:
                    return(None)



            for location in alignment['locations']:



                mask_bases = 'X' * ( location[1] - location[0])


                minus_strand = minus_strand[:location[0]] + mask_bases + minus_strand[location[1]:]

                hits += 1


                score = pssm(phred_to_prob(np.array(soft_clipped_read['qual'],dtype=np.float64)), encoded_nucs,
                             edlib_cigar_to_iterable(alignment['cigar']),minus_base_freqs,gap_open,gap_extend,nuc_and_ops,verbose)

                if score < min_score:
                    min_score = score


                top_hits[hits] = (location,alignment['cigar'],score,alignment['editDistance'])

        else:
            # the search was exaustive
            hits +=n_hits

    else:
        #min socre stops the search if the score is orders of magnitude smaller that the top score given the edit
        #distance
        while hits < n_hits and min_score >= -10:



            alignment = edlib.align(soft_clipped_read['seq'], plus_strand, mode='HW', task='path')
            #stop search if edit distance is to high
            if hits ==0:
                if alignment['editDistance'] > max_edit:
                    return (None)


            for location in alignment['locations']:

                mask_bases = 'X' * ( location[1] - location[0])

                plus_strand = plus_strand[:location[0]] + mask_bases + plus_strand[location[1]:]

                hits += 1

                score = pssm(phred_to_prob(np.array(soft_clipped_read['qual'],dtype=np.float64)), encoded_nucs,
                             edlib_cigar_to_iterable(alignment['cigar']), plus_base_freqs,gap_open,gap_extend,nuc_and_ops,verbose)

                if score < min_score:
                    min_score = score

                top_hits[hits] = (location,alignment['cigar'],score,alignment['editDistance'])

        else:

            hits +=n_hits




    return({'alignments':top_hits,'mapq_prior': soft_clipped_read['mapq']})


def edlib_cigar_to_iterable(edlib_cigar):
    """Function that takes as input the edlib cigar and parses it to get it in a iterable manner"""
    #encoding of DNA and operations A,T,C,G,=,X,ID
    #nuc_and_ops =  np.array([1,2,3,4,5,6,7])

    length = []
    operations = []

    for i in re.findall(r'\d+[IDX=]',edlib_cigar):
        length.append(int(i[0]))
        if i[1] == '=':
            operations.append(5)
        elif i[1] == 'X':
            operations.append(6)
        elif i[1] == 'I' or 'D':
            operations.append(7)


    return(np.array(length),np.array(operations))


@jit(nopython=True)
def pssm(seq_prob,seq_nucl,iterable_cigar,base_freqs,gap_open,gap_extend,nuc_and_ops,verbose):
    """Function that takes as input the sequencing probabilities and cigar string and returns the log2 pssm of the read"""





    #start positon to operate in the pssm. This is done to iterate over the operations in the cigar, and keep track of
    # were I am in the seq and quality values
    seq_pos = 0
    indel_penalty = 0



    #iterate trough CIGAR operations
    for index in range(0,len(iterable_cigar[0])):

        operation_length = iterable_cigar[0][index]
        end = operation_length + seq_pos



        operation = iterable_cigar[1][index]


        #match, 1 minus prob(base called wrong)
        if operation == nuc_and_ops[4]:

            for nucleotide in range(seq_pos,end):

                if seq_nucl[nucleotide] == nuc_and_ops[0]:

                    seq_prob[nucleotide] = np.log2((1 - (seq_prob[nucleotide]))/base_freqs[0])

                elif seq_nucl[nucleotide] == nuc_and_ops[1]:

                    seq_prob[nucleotide] = np.log2((1 - (seq_prob[nucleotide])) / base_freqs[1])

                elif seq_nucl[nucleotide] == nuc_and_ops[2]:

                    seq_prob[nucleotide] = np.log2((1 - (seq_prob[nucleotide])) / base_freqs[2])

                elif seq_nucl[nucleotide] == nuc_and_ops[3]:

                    seq_prob[nucleotide] = np.log2((1 - (seq_prob[nucleotide])) / base_freqs[3])



            seq_pos += operation_length



        elif operation == nuc_and_ops[5]:


            for nucleotide in range(seq_pos,end):

                if seq_nucl[nucleotide] == nuc_and_ops[0]:

                    seq_prob[nucleotide] = np.log2(
                        (seq_prob[nucleotide]/3)/base_freqs[0])


                elif seq_nucl[nucleotide] == nuc_and_ops[1]:

                    seq_prob[nucleotide] = np.log2(
                        (seq_prob[nucleotide]/3)/base_freqs[1])

                elif seq_nucl[nucleotide] == nuc_and_ops[2]:

                    seq_prob[nucleotide] = np.log2(
                        (seq_prob[nucleotide]/3)/base_freqs[2])


                elif seq_nucl[nucleotide] == nuc_and_ops[3]:

                    seq_prob[nucleotide] = np.log2(
                        (seq_prob[nucleotide]/3)/base_freqs[3])




                elif seq_nucl[nucleotide] == nuc_and_ops[6]:

                    if verbose < 2:

                        seq_prob[nucleotide] = 0
                        print("Warning:Ambiguous base found in nucleotide sequence. Assigning score of 0 in the log2 pssm")

            seq_pos += operation_length


        elif operation == nuc_and_ops[6]:

            #affine gap scoring model
            indel_penalty += gap_open + gap_extend*(operation_length-1)


    return(np.sum(seq_prob)-indel_penalty)


def realignment_probability(hit_dict,interval_length):
    """Function that takes as input the realignment dictionary and returns the alignment probability of the best hit"""


    best_hit = hit_dict['alignments'][1][2]

    #this might be included on the denominator

    posterior = 2**best_hit/(np.sum((2**value[2]) for key,value in hit_dict['alignments'].items()))

    return(posterior)







def fraction(start1,start2,end1,end2,read1,read2):
    """Function that performs a first round of merging. If the realigned intervals and SA intervals overlap, and are ca-
    lled within the same iteration (which means that it is the same circle probably) they will be merged"""

    #check that they come from the same read
    read_match = (read1 == read2)*1


    #calculate distance between the two intervals
    distance = (abs(start1-start2) + abs(end1-end2))

    #overlap of interval 1 on interval 2
    one_overlap_two = 1 - (distance/(end1-start1))
    #overlap of interval two on interval 1
    two_overlap_one =  1 - (distance/(end2-start2))

    return(one_overlap_two + two_overlap_one + read_match)



def merge_fraction(chrom1,x1,x2,chrom2,y1,y2):
    """compute overlap (reciprocal) of the interval y over interval x"""

    distance = (np.minimum(x2.values,y2.values) - np.maximum(x1.values,y1.values))



    one_overlap_two = distance/(y2.values-y1.values)

    two_overlap_one = distance/(x2.values-x1.values)


    # check if they are on the same chromosome and the amount of overlap if so
    return(pd.Series(chrom1 == chrom2) + pd.Series(two_overlap_one.clip(0)) + pd.Series(one_overlap_two.clip(0)))


def iteration_merge(only_discordants,results,fraction,splits,score,sc_len,bam,af,insert,std,n_discordant):
    """finction that merges the results of every iteration and filters the data by allele frequency"""

    norm_fraction = 3

    parsed_discordants = []
    for interval in only_discordants:
        interval.append(0)
        parsed_discordants.append(interval)


    discordant_bed = bt.BedTool(parsed_discordants)





    unparsed_pd = pd.DataFrame.from_records(results,
        columns=['chrom', 'start', 'end', 'read', 'iteration','score', 'discordants'])

    unparsed_pd = unparsed_pd.sort_values(['iteration','chrom','start','end']).reset_index()


    grouped = unparsed_pd.groupby(merge_fraction(unparsed_pd.iteration.shift(), unparsed_pd.start.shift(),
                                           unparsed_pd.end.shift(), unparsed_pd.iteration,
                                           unparsed_pd.start,
                                           unparsed_pd.end).lt(norm_fraction).cumsum()).agg(
        {'chrom': 'first', 'start': 'min', 'end': 'max', 'discordants': 'max', 'read': 'sum','score':'sum'})

    bedtool_output = bt.BedTool.from_dataframe(grouped)






    allele_free = bedtool_output.cat(discordant_bed, postmerge=False)
    write = []

    for interval in allele_free:
        try:
            if int(interval[4]) != 0:
                if (int(interval[4])) >= splits and float(interval[5]) > score:
                    start_cov = bam.count(contig=interval[0],
                                                   start=int(interval[1]), stop=int(interval[1])+1
                                                   ,read_callback='nofilter')

                    end_cov = bam.count(contig=interval[0],
                                                 start=int(interval[2])-1, stop=int(interval[2])
                                                 ,read_callback='nofilter')



                    circle_af = ((int(interval[4]) * 2)) / ((start_cov+end_cov+0.01)/2)
                    if circle_af >=af:
                        write.append(interval)
            else:
                if int(interval[3]) >= n_discordant:
                        start_cov = bam.count(contig=interval[0],start=int(interval[1]), stop=int(interval[1]) + 1,
                                                   read_callback='nofilter')

                        end_cov = bam.count(contig=interval[0],
                                                 start=int(interval[2]) - 1, stop=int(interval[2]),
                                                 read_callback='nofilter')

                        circle_af = (int(interval[3])) / ((start_cov+end_cov+0.01)/2)

                        if circle_af >= af:
                            write.append(interval)
        except BaseException as e:
            print(e)
            pass

    return(bt.BedTool(write))





def merge_final_output(bam,results,begin,splits,dir,fraction,pid):

    """Function that takes as input the final results, and merge reciprocal intervals (this is done to combine the output
    of different clusters)"""



    bam = ps.AlignmentFile(bam, "rb")
    os.chdir("temp_files_%s/" % pid)

    # multiply *2 for reciprocal overlap +1 to check chromosome
    norm_fraction = (fraction*2)+1

    unparsed_bed = bt.BedTool(results)




    print("Writting final output to disk")


    unparsed_pd = unparsed_bed.to_dataframe(
        names=['chrom', 'start', 'end', 'discordants', 'sc','score'])



    second_merging_round = unparsed_pd.sort_values(['chrom', 'start', 'end']).reset_index()
    #merge the output
    # merge_fraction calculates the degree of overlap between the two genomic intervals
    #lt(norm_freaction) looks the ones that surpass the merging threshold (returns 0 if true, 1 if not)
    # Cumsum calculates the cumulative sum over the output of lt. Which is then used for the grouping. 
    #If the cumulative sum is the same for two rows, they are merged
    final_output = second_merging_round.groupby(
        merge_fraction(second_merging_round.chrom.shift(), second_merging_round.start.shift(),
                     second_merging_round.end.shift(),second_merging_round.chrom,second_merging_round.start,second_merging_round.end).lt(norm_fraction).cumsum()).agg(
        {'chrom': 'first', 'start': 'min', 'end': 'max', 'discordants' : 'max', 'sc': 'sum','score':'sum'})

    unfiltered_output = bt.BedTool.from_dataframe(final_output)

    # filter splits

    filtered = []
    for interval in unfiltered_output:

        if (int(interval[4])+int(interval[3])) >= splits:
            if int(interval[1]) != 0:
                interval[1] = int(interval[1])+1
            filtered.append(interval)

    filtered_output = bt.BedTool(filtered)

    os.chdir("%s" % dir)


    print("Finished!")

    end = time.time()

    total_time = (end - begin) / 60


    print("\nCircle-Map realign finished indentifying circles in %s \n" % total_time)
    print("\nCircle-Map has identified %s circles\n" % len(filtered_output))



    return(filtered_output)


def write_to_disk(partial_bed,output,locker,dir,pid):

    """function that writes to disk the results of every worker thread"""


    locker.acquire()
    os.chdir("%s/temp_files_%s/" % (dir,pid))
    output_bed = bt.BedTool('%s' % output)
    writer_bed = output_bed.cat(partial_bed,postmerge=False)
    writer_bed.saveas('%s' % output)
    os.chdir("%s" % dir)
    locker.release()

def start_realign(circle_bam,output,threads,verbose,pid,clusters):
    """Function that start the realigner function
        - Splits the clusters to cores and removes the from disk the bedtools intermediates"""

    begin = time.time()

    print("\nRunning Circle-Map realign\n")

    print("Clustering structural variant reads\n")

    eccdna_bam = ps.AlignmentFile("%s" % circle_bam, "rb")

    sp.call("mkdir temp_files_%s" % pid, shell=True)


    sorted_bam,peaks = bam_circ_sv_peaks(eccdna_bam,circle_bam,threads,verbose,pid,clusters)

    chunks = int(math.ceil(len(peaks)/(threads*100)))


    splitted = [peaks[x:x + chunks] for x in range(0, len(peaks),chunks)]


    # split to cores

    print("\nSplitting coverage file to cores\n")
    os.chdir("temp_files_%s" % pid)
    sp.call("touch %s" % output, shell=True)
    os.chdir("../")

    #this releases from tmp file the unmerged and peak file
    bt.cleanup()

    return(splitted,sorted_bam,begin)

def start_simulate(pid):
    """Function for starting Circle-Map simulate"""

    print("\nRunning Circle-Map Simulate\n")


    sp.call("mkdir temp_files_%s" % pid, shell=True)


    return(pid)

def mutate(genome,pid,indel,snp,java_mem):
    """Function that takes as input the path of the genome,the indel ans substitution rate, and it will create a sinthetic
    genome introducing random mutations on the fasta sequence and providing a vcf"""

    print("Introducing mutations in the fasta genome")
    print("\t Indel rate: %s" % indel)
    print("\t Substitution rate: %s" % snp)
    sp.call("mutate.sh %s in=%s out=temp_files_%s/mutated.fa subrate=%s indelrate=%s" % (java_mem,genome,pid,snp,indel),shell=True)

    print("Simulating reads")

    return(None)






def check_size_and_write(results,only_discortants,output,lock,directory,fraction,pid):
    """Function that checks if the intervals in memory are to big. And writes them to disk to release memory."""


    if sys.getsizeof(results) < 100000000:
        return(False)


    else:

        partial_bed = iteration_merge(only_discortants, results,fraction)

        print("Writting %s circular intervals to disk" % len(partial_bed))
        write_to_disk(partial_bed,output,lock,directory,pid)

        return(True)

def merge_coverage_bed(results,frac,number):

    """Function that takes as bed file containing the coordinates of the double mapped reads and
    returns the merged bed file containing the information about the clusters"""

    fraction = (frac*2)+1

    unparsed_pd = pd.DataFrame.from_records(results,columns=['chrom', 'start', 'end','item'])



    sort = unparsed_pd.sort_values(by=['chrom', 'start', 'end']).reset_index(drop=True)

    merging_out  = sort.groupby(
        merge_fraction(sort.chrom, sort.start,
                     sort.end,sort.chrom.shift(),sort.start.shift(),sort.end.shift()).lt(fraction).cumsum()).agg(
        {'chrom': 'first', 'start': 'min', 'end': 'max','item': 'sum'})

    merging_out = merging_out.drop(merging_out[merging_out.item < number].index)

    merging_out = merging_out.sort_values(by=['chrom', 'start', 'end']).reset_index(drop=True)







    final_output = merging_out.groupby(
        merge_fraction(merging_out.chrom, merging_out.start,
                       merging_out.end, merging_out.chrom.shift(), merging_out.start.shift(),merging_out.end.shift()).lt(fraction).cumsum()).agg(
        {'chrom': 'first', 'start': 'min', 'end': 'last', 'item': 'sum'})


    bedtool_output = bt.BedTool.from_dataframe(final_output)

    return(bedtool_output)

def filter_by_ratio(eccdna_bed,cutoff):
    """Function that takes as input the eccDNA bed and returns the data filtered by tha change at the start and the end
    """

    #circle list is a shared memory object
    circle_list = []
    unparsed_pd = eccdna_bed.to_dataframe(
        names=['chrom', 'start', 'end', 'discordants', 'soft-clipped', 'score', 'mean','std','start_ratio','end_ratio','continuity'])
    for item, row in unparsed_pd.iterrows():


        if float(row[8]) > cutoff or float(row[9]) > cutoff:
            circle_list.append([row['chrom'],row['start'],row['end'],row['discordants'],row['soft-clipped'],
                             row['score'],row['mean'],row['std'],row['start_ratio'],row['end_ratio'],row['continuity']])

    output = pd.DataFrame.from_records(circle_list,columns=['chrom', 'start', 'end', 'discordants', 'soft-clipped', 'score', 'mean','std',
                          'start_ratio','end_ratio','continuity'])

    return(output)


def merge_bed(discordants_pd):
    """Function that takes as input a bed file and returns a pandas dataframe indicating if files should be merged. This
    function will merge everything that is overlapping by at least 1bp"""
    #check range overlap
    overlap = ((discordants_pd.start - discordants_pd.shift().end) - 1).lt(0)
    #check chr overlap
    chr_overlap = (discordants_pd.chrom == discordants_pd.shift().chrom)
    #if both bools are succesful returns a 2
    return ((overlap * 1 + chr_overlap * 1).lt(2).cumsum())

def assign_discordants(split_bed,discordant_bed,insert_mean,insert_std):
    """Function that takes as input the the discordant reads supporting an interval and assigns them to the
    interval if they are close by (using the insert size estimate)"""

    max_dist = (insert_mean / 2) + (3 * insert_std)

    splits = pd.DataFrame.from_records(split_bed, columns=['chrom', 'start', 'end', 'read', 'iteration',
                                                           'score']).sort_values(['chrom', 'start', 'end'])

    splits['score'] = splits['score'].astype(float)

    merged_splits = splits.groupby(['chrom', 'start', 'end', 'iteration']).agg(
        {'chrom': 'first', 'start': 'first', 'end': 'max', 'read': 'nunique', 'iteration': 'first', 'score': 'sum'})

    merged_splits['read'] = merged_splits['read'].astype(int)
    discordant_bed = pd.DataFrame.from_records(discordant_bed,columns=['chrom', 'start', 'end', 'read'])

    if len(discordant_bed) > 0:
        assigned_splits = []

        for index, row in merged_splits.iterrows():
            chrom_filt = discordant_bed[(discordant_bed['chrom'] == row['chrom'])]
            start_filt = chrom_filt[
                (chrom_filt['start'] > row['start']) & ((chrom_filt['start'] - row['start']) < max_dist)]
            end_filt = start_filt[(start_filt['end'] < row['end']) & ((row['end'] - start_filt['end']) < max_dist)]

            assigned_splits.append(
                [row['chrom'], row['start'], row['end'], row['read'], row['iteration'], float(row['score']), len(end_filt)])

        return (assigned_splits)

    else:
        assigned_splits = []
        for index,row in merged_splits.iterrows():
            assigned_splits.append(
                [row['chrom'], row['start'], row['end'], row['read'], row['iteration'], float(row['score']),
                 0])

        return(assigned_splits)

def adaptative_myers_k(sc_len,edit_frac):
    """Calculate the edit distance allowed as a function of the read length"""
    return(float(sc_len*edit_frac))

def non_colinearity(read,mate_interval):
    """Input a read and the mate interval in the graph. The function checks whether the alignment would be linear (splicing)
    or colinear. Will return false, in order to not attemp realignment. This is mainly thought for skipping deletions and
    RNA splicing"""

    #check left soft-clipped
    if read.cigar[0][0] == 4:
        #graph needs to be upstream or looping to itself
        if int(mate_interval.start) > read.pos:
            return (True)
        elif read.pos < int(mate_interval.end):
            #looping to itself
            return (True)
        else:
            return (False)
    #check right softclipped
    if read.cigar[-1][0] == 4:
    # graph needs to be downstream or looping to itself
        if int(mate_interval.end) < read.pos:
            return (True)
        elif read.pos > int(mate_interval.start):
            #looping to itself
            return (True)
        else:
            return (False)