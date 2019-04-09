import cythonize_realigner
import pysam as ps
from numba import jit
import numpy as np
import edlib
import re
import time



bam = ps.AlignmentFile("/home/iprada/faststorage/projects/6_aged_yeast/working_directory/aligned/B05/sorted_B05.bam","rb")

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
def phred_to_prob(values):
    """Function that takes as input a numpy array with phred base quality scores and returns an array with base probabi-
    lity scores"""
    return(10**((values*-1)/10))


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

def number_encoding(seq):
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

def realign(read,n_hits,plus_strand,minus_strand,plus_base_freqs,minus_base_freqs,gap_open,gap_extend,verbose,max_edit):


    """Function that takes as input a read, the number of hits to find and the plus and minus strand and will return
    the number of hits, the sequencing qualities for that read and the g+c content of the realignment interval"""


    #get soft-clipped read
    soft_clipped_read = get_longest_soft_clipped_bases(read)

    #encoding of DNA and operations A,T,C,G,=,X,DI
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


i = 0
st_seq = ''
seq = ''
st_read = ''
freqs = [0.25,0.25,0.25,0.25]
set_first = 0
for read in bam.fetch(until_eof=True):
    if set_first == 0:
        if is_soft_clipped(read) == True:
            set_first = 1

            if i == 0:
                st_seq = read.seq
                st_read = read
                seq = read.seq
                i +=1
        else:
            continue

    else:
        seq = seq + read.seq
        i+=1
    if i == 11:
        break
    else:
        continue

print(len(seq))
st_seq = st_seq[0:15]

begin = time.time()
for i in range(0,10):
    realign(st_read,10,seq,seq,freqs,freqs,5,1,1,0.05)
end = time.time()

print((end-begin)/60)





