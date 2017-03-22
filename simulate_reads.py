from Bio import SeqIO
import Bio as bio
from io import StringIO
import numpy as np
def sim_single_end(genome_fa,chr,chr_pos_start,chr_pos_end,read_length, unique_id):
    # create fastafile object
    fastafile = ps.FastaFile(genome_fa)
    #pick a random position in the circle
    start = np.random.randint(chr_pos_start,(chr_pos_end+1))
    end = start + read_length
    #if the end position if bigger than the chr_end_position, that indicates a circle
    if end > chr_pos_end:
        print("backspliced")
        #back spliced_read


        #number of right nucleotides to take
        right_backspliced = end - chr_pos_end
        #right nucleotides position
        right_nucleotides = chr_pos_start + right_backspliced
        # get the right split read
        right_split_read = fastafile.fetch(chr, chr_pos_start,right_nucleotides)
        #left split read nucleotide position
        left_nucleotides_start = chr_pos_end - (read_length - right_backspliced)
        #get left split read
        left_split_read = fastafile.fetch(chr, left_nucleotides_start, chr_pos_end)
        # reverse it, this is how circles are formed
        reversed_left_split_read = left_split_read[::-1]
        #put all together
        total_read = reversed_left_split_read + right_split_read
        unique_id = "split-" + unique_id


    else:
        #sample normal read
        print("normal")
        total_read = fastafile.fetch(chr, start,end)


    #get each entry of the fastq file
    seq_id = "id:%s-%s|left_pos:%s-%s|right:%s-%s " % (unique_id,chr, int(chr_pos_end - (read_length / 2)), int(chr_pos_end), int(chr_pos_start),int(chr_pos_start + (read_length / 2)))
    #right now, the quality is maximum
    quality = "I" * read_length
    #put all together
    fastq_string = "@%s\n%s\n+\n%s\n" % (seq_id, total_read, quality)
    new_record = SeqIO.read(StringIO(fastq_string), "fastq")
    return(new_record)



def sim_paired_end(insert_size,chr_pos_start,chr_pos_end,read_length, unique_id):
    """Function that simulates paired-end reads"""
    #left split read
    start = np.random.randint(chr_pos_start, (chr_pos_end + 1))
    end = start + insert_size
    if end > chr_pos_end:
        # left split read


