#!/data/xsh723/anaconda/bin/python3.6
#Author: Inigo Prada Luengo
#email: inigo.luengo@bio.ku.dk

def is_soft_clipped(read):
    """Function that checks the CIGAR string of the sam file and returns true if the read is soft-clipped"""

    # cigar 4 equals to H in pysam sam representation
    for cigar in read.cigar:
        if cigar[0] == 4:
            return (True)
        else:
            return (False)

def is_hard_clipped(read):
    """Function that checks the CIGAR string of the sam file and returns true if the read is hard-clipped"""

    #cigar 5 equals to H in pysam sam representation
    for cigar in read.cigar:
        if cigar[0] == 5:
            return (True)
        else:
            return (False)

def bam_circ_sv_peaks(bam):
    """Function that takes as input a bam file and returns a merged bed file of the genome covered by the bam"""

    peak_coverage = bam.genome_coverage(bg=True)

    #sort and merge as sanity
    sorted_peak_coverage = peak_coverage.sort()
    merged_peak_coverage = sorted_peak_coverage.merge()

    return(merged_peak_coverage)







