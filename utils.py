#!/data/xsh723/anaconda/bin/python3.6
#Author: Inigo Prada Luengo
#email: inigo.luengo@bio.ku.dk

def is_soft_clipped(read):
    """Function that checks the CIGAR string of the sam file and returns true if the read is soft-clipped"""
    for cigar in read.cigar:
        if cigar[0] == 4:
            return (True)
        else:
            return (False)



