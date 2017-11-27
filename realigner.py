#!/data/xsh723/anaconda/bin/python3.6

import subprocess as sb
import os
import pysam as ps
import pybedtools as bt
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import time
from utils import is_soft_clipped,is_hard_clipped,bam_circ_sv_peaks

class realigner:
    """Class for managing the realignment and eccDNA indetification of circle-map"""

    def __init__(self, input_bam,genome_fasta):
        self.input_bam = input_bam
        self.genome_fa = genome_fasta



    def realignment(self):
        """Function that will iterate trough the bam file containing reads indicating eccDNA structural variants and
        will output a bed file containing the soft-clipped reads, the discordant and the coverage within the interval"""

        eccdna_bam = ps.AlignmentFile("%s" % self.input_bam, "rb")

        print("aaaa")

        eccdna_bam.close()
