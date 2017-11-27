#!/data/xsh723/anaconda/bin/python3.6

import subprocess as sb
import os
import pysam as ps
import pybedtools as bt
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import time
from utils import is_soft_clipped

class realigner:
    """Class for managing the realignment and eccDNA indetification of circle-map"""

    def __init__(self, input_bam):
        self.input_bam = input_bam
