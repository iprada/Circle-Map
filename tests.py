from __future__ import division

from skbio.alignment import local_pairwise_align_ssw
from skbio import DNA
from Bio import pairwise2
import time

begin = time.time()



a,b,c,d = local_pairwise_align_ssw(DNA("ACTAAGGCTCTCTACCCCTCTCAGAGA"),DNA("ACTAAGGCTCCTAACCCCCTTTTCTCAGA"))



