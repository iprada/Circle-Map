from __future__ import division

from skbio.alignment import local_pairwise_align_ssw
from skbio import DNA
from Bio import pairwise2
import time

begin = time.time()


for i in range(0,1000):
    a = local_pairwise_align_ssw(DNA("ACTAAGGCTCTCTACCCCTCTCAGAGA"),DNA("ACTAAGGCTCCTAACCCCCTTTTCTCAGA"))
end = time.time()

print("Ellapsed time: using skbio",(end-begin)/60)

begin = time.time()
for j in range(0,1000):
    b = pairwise2.align.localxx("ACTAAGGCTCTCTACCCCTCTCAGAGA", "ACTAAGGCTCCTAACCCCCTTTTCTCAGA")


a = local_pairwise_align_ssw(DNA("ACTAAGGCTCTCTACCCCTCTCAGAGA"),DNA("ACTAAGGCTCCTAACCCCCTTTTCTCAGA"))
print(a)
b = pairwise2.align.localxx("ACTAAGGCTCTCTACCCCTCTCAGAGA", "ACTAAGGCTCCTAACCCCCTTTTCTCAGA")
print(b)
end = time.time()
print("Ellapsed time: using biopython",(end-begin)/60)