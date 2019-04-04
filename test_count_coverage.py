import pysam as ps
import time
import numpy as np

bam = ps.AlignmentFile("/home/iprada/faststorage/projects/6_aged_yeast/working_directory/aligned/B05/sorted_B05.bam")

begin = time.time()
for i in range(0,10):
    print(i)
    bam.count(contig="chrXII",start=460000,stop=460001)
end = time.time()
print((end-begin)/60)

begin = time.time()
print(bam.count(contig="chrXII",start=460000,stop=460003))
end = time.time()
print((end-begin)/60)
