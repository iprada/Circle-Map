from circlemap.utils import start_realign,insert_size_dist
import multiprocessing as mp
from circlemap.realigner import realignment
import os
from tqdm import *

input = "/home/iprada/faststorage/projects/6_aged_yeast/working_directory/aligned/BM3/sort_circle_qname_BM3.bam"
qbam = "/home/iprada/faststorage/projects/6_aged_yeast/working_directory/aligned/BM3/qname_BM3.bam"
sort_bam = "/home/iprada/faststorage/projects/6_aged_yeast/working_directory/aligned/BM3/sorted_BM3.bam"
fasta = "/home/iprada/faststorage/reference_Data/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/8_plasmids_genome/yeast_8_plasmids.fa"

splitted, sorted_bam, begin = start_realign(input,"profiling_output.bed", 1,3,1,500)

sorted_bam.close()
#get global insert size prior
metrics = insert_size_dist(100000,60,qbam)



lock = mp.Lock()

object = realignment(input, qbam, sort_bam, fasta,
                     os.getcwd(),
                     20,60, 4, 100000,5,1, 10, 0.99, 6,0.95,0.01,"profiling_output.bed",16,0.1, lock, 0,0.0,3,1,
                     0.05, False,False, 0,0.0, metrics, 3)



object.realign(splitted[0])
