import pysam as ps
import pandas as pd
import os



os.chdir("/home/inigo/msc_thesis/mapping_stats/")
os.system("samtools view -b -h -F 2 sorted_paired_end_sim_aln.bam.sv.bam > circ_calls.bam")
os.system("samtools index circ_calls.bam")

samfile = ps.Samfile("circ_calls.bam", "rb")




for read in samfile:
    print(read.query_name)
    break

