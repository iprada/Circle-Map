import pysam as ps
import pandas as pd
import os



os.chdir("/home/inigo/msc_thesis/mapping_stats/")
#os.system("samtools sort -n -o query_name_sorted_paired_end_sim_aln.bam.sv.bam sorted_paired_end_sim_aln.bam.sv.bam")
samfile = ps.Samfile("query_name_sorted_paired_end_sim_aln.bam.sv.bam", "rb")

for read in samfile:
    query = read.query_name
    query_info = query.split('|')

    break