import pysam as ps
import pandas as pd
import os



os.chdir("/home/jarvis/inigo/circDNA/")
os.system("samtools view -b -h -F 3 sorted_paired_end_sim_aln.bam.sv.bam > circ_structural_variants.bam")
os.system("samtools sort -n -o sorted_circ_structural_variants.bam circ_structural_variants.bam")
exit()
samfile = ps.Samfile("sorted_paired_end_sim_aln.bam.sv.bam", "rb")

<<<<<<< HEAD
=======
for read in samfile:
    print(read.query_name)
    print(read.flags)
    break
>>>>>>> 797a3c4b99dbb4655eb50904c9e05248b5bde38b
