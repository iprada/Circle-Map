

import pysam as ps
import time

#chr1 length

begin = time.time()
genome_fa = "/data/xsh723/scratch/hg38/canonical_hg38/hg38.fa"

fastafile = ps.FastaFile(genome_fa)
   # get the sequence
fasta = fastafile.fetch("chr1", 1320037,1320087)
fasta2 = fastafile.fetch("chr1", 1406347,1406397)

#CGGGGTAAGAAAAAAAAACCGGGGTACGTTAGATGCTAGGTTCATGAGGCAGGGAAGGGCAGGGGGCCAGCAGGAGTGCTGTGGCCGTCCAGACGAGGCC
print(fasta2.upper() + fasta.upper())
print("TCGGGGTAAGAAAAAAAAACCGGGGTACGTTAGATGCTAGGTTCATGAGGAGGGAAGGGCAGGGGGCCAGCAGGAGTGCTGTGGCCGTCCAGACGAGGCC")


exit()



end = time.time()

print((end-begin)/60)

a = fasta.upper()

end = time.time()

print((end-begin)/60)

