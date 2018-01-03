

import pysam as ps
import time

#chr1 length

begin = time.time()
genome_fa = "/data/xsh723/scratch/hg38/canonical_hg38/hg38.fa"

fastafile = ps.FastaFile(genome_fa)
   # get the sequence
fasta = fastafile.fetch("chr1", 0,248956422)



end = time.time()

print((end-begin)/60)

a = fasta.upper()

end = time.time()

print((end-begin)/60)

