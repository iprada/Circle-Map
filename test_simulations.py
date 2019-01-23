from simulations import sim_ecc_reads
import multiprocessing as mp

sim_ecc_reads("/home/iprada/small_size_projects/canonical_hg38/hg38.fa",100,
              directory="/home/iprada/faststorage/projects/simulations",reads=1000000,
              fastq="temp",insert_size=300,errors=False,mean_cov=10,locker=mp.Lock())
