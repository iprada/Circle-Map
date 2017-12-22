from simulations import sim_ecc_reads

sim_ecc_reads("/data/xsh723/scratch/hg38/canonical_hg38/hg38.fa","/data/xsh723/scratch/hg38/canonical_hg38/hg38.fa",100,paired_end = True,directory="/isdata/kroghgrp/xsh723/circle_map/test_data/simulations",reads=120000000,temp_fastq="temp.fastq",insert_size=300,errors=False)
