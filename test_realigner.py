#!/data/xsh723/anaconda/bin/python3.6
from realigner import realignment

create_object = realignment("test_chr1_small/chr1_629063_634996.bam","qname_t08.bam",
                            "/data/xsh723/scratch/hg38/canonical_hg38/hg38.fa","/isdata/kroghgrp/xsh723/circle_map/test_data/real_data_t01/",
                            mapq_cutoff=10,insert_size_mapq = 60,std_extension=4,insert_size_sample_size=100000,ncores=10,gap_open=17,
                            gap_ext=1,n_hits=10,prob_cutoff=0.1,min_soft_clipped_length=8,overlap_frac = 0.99,interval_p_cut=0.01)
create_object.realignment()


