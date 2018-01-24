#!/data/xsh723/anaconda/bin/python3.6
from realigner import realignment

create_object = realignment("coordinate_circle_qname_sorted_paired_end_sim_aln.bam","qname_sorted_paired_end_sim_aln.bam",
                            "/data/xsh723/scratch/hg38/canonical_hg38/hg38.fa","/isdata/kroghgrp/xsh723/circle_map/test_data/realigner/",
                            mapq_cutoff=10,insert_size_mapq = 60,std_extension=4,insert_size_sample_size=100000,ncores=10,gap_open=6,
                            gap_ext=1,n_hits=10,prob_cutoff=0.01)
create_object.realignment()


