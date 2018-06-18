#!/isdata/kroghgrp/xsh723/scratch/data_partition/anaconda/bin/python3.6
from realigner import realignment

create_object = realignment("sorted_circle_qname_B03_aln.bam","qname_B03_aln.bam",
"/isdata/kroghgrp/xsh723/scratch/data_partition/scratch/S288C_15_01_13/saccer_plasmids.fa","/isdata/kroghgrp/xsh723/projects/circle_map/test_data/yeast_data",
                            mapq_cutoff=10,insert_size_mapq = 60,std_extension=4,insert_size_sample_size=100000,ncores=10,
                            gap_open=17,gap_ext=5,n_hits=10,prob_cutoff=0.01,min_soft_clipped_length=19,
                            overlap_frac = 0.99,interval_p_cut=0.01,circle_peaks=)
create_object.realign()



