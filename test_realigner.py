#!/isdata/kroghgrp/xsh723/scratch/data_partition/anaconda/bin/python3.6
from realigner import realignment

create_object = realignment("coordinate_circle_qname_sorted_BWAB03.bam","qname_sorted_BWAB03.bam",
"/isdata/kroghgrp/xsh723/students/rasmus/reference_files/catfiles.fasta","/isdata/kroghgrp/xsh723/projects/circle_map/test_data/yeast_data/",
                            mapq_cutoff=10,insert_size_mapq = 60,std_extension=4,insert_size_sample_size=100000,ncores=10,
                            gap_open=17,gap_ext=1,n_hits=10,prob_cutoff=0.1,min_soft_clipped_length=8,
                            overlap_frac = 0.99,interval_p_cut=0.01)
create_object.realignment()



