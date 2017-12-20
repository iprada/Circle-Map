from realigner import realignment

create_object = realignment("coordinate_circle_qname_sorted_paired_end_sim_aln.bam",
                            "/data/xsh723/scratch/hg38/hg38.fa","/isdata/kroghgrp/xsh723/circle_map/test_data/realigner/",mapq_cutoff=10,ncores=1)
create_object.realignment()

