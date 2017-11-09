#!/data/xsh723/anaconda/bin/python3.6
from extract_circle_SV_reads import readExtractor

create_object = readExtractor("qname_sorted_paired_end_sim_aln.bam","circle_reads.bam","/isdata/kroghgrp/xsh723/circle_map/test_data/read_extractor/Read_Extractor/")
create_object.extract_sv_circleReads()
