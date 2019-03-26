import pysam as ps

circle = ["chrV",398445,409306,177,605]

bam = ps.AlignmentFile("/isdata/kroghgrp/xsh723/projects/6_aged_yeast/working_directory/aligned/B05/coordinate_circle_qname_sortB05.bam")
array= bam.count_coverage(contig='chrV',start=409306-10,stop=409306,read_callback='nofilter')
print(array)