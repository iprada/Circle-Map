import pybedtools as bt
import os
import pandas as pd

os.chdir("/isdata/kroghgrp/xsh723/projects/circle_map/test_data/realigner/test_output_parsing")

my_bed = bt.BedTool("circle_map_results_test_speed.bed")

bed_dataframe = my_bed.to_dataframe(names=['chrom','start','end','read','read_type','iteration'])

bed_dataframe = bed_dataframe.sort_values(by=['iteration','chrom', 'start', 'end'],ascending=[True,True, True, True])

print(bed_dataframe)

sorted_grouped = bed_dataframe.groupby(['chrom', 'start', 'end','iteration']).agg({'read':'nunique'}).reset_index()

print(sorted_grouped)