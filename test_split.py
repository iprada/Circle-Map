import pybedtools as bt
import os

os.chdir("/isdata/kroghgrp/xsh723/projects/circle_map/test_data/realigner/test_command_interfance")

print(os.listdir())

my_bed = bt.BedTool("circle_map_results_optimized_solved_groupby.bed")

#print(my_bed)


for i in range(0,100):
    intervals = []
    for interval in my_bed:
        intervals.append(interval)


