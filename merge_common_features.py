import pybedtools as bt
import os
import time

begin = time.time()



os.chdir("/isdata/kroghgrp/xsh723/circle_map/test_data/old_vs_new/")

my_bed = bt.BedTool("circle_map_results_optimized_2.bed")

sorted = my_bed.sort()


previous_interval = []
for interval in my_bed:

    #previous elements
    previous_intervals = bt.BedTool(previous_interval)


    #interval iterating on
    current_interval = [interval]
    current = bt.BedTool(current_interval)

    if current == previous_intervals:
        previous_interval = current_interval

    else:
        intersection = previous_intervals.intersect(current,wa=True,r=True,f=0.99)
        previous_interval = current_interval
        if len(intersection) > 0:
            print(intersection)




end = time.time()

print((end-begin)/60)