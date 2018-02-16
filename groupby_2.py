import pybedtools as bt
import pandas as pd
import os

os.chdir("/isdata/kroghgrp/xsh723/projects/circle_map/test_data/realigner/test_groupby/")

print(os.listdir())

my_myoutput = bt.BedTool("circle_map_results_working_branch.bed")

#print(my_myoutput)

def second_merge(start1,start2,end1,end2,chrom1,chrom2):

    chrom_match = (chrom1 == chrom2) * 1

    distance = (abs(start1 - start2) + abs(end1 - end2))

    one_overlap_two = 1 - (distance / (end1 - start1))
    two_overlap_one = 1 - (distance / (end2 - start2))
    return (one_overlap_two + two_overlap_one +chrom_match)

my_myoutput_pd = my_myoutput.to_dataframe(names=['chrom','start','end','sc','dr'])

my_myoutput_pd = my_myoutput_pd.sort_values(by=['chrom','start','end'])

print(my_myoutput_pd.head(22))

print(second_merge(my_myoutput_pd.start,my_myoutput_pd.start.shift(),
                                                    my_myoutput_pd.end,my_myoutput_pd.end.shift(),
                                                    my_myoutput_pd.chrom,my_myoutput_pd.chrom.shift()).lt(2.99).cumsum().head(22))

print(second_merge(my_myoutput_pd.start,my_myoutput_pd.start.shift(),
                                                    my_myoutput_pd.end,my_myoutput_pd.end.shift(),
                                                    my_myoutput_pd.chrom,my_myoutput_pd.chrom.shift()).lt(2.99).head(22))

print(second_merge(my_myoutput_pd.start,my_myoutput_pd.start.shift(),
                                                    my_myoutput_pd.end,my_myoutput_pd.end.shift(),
                                                    my_myoutput_pd.chrom,my_myoutput_pd.chrom.shift()).head(22))

#print(my_myoutput_pd)

final_output = my_myoutput_pd .groupby(second_merge(my_myoutput_pd.start,my_myoutput_pd.start.shift(),
                                                    my_myoutput_pd.end,my_myoutput_pd.end.shift(),
                                                    my_myoutput_pd.chrom,my_myoutput_pd.chrom.shift()).lt(2.9999).cumsum()).agg({'chrom':'first','start':'first','end':'last','sc':'sum','dr':'sum'})



bedtool_output = bt.BedTool.from_dataframe(final_output)

bedtool_output.saveas("final_output_working_branch.bed")

