import pandas as pd
import pybedtools as bt
import os
import numpy as np

os.chdir("/isdata/kroghgrp/xsh723/projects/6_aged_yeast/working_directory/aligned")

bed=bt.BedTool("test_B02.bed")

results = []
for interval in bed:
    results.append([interval.chrom,int(interval.start),int(interval.end),interval[3]])



def merge_fraction(chrom1,x1,x2,chrom2,y1,y2):
    """compute overlap of the interval y over interval x"""

    distance = (np.minimum(x2.values,y2.values) - np.maximum(x1.values,y1.values))



    one_overlap_two = distance/(y2.values-y1.values)

    two_overlap_one = distance/(x2.values-x1.values)


    # check if they are on the same chromosome and the amount of overlap if so
    return(pd.Series(chrom1 == chrom2) + pd.Series(two_overlap_one.clip(0)) + pd.Series(one_overlap_two.clip(0)))

frac=0.8
fraction = (frac*2)+1

unparsed_pd = pd.DataFrame.from_records(results,columns=['chrom', 'start', 'end','item'])



sort = unparsed_pd.sort_values(by=['chrom', 'start','end'],ascending=True).reset_index(drop=True)

merging_out  = sort.groupby(
        merge_fraction(sort.chrom, sort.start,
                     sort.end,sort.chrom.shift(),sort.start.shift(),sort.end.shift()).lt(fraction).cumsum()).agg(
        {'chrom': 'first', 'start': 'min', 'end': 'max','item': 'sum'})


merging_out = merging_out.sort_values(by=['chrom', 'start', 'end']).reset_index(drop=True)

bedtool_output = bt.BedTool.from_dataframe(merging_out)
print(bedtool_output)
exit()

final_output = merging_out.groupby(
        merge_fraction(merging_out.chrom, merging_out.start,
                       merging_out.end, merging_out.chrom.shift(), merging_out.start.shift(),merging_out.end.shift()).lt(fraction).cumsum()).agg(
        {'chrom': 'first', 'start': 'first', 'end': 'last', 'item': 'sum'})

print(final_output)
bedtool_output = bt.BedTool.from_dataframe(final_output)


#print(bedtool_output)