import os
import pandas as pd
import pybedtools as bt
import numpy as np

def iteration_merge(only_discordants,results):
    """a"""


    discordant_bed = bt.BedTool(only_discordants)
    unparsed_bed = bt.BedTool(results)


    unparsed_pd = unparsed_bed.to_dataframe(
        names=['chrom', 'start', 'end', 'read', 'iteration', 'discordants'])

    unparsed_pd = unparsed_pd.sort_values(['iteration','chrom','start'])


    grouped = unparsed_pd.groupby(merge_fraction(unparsed_pd.iteration, unparsed_pd.start,
                                           unparsed_pd.end, unparsed_pd.iteration.shift(),
                                           unparsed_pd.start.shift(),
                                           unparsed_pd.end.shift()).lt(2.98).cumsum()).agg(
        {'chrom': 'first', 'start': np.min, 'end': np.max, 'discordants': 'max', 'read': 'nunique'})

    bedtool_output = bt.BedTool.from_dataframe(grouped)
    print(bedtool_output)


    bed_to_write = bedtool_output.cat(discordant_bed, postmerge=False)

    return(bed_to_write)

def merge_fraction(chrom1,x1,x2,chrom2,y1,y2):
    """compute overlap of the interval y over interval x"""

    distance = (np.minimum(x2.values,y2.values) - np.maximum(x1.values,y1.values))




    one_overlap_two = distance/(y2.values-y1.values)

    two_overlap_one = distance/(x2.values-x1.values)


    # check if they are on the same chromosome and the amount of overlap if so
    return(pd.Series(chrom1 == chrom2) + pd.Series(two_overlap_one.clip(0)) + pd.Series(one_overlap_two.clip(0)))


os.chdir('/home/iprada/faststorage/projects/realigner')
bed = bt.BedTool('testing_p.bed')

print(bed.sort())

bed_pd = bed.to_dataframe(names=['chrom','start','end','read','iteration','discordants'])
bed_pd = bed_pd.sort_values(['iteration','chrom','start']).reset_index()
print(bed_pd)

grouped = bed_pd.groupby(merge_fraction(bed_pd.iteration.shift(), bed_pd.start.shift(),
                                           bed_pd.end.shift(), bed_pd.iteration,
                                           bed_pd.start,
                                           bed_pd.end).lt(2.98).cumsum()).agg(
        {'chrom': 'first', 'start': 'min', 'end': 'max', 'discordants': 'max', 'read': 'nunique'}).reset_index()

grouped = grouped.sort_values(['chrom','start','end']).reset_index()




grouped_2 = grouped.groupby(merge_fraction(grouped.chrom.shift(), grouped.start.shift(),
                                           grouped.end.shift(), grouped.chrom,
                                           grouped.start,
                                           grouped.end).lt(2.98).cumsum()).agg(
        {'chrom': 'first', 'start': 'min', 'end': 'max', 'discordants': 'max', 'read': 'sum'}).sort_values(['chrom','start','end'])

print(grouped)

print(merge_fraction(grouped.chrom.shift(), grouped.start.shift(),
                                           grouped.end.shift(), grouped.chrom,
                                           grouped.start,
                                           grouped.end).lt(2.98).cumsum())

print(grouped_2)





