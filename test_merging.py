import pandas as pd
import os
import pybedtools as bt

os.chdir('/home/iprada/faststorage/projects/realigner')

bed = bt.BedTool('testing.bed')



pandas_bed = bed.to_dataframe(names=['chrom','start','end','discordant','splits'])


print(pandas_bed)



#merged = discordant_bed.groupby(
#                                (discordant_bed.end - discordant_bed.start.shift()-1).lt(
#                                    0).cumsum()).agg(
#                                {'chrom': 'first', 'start': 'first', 'end': 'last', 'read': 'count'})


def merge_function(pandas_data_frame):
    overlap = ((pandas_data_frame.start-pandas_data_frame.shift().end)-1).lt(0)
    chr_overlap = (pandas_data_frame.chrom == pandas_data_frame.shift().chrom)
    return((overlap*1 + chr_overlap*1).lt(2).cumsum())




print(pandas_bed.groupby(merge_function(pandas_bed)).agg({'chrom':'first','start':'first','end':'last',
                                                          'discordant':'max','splits':'sum'}))

