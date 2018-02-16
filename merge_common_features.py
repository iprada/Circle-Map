import pybedtools as bt
import os
import time
import pandas as pd


os.chdir("/isdata/kroghgrp/xsh723/projects/circle_map/test_data/realigner/test_groupby/")
begin = time.time()

def fraction(start1,start2,end1,end2,read1,read2,chrom1,chrom2):

    read_match = (read1 == read2)*1
    chrom_match = (chrom1 == chrom2) * 1



    distance = (abs(start1-start2) + abs(end1-end2))

    one_overlap_two = 1 - (distance/(end1-start1))
    two_overlap_one =  1 - (distance/(end2-start1))




    return(one_overlap_two + two_overlap_one + read_match +chrom_match)


def second_merge(start1,start2,end1,end2,chrom1,chrom2):

    chrom_match = (chrom1 == chrom2) * 1

    distance = (abs(start1 - start2) + abs(end1 - end2))

    one_overlap_two = 1 - (distance / (end1 - start1))
    two_overlap_one = 1 - (distance / (end2 - start1))
    return (one_overlap_two + two_overlap_one +chrom_match)













names = ['chrom', 'start', 'end','probability','read']

my_list = [['chr1',1,10,0.99,'read1'],
            ['chr1',1,12,0.99,'read2'],
            ['chr1',1,12,0.99,'read2']
            ]
my_data_pd = pd.DataFrame.from_records(my_list, columns=names)

#print(my_data_pd)

#grouped = my_data_pd.groupby((my_data_pd.end.shift()-my_data_pd.start).lt(0).cumsum()).agg({'chrom':'first','start':'first','end':'last','probability':'sum','read':'nunique'})

#print(grouped)



#print(fraction(my_data_pd.start,my_data_pd.start.shift(),my_data_pd.end,my_data_pd.end.shift(),my_data_pd.read,my_data_pd.read.shift()).gt(2.98).cumsum())


simulated_data_bed = bt.BedTool("unparsed.bed")


simulated_data_pd = simulated_data_bed.to_dataframe(names=['chrom','start','end','read','iteration','discordants'])

simulated_sorted = simulated_data_pd.sort_values(by=['iteration','chrom','start','end'])

#print(simulated_sorted)
grouped = simulated_sorted.groupby(fraction(simulated_sorted.start,simulated_sorted.start.shift(),
                                  simulated_sorted.end,simulated_data_pd.end.shift(),
                                  simulated_sorted.iteration,
                                  simulated_sorted.iteration.shift(),simulated_sorted.chrom,simulated_sorted.chrom.shift()).lt(3.98).cumsum()).agg({'chrom':'first','start':'first','end':'last','discordants':'max','read':'nunique'})

print(grouped)


exit()
bedtool_output = bt.BedTool.from_dataframe(grouped)
bedtool_output.saveas("final_output_opt_my_test.bed")
exit()
second_merging_round = grouped.sort_values(by=['chrom','start','end'])

final_output = second_merging_round.groupby(second_merge(second_merging_round.start,second_merging_round.start.shift(),
                                                     second_merging_round.end,second_merging_round.end.shift(),
                                                         second_merging_round.chrom,second_merging_round.chrom.shift()).lt(2.99).cumsum()).agg({'chrom':'first','start':'first','end':'last','read':'sum','discordants':'sum'})



bedtool_output = bt.BedTool.from_dataframe(final_output)

bedtool_output.saveas("final_output_opt_my_test.bed")

