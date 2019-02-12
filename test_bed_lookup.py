import pybedtools as bt

file = bt.BedTool("/home/iprada/small_size_projects/circle_map_paper/raw_data/gap_regions_hg38_01_02_19.bed").sort()
print(file)



for i in range(0,10000):
    print(file.any_hits(bt.Interval('chr1',1000,2000)))