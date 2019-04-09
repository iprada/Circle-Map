from realigner import realignment
import multiprocessing as mp
import pandas as pd
import pybedtools as bt

create_object = realignment(input_bam="/home/iprada/faststorage/projects/6_aged_yeast/working_directory/aligned/BM3/sort_circle_qname_BM3.bam",qname_bam="qname_BM3.bam",sorted_bam="/home/iprada/faststorage/projects/6_aged_yeast/working_directory/aligned/BM3/sorted_BM3.bam",
genome_fasta="/home/iprada/faststorage/reference_Data/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/8_plasmids_genome/yeast_8_plasmids.fa",
                            directory="/home/iprada/faststorage/projects/6_aged_yeast/working_directory/aligned/BM3/",
                            mapq_cutoff=10,insert_size_mapq = 60,std_extension=4,insert_size_sample_size=100000,ncores=32,
                            gap_open=5,gap_ext=1,n_hits=10,prob_cutoff=0.01,min_soft_clipped_length=8,
                            overlap_frac = 0.95,interval_p_cut=0,af=0.1,edit_distance_frac=0.05,insert_size=[300,25],locker=mp.Lock(),
                            only_discordants=True,pid=200,ratio=0,remap_splits=False,score=8,splits=2,split=0,output_name="test",verbose=1,discordant_filter=3)
#print(pd.read_table("/home/iprada/faststorage/projects/6_aged_yeast/working_directory/aligned/BM3/merged.txt",names=["chrom","start","end"]))
parsed = []

for interval in bt.BedTool("/home/iprada/faststorage/projects/6_aged_yeast/working_directory/aligned/BM3/merged.txt"):
    parsed.append([interval.chrom,interval.start,interval.end])
create_object.realign(parsed)