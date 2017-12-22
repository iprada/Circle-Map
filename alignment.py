#!/data/xsh723/anaconda/bin/python3.6
import os
import pysam as ps
import pybedtools as bt
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import time


begin = time.time()
#fastafile = ps.FastaFile("/home/inigo/msc_thesis/genome_data/hg38.fa")

os.chdir("/isdata/kroghgrp/xsh723/circle_map/test_data/aligner")
read_length = 100


os.system("samtools index coordinate_circle_qname_sorted_paired_end_sim_aln.bam")
os.system("bedtools genomecov -bg -ibam coordinate_circle_qname_sorted_paired_end_sim_aln.bam | mergeBed > circ_read_coverage.bed")


def is_soft_clipped(read):
    for cigar in read.cigar:
        if cigar[0] == 4:
            return (True)
        else:
            return (False)

circ_bam = ps.AlignmentFile("coordinate_circle_qname_sorted_paired_end_sim_aln.bam", "rb")

circ_intervals = bt.BedTool("circ_read_coverage.bed")
mapq_filtered = []



number = 0
for interval in circ_intervals:
    high_mapq_support = 0
    for read in circ_bam.fetch(interval.chrom, interval.start, interval.end):
        if high_mapq_support < 1:
            if read.mapq >= 10:
                high_mapq_support +=1
            else:
                continue
    if high_mapq_support > 0:
        number += 1
        interval.append(number)
        mapq_filtered.append(interval)





mapq_filtered = bt.BedTool(mapq_filtered)
mapq_filtered.saveas("mapq_filter.bed")



mapq_filtered = bt.BedTool("mapq_filter.bed")

print(len(mapq_filtered))
print("I am here in the code")
exit()

results = []
i = 0
#loop trough genomecov intervals
for interval in mapq_filtered:
    print(i)
    i +=1
    #posible mate boundaries
    other_boundaries = []
    this_boundary_intervals = []
    for read in circ_bam.fetch(interval.chrom, interval.start, interval.end):
        if read.mapq > 10:
            #possible mate intervals
            mate_interval = [interval.chrom,read.next_reference_start - read_length,read.next_reference_start + read_length ]
            this_boundary_intervals.append(bt.create_interval_from_list(mate_interval))

    #sort_and_merge
    this_boundary_bedtool = bt.BedTool(this_boundary_intervals)

    this_boundary_bedtool = this_boundary_bedtool.sort()
    this_boundary_bedtool = this_boundary_bedtool.merge()

    #look for circle-SV
    #create hit matrix
    this_interval_output = []
    # add 0 values corresponding to the realignment/SA and discordant score
    for mate_interval in this_boundary_bedtool:
        mate_interval.append(0)
        mate_interval.append(0)
        this_interval_output.append(mate_interval)



    realignment_intervals = bt.BedTool(this_interval_output)


    #loop trough the possible mates
    for mate_interval in realignment_intervals:
        #loop trough reads
        for read in circ_bam.fetch(interval.chrom, interval.start, interval.end):

            if read.mapq > 10:
                if is_soft_clipped(read) == True:
                    if read.has_tag('SA'):
                        # if there is SA (split alignment) we do not need to realign
                        suplementary = read.get_tag('SA')
                        # this list will have the following information [chr,left_most start,"strand,CIGAR,mapq, edit_distance]
                        supl_info = [x.strip() for x in suplementary.split(',')]

                        # check the chromosome
                        if int(supl_info[4]) > 10:
                            if supl_info[0] == interval.chrom:
                                # aligned in the left boundary, supplementary on right
                                if mate_interval.start < int(supl_info[1]) < mate_interval.end:

                                    mate_interval[3] = int(mate_interval[3]) +1

                    else:
                        #realignment
                        if read.cigar[0][0] == 4:
                            # get the nucleotides of the beginning
                            nucleotides = read.cigar[0][1]
                            soft_clip_seq = read.seq[0:nucleotides]

                            if 'N' in soft_clip_seq:
                                continue
                            else:

                                if read.is_reverse == True:
                                    # get support interval seq
                                    seq = fastafile.fetch(mate_interval.chrom, mate_interval.start, mate_interval.end)
                                    # reverse it
                                    seq = Seq(seq, generic_dna)
                                    seq = seq.reverse_complement()

                                    # smith watermann
                                    pairwise_al = pairwise2.align.localms(seq.upper(), soft_clip_seq, 1, -1, -.5, -.1)

                                    # if the score is equal to read length
                                    try:
                                        if int(pairwise_al[0][2]) == len(soft_clip_seq):
                                            mate_interval[3] = int(mate_interval[3]) + 1
                                    except:
                                        continue


                                else:
                                    # read mapped to + strand
                                    seq = fastafile.fetch(mate_interval.chrom, mate_interval.start,mate_interval.end)
                                    # smith watermann
                                    pairwise_al = pairwise2.align.localms(seq.upper(), soft_clip_seq, 1, -1, -.5, -.1)

                                    try:
                                        if int(pairwise_al[0][2]) == len(soft_clip_seq):
                                            mate_interval[3] = int(mate_interval[3]) + 1
                                    except:
                                        continue


                        elif read.cigar[-1][0] == 4:
                            # get nucleotides in the end
                            nucleotides = read.cigar[-1][1]
                            soft_clip_seq = read.seq[nucleotides:]

                            if 'N' in soft_clip_seq:
                                continue

                            else:

                                if read.is_reverse == True:
                                    # get support interval seq
                                    seq = fastafile.fetch(mate_interval.chrom, mate_interval.start,
                                                          mate_interval.end)
                                    # reverse it
                                    seq = Seq(seq, generic_dna)
                                    seq = seq.reverse_complement()
                                    pairwise_al = pairwise2.align.localms(seq.upper(), soft_clip_seq, 1, -1, -.5, -.1)
                                    try:
                                        if int(pairwise_al[0][2]) == len(soft_clip_seq):
                                            mate_interval[3] = int(mate_interval[3]) + 1
                                    except:
                                        continue


                                else:
                                    seq = fastafile.fetch(mate_interval.chrom, mate_interval.start,
                                                          mate_interval.end)

                                    pairwise_al = pairwise2.align.localms(seq.upper(), soft_clip_seq, 1, -1, -.5, -.1)

                                    try:
                                        if int(pairwise_al[0][2]) == len(soft_clip_seq):
                                            mate_interval[3] = int(mate_interval[3]) + 1
                                    except:
                                        pass
                else:
                    if mate_interval.start < read.next_reference_start < mate_interval.end:
                        mate_interval[4] = int(mate_interval[4]) + 1

        #output values
        values = [int(mate_interval.start),int(mate_interval.end),int(interval.start),int(interval.end)]
        sorted_values =sorted(values)

        output = [interval.chrom,sorted_values[0],sorted_values[3],mate_interval[3],mate_interval[4]]

        intersection = mapq_filtered.intersect([mate_interval])


        #if there are more than two SV
        if int(mate_interval[3]) + int(mate_interval[4]) >2:
            if len(intersection) > 0:
                output.append(interval[3])
                output.append(intersection[0][3])
                results.append(output)



            print(output)




results = bt.BedTool(results)
results.saveas("results.bed")












os.system("cat results.bed  | awk '$7!=$6' | awk -v OFS='\t' '{if ($7>$6) {print $0} else {print $1,$2,$3,$4,$5,$7,$6}}' | sortBed | bedtools groupby -g 1,6,7 -c 2,3,4,5 -o min,max,sum,sum | awk -v OFS='\t' '{print $1,$4,$5,$6,$7}' > possible_results.bed")





end = time.time()
time = end -begin
print(time)


