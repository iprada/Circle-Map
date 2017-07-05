import pysam
import time

raw_bam = pysam.AlignmentFile("/home/inigo/test_variant/sp1_query_name_sorted.bam", "rb")
#without_plasmids = pysam.AlignmentFile("/home/inigo/test_variant/no_plasmids.bam", "wb", template=raw_bam)
circle_reads = pysam.AlignmentFile("/home/inigo/test_variant/circle_reads_v2.bam", "wb", template=raw_bam)


def is_soft_clipped(read):
    for cigar in read.cigar:
        if cigar[0] == 4:
            return (True)
        else:
            return (False)

#remove reads mapped to plasmids

#for read in raw_bam:
#    try:
#        contig =raw_bam.get_reference_name(read.reference_id)
#        if contig[0:3] == "chr":
#            without_plasmids.write(read)
#    except:
#        pass

#without_plasmids.close()


#without_plasmids = pysam.AlignmentFile("/home/inigo/test_variant/no_plasmids.bam", "rb")
without_plasmids = pysam.AlignmentFile("/home/inigo/test_variant/sorted_no_plasmids.bam", "rb")



start = time.time()
i = 0




for read in without_plasmids.fetch('chr10',102032827,102048877):
    print(i)
    i +=1
    if read.is_paired and read.is_read1:
        if read.mate_is_reverse and read.is_reverse == False:
            if read.reference_start > read.next_reference_start:
                try:
                    mate = without_plasmids.mate(read)
                    circle_reads.write(mate)
                    circle_reads.write(read)
                    print("done perfect")

                except:
                    print("unable to parse and get mate")
                    continue


    else:
        if is_soft_clipped(read) == True:
            circle_reads.write(read)
        else:
            try:
                mate = without_plasmids.mate(read)
                if is_soft_clipped(mate) == True:
                    circle_reads.write(mate)
            except:
                continue
end = time.time()
print(end-start)
print(i)
