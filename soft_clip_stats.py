import pysam as ps
import os
import matplotlib.pyplot as plt
os.chdir("/home/inigo/msc_thesis/mapping_stats/")



#sort the bam file by query name, the reads will be sorted in the order they were generated
os.system("samtools sort -n -o query_name_sorted_paired_end_sim_aln.bam.sv.bam sorted_paired_end_sim_aln.bam.sv.bam")


#open the bam file
samfile = ps.Samfile("query_name_sorted_paired_end_sim_aln.bam.sv.bam", "rb")


#file name of all soft-clipped reads
soft_cliped_bam = ps.AlignmentFile("soft_clipped.bam","wb",template=samfile)



for read in samfile:

    #append the cigar codes of each reads, the soft clipped reads have cigar code = 4
    cigar_codes = []
    # iterate trough cigar
    for cigar in read.cigar:
        cigar_codes.append(cigar[0])

    if (4 in cigar_codes) == True:
        #if soft-clipped, add to the soft-clipped file
        soft_cliped_bam.write(read)
    else:
        pass


soft_cliped_bam.close()

#open the soft-clipped bam to iterate
soft_clipped = ps.Samfile("soft_clipped.bam", "rb")

#file of the unique soft-clipped, meaning that they only have mapped to one part of the breakpoint
uniq_soft_cliped_bam = ps.AlignmentFile("unique_soft_clipped.bam","wb",template=samfile)


#build an index so that one can find the how many times in each read present

index = ps.IndexedReads(soft_clipped)
index.build()





#get query names of the reads that are unique
os.system("samtools view soft_clipped.bam | awk '{print $1}' | uniq -c | awk '{if ($1==1) {print $2}}' > uniq_soft_clip.txt")

soft_clip_queries = open("uniq_soft_clip.txt","r")
soft_clip_queries = soft_clip_queries.readlines()



for query in soft_clip_queries:
    iterator = index.find(query.strip())
    for read in iterator:
        uniq_soft_cliped_bam.write(read)


uniq_soft_cliped_bam.close()

#remove intermediate files
os.system("rm uniq_soft_clip.txt")
os.system("rm query_name_sorted_paired_end_sim_aln.bam.sv.bam")
os.system("rm soft_clipped.bam")



uniq_soft_cliped_bam = ps.AlignmentFile("unique_soft_clipped.bam","rb")


soft_clip_len = []

for read in uniq_soft_cliped_bam:
    if ("N" in read.seq) == True:
        continue
    else:


        for cigar in read.cigar:
            if cigar[0] == 4:
                soft_clip_len.append(cigar[1])
            else:
                continue


uniq_soft_cliped_bam.close()


plt.hist(soft_clip_len,bins=30)
plt.ylabel("Frequency")
plt.xlabel("Soft clip length")
plt.title("%s single soft-clip" % len(soft_clip_queries))
plt.savefig('soft_clip_len.png')
plt.close()

