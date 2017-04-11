import pysam as ps
import os
import matplotlib.pyplot as plt
os.chdir("/home/inigo/msc_thesis")



#sort the bam file by query name, the reads will be sorted in the order they were generated
#os.system("samtools sort -n -o query_name_sorted_paired_end_sim_aln.bam.sv.bam sorted_paired_end_sim_aln.bam.sv.bam")


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

first_read = uniq_soft_cliped_bam.head(1)





previous_read = None

for read in first_read:
    previous_read = read.query_name


for read in soft_clipped:
    if previous_read == read.query_name:
        continue
    else:
        uniq_soft_cliped_bam.write(read)

    previous_read = read.query_name


uniq_soft_cliped_bam.close()


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


print(len(soft_clip_len))

plt.hist(soft_clip_len,bins=70)
plt.ylabel("Frequency")
plt.xlabel("Soft clip length")
plt.show()


