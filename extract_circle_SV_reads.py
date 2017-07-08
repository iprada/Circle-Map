import pysam
from circ_dna.alignment import alignment
import os
class readExtractor:
    def __init__(self,sorted_bam,output_bam,working_dir):
        self.sorted_bam = sorted_bam
        self.output_bam = output_bam
        self.no_plasmids_bam = "without_plasmids_%s" % str(sorted_bam)
        self.working_dir = working_dir

    def is_soft_clipped(self, read):
        for cigar in read.cigar:
            if cigar[0] == 4:
                return (True)
            else:
                return (False)

    def extract_sv_circles(self):
        os.chdir("%s" % self.working_dir)
        """Function that extract R2F1 discordant reads and soft-clipped reads"""
        raw_bam = pysam.AlignmentFile(self.working_dir + self.sorted_bam, "rb")
        without_plasmids = pysam.AlignmentFile(self.working_dir + self.no_plasmids_bam,"wb", template=raw_bam)
        circle_reads = pysam.AlignmentFile(self.working_dir + self.output_bam, "wb", template=raw_bam)

        print("all files loaded")

        print("removing plasmids")
        #remove reads aligned to plasmids
        for read in raw_bam:
            try:
                contig = raw_bam.get_reference_name(read.reference_id)
                if contig[0:3] == "chr":
                    without_plasmids.write(read)
            except:
                pass

        print("plasmids removed")

        without_plasmids.close()

        #index the without plasmids file

        os.system("samtools index %s" % self.no_plasmids_bam)

        #reopen the without_plasmids file as read to extract SV reads
        without_plasmids = pysam.AlignmentFile(self.working_dir + self.no_plasmids_bam, "rb")


        print("getting SV reads")

        i = 0
        for read in without_plasmids:
            i+=1
            print(i)
            #only look at read1, (speed up), this part looks for discordants
            if read.is_paired and read.is_read1:
                #look for R1F2 orientation
                if read.mate_is_reverse and read.is_reverse == False:
                    #R2 start must be smaller than F1, meaning the following:
                    # <-(R2) ->(F1) a circle :)
                    if read.reference_start > read.next_reference_start:
                        #only extract the discordants (<- ->) if the mate is present
                        try:
                            mate = without_plasmids.mate(read)
                            circle_reads.write(mate)
                            circle_reads.write(read)

                        except:
                            continue


            else:
                #this part will extract soft-clipped reads
                if self.is_soft_clipped(read) == True:
                    circle_reads.write(read)
                else:
                    try:
                        mate = without_plasmids.mate(read)
                        if self.is_soft_clipped(mate) == True:
                            circle_reads.write(mate)
                    except:
                        continue
        circle_reads.close()
        print("finished extracting reads")











