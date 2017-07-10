import pysam
import os

class readExtractor:
    def __init__(self,sorted_bam,output_bam,working_dir):
        self.sorted_bam = sorted_bam
        self.output_bam = output_bam
        self.no_plasmids_bam = "without_plasmids_%s" % str(sorted_bam)
        self.working_dir = working_dir
        self.qname_sorted = "qname_sorted_%s" % str(self.no_plasmids_bam)





    def extract_sv_circles(self):
        """Function that extracts Structural Variant reads that indicate circular DNA,
        The programme with extract soft-clipped reads and R2F1 (<- ->) oriented reads"""

        def is_soft_clipped(read):
            for cigar in read.cigar:
                if cigar[0] == 4:
                    return (True)
                else:
                    return (False)


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

        os.system("samtools sort -n -o %s %s" % (self.qname_sorted,self.no_plasmids_bam))

        #reopen the without_plasmids file as read to extract SV reads
        without_plasmids = pysam.AlignmentFile(self.working_dir + self.qname_sorted, "rb")


        print("getting SV reads")

        #cache read1
        read1=''
        for read in without_plasmids:
            #check that is read 1 and cache it
            if read.is_read1:
                read1=read


            else:
                #now, assert that reads are pairs
                if read.is_read2 and read1.query_name == read.query_name:
                    read2 = read
                    #assert that it is paired
                    if read1.is_paired and read1.is_read1:
                        #assert that mate is R and read1 is F
                        if read1.mate_is_reverse and read1.is_reverse == False:
                            # check that R2 reference start is < than F1
                            if read1.reference_start > read1.next_reference_start:
                                #only save then if they are mapped to the same chromosome
                                if read1.reference_id == read2.reference_id:
                                    circle_reads.write(read1)
                                    circle_reads.write(read2)
                    else:
                        # if they are not paired, check it r1 and r2 are soft clipped
                        if is_soft_clipped(read1) == True:
                            circle_reads.write(read1)
                        if is_soft_clipped(read):
                            circle_reads.write(read)



        circle_reads.close()
        print("finished extracting reads")











