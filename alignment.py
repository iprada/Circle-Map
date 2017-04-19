import pysam as ps
import os
import pandas as pd
import pybedtools as bt
import multiprocessing as mp





class alignment:

    def __init__(self,bam,working_dir):
        self.bam = bam
        self.circ_bam = "circ_calls.bam"
        self.circ_bedgraph = "circ_dna_boundaries.bedGraph"
        self.circ_boundaries = "circ_dna_boundaries.bed"
        self.working_dir = working_dir
        self.all_bam = "sorted_paired_end_sim_aln.bam"

        self.number_of_cores = 3


    def remove_concordant_pairs(self):
        """Function that removes concordant pairs"""
        os.system("samtools view -hb -F 2 %s > %s" % (self.bam,self.circ_bam))




    def get_circ_dna_boundaries(self):
        """function that computes circDNA boundaries from a bamfile with soft-clipped and split reads"""
        print("computing bedGraph")

        alignment.remove_concordant_pairs(self)
        os.system("bedtools genomecov -ibam %s -bg  > %s" % (self.circ_bam,self.circ_bedgraph))
        print("Removing lonely reads")
        os.system("cat %s | awk -v OFS='\t' '{ if ($4 > 5) { print $1,$2,$3} }' | mergeBed | awk -v OFS='\t' '{ print $1,$2,$3,'0'}' > %s" % (self.circ_bedgraph,self.circ_boundaries))
        circ_boundaries= bt.BedTool(self.circ_boundaries)
        return(circ_boundaries)

    def query_name_sorted_circs(self):
        os.system("samtools sort -n -o query_name_sorted_circs.bam circ_calls.bam")
        os.system("samtools index circ_calls.bam ")
        sorted_circs = ps.AlignmentFile("query_name_sorted_circs.bam","rb")
        return(sorted_circs)

    def is_soft_clipped(self,read):
        for cigar in read.cigar:
            if cigar[0] == 4:
                return(True)
            else:
                return(False)





    def call_circles(self):
        """Function that aims to call circles based on the boundaries identified using the discordant reads and soft clipped reads"""
        print("Hunting circles")
        circ_boundaries = alignment.get_circ_dna_boundaries(self)
        circ_boundaries.columns = ['chr','start','end']
        #parameters to run in parallel
        number_of_rows = len(circ_boundaries)
        l = range(1,number_of_rows)
        medium_process = number_of_rows / self.number_of_cores
        medium_process = (int(medium_process))
        n = medium_process

        temp_file_names = []
        for i in range(0, self.number_of_cores):
            temp_file = "circles_" + str(i) + ".bed"
            temp_file_names.append(temp_file)

        spawn_reads = [l[i:i + n] for i in [0, len(l), n]]

        # iterate trough each entry in the bed file


        sorted_circs = alignment.query_name_sorted_circs(self)
        #index by read name

        index = ps.IndexedReads(sorted_circs)
        index.build()
        os.system("samtools index circ_calls.bam")
        circs = ps.AlignmentFile("circ_calls.bam","rb")

        i = 0

        for interval in circ_boundaries:
            print(i)
            i += 1
            if int(interval.name) == 0:


                mate_file = ps.AlignmentFile("mates.bam", "wb", template=circs)

                soft_clipped = ps.AlignmentFile("soft_clipped.bam","wb",template=circs)

                soft_clipped_mates = ps.AlignmentFile("soft_clipped_mates.bam","wb",template=circs)

                # loop trough reads in the boundary
                for read in circs.fetch(interval.chrom,interval.start,interval.end):

                    # check if it is soft-clipped
                    if alignment.is_soft_clipped(self,read) == True:
                        soft_clipped.write(read)
                        try:

                            mates = index.find(read.query_name)

                            for mate in mates:

                                if mate == read:
                                    continue
                                else:
                                    if alignment.is_soft_clipped(self, mate) == True:
                                        soft_clipped_mates.write(mate)

                                    mate_file.write(mate)


                        except:
                            continue
                    else:

                        try:
                            mates = index.find(read.query_name)

                            for mate in mates:
                                if mate == read:
                                    continue
                                else:
                                    if alignment.is_soft_clipped(self, mate) == True:
                                        soft_clipped_mates.write(mate)

                                    mate_file.write(mate)
                        except:
                            continue




                mate_file.close()
                soft_clipped.close()
                soft_clipped_mates.close()

                # get the mate boundaries
                mate_bam = bt.BedTool('mates.bam')
                os.system("samtools view mates.bam | awk '{print $1,$6}'")
                mate_bam.bam_to_bed(output='mate.bed')

                #os.system("cat mate.bed | uniq > mate.bed")

                mate_bed = bt.BedTool('mate.bed')
                mate_bed = mate_bed.sort()
                mate_bed = mate_bed.merge()

                exit()



                print(mate_bed)






            else:
                continue





            os.system("rm mates.bam")
            os.system("rm soft_clipped.bam")
            os.system("rm soft_clipped_mates.bam")




































