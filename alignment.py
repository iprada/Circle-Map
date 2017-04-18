import pysam as ps
import os
import pandas as pd
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
        os.system("cat %s | awk -v OFS='\t' '{ if ($4 > 2) { print $1,$2,$3} }' | mergeBed > %s" % (self.circ_bedgraph,self.circ_boundaries))
        circ_boundaries = pd.read_table("%s" % self.circ_boundaries,sep='\t')

        return(circ_boundaries)





    def call_circles(self):
        """Function that aims to call circles based on the boundaries identified using the discordant reads and soft clipped reads"""
        print("Hunting circles")
        circ_boundaries = alignment.get_circ_dna_boundaries(self)
        circ_boundaries.columns = ['chr','start','end']
        #parameters to run in parallel
        number_of_rows = len(circ_boundaries.index)
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

        for index, row in circ_boundaries.iterrows():

            chromosome = row['chr']
            start = row['start']

            end = row['end']


            os.system("samtools index %s" % self.circ_bam)
            circ_file = ps.AlignmentFile(self.circ_bam,"rb")
            boundary_file = ps.AlignmentFile("boundary_file.bam", "wb", template=circ_file)

            mate_file = ps.AlignmentFile("mates.bam", "wb", template=circ_file)

            # loop trough reads in the boundary
            for read in circ_file.fetch('chr1',95047748,95135811):
                try:
                    print("pair")
                    print(circ_file.mate(read))
                    print(read)
                except:
                    print("alone")
                    print(read)







            break
























