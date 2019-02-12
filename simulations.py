import numpy as np
import os
import pysam as ps
from Bio import SeqIO
from io import StringIO
import random as rd
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import sys
import pybedtools as bt
import time
import subprocess as sp

def sim_ecc_reads(genome_fasta,read_length,directory,reads,exclude_regions,fastq,insert_size,errors,mean_cov,locker,process,sim_circles,paired_end_fastq_1,paired_end_fastq_2,skipped,correct):
    """Function that takes as arguments a genome fasta file, weights each chromosome based on the length
    and simulates single end eccDNA reads
    """


    # Get the length of the chromosomes and store them in a sequence dictionary
    chromosomes = {}
    whole_genome_len = 0
    for rec in SeqIO.parse(genome_fasta, 'fasta'):
        name = rec.id
        seqLen = len(rec)
        whole_genome_len += seqLen
        chromosomes[name] = seqLen
    #chromosome sampling probability weighted based on its length
    weighted_chromosomes = {}
    for contigs in chromosomes:
        weighted_chromosomes[contigs] = chromosomes[contigs]/whole_genome_len


    contig_list = []
    weights = []
    for contigs, value in weighted_chromosomes.items():
        weights.append(value)
        contig_list.append(contigs)

    #Simulate the reads:



    os.chdir(directory)
    circle_bed = []




    set_of_reads = []
    set_of_left_reads = []
    set_of_right_reads = []

    circle_number = 0
    n_of_reads = 0

    begin = time.time()
    #simulated reads
    while n_of_reads < reads + 1:


        #sample weighted chromosome
        #set random seed, important for paralell
        np.random.seed()
        chr = np.random.choice(contig_list, p=weights)

        # decide ecDNA length

        #sample circle length
        circle_length = rd.randint(150,100000)


        # linear decrease in coverage based on circle length





        # compute circles sequencing coverage

        rounds_of_sim = (circle_length * mean_cov)/(read_length*2)







        # take in to account short length contigs
        #start position can't be bigger than (chr_length-circle_length)
        chr_pos_start =  rd.randint(0,(chromosomes[chr] - circle_length))
        #set end
        if chromosomes[chr] == (chromosomes[chr] - circle_length):
            chr_pos_end = chromosomes[chr]
        else:
            chr_pos_end = chr_pos_start + circle_length
            
        #if user of provides regions to exclude, check within it is on the region. and skip it
        if exclude_regions != None and bt.BedTool(exclude_regions).sort().any_hits(bt.Interval(chr,chr_pos_start,chr_pos_end)) != 0:
            #hit in a gap region
            # shared memory object between processes. It is use to track the number of skipped circles
            with skipped.get_lock():
                skipped.value+=1
            continue
        else:
            #shared memory object between processes. It is use to track the number of correctly simulated circles
            with correct.get_lock():
                correct.value+=1
        #save each circle positions, so that then I can check true circles


            first_line = [chr, chr_pos_start, chr_pos_end]

            #create class object outside the loop
            new_read = sim_paired_end(n_of_reads, insert_size, genome_fasta, chr, chr_pos_start,
                                      chr_pos_end, read_length, circle_number,process)

            #simulation rounds
            for each_sim in range(0,round(int(rounds_of_sim))):


                if errors == True:



                    if ((n_of_reads + 1) / 10000).is_integer() == False:
                        try:

                            # sim the read
                            get_seq = new_read.simulate_read()
                            # put it in fastq format
                            simulated_reads = sim_paired_end.simulate_read_with_errors(new_read, get_seq[0], get_seq[1],
                                                                                   get_seq[2])
                            # save the read
                            set_of_left_reads.append(simulated_reads[0])
                            set_of_right_reads.append(simulated_reads[1])
                            n_of_reads += 1
                            print(n_of_reads)

                        except BaseException as e:
                            print("I catched the following exception:",e)
                            pass

                    else:
                        try:
                            # simulate reads and save to disk
                            get_seq = new_read.simulate_read()
                            simulated_reads = sim_paired_end.simulate_read_with_errors(new_read, get_seq[0], get_seq[1],
                                                                                   get_seq[2])
                            set_of_left_reads.append(simulated_reads[0])
                            set_of_right_reads.append(simulated_reads[1])

                            # save to disk
                            locker.acquire()
                            print("Process %s: writting to disk 10000 reads" % process )
                            SeqIO.write(set_of_left_reads, paired_end_fastq_1, "fastq")
                            SeqIO.write(set_of_right_reads, paired_end_fastq_2, "fastq")
                            locker.release()

                            n_of_reads += 1
                            # sim the first read of the list
                            new_read = sim_paired_end(n_of_reads, insert_size, genome_fasta, chr, chr_pos_start,
                                                      chr_pos_end, read_length, circle_number,process)
                            get_seq = new_read.simulate_read()
                            simulated_reads = sim_paired_end.simulate_read_with_errors(new_read, get_seq[0], get_seq[1],
                                                                                   get_seq[2])

                            set_of_left_reads = [simulated_reads[0]]
                            set_of_right_reads = [simulated_reads[1]]
                            n_of_reads += 1
                        except BaseException as e:
                            print("I catched the following exception:",e)
                            pass

                else:

                    if ((n_of_reads + 1) / 10000).is_integer() == False:
                        try:
                            #sim the read
                            get_seq = new_read.simulate_read()
                            #put it in fastq format
                            simulated_reads = sim_paired_end.simulate_perfect_read(new_read,get_seq[0], get_seq[1], get_seq[2])
                            #save the read
                            set_of_left_reads.append(simulated_reads[0])
                            set_of_right_reads.append(simulated_reads[1])
                            n_of_reads +=1

                        except BaseException as e:
                            print("I catched the following exception:",e)

                            pass

                    else:
                        try:
                            #simulate reads and save to disk
                            get_seq = new_read.simulate_read()
                            simulated_reads = sim_paired_end.simulate_perfect_read(new_read, get_seq[0], get_seq[1],
                                                                                   get_seq[2])
                            set_of_left_reads.append(simulated_reads[0])
                            set_of_right_reads.append(simulated_reads[1])

                            #save to disk
                            locker.acquire()
                            print("Process %s: writting to disk 10000 reads" % process)
                            SeqIO.write(set_of_left_reads,paired_end_fastq_1 , "fastq")
                            SeqIO.write(set_of_right_reads, paired_end_fastq_2, "fastq")
                            locker.release()

                            n_of_reads += 1
                            #sim the first read of the list
                            new_read = sim_paired_end(n_of_reads, insert_size, genome_fasta, chr, chr_pos_start,
                                                      chr_pos_end, read_length, circle_number,process)
                            get_seq = new_read.simulate_read()
                            simulated_reads = sim_paired_end.simulate_perfect_read(new_read, get_seq[0], get_seq[1],
                                                                                   get_seq[2])

                            set_of_left_reads = [simulated_reads[0]]
                            set_of_right_reads = [simulated_reads[1]]
                            n_of_reads +=1
                        except BaseException as e:
                            print("I catched the following exception:",e)
                            pass





            circle_bed.append(first_line)

    #shared memory between the processes.This is a list that every process will rate the simulated circles
    for element in circle_bed:
        sim_circles.append(element)




class sim_paired_end:

    #init the class
    def __init__(self,read_number,insert_size,genome_fa,chr,chr_pos_start,chr_pos_end,read_length,circle_id,process):
        self.read_number = read_number
        self.insert_size = insert_size
        self.genome_fa = genome_fa
        self.chr = chr
        self.chr_pos_start = chr_pos_start
        self.chr_pos_end = chr_pos_end
        self.read_length = read_length
        self.circle_id = circle_id
        self.process = process

    def simulate_read(self):
        """Function that simulates perfect paired-end reads"""

        fastafile = ps.FastaFile(self.genome_fa)
        # left split read

        insert = int(np.random.normal(self.insert_size, (self.insert_size / 12), 1))
        start = int(np.random.randint(self.chr_pos_start, (self.chr_pos_end + 1)))
        left_end = start + self.read_length
        total_end = start + int(np.round(insert))
        right_start = total_end - self.read_length
        if total_end > self.chr_pos_end:
            # split read scenario or insert spanning split read scenario
            if left_end > self.chr_pos_end:
                # left read spanning split read scenario
                # left_read
                left_dntps = self.chr_pos_end - start
                right_dntps = self.read_length - left_dntps

                # the error could be here
                left_split_read = fastafile.fetch(self.chr, start, self.chr_pos_end)
                right_split_read = fastafile.fetch(self.chr, self.chr_pos_start, (self.chr_pos_start + right_dntps))
                left_read = left_split_read + right_split_read

                # right_read
                right_start = self.chr_pos_start + int(round(self.insert_size - left_dntps - self.read_length))
                right_read = fastafile.fetch(self.chr, right_start, (right_start + self.read_length))

                # assertion to check the error here

                common_id = "%s|%s|%s:%s-%s:%s|%s:%s|1|%s" % (
                self.read_number, self.chr, start, self.chr_pos_end, self.chr_pos_start, (self.chr_pos_start + right_dntps), right_start,
                (right_start + self.read_length), self.circle_id)




            else:
                if right_start > self.chr_pos_end:
                    # insert spanning split read scenario
                    left_read = fastafile.fetch(self.chr, start, (start + self.read_length))
                    right_start = self.chr_pos_start + (right_start - self.chr_pos_end)
                    right_read = fastafile.fetch(self.chr, right_start, (right_start + self.read_length))
                    common_id = "%s|%s|%s:%s|%s:%s|3|%s" % (
                        self.read_number, self.chr, start, (start + self.read_length), right_start, (right_start + self.read_length), self.circle_id)
                else:
                    # right split read scenario
                    assert right_start <= self.chr_pos_end
                    assert (right_start + self.read_length) > self.chr_pos_end
                    left_read = fastafile.fetch(self.chr, start, (start + self.read_length))

                    # compute right dntps
                    left_dntps = self.chr_pos_end - right_start
                    right_dntps = self.read_length - left_dntps
                    left_split_read = fastafile.fetch(self.chr, right_start, self.chr_pos_end)
                    right_split_read = fastafile.fetch(self.chr, self.chr_pos_start, (self.chr_pos_start + right_dntps))
                    right_read = left_split_read + right_split_read
                    common_id = "%s|%s|%s:%s|%s:%s-%s:%s|2|%s" % (
                        self.read_number,self.chr, start, (start + self.read_length), right_start, self.chr_pos_end, self.chr_pos_start,
                        (self.chr_pos_start +
                         right_dntps), self.circle_id)


        else:
            # non split read scenario
            left_read = fastafile.fetch(self.chr, start, (start + self.read_length))
            # correct right read start
            right_read = fastafile.fetch(self.chr, right_start, (right_start + self.read_length))
            common_id = "%s|%s|%s:%s|%s:%s|0|%s" % (
                self.read_number, self.chr, start, (start + self.read_length), right_start, (right_start + self.read_length), self.circle_id)

        return(right_read,left_read,common_id)



    def simulate_perfect_read(self,right_read,left_read,common_id):
        # put all together
        # unique identifiers for right and left reads
        right_read_id = "2:N:0:CGCTGTG"
        right_id = common_id + "  " + right_read_id
        left_read_id = "1:N:0:CGCTGTG"
        left_id = common_id + "  " + left_read_id
        quality = "I" * self.read_length
        # get the reverse complement of the right read
        right_read = Seq(right_read, generic_dna)
        right_read = right_read.reverse_complement()

        left_read = left_read.upper()

        right_read = right_read.upper()


        fastq_left = "@%s\n%s\n+\n%s\n" % (left_id, left_read, quality)
        fastq_right = "@%s\n%s\n+\n%s\n" % (right_id, right_read, quality)

        right_record = SeqIO.read(StringIO(fastq_right), "fastq")
        left_record = SeqIO.read(StringIO(fastq_left), "fastq")
        return (left_record, right_record)


    def simulate_read_with_errors(self,right_read, left_read, common_id):
        # put all together
        # unique identifiers for right and left reads
        right_read_id = "2:N:0:CGCTGTG"
        right_id = common_id + "space" + right_read_id
        left_read_id = "1:N:0:CGCTGTG"
        left_id = common_id + "space" + left_read_id

        # attemp to use art to simulate the quality scores and the error rate
        #create a one read genome
        left_fasta = open("left_read_%s.fa" % self.process, "w")
        left_fasta.write(">" + left_id + "\n" + str(left_read) + "\n")
        # sim the read with art
        left_fasta.close()
        sp.call("art_illumina -q -na -ss HS25 -nf 0 -i left_read_%s.fa -l %s -f 1 -o left%s" %
                (self.process,self.read_length,self.process),shell=True,stdout=sp.DEVNULL, stderr=sp.STDOUT)
        with open("left%s.fq" % self.process, 'r') as left:
            left_read = left.read().replace('space', '   ').replace('1:N:0:CGCTGTG-1', '1:N:0:CGCTGTG')


        left_record = SeqIO.read(StringIO(left_read), "fastq")
        # get the reverse complement of the right read
        right_read = Seq(right_read, generic_dna)
        right_read = right_read.reverse_complement()

        right_fasta = open("right_read_%s.fa" % self.process, "w")
        right_fasta.write(">" + right_id + "\n" + str(right_read) + "\n")
        right_fasta.close()
        # sim the read with art
        sp.call("art_illumina -na -q -ss HS25  -nf 0 -i right_read_%s.fa -l %s -f 1 -o right%s" %
                (self.process,self.read_length,self.process),shell=True,stdout=sp.DEVNULL, stderr=sp.STDOUT)
        with open("right%s.fq" % self.process, 'r') as right:
            right_read = right.read().replace('space', '   ').replace('1:N:0:CGCTGTG-1', '2:N:0:CGCTGTG')

        right_record = SeqIO.read(StringIO(right_read), "fastq")
        return (left_record, right_record)









