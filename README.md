# Author Iñigo Prada

Format of the fastq header(@line)

For paired-end datasets, the script is designed as follows, all of them will have six fields in the header:

@read_number|chromosome|right_read_coordinates|left_read_coordinates|type|circle_number

When the a read covers the split read, that will be denoted as start:circle_end-circle_start:end

0 = normal read
1 = left split read
2 = Right split read
3 = split read scenario



For single-end datasets the header is designed as follows:

For the split reads the information is encoded as follows:

In the split read case, type will be equal to 1, and the header will be:

@read_number|chromosome|left_nucleotides_start:circle_end-circle_start:right_nucleotides_end|type|circle_number



For the reads that do not span the split read:

@read_number|chromosome|read_start:read_end|type|circle_number

