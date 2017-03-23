Format of the fastq header(@line)

For paired-end datasets:

@0 = Normal read
@1 = Left split read
@2 = Right split read
@3 = Split scenario

Then, the next lines in the header indicate the position:
@scenario|chr|non_split_start:non_split_end|start_left_split:end_left_split-start_right_split:end_right_split
