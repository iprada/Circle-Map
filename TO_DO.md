##Bugs to fix


## To MSc thesis

-Run all the simulations again with the fixed issue of the repeat read when saving to disk
-Performance on real data
-Performance on simulated data
- Warnings in the loops of alignment.py

## To publish

- Tests



- Make friendly command line interface with argparse
- Implement coverage
- Change some containers to dictionary for easy code reading


## alignment
- I am going to create the coverage file that serves as input for generating the realignment priors. This helps in the development,
remember to remove this when the pipeline is finished.

- During the realignment step check that that the soft-clipped reads have their mate inside the putative eccDNA bounda-
ries

- When programming the insert size distribution estimation code warnings if the input bam is not queryname sorted

- To check out: Take out all possible calculations from the down loops. The interval extension with the insert size for example. UPDATE=25/01/2018 This does not seem to be an issue

- Make a rescue round once Circle-Map has finished with the intervals were there are only soft-clipped reads. For read in
the sc interval check that they do not have supplementary alignments and align them with bwa mem to the whole chromosome? UPDATE=25/01/2018: Use Circle length distribtion as a prior
or skip lonely soft-clipped reads

- Regarding the above comment. In the function get_realignment_interval there is a if-elif-else loop. In the else part of the
loop there should be a saved bed file for the rescue of the lonely soft-clips

## May be

- Try efficient looping with itertools

## Limitations

- Circle-Map realignment of the soft-clipped reads does not support reads with ONLY N. UPDATE=25/01/2018: This is not a limitation. SUpported
