##Bugs to fix


## To MSc thesis

-Run all the simulations again with the fixed issue of the repeat read when saving to disk
-In the realignment of the soft-clip reads, compare to the origin
-Performance on real data
-Performance on simulated data
- Remove os.system
- Use pybedtools
- Warnings in the loops of alignment.py

## To publish

- Tests
- Check that read extractor is doing well
- Check the orientation of everything
- Check 0-1 based indexes
- Make friendly command line interface with argparse
- Implement coverage

## ReadExtractor
- When I am extracting the discordant reads. Could there be discordant reads aligned to the breakpoint?

## alignment
- I am going to create the coverage file that serves as input for generating the realignment priors. This helps in the development,
remember to remove this when the pipeline is finished.

- During the realignment step check that that the soft-clipped reads have their mate inside the putative eccDNA bounda-
ries

- If there are only discordant readsin the realignment prior interval, sample 100 F1R2 reads and compute the standard deviation of the insert size. Then, create a realignment interval of 4 standard deviations of the insert size

- When programming the insert size distribution estimation code warnings if the input bam is not queryname sorted

- Take out some calculations from the down loops. The interval extension with the insert size for example.

- Check why some realignment intervals are extremenly big

- Load contig fasta in memory
- Make a rescue round once Circle-Map has finished with the intervals were there are only soft-clipped reads. For read in
the sc interval check that they do not have supplementary alignments and align them with bwa mem to the whole chromosome?
or skip