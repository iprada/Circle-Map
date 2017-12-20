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
- During the realignment step check that that the soft-clipped reads have their mate inside the putative eccDNA bounda-
ries

- I am doing the fetch of the bam with until_eof= True. Which does not require file indexing. If the code is to slow,
remember that it could be due to that. If it is slow, check and change accordingly


