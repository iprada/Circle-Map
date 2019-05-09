# Circle-Map

Welcome to Circle-Map official repository!

Circle-Map is an easy to install, python package that implements all the steps to easily detect extrachromosomal DNA circles. The package  contains easy to run algorithms for accurately detect circular DNA formed from mappable and non mappable regions of a genome.

    
## Why should I use Circle-Map?

Circle-Map takes as input an alignment of reads to a reference genome (e.g. a *BWA-MEM* generated *BAM* file) and like other methods, it will use those alignments to detect cases were the read has been splitted in to two segments (e.g. split reads) to detect genomic rearrangements supporting a circular DNA structure.

However, this approach results in many split read alignments being missed because the aligner is not able to map both split segments, reporting read alignments with some of the bases unmapped (e.g soft-clipped reads).



- Circle-Map is a split read detection method that like other methods is build on top of bwa
- Like other methods Circle-Map detects split reads by identifying reads that have been remaped by the aligner
- However, in many cases this reads are to short to be aligned
- Unlike other methods Circle-Map goes one step further. Circle-Map tries to obtain the location of ultra short read (down 5bp) by building a probabilistic model of the breakpoint of the reads
- This allows to narrow down the search

## Getting started

## Contact

## Citing

## Third party software

## Acknowledgements

Circle-Map is being developed by Iñigo Prada-Luengo, Lasse Maretty and Birgitte Regenberg at the University of Copenhagen
