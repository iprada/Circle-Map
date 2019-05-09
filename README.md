# Circle-Map

Welcome to Circle-Map official repository!

Circle-Map is an easy to install, python package that implements all the steps required to detect extrachromosomal DNA circles. The package  contains easy to run algorithms for accurately detect circular DNA formed from mappable and non mappable regions of a genome.

    
## Why should I use Circle-Map?

Circle-Map takes as input an alignment of reads to a reference genome (e.g. a *BWA-MEM* generated *BAM* file) and like other methods, it will use those alignments to detect cases were the read has been splitted into two segments (e.g. split reads) to detect genomic rearrangements supporting a circular DNA structure.

However, this approach results in many split read alignments being missed because the aligner is not able to map both split segments. In this cases, the aligner will report a read alignment containing some of the bases unmapped (e.g soft-clipped reads). Unlike other methods Circle-Map is able to map both segments of this reads by realigning them probabilistically to a graph representation of the circular DNA breakpoints, which in turn allows for a more accurate detection of the circular DNA breakpoints.


## Getting started

## Contact

## Citing

## Third party software

## Acknowledgements

Circle-Map is being developed by Iñigo Prada-Luengo, Lasse Maretty and Birgitte Regenberg at the University of Copenhagen
