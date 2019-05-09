# Welcome to Circle-Map official repository!

Circle-Map is an easy to install, python package that implements all the steps required to detect extrachromosomal DNA circles. The package  contains easy to run algorithms for accurately detect circular DNA formed from mappable and non mappable regions of a genome.

    
## Why should I use Circle-Map?

Circle-Map takes as input an alignment of reads to a reference genome (e.g. a *BWA-MEM* generated *BAM* file) and like other methods, it will use those alignments to detect cases were the read has been split into two segments (e.g. split reads) to detect genomic rearrangements supporting a circular DNA structure.

However, this approach results in many split read alignments being missed because the aligner is not able to map both split segments of the read, either because they are too short or because they align to many places. In this cases, the aligner will report a read alignment containing some of the bases unmapped (e.g soft-clipped reads). 

Unlike other methods, Circle-Map is able to map both segments of the soft-clipped reads by realigning the unmapped parts probabilistically to a graph representation of the circular DNA breakpoints, which in turn allows for a more accurate detection of the circular DNA breakpoints. In our recent paper/preprint (**TO BE ADDED**) we show how this approach dramatically increases sensitivity while retaining high precision.


## Getting started

### Installation

Circle-Map, can be easily installed through **conda** using the following command:

     conda install -c bioconda I_NEED_TO_CREATE_THE_RECIPE

#### Dependencies

## Contact

## Citing

## Third party software

## Acknowledgements

Circle-Map is being developed by Iñigo Prada-Luengo, Lasse Maretty and Birgitte Regenberg at the University of Copenhagen
