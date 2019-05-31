# Welcome to Circle-Map official repository!

Circle-Map is an easy to install, python package that implements all the steps required to detect extrachromosomal DNA circles. The package  contains easy to run algorithms for accurately detect circular DNA formed from mappable and non mappable regions of a genome.

    
## Why should I use Circle-Map?

Circle-Map takes as input an alignment of reads to a reference genome (e.g. a *BWA-MEM* generated *BAM* file) and like other methods, it will use those alignments to detect cases were the read has been split into two segments (e.g. split reads) to detect genomic rearrangements supporting a circular DNA structure.

However, this approach results in many split read alignments being missed because the aligner is not able to map both split segments of the read, either because they are too short or because they align to too many places. In this cases, the aligner will report a read alignment containing some of the bases unmapped (e.g soft-clipped reads). 

Unlike other methods, Circle-Map is able to map both segments of the soft-clipped reads by realigning the unmapped parts probabilistically to a graph representation of the circular DNA breakpoints. Which allows for a more accurate detection of the circular DNA breakpoints. In our recent preprint (https://www.biorxiv.org/content/10.1101/654194v1) we show how this approach dramatically increases sensitivity while retaining high precision.


## Getting started

### Installation

Circle-Map runs on GNU/Linux operating systems, and can be installed and set-up using the following ways:

This will install Circle-Map, and all the external packages required run every part of Circle-Map software.

Installation using **pip**:

     python -m pip install Circle-Map
     
Installation through conda/bioconda will be supported soon.
     
**Note**: If you have use pip and you want to simulate circular DNA short reads, you will need to install [BBMap](https://sourceforge.net/projects/bbmap/) and [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) on your system.

### Using Circle-Map

Now you are ready to get started detecting circular DNA. We have created a [Circle-Map wiki](https://github.com/iprada/Circle-Map/wiki) that explains step by step how you can go from your raw sequencing reads to interpretable results. In the wiki, you can try and learn using Circle-Map with the following tutorials:

   * [Tutorial: identification of mappable circular DNA using Circle Map Realign](https://github.com/iprada/Circle-Map/wiki/Tutorial:-identification-of-mappable-circular-DNA-using-Circle-Map-Realign) 
    

## Getting help

The best place for getting help, feedback,report bugs or request new features is to post an [issue](https://github.com/iprada/Circle-Map/issues). You can also reach me at xsh723 at dot binf dot ku dot dk 
## Citing

*Sensitive detection of circular DNA at single-nucleotide resolution using guided realignment of partially aligned reads*

Inigo Prada-Luengo, Anders Krogh, Lasse Maretty, Birgitte Regenberg
**bioRxiv** 654194; doi: https://doi.org/10.1101/654194

A manuscript describing the methods is currently under revision.

## License

Circle-Map is freely available under the [MIT license](https://opensource.org/licenses/MIT)

## Acknowledgements

Circle-Map is being developed by Iñigo Prada-Luengo, Anders Krogh, Lasse Maretty and Birgitte Regenberg at the University of Copenhagen
