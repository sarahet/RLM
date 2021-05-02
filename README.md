# RLM 
Read level DNA methylation analysis of bisulfite converted sequencing data

![](https://github.com/sarahet/RLM/actions/workflows/ci.yml/badge.svg)

## Dependencies
* GCC   (minimum required version: 7)
* CMake (minimum required version: 3.8)

## Installation
To install RLM from github, run:
```
git clone --recurse-submodules https://github.com/sarahet/RLM.git

mkdir build && cd build
cmake ../RLM
make
```

## Usage
To use RLM from your build directory, run:
```
bin/RLM -b /path/to/bam/<sample>.bam \
        -r /path/to/reference/<reference>.fa \
        -m <sequencing_mode> \
        -s <desired_score> \
        -a <aligner>
```

To see the help page and learn more about available basic options, run:
```
bin/RLM -h
```

To see available advanced options, run:
```
bin/RLM -hh
```

All options available in RLM:
```
-h, --help                Prints the help page.

-hh, --advanced-help      Prints the help page including advanced options.

--version                 Prints the version information.

--copyright               Prints the copyright/license information.

--export-help             Export the help page information. Value must be one of [html, man].

--version-check           Whether to check for the newest app version. Default: true.

-b, --bam                 BAM or SAM file (best sorted by position and after deduplication). The
                          input file must exist and read permissions must be granted. Valid file
                          extensions are: [sam, bam].

-r, --reference           Reference genome used to align the BAM file. The input file must exist
                          and read permissions must be granted. Valid file extensions are:
                          [fa, fasta].

-m, --mode                Sequencing mode. Value must be one of [SE,PE].

-s, --score               The score(s) to compute. For 'entropy', 'pdr' and 'all' the single
                          read output is always computed. Value must be one of
                          [single_read,entropy,pdr,all].

-a, --aligner             The alignment tool used to create the BAM file. Default: bsmap. Value
                          must be one of [bsmap,bismark,segemehl].

-c, --coverage            Minimum number of reads required to report a CpG or kmer for 'pdr'
                          and 'entropy' mode. Default: 10. Value must be in range [1,1000].

-q, --mapping_quality     Minimum mapping quality required to consider a read. Default: 30.
                          Value must be in range [0,255].

-d, --rrbs                Set if BAM file contains reads from an RRBS experiment.

-o, --output_single_read  Output file with DNA methylation information for every single read
                          with at least 3 CpGs. Default: "output_single_read_info.bed". Write
                          permissions must be granted.
                          Valid file extensions are: [bed, tsv, txt].

-e, --output_entropy      Output file with entropy, epipolymorphism and epiallele information
                          for every 4-mer spanned by complete reads.
                          Default: "output_entropy.bed". Write permissions must be granted.
                          Valid file extensions are: [bed, tsv, txt].

-p, --output_pdr          Output file with read-transition score and percent discordant reads
                          for every CpG spanned by complete reads. Only reads that cover at
                          least 3 CpGs are considered. Default: "output_pdr.bed". Write
                          permissions must be granted.
                          Valid file extensions are: [bed, tsv, txt].
```
