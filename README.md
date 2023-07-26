# RLM
Read level DNA methylation analysis of bisulfite converted sequencing data

![](https://github.com/sarahet/RLM/actions/workflows/ci_linux.yml/badge.svg)
![](https://github.com/sarahet/RLM/actions/workflows/ci_macos.yml/badge.svg)

For detailed documentation and usage examples, please visit the [Wiki](https://github.com/sarahet/RLM/wiki).

## Dependencies
* GCC   (minimum required version: 7)
* CMake (minimum required version: 3.8)
* zlib (minimum required version: 1.2)

**Attention:** Due to the requirements of the SeqAn3 library, only the latest minor GCC releases are supported for each major version. 

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
                          must be one of [bsmap,bismark,segemehl,gem].

-c, --coverage            Minimum number of reads required to report a CpG or kmer for 'pdr'
                          and 'entropy' mode. Default: 10. Value must be in range [1,1000].

-q, --mapping_quality     Minimum mapping quality required to consider a read. Default: 30.
                          Value must be in range [0,255].

-d, --rrbs                If BAM file contains reads from an RRBS experiment and reads should
                          be trimmed in order to avoid bias of artifical CpGs. Do NOT use if
                          you already accounted for this problem during trimming.

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

## Visualization with R

RLM is a standalone C++ application. However, to provide summary statistics and figures as well as ideas for post-processing, we provide an [R Markdown script](https://github.com/sarahet/RLM/blob/main/post_processing/summarize_read_level_stats.Rmd) in the ```post_processsing``` folder. In order to use the script, R needs to be installed including the following packages:

* ```knitr```
* ```data.table```
* ```GenomicRanges```
* ```RColorBrewer```
* ```vioplot```
* ```ggplot2```
* ```ggpubr```

The script can be called the following way:

```
Rscript -e "rmarkdown::render('summarize_read_level_stats.Rmd', 
params=list(
single_read_input_file = '/path/to/output_single_read_info.bed',
pdr_input_file = '/path/to/output_pdr.bed',
entropy_input_file = '/path/to/output_entropy.bed',
sample_name = 'my_sample',
feature_input_file = '/path/to/features.bed'), 
output_file = 'my_output.pdf')"
```

The parameter ```feature_input_file``` is optional and if not provided, no feature-wise figures will be reported. If provided, it should be a bedgraph file of the following format: ```<chr> <start> <end> <feature_name>``` where ```<feature_name>``` should be the name of the feature type the region belongs to. For more information and an example run with figures please visit our [Wiki](https://github.com/sarahet/RLM/wiki/Post-processing-and-use-cases).
