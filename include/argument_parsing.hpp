// ==========================================================================
//                                  RLM
// ==========================================================================
// Copyright (c) 2021, Sara Hetzel <hetzel @ molgen.mpg.de>
// Copyright (c) 2021, Max-Planck-Institut für Molekulare Genetik
// All rights reserved.
//
// This file is part of RLM.
//
// RLM is Free Software: you can redistribute it and/or modify it
// under the terms found in the LICENSE[.md|.rst] file distributed
// together with this file.
//
// RLM is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//
// ==========================================================================
// Functions and data structures for argument parsing
// ==========================================================================

#pragma once

// Struct that stores command line arguments
struct cmd_arguments
{
    std::filesystem::path bam_file{};
    std::filesystem::path fasta_file{};

    std::filesystem::path output_file_single_reads{"output_single_read_info.bed"};
    std::filesystem::path output_file_entropy{"output_entropy.bed"};
    std::filesystem::path output_file_pdr{"output_pdr.bed"};

    uint32_t verbosity = 0;
    uint32_t mapq_filter = 30;
    uint32_t coverage_filter = 10;

    bool rrbs = false;

    std::string mode;
    std::string score = "single_read";
    std::string aligner = "bsmap";
};

// Function to initialize the argument parser
void initialise_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Sara Hetzel";
    parser.info.version = "1.0.0";
    parser.info.short_description = "Read level DNA methylation analysis of bisulfite converted sequencing data.";

    parser.add_option(args.bam_file, 'b', "bam",
                      "BAM or SAM file (best sorted by position and after deduplication).",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{{"sam", "bam"}});

    parser.add_option(args.fasta_file, 'r', "reference",
                      "Reference genome used to align the BAM file.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{{"fa", "fasta"}});

    parser.add_option(args.mode, 'm', "mode",
                      "Sequencing mode.",
                      seqan3::option_spec::required,
                      seqan3::value_list_validator{"SE", "PE"});

    parser.add_option(args.score, 's', "score",
                      "The score(s) to compute. For 'entropy', 'pdr' and 'all' the single read output is always computed.",
                      seqan3::option_spec::required,
                      seqan3::value_list_validator{"single_read", "entropy", "pdr", "all"});

    parser.add_option(args.aligner, 'a', "aligner",
                      "The alignment tool used to create the BAM file.",
                      seqan3::option_spec::standard,
                      seqan3::value_list_validator{"bsmap", "bismark", "segemehl"});

    parser.add_option(args.coverage_filter, 'c', "coverage",
                      "Minimum number of reads required to report a CpG or kmer for 'pdr' and 'entropy' mode.",
                      seqan3::option_spec::advanced,
                      seqan3::arithmetic_range_validator{1, 1000});

    parser.add_option(args.mapq_filter, 'q', "mapping_quality",
                      "Minimum mapping quality required to consider a read.",
                      seqan3::option_spec::advanced,
                      seqan3::arithmetic_range_validator{0, 255});

    parser.add_flag(args.rrbs, 'd', "rrbs",
                    "If BAM file contains reads from an RRBS experiment and reads should be trimmed in order to avoid bias of artifical CpGs. "
                    "Do NOT use if you already accounted for this problem during trimming.",
                    seqan3::option_spec::advanced);

    parser.add_option(args.output_file_single_reads, 'o', "output_single_read",
                      "Output file with DNA methylation information for every single read with at least 3 CpGs.",
                      seqan3::option_spec::standard,
                      seqan3::output_file_validator{seqan3::output_file_open_options::open_or_create, {"bed", "tsv", "txt"}});

    parser.add_option(args.output_file_entropy, 'e', "output_entropy",
                      "Output file with entropy, epipolymorphism and epiallele information for every 4-mer spanned by complete reads.",
                      seqan3::option_spec::standard,
                      seqan3::output_file_validator{seqan3::output_file_open_options::open_or_create, {"bed", "tsv", "txt"}});

    parser.add_option(args.output_file_pdr, 'p', "output_pdr",
                      "Output file with read-transition score and percent discordant reads for every CpG spanned by complete reads. "
                      "Only reads that cover at least 3 CpGs are considered.",
                      seqan3::option_spec::standard,
                      seqan3::output_file_validator{seqan3::output_file_open_options::open_or_create, {"bed", "tsv", "txt"}});
}
