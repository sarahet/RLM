// ==========================================================================
//                                  RLM
// ==========================================================================
// Copyright (c) 2021-2025, Sara Hetzel <hetzel @ molgen.mpg.de>
// Copyright (c) 2021-2025, Max-Planck-Institut für Molekulare Genetik
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
void initialise_argument_parser(sharg::parser & parser, cmd_arguments & args)
{
    parser.info.author = "Sara Hetzel";
    parser.info.short_description = "Read level DNA methylation analysis of bisulfite converted sequencing data.";
    parser.info.version = "1.2.0";

    parser.add_option(args.bam_file,
                      sharg::config{.short_id    = 'b',
                                    .long_id     = "bam",
                                    .description = "BAM or SAM file (best sorted by position and after deduplication).",
                                    .required    = true,
                                    .validator   = sharg::input_file_validator{{"sam", "bam"}}});

    parser.add_option(args.fasta_file,
                      sharg::config{.short_id    = 'r',
                                    .long_id     = "reference",
                                    .description = "Reference genome used to align the BAM file.",
                                    .required    = true,
                                    .validator   = sharg::input_file_validator{{"fa", "fasta"}}});

    parser.add_option(args.mode,
                      sharg::config{.short_id    = 'm',
                                    .long_id     = "mode",
                                    .description = "Sequencing mode.",
                                    .required    = true,
                                    .validator   = sharg::value_list_validator{"SE", "PE"}});

    parser.add_option(args.score,
                      sharg::config{.short_id    = 's',
                                    .long_id     = "score",
                                    .description = "The score(s) to compute. For 'entropy', 'pdr' and 'all' the single read output is always computed.",
                                    .required    = true,
                                    .validator   = sharg::value_list_validator{"single_read", "entropy", "pdr", "all"}});

    parser.add_option(args.aligner,
                      sharg::config{.short_id    = 'a',
                                    .long_id     = "aligner",
                                    .description = "The alignment tool used to create the BAM file.",
                                    .validator   = sharg::value_list_validator{"bsmap", "bismark", "segemehl", "gem"}});

    parser.add_option(args.coverage_filter,
                      sharg::config{.short_id    = 'c',
                                    .long_id     = "coverage",
                                    .description = "Minimum number of reads required to report a CpG or kmer for 'pdr' and 'entropy' mode.",
                                    .validator   = sharg::arithmetic_range_validator{1, 1000}});

    parser.add_option(args.mapq_filter,
                      sharg::config{.short_id    = 'q',
                                    .long_id     = "mapping_quality",
                                    .description = "Minimum mapping quality required to consider a read.",
                                    .validator   = sharg::arithmetic_range_validator{0, 255}});

    parser.add_flag(args.rrbs,
                    sharg::config{.short_id    = 'd',
                                  .long_id     = "rrbs",
                                  .description =
                                  "If BAM file contains reads from an RRBS experiment and reads should be trimmed in order to avoid bias of artifical CpGs. "
                                  "Do NOT use if you already accounted for this problem during trimming."});

    parser.add_option(args.output_file_single_reads,
                      sharg::config{.short_id    = 'o',
                                    .long_id     = "output_single_read",
                                    .description = "Output file with DNA methylation information for every single read with at least 3 CpGs.",
                                    .validator   = sharg::output_file_validator{sharg::output_file_open_options::open_or_create, {"bed", "tsv", "txt"}}});

    parser.add_option(args.output_file_entropy,
                      sharg::config{.short_id    = 'e',
                                    .long_id     = "output_entropy",
                                    .description = "Output file with entropy, epipolymorphism and epiallele information for every 4-mer spanned by complete reads.",
                                    .validator   = sharg::output_file_validator{sharg::output_file_open_options::open_or_create, {"bed", "tsv", "txt"}}});

    parser.add_option(args.output_file_pdr,
                      sharg::config{.short_id    = 'p',
                                    .long_id     = "output_pdr",
                                    .description =
                                    "Output file with read-transition score and percent discordant reads for every CpG spanned by complete reads. "
                                    "Only reads that cover at least 3 CpGs are considered.",
                                    .validator   = sharg::output_file_validator{sharg::output_file_open_options::open_or_create, {"bed", "tsv", "txt"}}});
}
