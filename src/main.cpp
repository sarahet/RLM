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
// Main program
// ==========================================================================

#include <algorithm>
#include <fstream>
#include <numeric>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include <seqan3/range/views/slice.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/std/ranges>

#include "../include/argument_parsing.hpp"
#include "../include/data_structures.hpp"
#include "../include/methylation_scores.hpp"
#include "../include/output.hpp"
#include "../include/process_record.hpp"

using seqan3::operator""_tag;
using seqan3::operator""_dna5;
using seqan3::operator""_cigar_operation;

// Forward declaration
template <bool calc_pdr_score,  bool calc_entropy_score>
int arg_conv1(cmd_arguments & args);

template <bool calc_pdr_score,  bool calc_entropy_score, bool rrbs>
int arg_conv2(cmd_arguments & args);

template <bool calc_pdr_score,  bool calc_entropy_score, bool rrbs, align_type aligner>
int real_main(cmd_arguments & args);

// Main function to parse arguments and set template arguments depending on input score selected
int main(int argc, char ** argv)
{
    // The argument parser
    seqan3::argument_parser parser{"RLM", argc, argv};
    cmd_arguments args{};

    initialise_argument_parser(parser, args);

    try
    {
         parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "Parsing error. " << ext.what() << "\n";
        return -1;
    }

    // Get score enum
    score_type score = _score_name_to_enum(args.score);

    try
    {
        switch (score)
        {
            case score_type::SINGLE_READ:  return arg_conv1<false, false>(args);
            case score_type::PDR:          return arg_conv1<true, false>(args);
            case score_type::ENTROPY:      return arg_conv1<false, true>(args);
            case score_type::ALL:          return arg_conv1<true, true>(args);
            default: throw "Undefined score requested.";
        }
    }
    catch (const char * e)
    {
        std::cerr << "Error: " << e << std::endl;
        return -1;
    }
}

template <bool calc_pdr_score,  bool calc_entropy_score>
int arg_conv1(cmd_arguments & args)
{
    sequencing_type type = _sequencing_type_to_enum(args.rrbs);
    switch (type)
    {
        case sequencing_type::RRBS:  return arg_conv2<calc_pdr_score, calc_entropy_score, true>(args);
        case sequencing_type::WGBS:  return arg_conv2<calc_pdr_score, calc_entropy_score, false>(args);
    }
}

template <bool calc_pdr_score,  bool calc_entropy_score, bool rrbs>
int arg_conv2(cmd_arguments & args)
{
    align_type type = _aligner_name_to_enum(args.aligner);

    try
    {
        switch (type)
        {
            case align_type::BSMAP:     return real_main<calc_pdr_score, calc_entropy_score, rrbs, align_type::BSMAP>(args);
            case align_type::BISMARK:   return real_main<calc_pdr_score, calc_entropy_score, rrbs, align_type::BISMARK>(args);
            case align_type::SEGEMEHL:  return real_main<calc_pdr_score, calc_entropy_score, rrbs, align_type::SEGEMEHL>(args);
            default: throw "Undefined alignment tool requested.";
        }
    }
    catch (const char * e)
    {
        std::cerr << "Error: " << e << std::endl;
        return -1;
    }
}

// Real main function containing the program
template <bool calc_pdr_score, bool calc_entropy_score, bool rrbs, align_type aligner>
int real_main(cmd_arguments & args)
{
    std::cout << "Starting RLM" << std::endl;

    // Set threads for BAM decompression
    seqan3::contrib::bgzf_thread_count = args.threads;

    // Load genome reference file
    std::cout << "Reading the reference genome" << std::endl;

    seqan3::sequence_file_input reference_file{args.fasta_file};
    std::vector<std::string> genome_seqs_ids{};
    std::vector<seqan3::dna5_vector> genome_seqs{};

    for (auto && record : reference_file)
    {
        genome_seqs_ids.push_back(std::move(record.id()));
        genome_seqs.push_back(std::move(record.sequence()));
    }

    // Initialize BAM file stream
    std::cout << "Opening the bam file" << std::endl;

    using field_type = seqan3::fields<seqan3::field::id,
                                      seqan3::field::flag,
                                      seqan3::field::ref_id,
                                      seqan3::field::ref_offset,
                                      seqan3::field::mapq,
                                      seqan3::field::seq,
                                      seqan3::field::cigar,
                                      seqan3::field::tags>;
    seqan3::sam_file_input mapping_file{args.bam_file, field_type{}};

    // Validate whether reference genome and BAM file have the same order of chromosomes
    try
    {
        if (mapping_file.header().ref_ids().size() != genome_seqs_ids.size())
            throw "Different number of sequences in references and BAM file.";

        for (size_t i = 0; i < mapping_file.header().ref_ids().size(); i++)
        {
            if (mapping_file.header().ref_ids()[i] != genome_seqs_ids[i])
                throw "Different order of sequences in references and BAM file.";
        }
    }
    catch (const char * e)
    {
        std::cerr << "Error: " << e << std::endl;
        return -1;
    }

    // Map used to store reads until mate is read
    using record_t = std::ranges::range_value_t<decltype(mapping_file)>;
    std::map<std::string, record_t> records;

    // Map to store CpGs
    using num_reads_t = uint32_t;
    using num_discordant_reads_t = uint32_t;
    using sum_transitions_t = double;

    std::map<GenomePosition, std::tuple<num_reads_t, num_discordant_reads_t, sum_transitions_t> > all_CpGs;

    // Map to store 4-mers with epialleles
    std::map<GenomePosition, std::vector<uint32_t> > all_kmers;

    // Set mode for calculations
    using score_tag = score_tag<calc_pdr_score, calc_entropy_score>;

    std::ofstream output_stream;
    output_stream.open(args.output_file_single_reads);
    write_header_read_info(output_stream);

    std::cout << "Starting BAM file processing" << std::endl;

    // for (auto & rec : mapping_file | mapq_filter)
    for (auto & rec : mapping_file)
    {
        // Check if read is properly mapped and paired (if PE mode), not vendor-failed, not supplementary
        // or secondary alignment and not PCR duplicate
        if ((!static_cast<bool>(rec.flag() & seqan3::sam_flag::paired) && args.mode == "PE") ||
            static_cast<bool>(rec.flag() & seqan3::sam_flag::unmapped) ||
            static_cast<bool>(rec.flag() & seqan3::sam_flag::secondary_alignment) ||
            static_cast<bool>(rec.flag() & seqan3::sam_flag::failed_filter) ||
            static_cast<bool>(rec.flag() & seqan3::sam_flag::duplicate) ||
            static_cast<bool>(rec.flag() & seqan3::sam_flag::supplementary_alignment) ||
            rec.mapping_quality() < args.mapq_filter)
            continue;

        // Check if alignment contains indels: If yes, skip record.
        using seqan3::get;
        bool indel = false;
        for (auto c : rec.cigar_sequence())
        {
            if (get<seqan3::cigar::operation>(c) != 'M'_cigar_operation)
            {
                indel = true;
                break;
            }
        }

        if (indel)
            continue;

        // If RRBS mode, change potentially artifical bases to 'N'
        if constexpr (rrbs)
        {
            if (args.mode == "SE" ||
                (args.mode == "PE" && static_cast<bool>(rec.flag() & seqan3::sam_flag::first_in_pair)))
            {
                rec.sequence()[0] = 'N'_dna5;
                rec.sequence()[1] = 'N'_dna5;
            }
            else
            {
                rec.sequence()[rec.sequence().size() - 1] = 'N'_dna5;
                rec.sequence()[rec.sequence().size() - 2] = 'N'_dna5;
            }
        }

        // Set tag depending on aligner
        read_type rec_type;

        if constexpr (aligner == align_type::BSMAP)
        {
            rec_type = _read_tag_bsmap_to_enum(rec.tags().get<"ZS"_tag>());
        }
        else if (aligner == align_type::BISMARK)
        {
            rec_type = _read_tag_bismark_to_enum(rec.tags().get<"XG"_tag>());
        }
        else
        {
            rec_type = _read_tag_segemehl_to_enum(rec.tags().get<"XB"_tag>());
        }

        if (args.mode == "SE")
        {
            process_bam_record(output_stream,
                               rec_type,
                               rec.reference_id().value(),
                               rec.reference_position().value(),
                               rec.sequence(),
                               rec.id(),
                               mapping_file.header().ref_ids(),
                               genome_seqs,
                               all_CpGs,
                               all_kmers,
                               score_tag{});
        }
        else
        {
            // Check if mate has already been read
            if (!records.insert(std::make_pair(rec.id(), rec)).second)
            {
                // Check if reads are overlapping
                int overlap = std::min(records.at(rec.id()).sequence().size() + records.at(rec.id()).reference_position().value(),
                                       rec.sequence().size() + rec.reference_position().value()) -
                              std::max(records.at(rec.id()).reference_position().value(), rec.reference_position().value());

                if (overlap > 0)
                {
                    // First check if one read is included in the other - just process the longer one in this case
                    if (overlap == records.at(rec.id()).sequence().size())
                    {
                        // First read included in second read
                        process_bam_record(output_stream,
                                           rec_type,
                                           rec.reference_id().value(),
                                           rec.reference_position().value(),
                                           rec.sequence(),
                                           rec.id(),
                                           mapping_file.header().ref_ids(),
                                           genome_seqs,
                                           all_CpGs,
                                           all_kmers,
                                           score_tag{});

                    }
                    else if (overlap == rec.sequence().size())
                    {
                        // Second read included in first read
                        process_bam_record(output_stream,
                                           rec_type,
                                           records.at(rec.id()).reference_id().value(),
                                           records.at(rec.id()).reference_position().value(),
                                           records.at(rec.id()).sequence(),
                                           records.at(rec.id()).id(),
                                           mapping_file.header().ref_ids(),
                                           genome_seqs,
                                           all_CpGs,
                                           all_kmers,
                                           score_tag{});
                    }
                    else
                    {
                        // Merge reads
                        // Determine which read comes first
                        bool is_first = records.at(rec.id()).reference_position().value() <= rec.reference_position().value();
                        auto & rec1 = is_first ? records.at(rec.id()) : rec;
                        auto & rec2 = is_first ? rec : records.at(rec.id());

                        seqan3::dna5_vector second_seq_part = rec2.sequence()
                                                            | seqan3::views::slice(overlap, rec2.sequence().size())
                                                            | seqan3::views::to<seqan3::dna5_vector>;
                        rec1.sequence().insert(rec1.sequence().end(), second_seq_part.begin(), second_seq_part.end());

                        process_bam_record(output_stream,
                                           rec_type,
                                           rec1.reference_id().value(),
                                           rec1.reference_position().value(),
                                           rec1.sequence(),
                                           rec1.id(),
                                           mapping_file.header().ref_ids(),
                                           genome_seqs,
                                           all_CpGs,
                                           all_kmers,
                                           score_tag{});
                    }
                }
                else
                {
                    // Process both mate records
                    // Record already stored in map
                    process_bam_record(output_stream,
                                       rec_type,
                                       records.at(rec.id()).reference_id().value(),
                                       records.at(rec.id()).reference_position().value(),
                                       records.at(rec.id()).sequence(),
                                       records.at(rec.id()).id(),
                                       mapping_file.header().ref_ids(),
                                       genome_seqs,
                                       all_CpGs,
                                       all_kmers,
                                       score_tag{});
                    // Current record
                    process_bam_record(output_stream,
                                       rec_type,
                                       rec.reference_id().value(),
                                       rec.reference_position().value(),
                                       rec.sequence(),
                                       rec.id(),
                                       mapping_file.header().ref_ids(),
                                       genome_seqs,
                                       all_CpGs,
                                       all_kmers,
                                       score_tag{});
                }

                // Remove from map, now not used anymore
                records.erase(rec.id());
            }
            else
            {
                continue;
            }
        }
    }

    output_stream.close();
    std::cout << "Finished BAM file processing" << std::endl;
    std::cout << "Finished writing 'single_read' output" << std::endl;

    if constexpr (calc_pdr_score)
    {
        std::cout << "Starting PDR and RTS calculations" << std::endl;

        std::ofstream output_stream_pdr;
        output_stream_pdr.open(args.output_file_pdr);
        write_header_pdr(output_stream_pdr);

        for (auto it = all_CpGs.begin(); it != all_CpGs.end(); it++)
        {
            write_record_pdr(output_stream_pdr, mapping_file.header().ref_ids(), it->first, it->second, args.coverage_filter);
        }

        output_stream_pdr.close();

        std::cout << "Finished writing 'pdr' output" << std::endl;
    }

    if constexpr (calc_entropy_score)
    {
        std::cout << "Starting entropy and epipolymorphism calculations" << std::endl;

        std::ofstream output_stream_entropy;
        output_stream_entropy.open(args.output_file_entropy);
        write_header_entropy(output_stream_entropy);

        for (auto it = all_kmers.begin(); it != all_kmers.end(); it++)
        {
            write_record_entropy(output_stream_entropy, mapping_file.header().ref_ids(), it->first, it->second, args.coverage_filter);
        }

        output_stream_entropy.close();

        std::cout << "Finished writing 'entropy' output" << std::endl;
    }

    std::cout << "Terminating RLM" << std::endl;

    return 0;
}
