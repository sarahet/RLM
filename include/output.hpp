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
// Functions for the output
// ==========================================================================

#pragma once

#include "methylation_scores.hpp"

using num_reads_t = uint32_t;
using num_discordant_reads_t = uint32_t;
using sum_transitions_t = double;

// Write header for 'single_read' mode
void write_header_read_info(std::ofstream & output_stream)
{
    if (output_stream.is_open())
    {
        output_stream << "#chr\t"
                      << "start\t"
                      << "end\t"
                      << "read_name\t"
                      << "CpG_pattern\t"
                      << "n_CpGs\t"
                      << "n_CpGs_methyl\t"
                      << "discordance_score\t"
                      << "transitions_score\t"
                      << "mean_methylation\n";
    }
    else
    {
        throw std::runtime_error("ERROR: Could not open single read information output file.");
    }
}

// Write header for 'entropy' mode
void write_header_entropy(std::ofstream & output_stream)
{
    if (output_stream.is_open())
    {
        output_stream << "#chr\t"
                      << "start\t"
                      << "end\t"
                      << "entropy\t"
                      << "epipolymorphism\t"
                      << "gggg\t"
                      << "gggG\t"
                      << "ggGg\t"
                      << "ggGG\t"
                      << "gGgg\t"
                      << "gGgG\t"
                      << "gGGg\t"
                      << "gGGG\t"
                      << "Gggg\t"
                      << "GggG\t"
                      << "GgGg\t"
                      << "GgGG\t"
                      << "GGgg\t"
                      << "GGgG\t"
                      << "GGGg\t"
                      << "GGGG\t"
                      << "coverage\n";
    }
    else
    {
        throw std::runtime_error("ERROR: Could not open entropy output file.");
    }
}

// Write header for 'pdr' mode
void write_header_pdr(std::ofstream & output_stream)
{
    if (output_stream.is_open())
    {
        output_stream << "#chr\t"
                      << "start\t"
                      << "end\t"
                      << "PDR\t"
                      << "RTS\t"
                      << "coverage\n";
    }
    else
    {
        throw std::runtime_error("ERROR: Could not open PDR output file.");
    }
}

// Write record for 'entropy' mode
void write_record_entropy(std::ofstream & output_stream,
                          std::deque<std::string> const & ref_ids,
                          GenomePosition const & pos,
                          std::vector<uint32_t> const & epialleles,
                          uint32_t const & coverage_filter)
{
    uint32_t coverage = std::accumulate(epialleles.begin(), epialleles.end(), 0);
    if (coverage < coverage_filter)
        return;

    output_stream << ref_ids[pos.ref_id] << "\t"
                  << pos.start << "\t"
                  << pos.start + 2 << "\t"
                  << calculate_entropy_across_reads(epialleles) << "\t"
                  << calculate_epipolymorphism_across_reads(epialleles) << "\t"
                  << epialleles[0] << "\t"
                  << epialleles[1] << "\t"
                  << epialleles[2] << "\t"
                  << epialleles[3] << "\t"
                  << epialleles[4] << "\t"
                  << epialleles[5] << "\t"
                  << epialleles[6] << "\t"
                  << epialleles[7] << "\t"
                  << epialleles[8] << "\t"
                  << epialleles[9] << "\t"
                  << epialleles[10] << "\t"
                  << epialleles[11] << "\t"
                  << epialleles[12] << "\t"
                  << epialleles[13] << "\t"
                  << epialleles[14] << "\t"
                  << epialleles[15] << "\t"
                  << coverage << "\n";
}

// Write record for 'pdr' mode
void write_record_pdr(std::ofstream & output_stream,
                      std::deque<std::string> const & ref_ids,
                      GenomePosition const & pos,
                      std::tuple<num_reads_t, num_discordant_reads_t, sum_transitions_t> const & position_counts,
                      uint32_t const & coverage_filter)
{
    if (std::get<0>(position_counts) < coverage_filter)
        return;

    output_stream << ref_ids[pos.ref_id] << "\t"
                  << pos.start << "\t"
                  << pos.start + 2 << "\t"
                  << calculate_avg_discordance_across_reads(position_counts) << "\t"
                  << calculate_avg_transitions_across_reads(position_counts) << "\t"
                  << std::get<0>(position_counts) << "\n";
}
