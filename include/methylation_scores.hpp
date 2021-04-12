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
// Scores for read level methylation analysis
// ==========================================================================

#pragma once

#include <cmath>

#include "data_structures.hpp"

using num_reads_t = uint32_t;
using num_discordant_reads_t = uint32_t;
using sum_transitions_t = double;

// Calculate transirion score of a single read
double calculate_transitions_per_read(std::vector<uint16_t> const & cpg_config)
{
    uint16_t transitions = 0;
    for (size_t i = 0; i < cpg_config.size() - 1; i++)
    {
        if (cpg_config[i] != cpg_config[i+1])
            transitions++;
    }
    return static_cast<double>(transitions) / (cpg_config.size() - 1);
}

// Calculate discordance of a single read
uint16_t calculate_discordance_per_read(std::vector<uint16_t> const & cpg_config)
{
    uint16_t transitions = 0;
    for (size_t i = 0; i < cpg_config.size() - 1; i++)
    {
        if (cpg_config[i] != cpg_config[i+1])
            transitions++;
    }
    return transitions == 0 ? 0 : 1;
}

// Calculate average RTS for a CpG
double calculate_avg_transitions_across_reads(std::tuple<num_reads_t, num_discordant_reads_t, sum_transitions_t> const & position_counts)
{
    return std::get<2>(position_counts) / std::get<0>(position_counts);
}

// Calculate average discordance for a CpG
double calculate_avg_discordance_across_reads(std::tuple<num_reads_t, num_discordant_reads_t, sum_transitions_t> const & position_counts)
{
    return static_cast<double>(std::get<1>(position_counts)) / std::get<0>(position_counts);
}

// Calculate entropy for a 4-mer
double calculate_entropy_across_reads(std::vector<uint32_t> const & epialleles)
{
    uint32_t num_reads = std::accumulate(epialleles.begin(), epialleles.end(), 0);
    double entropy = 0;

    for (size_t i = 0; i < epialleles.size(); i++)
    {
        if (epialleles[i] != 0)
            entropy += (-(static_cast<double>(epialleles[i]) / num_reads) * std::log2(static_cast<double>(epialleles[i]) / num_reads));
    }

    entropy = entropy / 4;
    return entropy;
}

// Calculate epipolymorphism for a 4-mer
double calculate_epipolymorphism_across_reads(std::vector<uint32_t> const & epialleles)
{
    uint32_t num_reads = std::accumulate(epialleles.begin(), epialleles.end(), 0);
    double epipolymorphism = 0;

    for (size_t i = 0; i < epialleles.size(); i++)
    {
        epipolymorphism += std::pow(static_cast<double>(epialleles[i]) / num_reads, 2);
    }

    epipolymorphism = 1 - epipolymorphism;
    return epipolymorphism;
}
