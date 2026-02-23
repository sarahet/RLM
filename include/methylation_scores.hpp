// ==========================================================================
//                                  RLM
// ==========================================================================
// Copyright (c) 2021-2026, Sara Hetzel <hetzel @ molgen.mpg.de>
// Copyright (c) 2021-2026, Max-Planck-Institut für Molekulare Genetik
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

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <tuple>
#include <vector>


#include "data_structures.hpp"

using num_reads_t = uint32_t;
using num_discordant_reads_t = uint32_t;
using sum_transitions_t = double;
using num_methyl_cpgs_t = uint32_t;

// Calculate transition score of a single read
inline double calculate_transitions_per_read(std::vector<uint16_t> const & cpg_config)
{
    uint16_t transitions = 0;

    if (cpg_config.size() > 1)
    {
        for (size_t i = 0; i < (cpg_config.size() - 1); i++)
        {
            if (cpg_config[i] != cpg_config[i+1])
                transitions++;
        }
        return static_cast<double>(transitions) / (cpg_config.size() - 1);
    }
    else return 0;
}

// Calculate discordance of a single read
inline uint16_t calculate_discordance_per_read(std::vector<uint16_t> const & cpg_config)
{
    uint16_t transitions = 0;
    
    if (cpg_config.size() > 1)
    {
        for (size_t i = 0; i < (cpg_config.size() - 1); i++)
        {
            if (cpg_config[i] != cpg_config[i+1])
                transitions++;
        }
        return transitions == 0 ? 0 : 1;
    }
    else return 0;
}

// Calculate average RTS for a CpG
inline double calculate_avg_transitions_across_reads(std::tuple<num_reads_t, num_discordant_reads_t, sum_transitions_t, num_methyl_cpgs_t> const & position_counts)
{
    assert(std::get<0>(position_counts) > 0);
    return std::get<2>(position_counts) / std::get<0>(position_counts);
}

// Calculate average discordance for a CpG
inline double calculate_avg_discordance_across_reads(std::tuple<num_reads_t, num_discordant_reads_t, sum_transitions_t, num_methyl_cpgs_t> const & position_counts)
{
    assert(std::get<0>(position_counts) > 0);
    return static_cast<double>(std::get<1>(position_counts)) / std::get<0>(position_counts);
}

// Calculate average methylation for a CpG
inline double calculate_avg_methylation_across_reads(std::tuple<num_reads_t, num_discordant_reads_t, sum_transitions_t, num_methyl_cpgs_t> const & position_counts)
{
    assert(std::get<0>(position_counts) > 0);
    return static_cast<double>(std::get<3>(position_counts)) / std::get<0>(position_counts);
}

// Calculate entropy for a 4-mer
inline double calculate_entropy_across_reads(std::vector<uint32_t> const & epialleles, uint32_t num_reads)
{
    assert(num_reads > 0);

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
inline double calculate_epipolymorphism_across_reads(std::vector<uint32_t> const & epialleles, uint32_t num_reads)
{
    assert(num_reads > 0);

    double epipolymorphism = 0;

    for (size_t i = 0; i < epialleles.size(); i++)
    {
        epipolymorphism += std::pow(static_cast<double>(epialleles[i]) / num_reads, 2);
    }

    epipolymorphism = 1 - epipolymorphism;
    return epipolymorphism;
}

// Calculate average methylation for a 4-mer
inline double calculate_avg_kmer_methylation_across_reads(std::vector<uint32_t> const & epialleles, uint32_t num_reads)
{
    assert(num_reads > 0);
    assert(epialleles.size() == 16);

    uint32_t methylated_cpgs =
    epialleles[1] * 1 +
    epialleles[2] * 1 +
    epialleles[3] * 2 +
    epialleles[4] * 1 +
    epialleles[5] * 2 +
    epialleles[6] * 2 +
    epialleles[7] * 3 +
    epialleles[8] * 1 +
    epialleles[9] * 2 +
    epialleles[10] * 2 +
    epialleles[11] * 3 +
    epialleles[12] * 2 +
    epialleles[13] * 3 +
    epialleles[14] * 3 +
    epialleles[15] * 4;

    return static_cast<double>(methylated_cpgs) / (num_reads * 4);
}
