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
// Functions needed to process a single BAM record
// ==========================================================================

#pragma once

#include "methylation_scores.hpp"

using seqan3::operator""_dna5;
using num_reads_t = uint32_t;
using num_discordant_reads_t = uint32_t;
using sum_transitions_t = double;

// Find positions of all CpGs in a sequence
template <typename ref_type>
std::vector<uint16_t> find_cpg_pos(ref_type const & reference)
{
    std::vector<uint16_t> occurrencs;
    seqan3::dna5_vector pattern = "CG"_dna5;
    auto cpg_pos = std::search(reference.begin(), reference.end(), pattern.begin(), pattern.end());
    while (cpg_pos != reference.end())
    {
        occurrencs.push_back(cpg_pos - reference.begin());
        ++cpg_pos;
        cpg_pos = std::search(cpg_pos, reference.end(), pattern.begin(), pattern.end());
    }
    return occurrencs;
}

// Insert CpG into map to store it until all BAM records are read
void insert_CpG(size_t const & reference_id,
                size_t const & reference_position,
                size_t const & offset,
                std::map<GenomePosition, std::tuple<num_reads_t, num_discordant_reads_t, sum_transitions_t> > & all_CpGs,
                std::vector<uint16_t> const & cpg_pos,
                std::vector<uint16_t> const & cpg_config)
{
    for (size_t i; i < cpg_pos.size(); i++)
    {
        GenomePosition pos;
        pos.ref_id = reference_id;
        pos.start = reference_position + offset + cpg_pos[i];

        auto it = all_CpGs.find(pos);

        if (it != all_CpGs.end())
        {
            std::get<0>(it->second)++;
            std::get<1>(it->second) += cpg_config[i];
            std::get<2>(it->second) += calculate_transitions_per_read(cpg_config);
        }
        else
        {
            all_CpGs.insert(std::make_pair(pos, std::make_tuple(1, cpg_config[i], calculate_transitions_per_read(cpg_config))));
        }
    }
}

// Insert kmer into map to store it until all BAM records are read
void insert_kmer(size_t const & reference_id,
                 size_t const & reference_position,
                 size_t const & offset,
                 std::map<GenomePosition, std::vector<uint32_t> > & all_kmers,
                 std::vector<uint16_t> const & cpg_pos,
                 std::vector<uint16_t> const & cpg_config)
{
    if (cpg_pos.size() < 4)
        return;

    for (size_t i; i < (cpg_pos.size() - 3); i++)
    {
        GenomePosition pos;
        pos.ref_id = reference_id;
        pos.start = reference_position + offset + cpg_pos[i];

        auto it = all_kmers.find(pos);

        if (it != all_kmers.end())
        {
            (it->second)[epiallele_pos_array[cpg_config[i]][cpg_config[i+1]][cpg_config[i+2]][cpg_config[i+3]]]++;
        }
        else
        {
            std::vector<uint32_t> epialleles(16, 0);
            epialleles[epiallele_pos_array[cpg_config[i]][cpg_config[i+1]][cpg_config[i+2]][cpg_config[i+3]]]++;
            all_kmers.insert(std::make_pair(pos, epialleles));
        }
    }
}

// Internal function to process a single BAM record
void process_bam_record_impl(std::ofstream & output_stream,
                             read_type const & tag,
                             size_t const & reference_id,
                             size_t const & reference_position,
                             seqan3::dna5_vector const & sequence,
                             std::string const & id,
                             std::deque<std::string> const & ref_ids,
                             std::vector<seqan3::dna5_vector> const & genome_seqs,
                             std::vector<uint16_t> & cpg_pos,
                             std::vector<uint16_t> & cpg_config)
{
    // Define look-up for methylated or unmethylated CpGs (depending on base that needs to be evaluated).
    static constexpr std::array<uint16_t, 5> methyl_context = {0, 1, 1, 0, 0};
    static constexpr std::array<char, 2> methyl_context_char = {'g', 'G'};

    // Extract reference sequence matching the reads.
    // For reads coming from the forward strand, the reference position is shifted up by one.
    // For reads coming from the reverse strand, the reference position is shifted down by one.
    seqan3::dna5_vector ref_sequence;
    if (tag == read_type::REV)
    {
        ref_sequence = genome_seqs[reference_id]
                     | seqan3::views::slice(reference_position - 1, reference_position + sequence.size() - 1)
                     | seqan3::views::to<seqan3::dna5_vector>;
    }
    else
    {
        ref_sequence = genome_seqs[reference_id]
                     | seqan3::views::slice(reference_position + 1, reference_position + sequence.size() + 1)
                     | seqan3::views::to<seqan3::dna5_vector>;
    }

    // Find all CpG positions
    cpg_pos = find_cpg_pos(ref_sequence);

    if (cpg_pos.size() < 3)
        return;

    // For every CpG determine unmethylated/methylated status.
    // For reads coming from the forward strand, the position of the 'C' needs to be evaluated.
    // For reads coming from the reverse strand, the position of the 'G' needs to be evaluated.
    bool skip = false;
    for (size_t i = 0; i < cpg_pos.size(); i++)
    {
        if (tag == read_type::REV)
        {
            if (sequence[cpg_pos[i]] == 'A'_dna5 || sequence[cpg_pos[i]] == 'G'_dna5)
            {
                cpg_config.push_back(methyl_context[sequence[cpg_pos[i]].to_rank()]);
            }
            else
            {
                skip = true;
                break;
            }
        }
        else
        {
            if (sequence[cpg_pos[i] + 1] == 'C'_dna5 || sequence[cpg_pos[i] + 1] == 'T'_dna5)
            {
                cpg_config.push_back(methyl_context[sequence[cpg_pos[i] + 1].to_rank()]);
            }
            else
            {
                skip = true;
                break;
            }
        }
    }
    if (skip)
        return;

    // Prepare output
    uint16_t num_methyl_cpgs = std::accumulate(cpg_config.begin(), cpg_config.end(), 0);

    output_stream << ref_ids[reference_id] << "\t"
                  << reference_position << "\t"
                  << reference_position + sequence.size() << "\t"
                  << id << "\t";

    for (auto i : cpg_config)
        output_stream << methyl_context_char[i];

    output_stream << "\t"
                  << cpg_config.size() << "\t"
                  << num_methyl_cpgs << "\t"
                  << calculate_discordance_per_read(cpg_config) << "\t"
                  << calculate_transitions_per_read(cpg_config) << "\t"
                  << static_cast<double>(num_methyl_cpgs) / cpg_config.size() << "\n";
}

// Outer wrapper function overload for single read score only
void process_bam_record(std::ofstream & output_stream,
                        read_type const & tag,
                        size_t const & reference_id,
                        size_t const & reference_position,
                        seqan3::dna5_vector const & sequence,
                        std::string const & id,
                        std::deque<std::string> const & ref_ids,
                        std::vector<seqan3::dna5_vector> const & genome_seqs,
                        std::map<GenomePosition, std::tuple<num_reads_t, num_discordant_reads_t, sum_transitions_t> > & all_CpGs,
                        std::map<GenomePosition, std::vector<uint32_t> > & all_kmers,
                        score_tag<false, false>)
{
    std::vector<uint16_t> cpg_pos;
    std::vector<uint16_t> cpg_config;

    process_bam_record_impl(output_stream,
                            tag,
                            reference_id,
                            reference_position,
                            sequence,
                            id,
                            ref_ids,
                            genome_seqs,
                            cpg_pos,
                            cpg_config);
}

// Outer wrapper function overload for PDR/RTS scores
void process_bam_record(std::ofstream & output_stream,
                        read_type const & tag,
                        size_t const & reference_id,
                        size_t const & reference_position,
                        seqan3::dna5_vector const & sequence,
                        std::string const & id,
                        std::deque<std::string> const & ref_ids,
                        std::vector<seqan3::dna5_vector> const & genome_seqs,
                        std::map<GenomePosition, std::tuple<num_reads_t, num_discordant_reads_t, sum_transitions_t> > & all_CpGs,
                        std::map<GenomePosition, std::vector<uint32_t> > & all_kmers,
                        score_tag<true, false>)
{
    std::vector<uint16_t> cpg_pos;
    std::vector<uint16_t> cpg_config;

    process_bam_record_impl(output_stream,
                            tag,
                            reference_id,
                            reference_position,
                            sequence,
                            id,
                            ref_ids,
                            genome_seqs,
                            cpg_pos,
                            cpg_config);

    int16_t offset = tag == read_type::REV ? -1 : 1;
    insert_CpG(reference_id, reference_position, offset, all_CpGs, cpg_pos, cpg_config);
}

// Outer wrapper function overload for entropy/epipolymorphism scores
void process_bam_record(std::ofstream & output_stream,
                        read_type const & tag,
                        size_t const & reference_id,
                        size_t const & reference_position,
                        seqan3::dna5_vector const & sequence,
                        std::string const & id,
                        std::deque<std::string> const & ref_ids,
                        std::vector<seqan3::dna5_vector> const & genome_seqs,
                        std::map<GenomePosition, std::tuple<num_reads_t, num_discordant_reads_t, sum_transitions_t> > & all_CpGs,
                        std::map<GenomePosition, std::vector<uint32_t> > & all_kmers,
                        score_tag<false, true>)
{
    std::vector<uint16_t> cpg_pos;
    std::vector<uint16_t> cpg_config;

    process_bam_record_impl(output_stream,
                            tag,
                            reference_id,
                            reference_position,
                            sequence,
                            id,
                            ref_ids,
                            genome_seqs,
                            cpg_pos,
                            cpg_config);

    int16_t offset = tag == read_type::REV ? -1 : 1;
    insert_kmer(reference_id, reference_position, offset, all_kmers, cpg_pos, cpg_config);
}

// Outer wrapper function overload for all scores
void process_bam_record(std::ofstream & output_stream,
                        read_type const & tag,
                        size_t const & reference_id,
                        size_t const & reference_position,
                        seqan3::dna5_vector const & sequence,
                        std::string const & id,
                        std::deque<std::string> const & ref_ids,
                        std::vector<seqan3::dna5_vector> const & genome_seqs,
                        std::map<GenomePosition, std::tuple<num_reads_t, num_discordant_reads_t, sum_transitions_t> > & all_CpGs,
                        std::map<GenomePosition, std::vector<uint32_t> > & all_kmers,
                        score_tag<true, true>)
{
    std::vector<uint16_t> cpg_pos;
    std::vector<uint16_t> cpg_config;

    process_bam_record_impl(output_stream,
                            tag,
                            reference_id,
                            reference_position,
                            sequence,
                            id,
                            ref_ids,
                            genome_seqs,
                            cpg_pos,
                            cpg_config);

    int16_t offset = tag == read_type::REV ? -1 : 1;
    insert_CpG(reference_id, reference_position, offset, all_CpGs, cpg_pos, cpg_config);
    insert_kmer(reference_id, reference_position, offset, all_kmers, cpg_pos, cpg_config);
}
