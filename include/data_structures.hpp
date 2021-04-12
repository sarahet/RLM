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
// Data structures
// ==========================================================================

#pragma once

#include <vector>

#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>

using seqan3::operator""_tag;

// Define score types
enum class score_type : uint8_t
{
    SINGLE_READ,
    PDR,
    ENTROPY,
    ALL
};

inline score_type
_score_name_to_enum(std::string const & str)
{
    if (str == "single_read")
        return score_type::SINGLE_READ;
    else if (str == "entropy")
        return score_type::ENTROPY;
    else if (str == "pdr")
        return score_type::PDR;
    else if (str == "all")
        return score_type::ALL;

    return score_type::SINGLE_READ;
}

// Define aligner types
enum class align_type : uint8_t
{
    BSMAP,
    BISMARK,
    SEGEMEHL
};

inline align_type
_aligner_name_to_enum(std::string const & str)
{
    if (str == "bsmap")
        return align_type::BSMAP;
    else if (str == "bismark")
        return align_type::BISMARK;
    else if (str == "segemehl")
        return align_type::SEGEMEHL;

    return align_type::BSMAP;
}

// Define sequencing types
enum class sequencing_type : uint8_t
{
    WGBS,
    RRBS
};

inline sequencing_type
_sequencing_type_to_enum(bool const & rrbs)
{
    if (rrbs)
        return sequencing_type::RRBS;

    return sequencing_type::WGBS;
}

// Define read type
enum class read_type : uint8_t
{
    FWD,   // original forward strand
    REV    // original reverse strand
};

// Get original strand from BSMAP alignment
inline read_type
_read_tag_bsmap_to_enum(std::string const & str)
{
    if (str == "++")
        return read_type::FWD;
    else if (str == "+-")
        return read_type::FWD;
    else if (str == "-+")
        return read_type::REV;
    else if (str == "--")
        return read_type::REV;

    return read_type::FWD;
}

// Get original strand from BISMARK alignment
inline read_type
_read_tag_bismark_to_enum(std::string const & str)
{
    if (str == "CT")
        return read_type::FWD;
    else if (str == "GA")
        return read_type::REV;

    return read_type::FWD;
}

// Get original strand from segemehl alignment
inline read_type
_read_tag_segemehl_to_enum(std::string const & str)
{
    if (str == "F1/CT")
        return read_type::FWD;
    else if (str == "F1/GA")
        return read_type::REV;

    return read_type::FWD;
}

// Tag used to separate overloads for different score calculations
template <bool calc_pdr_score, bool calc_entropy_score>
struct score_tag{};

// Overload for ZS tag
template <>
struct seqan3::sam_tag_type<"ZS"_tag>
{
    using type = std::string;
};

// Overload for XG tag
template <>
struct seqan3::sam_tag_type<"XG"_tag>
{
    using type = std::string;
};

// Overload for XB tag
template <>
struct seqan3::sam_tag_type<"XB"_tag>
{
    using type = std::string;
};

// Store a CpG with methylation values from all reads that cover it
struct GenomePosition
{
    uint16_t ref_id;
    uint64_t start;

    inline bool operator== (GenomePosition const & c2) const
    {
         return std::tie(ref_id, start)
             == std::tie(c2.ref_id, c2.start);
    }

    inline bool operator< (GenomePosition const & c2) const
    {
        return ref_id < c2.ref_id || ref_id == c2.ref_id && start < c2.start;
    }
};

// Epiallele configuration used for entropy/epipolymorphism calculations
static constexpr std::array<std::array<std::array<std::array<uint16_t, 2>, 2>, 2>, 2> epiallele_pos_array
{
    [] () constexpr
    {
        std::array<std::array<std::array<std::array<uint16_t, 2>, 2>, 2>, 2> ret{};

        ret[0][0][0][0] = 0;
        ret[0][0][0][1] = 1;
        ret[0][0][1][0] = 2;
        ret[0][0][1][1] = 3;
        ret[0][1][0][0] = 4;
        ret[0][1][0][1] = 5;
        ret[0][1][1][0] = 6;
        ret[0][1][1][1] = 7;
        ret[1][0][0][0] = 8;
        ret[1][0][0][1] = 9;
        ret[1][0][1][0] = 10;
        ret[1][0][1][1] = 11;
        ret[1][1][0][0] = 12;
        ret[1][1][0][1] = 13;
        ret[1][1][1][0] = 14;
        ret[1][1][1][1] = 15;

        return ret;
    }()
};
