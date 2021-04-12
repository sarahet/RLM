#include <vector>

#include <gtest/gtest.h>

#include "../../include/methylation_scores.hpp"

TEST(scores, entropy)
{
    std::vector<uint32_t> vec1{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double entropy1 = calculate_entropy_across_reads(vec1);
    EXPECT_EQ(static_cast<double>(0), entropy1);

    std::vector<uint32_t> vec2{0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double entropy2 = calculate_entropy_across_reads(vec2);
    EXPECT_EQ(static_cast<double>(0), entropy2);

    std::vector<uint32_t> vec3(16, 1);
    double entropy3 = calculate_entropy_across_reads(vec3);
    EXPECT_EQ(static_cast<double>(1), entropy3);

    std::vector<uint32_t> vec4{0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    double entropy4 = calculate_entropy_across_reads(vec4);
    EXPECT_EQ(static_cast<double>(0.25), entropy4);
}

TEST(scores, epipolymorphism)
{
    std::vector<uint32_t> vec1{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double epipolymorphism1 = calculate_epipolymorphism_across_reads(vec1);
    EXPECT_EQ(static_cast<double>(0), epipolymorphism1);

    std::vector<uint32_t> vec2{0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double epipolymorphism2 = calculate_epipolymorphism_across_reads(vec2);
    EXPECT_EQ(static_cast<double>(0), epipolymorphism2);

    std::vector<uint32_t> vec3(16, 1);
    double epipolymorphism3 = calculate_epipolymorphism_across_reads(vec3);
    EXPECT_EQ(static_cast<double>(0.9375), epipolymorphism3);

    std::vector<uint32_t> vec4{0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    double epipolymorphism4 = calculate_epipolymorphism_across_reads(vec4);
    EXPECT_EQ(static_cast<double>(0.5), epipolymorphism4);
}

TEST(scores, read_transitions)
{
    std::vector<uint16_t> vec1{1, 0, 0, 0, 0};
    double transitions1 = calculate_transitions_per_read(vec1);
    EXPECT_EQ(static_cast<double>(0.25), transitions1);

    std::vector<uint16_t> vec2{0, 0, 1, 0, 0};
    double transitions2 = calculate_transitions_per_read(vec2);
    EXPECT_EQ(static_cast<double>(0.5), transitions2);

    std::vector<uint16_t> vec3(5, 0);
    double transitions3 = calculate_transitions_per_read(vec3);
    EXPECT_EQ(static_cast<double>(0), transitions3);

    std::vector<uint16_t> vec4(5, 1);
    double transitions4 = calculate_transitions_per_read(vec4);
    EXPECT_EQ(static_cast<double>(0), transitions4);

    std::vector<uint16_t> vec5{1, 0, 1, 0, 1};
    double transitions5 = calculate_transitions_per_read(vec5);
    EXPECT_EQ(static_cast<double>(1), transitions5);
}

TEST(scores, read_discordance)
{
    std::vector<uint16_t> vec1{1, 0, 0, 0, 0};
    uint16_t discordance1 = calculate_discordance_per_read(vec1);
    EXPECT_EQ(static_cast<uint16_t>(1), discordance1);

    std::vector<uint16_t> vec2{0, 0, 1, 0, 0};
    uint16_t discordance2 = calculate_discordance_per_read(vec2);
    EXPECT_EQ(static_cast<uint16_t>(1), discordance2);

    std::vector<uint16_t> vec3(5, 0);
    uint16_t discordance3 = calculate_discordance_per_read(vec3);
    EXPECT_EQ(static_cast<uint16_t>(0), discordance3);

    std::vector<uint16_t> vec4(5, 1);
    uint16_t discordance4 = calculate_discordance_per_read(vec4);
    EXPECT_EQ(static_cast<uint16_t>(0), discordance4);

    std::vector<uint16_t> vec5{1, 0, 1, 0, 1};
    uint16_t discordance5 = calculate_discordance_per_read(vec5);
    EXPECT_EQ(static_cast<uint16_t>(1), discordance5);
}
