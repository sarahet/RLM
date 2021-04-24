#include <string>
#include <fstream>

#include "cli_test.hpp"

TEST_F(RLM, skipped_reads)
{
    cli_test_result result = execute_app("RLM", "-b", data("test_skipped_reads.sam"), "-r", data("chrM.fa"), "-m", "PE", "-s", "single_read", "-a", "bsmap");

    std::ifstream output ("output_single_read_info.bed");
    std::ifstream control (data("control_skipped_reads.bed"));

    std::string line;
    std::vector<std::string> output_vec;
    std::vector<std::string> control_vec;

    while (std::getline(output, line))
    {
        output_vec.push_back(line);
    }
    output.close();

    while (std::getline(control, line))
    {
        control_vec.push_back(line);
    }
    control.close();

    EXPECT_RANGE_EQ(output_vec, control_vec);

    for (size_t i = 0; i < output_vec.size(); i++)
        EXPECT_EQ(output_vec[i], control_vec[i]);
}

TEST_F(RLM, overlap_reads)
{
    cli_test_result result = execute_app("RLM", "-b", data("test_overlap_reads.bam"), "-r", data("chrM.fa"), "-m", "PE", "-s", "single_read", "-a", "bsmap");

    std::ifstream output ("output_single_read_info.bed");
    std::ifstream control (data("control_overlap_reads.bed"));

    std::string line;
    std::vector<std::string> output_vec;
    std::vector<std::string> control_vec;

    while (std::getline(output, line))
    {
        output_vec.push_back(line);
    }
    output.close();

    while (std::getline(control, line))
    {
        control_vec.push_back(line);
    }
    control.close();

    EXPECT_RANGE_EQ(output_vec, control_vec);

    for (size_t i = 0; i < output_vec.size(); i++)
        EXPECT_EQ(output_vec[i], control_vec[i]);
}

TEST_F(RLM, single_reads)
{
    cli_test_result result = execute_app("RLM", "-b", data("test_single_reads.bam"), "-r", data("test_ref.fa"), "-m", "SE", "-s", "single_read", "-a", "bsmap");

    std::ifstream output ("output_single_read_info.bed");
    std::ifstream control (data("control_single_reads.bed"));

    std::string line;
    std::vector<std::string> output_vec;
    std::vector<std::string> control_vec;

    while (std::getline(output, line))
    {
        output_vec.push_back(line);
    }
    output.close();

    while (std::getline(control, line))
    {
        control_vec.push_back(line);
    }
    control.close();

    EXPECT_RANGE_EQ(output_vec, control_vec);

    for (size_t i = 0; i < output_vec.size(); i++)
        EXPECT_EQ(output_vec[i], control_vec[i]);
}
