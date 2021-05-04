#include <string>
#include <fstream>

#include "cli_test.hpp"

TEST_F(RLM, compare_bsmap_bismark)
{
    cli_test_result result_bsmap = execute_app("RLM", "-b", data("test_bsmap.bam"), "-r", data("test_ref.fa"), "-m", "SE", "-s", "single_read", "-a", "bsmap", "-o", "output_bsmap.bed");
    cli_test_result result_bismark = execute_app("RLM", "-b", data("test_bismark.bam"), "-r", data("test_ref.fa"), "-m", "SE", "-s", "single_read", "-a", "bismark", "-o", "output_bismark.bed");

    std::ifstream output_bsmap ("output_bsmap.bed");
    std::ifstream output_bismark ("output_bismark.bed");

    std::string line;
    std::string field;

    std::vector<std::vector<std::string> > output_vec_bsmap;
    std::vector<std::vector<std::string> > output_vec_bismark;

    while (std::getline(output_bsmap, line))
    {
        std::vector<std::string> current_line;
        std::istringstream iss(line);
        while(std::getline(iss, field, '\t'))
            current_line.push_back(field);
        output_vec_bsmap.push_back(current_line);
    }
    output_bsmap.close();

    while (std::getline(output_bismark, line))
    {
        std::vector<std::string> current_line;
        std::istringstream iss(line);
        while(std::getline(iss, field, '\t'))
            current_line.push_back(field);
        output_vec_bismark.push_back(current_line);
    }
    output_bismark.close();

    for (size_t i = 0; i < output_vec_bsmap.size(); i++)
    {
        for (size_t j = 0; j < output_vec_bismark.size(); j++)
        {
            if ((output_vec_bsmap[i][0] == output_vec_bismark[j][0]) &&
                (output_vec_bsmap[i][1] == output_vec_bismark[j][1]) &&
                (output_vec_bsmap[i][2] == output_vec_bismark[j][2]) &&
                (output_vec_bsmap[i][3] == output_vec_bismark[j][3]))
            {
                EXPECT_EQ(output_vec_bsmap[i][4], output_vec_bismark[j][4]);
                EXPECT_EQ(output_vec_bsmap[i][5], output_vec_bismark[j][5]);
                EXPECT_EQ(output_vec_bsmap[i][6], output_vec_bismark[j][6]);
                EXPECT_EQ(output_vec_bsmap[i][7], output_vec_bismark[j][7]);
                EXPECT_EQ(output_vec_bsmap[i][8], output_vec_bismark[j][8]);
                EXPECT_EQ(output_vec_bsmap[i][9], output_vec_bismark[j][9]);
            }
        }
    }

}

TEST_F(RLM, compare_bsmap_segemehl)
{
    cli_test_result result_bsmap = execute_app("RLM", "-b", data("test_bsmap.bam"), "-r", data("test_ref.fa"), "-m", "SE", "-s", "single_read", "-a", "bsmap", "-o", "output_bsmap.bed");
    cli_test_result result_segemehl = execute_app("RLM", "-b", data("test_segemehl.bam"), "-r", data("test_ref.fa"), "-m", "SE", "-s", "single_read", "-a", "segemehl", "-o", "output_segemehl.bed");

    std::ifstream output_bsmap ("output_bsmap.bed");
    std::ifstream output_segemehl ("output_segemehl.bed");

    std::string line;
    std::string field;

    std::vector<std::vector<std::string> > output_vec_bsmap;
    std::vector<std::vector<std::string> > output_vec_segemehl;

    while (std::getline(output_bsmap, line))
    {
        std::vector<std::string> current_line;
        std::istringstream iss(line);
        while(std::getline(iss, field, '\t'))
            current_line.push_back(field);
        output_vec_bsmap.push_back(current_line);
    }
    output_bsmap.close();

    while (std::getline(output_segemehl, line))
    {
        std::vector<std::string> current_line;
        std::istringstream iss(line);
        while(std::getline(iss, field, '\t'))
            current_line.push_back(field);
        output_vec_segemehl.push_back(current_line);
    }
    output_segemehl.close();

    for (size_t i = 0; i < output_vec_bsmap.size(); i++)
    {
        for (size_t j = 0; j < output_vec_segemehl.size(); j++)
        {
            if ((output_vec_bsmap[i][0] == output_vec_segemehl[j][0]) &&
                (output_vec_bsmap[i][1] == output_vec_segemehl[j][1]) &&
                (output_vec_bsmap[i][2] == output_vec_segemehl[j][2]) &&
                (output_vec_bsmap[i][3] == output_vec_segemehl[j][3]))
            {
                EXPECT_EQ(output_vec_bsmap[i][4], output_vec_segemehl[j][4]);
                EXPECT_EQ(output_vec_bsmap[i][5], output_vec_segemehl[j][5]);
                EXPECT_EQ(output_vec_bsmap[i][6], output_vec_segemehl[j][6]);
                EXPECT_EQ(output_vec_bsmap[i][7], output_vec_segemehl[j][7]);
                EXPECT_EQ(output_vec_bsmap[i][8], output_vec_segemehl[j][8]);
                EXPECT_EQ(output_vec_bsmap[i][9], output_vec_segemehl[j][9]);
            }
        }
    }

}

TEST_F(RLM, compare_bismark_segemehl)
{
    cli_test_result result_bsmap = execute_app("RLM", "-b", data("test_bismark.bam"), "-r", data("test_ref.fa"), "-m", "SE", "-s", "single_read", "-a", "bismark", "-o", "output_bismark.bed");
    cli_test_result result_segemehl = execute_app("RLM", "-b", data("test_segemehl.bam"), "-r", data("test_ref.fa"), "-m", "SE", "-s", "single_read", "-a", "segemehl", "-o", "output_segemehl.bed");

    std::ifstream output_bismark ("output_bismark.bed");
    std::ifstream output_segemehl ("output_segemehl.bed");

    std::string line;
    std::string field;

    std::vector<std::vector<std::string> > output_vec_bismark;
    std::vector<std::vector<std::string> > output_vec_segemehl;

    while (std::getline(output_bismark, line))
    {
        std::vector<std::string> current_line;
        std::istringstream iss(line);
        while(std::getline(iss, field, '\t'))
            current_line.push_back(field);
        output_vec_bismark.push_back(current_line);
    }
    output_bismark.close();

    while (std::getline(output_segemehl, line))
    {
        std::vector<std::string> current_line;
        std::istringstream iss(line);
        while(std::getline(iss, field, '\t'))
            current_line.push_back(field);
        output_vec_segemehl.push_back(current_line);
    }
    output_segemehl.close();

    for (size_t i = 0; i < output_vec_bismark.size(); i++)
    {
        for (size_t j = 0; j < output_vec_segemehl.size(); j++)
        {
            if ((output_vec_bismark[i][0] == output_vec_segemehl[j][0]) &&
                (output_vec_bismark[i][1] == output_vec_segemehl[j][1]) &&
                (output_vec_bismark[i][2] == output_vec_segemehl[j][2]) &&
                (output_vec_bismark[i][3] == output_vec_segemehl[j][3]))
            {
                EXPECT_EQ(output_vec_bismark[i][4], output_vec_segemehl[j][4]);
                EXPECT_EQ(output_vec_bismark[i][5], output_vec_segemehl[j][5]);
                EXPECT_EQ(output_vec_bismark[i][6], output_vec_segemehl[j][6]);
                EXPECT_EQ(output_vec_bismark[i][7], output_vec_segemehl[j][7]);
                EXPECT_EQ(output_vec_bismark[i][8], output_vec_segemehl[j][8]);
                EXPECT_EQ(output_vec_bismark[i][9], output_vec_segemehl[j][9]);
            }
        }
    }

}
