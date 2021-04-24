#include <string>

#include "cli_test.hpp"

TEST_F(RLM, no_options)
{
    cli_test_result result = execute_app("RLM");
    std::string expected
    {
        "RLM - Read level DNA methylation analysis of bisulfite converted sequencing data.\n"
        "=================================================================================\n"
        "    Try -h or --help for more information.\n"
    };

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(RLM, missing_bam_file)
{
    cli_test_result result = execute_app("RLM", "-r", data("chrM.fa"), "-m", "PE", "-s", "single_read", "-a", "bsmap");

    std::string expected
    {
        "Parsing error. Option -b/--bam is required but not set.\n"
    };

    EXPECT_EQ(result.exit_code, 0xFF00);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}
