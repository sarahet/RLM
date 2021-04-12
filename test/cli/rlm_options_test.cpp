#include <seqan3/std/ranges>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>

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
