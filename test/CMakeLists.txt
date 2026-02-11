# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

cmake_minimum_required (VERSION 3.20)

# Set directories for test output files, input data and binaries.
file (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)
add_definitions (-DOUTPUTDIR=\"${CMAKE_CURRENT_BINARY_DIR}/output/\")
add_definitions (-DDATADIR=\"${CMAKE_CURRENT_BINARY_DIR}/data/\")
add_definitions (-DBINDIR=\"${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/\")

# This includes `cmake/test/config.cmake` which takes care of setting up the test infrastructure. It also provides
# the `add_app_test` macro, which is used to add tests to the CMake test suite.
include("${PROJECT_SOURCE_DIR}/cmake/test/config.cmake")

# This makes test data available
include (data/datasources.cmake)

# If you previously had subfolders that call add_api_test/add_cli_test, keep them:
add_subdirectory (api)
add_subdirectory (cli)

message (STATUS "${FontBold}You can run `make check` to build and run tests.${FontReset}")
