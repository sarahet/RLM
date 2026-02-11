cmake_minimum_required (VERSION 3.20)

## CUSTOMISE

# Define the application name and version.
project (RLM VERSION 1.2.0)

if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR
        "RLM currently supports GCC only. Detected: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
endif ()

# C++20 (required by SeqAn3)
set (CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

## BUILD

# Make Release default build type
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release RelWithDebInfo"
         FORCE)
endif ()

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Messages
string (ASCII 27 Esc)
set (FontBold "${Esc}[1m")
set (FontReset "${Esc}[m")

# Dependencies via CPM
include (${CMAKE_CURRENT_SOURCE_DIR}/cmake/CPM.cmake)

# Pin versions
CPMAddPackage(
  NAME seqan3
  GITHUB_REPOSITORY seqan/seqan3
  GIT_TAG 3.4.0
)

CPMAddPackage(
  NAME sharg-parser
  GITHUB_REPOSITORY seqan/sharg-parser
  GIT_TAG 1.2.1
)

CPMAddPackage(
  NAME googletest
  GITHUB_REPOSITORY google/googletest
  GIT_TAG v1.15.2
  OPTIONS "INSTALL_GTEST OFF"
)

# GCC12 and above: Disable warning about std::hardware_destructive_interference_size not being ABI-stable.
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 12)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-interference-size")
    endif ()
endif ()

add_subdirectory (src)
message (STATUS "${FontBold}You can run `make` to build the application.${FontReset}")

## TEST

enable_testing ()
add_subdirectory (test EXCLUDE_FROM_ALL)
