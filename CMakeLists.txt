cmake_minimum_required (VERSION 3.8)

## CUSTOMISE

# Define the application name and version.
project (RLM VERSION 1.2.0)

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

# Dependency: SeqAn3, sharg-parser.
find_package (SeqAn3 QUIET REQUIRED HINTS lib/seqan3/build_system)
find_package (sharg QUIET REQUIRED HINTS lib/sharg-parser/build_system)

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
