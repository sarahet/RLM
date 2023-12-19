cmake_minimum_required (VERSION 3.8)

## CUSTOMISE

# Define the application name and version.
project (RLM CXX)

# ----------------------------------------------------------------------------
# App version
# ----------------------------------------------------------------------------

# change this after every release
set (SEQAN_APP_VERSION_MAJOR "1")
set (SEQAN_APP_VERSION_MINOR "2")
set (SEQAN_APP_VERSION_PATCH "0")

# don't change the following
set (SEQAN_APP_VERSION "${SEQAN_APP_VERSION_MAJOR}.${SEQAN_APP_VERSION_MINOR}.${SEQAN_APP_VERSION_PATCH}")

message (STATUS "RLM version is: ${SEQAN_APP_VERSION}")

# ----------------------------------------------------------------------------
# Make "Release" the default cmake build type
# ----------------------------------------------------------------------------

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release RelWithDebInfo"
         FORCE)
endif ()

# ----------------------------------------------------------------------------
# Option static build
# ----------------------------------------------------------------------------

option (RLM_STATIC_BUILD "Include all libraries in the binaries." OFF)

if (RLM_STATIC_BUILD)
    add_definitions (-DRLM_STATIC_BUILD=1)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    # apple does not support fully static builds, but at least libgcc and libstdc++
    if (APPLE)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
        message (WARNING "WARNING: Builds on Mac are never fully static.")
    else (APPLE)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    endif (APPLE)
    # on linux cmake adds -rdynamic automatically which clang can't handle in static builds
    if (CMAKE_SYSTEM_NAME MATCHES "Linux")
        SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
    endif (CMAKE_SYSTEM_NAME MATCHES "Linux")
endif (RLM_STATIC_BUILD)

message(STATUS "The following options are selected for the build:")
message(   "     RLM_STATIC_BUILD             ${RLM_STATIC_BUILD}")
message(STATUS "Run 'cmake -LH' to get a comment on each option.")
message(STATUS "Remove CMakeCache.txt and re-run cmake with -DOPTIONNAME=ON|OFF to change an option.")

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Messages
string (ASCII 27 Esc)
set (FontBold "${Esc}[1m")
set (FontReset "${Esc}[m")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

find_package (SeqAn3 QUIET REQUIRED HINTS lib/seqan3/build_system)
find_package (sharg QUIET REQUIRED HINTS lib/sharg-parser/build_system)
find_package (OpenMP QUIET)
find_package (ZLIB   QUIET)

message(STATUS "These dependencies were found:")
message(   "     OPENMP       ${OPENMP_FOUND}      ${OpenMP_CXX_FLAGS}")
message(   "     SEQAN3       ${SEQAN3_FOUND}      ${SEQAN3_VERSION_STRING}")
message(   "     SHARG        ${SHARG_FOUND}      ${SHARG_VERSION}")
message(   "     ZLIB         ${ZLIB_FOUND}      ${ZLIB_VERSION_STRING}")

# Warn if Zlib was not found.
if (NOT ZLIB_FOUND)
  message (WARNING "WARNING: Zlib not found. Building RLM without support for gzipped input and output (this includes support for .bam).")
endif (NOT ZLIB_FOUND)

message(STATUS "The requirements were met.")

# ----------------------------------------------------------------------------
# Compiler specifics
# ----------------------------------------------------------------------------

# GCC12 and above: Disable warning about std::hardware_destructive_interference_size not being ABI-stable.
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 12)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-interference-size")
    endif ()
endif ()

# ----------------------------------------------------------------------------
# Add RLM targets
# ----------------------------------------------------------------------------

add_subdirectory (src)
message (STATUS "${FontBold}You can run `make` to build the application.${FontReset}")

# ----------------------------------------------------------------------------
# Man-pages
# ----------------------------------------------------------------------------

# Umbrella man-page
add_custom_command (OUTPUT RLM.1
                    COMMAND RLM --export-help man > RLM.1
                    DEPENDS RLM)

add_custom_target (manual ALL DEPENDS RLM.1)

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Adapt to system paths
include (GNUInstallDirs)

# Install RLM binaries into LIBEXECDIR
install (TARGETS RLM
         DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install non-binary files for the package to DOCDIR, usually ${PREFIX}/share/doc/RLM
install (FILES LICENSE.md
               README.md
         DESTINATION ${CMAKE_INSTALL_DOCDIR})

# Man pages into MANDIR, usually ${PREFIX}/share/man/man1 (or without share)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/RLM.1
         DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

# Information
set (CPACK_PACKAGE_NAME "RLM")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "RLM - Read level DNA methylation analysis of bisulfite converted sequencing data.")
set (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.md")
set (CPACK_PACKAGE_VENDOR "Sara Hetzel <hetzel@molgen.mpg.de>")
set (CPACK_PACKAGE_CONTACT "${CPACK_PACKAGE_VENDOR}")
set (CPACK_PACKAGE_VERSION_MAJOR "${SEQAN_APP_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${SEQAN_APP_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${SEQAN_APP_VERSION_PATCH}")
set (CPACK_PACKAGE_VERSION "${SEQAN_APP_VERSION}")
set (CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_NAME} ${CPACK_PACKAGE_VERSION}")

# Package format(s)
if (CMAKE_SYSTEM_NAME MATCHES "Windows")
    set(CPACK_GENERATOR "ZIP;NSIS")
elseif (CMAKE_SYSTEM_NAME MATCHES "Darwin")
    set(CPACK_GENERATOR "ZIP;DragNDrop")
elseif (CMAKE_VERSION VERSION_LESS "3.1") # TXZ support since 3.1
    set(CPACK_GENERATOR "TBZ2")
else()
    set(CPACK_GENERATOR "TXZ")
endif ()

if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    set(CPACK_GENERATOR "${CPACK_GENERATOR};DEB;RPM")
endif ()

# Package architecture
if (CMAKE_SYSTEM_PROCESSOR STREQUAL "amd64")
    set(CMAKE_SYSTEM_PROCESSOR "x86_64")
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
endif ()

# Include architecture in package name
if (NOT DEFINED CPACK_SYSTEM_NAME)
  set(CPACK_SYSTEM_NAME "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
endif (NOT DEFINED CPACK_SYSTEM_NAME)

include (CPack)

# ----------------------------------------------------------------------------
# Add tests
# ----------------------------------------------------------------------------

enable_testing ()
add_subdirectory (test EXCLUDE_FROM_ALL)
