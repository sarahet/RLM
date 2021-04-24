cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

declare_datasource (FILE chrM.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/chrM.fa
                    URL_HASH SHA256=d9fd5e2d3240e9ea0b14cde42f89598bffe7448a65ad9329f5910dba1c29d577)

declare_datasource (FILE test_skipped_reads.sam
                    URL ${CMAKE_SOURCE_DIR}/test/data/test_skipped_reads.sam
                    URL_HASH SHA256=46a4a96b4242a8494ddf29efc209533f14a6860d77d4883acc3921b2cf14b487)

declare_datasource (FILE control_skipped_reads.bed
                    URL ${CMAKE_SOURCE_DIR}/test/data/control_skipped_reads.bed
                    URL_HASH SHA256=6cfd3fce1faf3d48466a86763e831576c194bf278dc022bc56948378d33da942)

declare_datasource (FILE test_overlap_reads.bam
                    URL ${CMAKE_SOURCE_DIR}/test/data/test_overlap_reads.bam
                    URL_HASH SHA256=890113bdd885aea3a9d4b5fd3379f4455d7382543aa23c265fd7c019cee78660)

declare_datasource (FILE control_overlap_reads.bed
                    URL ${CMAKE_SOURCE_DIR}/test/data/control_overlap_reads.bed
                    URL_HASH SHA256=e13dbeb88efde592b44fb5657c4e7bbb63111fce78789900b2822307ec1fde4b)

declare_datasource (FILE test_ref.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/test_ref.fa
                    URL_HASH SHA256=0d0e4d6e8e7a429e0517e63655fb570ecab40b7c088d824f08c1617c00161fda)

declare_datasource (FILE test_single_reads.bam
                    URL ${CMAKE_SOURCE_DIR}/test/data/test_single_reads.bam
                    URL_HASH SHA256=8e9627d7ffa26d848cf00045a41d9591815f6e71dd80cde496fe9a964b8b1109)

declare_datasource (FILE control_single_reads.bed
                    URL ${CMAKE_SOURCE_DIR}/test/data/control_single_reads.bed
                    URL_HASH SHA256=1199da4cd754cee67168e6578709379fdc9405c0e55834d8eebbc8d199672c9d)
