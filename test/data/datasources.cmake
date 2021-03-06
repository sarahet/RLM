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

declare_datasource (FILE test_clipped_reads.bam
                    URL ${CMAKE_SOURCE_DIR}/test/data/test_clipped_reads.bam
                    URL_HASH SHA256=0d4b3d695555bec54d35cea011fa977e4e678769501782375e960fd9344c5caa)

declare_datasource (FILE control_clipped_reads.bed
                    URL ${CMAKE_SOURCE_DIR}/test/data/control_clipped_reads.bed
                    URL_HASH SHA256=b049b9de05bfe92eb6473314bc7fa2b471bc09c1876fe47d68045997d9dadfa5)

declare_datasource (FILE test_invalid_mate.bam
                    URL ${CMAKE_SOURCE_DIR}/test/data/test_invalid_mate.bam
                    URL_HASH SHA256=448032cfede9644f66993a1be4f0de7760b33a521a44b63dab8cd28c58758632)

declare_datasource (FILE control_invalid_mate.bed
                    URL ${CMAKE_SOURCE_DIR}/test/data/control_invalid_mate.bed
                    URL_HASH SHA256=4d3e2a922a0982bb722973d49802f462377bf41e14962b5b7fccb6d09dfb69ee)

declare_datasource (FILE test_ref.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/test_ref.fa
                    URL_HASH SHA256=0d0e4d6e8e7a429e0517e63655fb570ecab40b7c088d824f08c1617c00161fda)

declare_datasource (FILE test_single_reads.bam
                    URL ${CMAKE_SOURCE_DIR}/test/data/test_single_reads.bam
                    URL_HASH SHA256=8e9627d7ffa26d848cf00045a41d9591815f6e71dd80cde496fe9a964b8b1109)

declare_datasource (FILE test_single_reads_name_sorted.bam
                    URL ${CMAKE_SOURCE_DIR}/test/data/test_single_reads_name_sorted.bam
                    URL_HASH SHA256=2ccbac9ae58399d93bb7b2546c44d8adc0ef9f0d6421ac19504b1aaed4492cd8)

declare_datasource (FILE control_single_reads.bed
                    URL ${CMAKE_SOURCE_DIR}/test/data/control_single_reads.bed
                    URL_HASH SHA256=1199da4cd754cee67168e6578709379fdc9405c0e55834d8eebbc8d199672c9d)

declare_datasource (FILE control_single_reads_rrbs.bed
                    URL ${CMAKE_SOURCE_DIR}/test/data/control_single_reads_rrbs.bed
                    URL_HASH SHA256=0bc0cb0079db972a3650ddb2c1e79f7c213813aa9bd8b8c4726a2edf93bcb157)

declare_datasource (FILE test_bsmap.bam
                    URL ${CMAKE_SOURCE_DIR}/test/data/test_bsmap.bam
                    URL_HASH SHA256=d9817ccb301fccf4dcb83f3fe715daa68a189a337ea504d6145c88310db4fe8e)

declare_datasource (FILE test_segemehl.bam
                    URL ${CMAKE_SOURCE_DIR}/test/data/test_segemehl.bam
                    URL_HASH SHA256=21b1aa9b1f7ca6cd701cd118908dea7400325d6977cba6b5a43fb6fa1befc9fd)

declare_datasource (FILE test_bismark.bam
                    URL ${CMAKE_SOURCE_DIR}/test/data/test_bismark.bam
                    URL_HASH SHA256=ec426d6d992d054b5073c5d3ea7bb48cb1c702b64d73d767a02a8bd23edc45c3)

declare_datasource (FILE test_gem.bam
                    URL ${CMAKE_SOURCE_DIR}/test/data/test_gem.bam
                    URL_HASH SHA256=41f1fcf96096db12361962c67e6caa1f998ef16efe6f110da92d700a613f7e4e)
