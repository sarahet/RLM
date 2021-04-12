cmake_minimum_required (VERSION 3.8)

# An object library (without main) to be used in multiple targets.
# add_library ("${PROJECT_NAME}_lib" STATIC fastq_conversion.cpp)
# target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC seqan3::seqan3)
# target_include_directories ("${PROJECT_NAME}_lib" PUBLIC ../include)

add_executable ("${PROJECT_NAME}" main.cpp)
target_link_libraries ("${PROJECT_NAME}" PUBLIC seqan3::seqan3)
