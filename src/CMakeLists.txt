cmake_minimum_required (VERSION 3.8)

add_executable ("${PROJECT_NAME}" RLM.cpp)
target_link_libraries ("${PROJECT_NAME}" PUBLIC seqan3::seqan3 sharg::sharg)
