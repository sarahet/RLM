cmake_minimum_required (VERSION 3.20)

add_library (${PROJECT_NAME}_lib INTERFACE)
target_include_directories (${PROJECT_NAME}_lib INTERFACE "${${PROJECT_NAME}_SOURCE_DIR}/include")
target_link_libraries (${PROJECT_NAME}_lib INTERFACE seqan3::seqan3 sharg::sharg)

add_executable ("${PROJECT_NAME}" RLM.cpp)
target_link_libraries (${PROJECT_NAME} PRIVATE ${PROJECT_NAME}_lib)
