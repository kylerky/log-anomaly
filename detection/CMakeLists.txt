cmake_minimum_required (VERSION 2.8.2)
project (detection)

add_subdirectory(stemmer)

add_executable(construct construct.cxx)
target_link_libraries(construct infomap stemmer journal Eigen3::Eigen)
