cmake_minimum_required (VERSION 2.8.2)
project (detection)

add_subdirectory(stemmer)

include_directories(${PROJECT_SOURCE_DIR}/io)

add_executable(construct construct.cxx)
target_link_libraries(construct infomap stemmer journal Eigen3::Eigen)

find_package(Curses REQUIRED)

include_directories(${CURSES_INCLUDE_DIRS})

add_executable(check check.cxx)
target_link_libraries(check ${CURSES_LIBRARIES} journal)
