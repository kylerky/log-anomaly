cmake_minimum_required (VERSION 2.8.2)
project (journal_test)

include_directories(${SYSTEMD_INCLUDE_DIRS})
add_definitions(-DHAVE_LIBSYSTEMD)

add_executable(journaltest sd_journal_test.cxx)
add_library(journal sd_journal.cxx)

# target_link_libraries(journaltest journal gtest_main ${SYSTEMD_LIBRARIES})
target_link_libraries(journal ${SYSTEMD_LIBRARIES})
target_link_libraries(journaltest journal)
