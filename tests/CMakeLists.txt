

add_executable(
  dormancy_test
  dormancy_test.cpp
)
target_link_libraries(
  dormancy_test
  tskit
  gtest_main
)

include(GoogleTest)
gtest_discover_tests(dormancy_test)

