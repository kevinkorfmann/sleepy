find_package(Boost 1.40 COMPONENTS program_options REQUIRED)

add_executable(sleepy sleepy.cpp)
target_link_libraries(sleepy tskit Boost::program_options) 