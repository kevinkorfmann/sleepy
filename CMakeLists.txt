project(sleepy)
cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_FLAGS "-O3")

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin/")
set(tsk ./tskit/c/tskit)

include_directories(./kastore/c)

add_library(kastore SHARED ./kastore/c/kastore.c)

include_directories(./tskit/c)
add_library(tskit SHARED 
    ${tsk}/convert.c 
    ${tsk}/core.c 
    ${tsk}/genotypes.c
    ${tsk}/stats.c
    ${tsk}/tables.c
    ${tsk}/trees.c)

target_link_libraries(tskit kastore)

add_subdirectory(src)