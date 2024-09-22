cmake_minimum_required(VERSION 3.29.0)

project(benchmark CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../dct ${CMAKE_CURRENT_BINARY_DIR}/dct)

set(executable benchmark)
add_executable(${executable} main.cpp)
target_include_directories(${executable} PRIVATE ${dct_SOURCE_DIR}/include)
target_link_libraries(${executable} PRIVATE dct)
