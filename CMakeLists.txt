cmake_minimum_required(VERSION 3.16)
project(syncmer-toolkit C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(cmake/syng)
add_subdirectory(src)
