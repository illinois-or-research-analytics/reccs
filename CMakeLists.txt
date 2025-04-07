cmake_minimum_required(VERSION 3.9)
project(reccs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the library directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)

# Add the executable
add_executable(reccs app/reccs.cpp)
