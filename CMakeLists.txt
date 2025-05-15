cmake_minimum_required(VERSION 3.9)
project(reccs)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find OpenMP package
find_package(OpenMP REQUIRED)

# Add the library directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)

# Add the executable
add_executable(reccs app/reccs.cpp)

# Add the loader test executable
# add_executable(loader_test tests/loader_test.cpp)

# Link with OpenMP
target_link_libraries(reccs PUBLIC OpenMP::OpenMP_CXX)
# target_link_libraries(loader_test PUBLIC OpenMP::OpenMP_CXX)

# Copy extlib directory to build directory
file(COPY ${CMAKE_SOURCE_DIR}/extlib DESTINATION ${CMAKE_BINARY_DIR})
