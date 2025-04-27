cmake_minimum_required(VERSION 3.9)
project(reccs)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Point to your local installation of igraph
if(DEFINED ENV{igraph_DIR})
  set(igraph_DIR $ENV{igraph_DIR})
else()
  set(igraph_DIR $ENV{HOME}/local/lib/cmake/igraph)
endif()

# Add the library directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)

# Find the installed igraph package
find_package(igraph 0.10.0 REQUIRED)

# Add igraph-cpp
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extlib/igraph-cpp/include)

# Add NetworKit as a submodule
add_subdirectory(extlib/networkit)

# Add NetworKit include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extlib/networkit/include)

# Add the executable
add_executable(reccs app/reccs.cpp)

# Link against igraph and NetworKit
target_link_libraries(reccs PRIVATE 
    igraph::igraph
    networkit
)

# Optional: Disable NetworKit tests and benchmarks to speed up build
set(NETWORKIT_BUILD_TESTS OFF CACHE BOOL "Build NetworKit tests")
set(NETWORKIT_BUILD_BENCHMARKS OFF CACHE BOOL "Build NetworKit benchmarks")
