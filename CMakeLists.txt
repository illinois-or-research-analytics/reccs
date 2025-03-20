cmake_minimum_required(VERSION 3.9)
project(reccs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the library directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)

# If you have source files in the lib directory that need to be compiled
# file(GLOB LIB_SOURCES "lib/*.cpp")
# add_library(reccs_lib ${LIB_SOURCES})
# target_link_libraries(reccs reccs_lib)

# Add the executable
add_executable(reccs app/reccs.cpp)

# Include directories if needed
# include_directories(include)

# Link libraries if needed
# target_link_libraries(reccs some_library)

# Optional: Set output directory
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Optional: Install configuration
# install(TARGETS reccs DESTINATION bin)

# Optional: Enable testing
# enable_testing()
# add_test(NAME reccs_test COMMAND reccs)

# Find required packages
find_package(Boost REQUIRED COMPONENTS graph random)
find_package(PkgConfig REQUIRED)

# Find graph-tool using pkg-config
pkg_check_modules(GRAPHTOOL REQUIRED graph-tool-2.0)

# Include directories
include_directories(
    ${Boost_INCLUDE_DIRS}
    ${GRAPHTOOL_INCLUDE_DIRS}
)

# Add executable
add_executable(sbm_generator sbm.cpp)

# Link libraries
target_link_libraries(sbm_generator
    ${Boost_LIBRARIES}
    ${GRAPHTOOL_LIBRARIES}
)

# Set compiler flags
target_compile_options(sbm_generator PRIVATE
    -O3
    -Wall
    -Wextra
    ${GRAPHTOOL_CFLAGS_OTHER}
)
