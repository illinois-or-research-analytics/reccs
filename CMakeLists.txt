cmake_minimum_required(VERSION 3.9)
project(reccs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
