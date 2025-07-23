cmake_minimum_required(VERSION 3.9)
project(reccs)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Ensure std::filesystem is available
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")

# Find OpenMP package
find_package(OpenMP REQUIRED)

# Add the library directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)

# Build TLX library separately to avoid circular dependencies
add_subdirectory(extlib/VieCut/extlib/tlx)

include_directories(
    extlib/VieCut/app 
    extlib/VieCut/lib             
    extlib/VieCut/extlib/growt
    extlib/VieCut/extlib/tlx
)

# Add nlohmann/json
add_subdirectory(extlib/json)

# Create shared library components for common functionality
# This helps avoid code duplication and ensures both executables use the same libraries
set(RECCS_LIBRARIES
    OpenMP::OpenMP_CXX 
    nlohmann_json::nlohmann_json
    tlx  # Link against TLX library
    # stdc++fs  # Add filesystem library
)

# Add the original reccs executable
add_executable(reccs app/reccs.cpp)
target_link_libraries(reccs PUBLIC ${RECCS_LIBRARIES})

# Add the new reccspp executable
add_executable(reccspp app/reccspp.cpp)
target_link_libraries(reccspp PUBLIC ${RECCS_LIBRARIES})

# Copy extlib directory to build directory
file(COPY ${CMAKE_SOURCE_DIR}/extlib DESTINATION ${CMAKE_BINARY_DIR})
