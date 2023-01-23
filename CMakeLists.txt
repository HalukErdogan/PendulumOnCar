cmake_minimum_required(VERSION 3.10)
project(Simulation)

# Set C++ standard to C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find Eigen
find_package(Eigen3 REQUIRED)

# Add include directory
include_directories(include ${EIGEN3_INCLUDE_DIR})

# Add sources
set(SOURCES src/main.cpp)

# Add executables
add_executable(simulation ${SOURCES})

# link the library 
target_link_libraries(simulation -lm)
