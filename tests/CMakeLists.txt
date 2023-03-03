cmake_minimum_required(VERSION 3.10)

# Set C++ standard to C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add testing framework
enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Add sources for tests
set(TEST_SOURCES 
    ./test_dynamic_model_base.cpp 
    # add more test source files as needed
)

# Add test executable
add_executable(test_simulation ${TEST_SOURCES})

# Link with the testing framework and Eigen
target_link_libraries(test_simulation ${GTEST_LIBRARIES} Eigen3::Eigen pthread)

# Add test
add_test(NAME test_simulation COMMAND test_simulation)
