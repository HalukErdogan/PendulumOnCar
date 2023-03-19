cmake_minimum_required(VERSION 3.10)
project(Simulation)

# Set the build type
set(CMAKE_BUILD_TYPE Release)

# Set C++ standard to C++14
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find Eigen
find_package(Eigen3 REQUIRED)

# Add include directory
include_directories(include ${EIGEN3_INCLUDE_DIR})

# Add libraries
add_library(DynamicModel SHARED 
    include/DynamicModel/DynamicModelBase.h
    include/DynamicModel/PendulumOnCar/PendulumOnCar.h
)

set_target_properties(DynamicModel PROPERTIES LINKER_LANGUAGE CXX)

add_library(IntegrationMethod SHARED 
    include/IntegrationMethod/EulerIntegration/EulerIntegration.h
    include/IntegrationMethod/IntegrationMethodBase.h
    include/IntegrationMethod/RKDPIntegration/RKDPIntegration.h
    include/IntegrationMethod/RungeKuttaIntegration/RungeKuttaIntegration.h
)

set_target_properties(IntegrationMethod PROPERTIES LINKER_LANGUAGE CXX)

add_library(ODESolver SHARED 
    include/ODESolver/ODESolverBase.h
    include/ODESolver/SimpleODESolver/SimpleODESolver.h
)

set_target_properties(ODESolver PROPERTIES LINKER_LANGUAGE CXX)

add_library(TrajectoryOptimizer SHARED 
    include/TrajectoryOptimizer/DDP/DDP.h
    include/TrajectoryOptimizer/TrajectoryOptimizerBase.h
)

set_target_properties(TrajectoryOptimizer PROPERTIES LINKER_LANGUAGE CXX)

# Add sources
set(SOURCES src/main.cpp)

# Add executables
add_executable(simulation ${SOURCES})

# link the libraries
target_link_libraries(simulation 
    DynamicModel 
    IntegrationMethod 
    ODESolver 
    TrajectoryOptimizer 
    -lm
)

# Add tests
enable_testing()
add_subdirectory(tests)
