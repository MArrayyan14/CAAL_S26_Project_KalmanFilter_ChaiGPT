cmake_minimum_required(VERSION 3.10)
project(KalmanFilter_Milestone2 VERSION 1.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen3
find_package(Eigen3 REQUIRED NO_MODULE)

# Linear Kalman Filter executable
add_executable(lkf lkf.cpp)
target_link_libraries(lkf Eigen3::Eigen)

# Extended Kalman Filter executable
add_executable(ekf ekf.cpp)
target_link_libraries(ekf Eigen3::Eigen)

# Compiler optimizations
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(lkf PRIVATE -Wall -Wextra -O3 -march=native)
    target_compile_options(ekf PRIVATE -Wall -Wextra -O3 -march=native)
endif()

# Display information
message(STATUS "==================================")
message(STATUS "Kalman Filter - Milestone 2")
message(STATUS "Team ChaiGPT")
message(STATUS "==================================")
message(STATUS "Eigen3 Version: ${EIGEN3_VERSION}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "==================================")
