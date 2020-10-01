cmake_minimum_required(VERSION 3.17)
project(cholesky)

set(CMAKE_CXX_STANDARD 14)

add_executable(cholesky main.cpp matrix.h matrix.cpp)