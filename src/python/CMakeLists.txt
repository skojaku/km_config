cmake_minimum_required(VERSION 3.0)

project(km_config)
set(PYBIND11_CPP_STANDARD -std=c++11)
set(CMAKE_CXX_FLAGS "-O3 -fopenmp")
add_subdirectory(pybind11)

pybind11_add_module(km_config SHARED km_config.cpp)  
