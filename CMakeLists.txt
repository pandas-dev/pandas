cmake_minimum_required(VERSION 3.21)

project(pandas VERSION 1.4)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

set(CMAKE_SHARED_LIBRARY_PREFIX "")
set(CMAKE_SHARED_LIBRARY_SUFFIX .cpython-38-darwin.so)

find_package(Python REQUIRED COMPONENTS Development NumPy)
add_subdirectory("pandas/_libs")
add_subdirectory("pandas/io/sas")
