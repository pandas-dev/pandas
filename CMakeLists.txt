cmake_minimum_required(VERSION 3.21)

project(pandas VERSION 1.4)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

find_package(Python REQUIRED COMPONENTS Development NumPy)
add_subdirectory("pandas/_libs")
add_subdirectory("pandas/io/sas")
