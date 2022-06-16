cmake_minimum_required(VERSION 3.21)

project(pandas VERSION 1.4)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_SHARED_LIBRARY_PREFIX "")
set(CMAKE_SHARED_LIBRARY_SUFFIX .cpython-38-x86_64-linux-gnu.so)
set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9")

find_package(Python REQUIRED COMPONENTS Development NumPy)

add_compile_definitions(NPY_NO_DEPRECATED_API=0)
add_subdirectory("pandas/_libs")
add_subdirectory("pandas/io/sas")
