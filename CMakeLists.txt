cmake_minimum_required(VERSION 3.17)

project(pandas)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_SOURCE_DIR .)
set(CMAKE_BINARY_DIR .)

set(CMAKE_SHARED_LIBRARY_PREFIX "")
# for macOS
if (APPLE)
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")
endif()

find_package(Python REQUIRED COMPONENTS Interpreter Development NumPy)
message("${Python_LINK_OPTIONS}")
set(CMAKE_SHARED_LIBRARY_SUFFIX ".${Python_SOABI}.so")
add_compile_definitions(NPY_NO_DEPRECATED_API=0)
add_subdirectory("pandas/_libs")
add_subdirectory("pandas/io/sas")
