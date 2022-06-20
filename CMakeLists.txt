cmake_minimum_required(VERSION 3.17)

set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9")
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(Python REQUIRED COMPONENTS Interpreter Development NumPy)

execute_process(
  COMMAND ${Python_EXECUTABLE} -c "import versioneer as v; print(str('.'.join(v.get_version().split('.')[:3])))"
  OUTPUT_VARIABLE PANDAS_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
project(pandas VERSION ${PANDAS_VERSION})

add_compile_definitions(NPY_NO_DEPRECATED_API=0)
add_subdirectory("pandas/_libs")
add_subdirectory("pandas/io/sas")
