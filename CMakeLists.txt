cmake_minimum_required(VERSION 3.18)

set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9")
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(WIN32)
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "")
  link_directories(${Python_LIBRARY_DIRS})
else()
  # we only choose Development.Module to support virtual environments where
  # libpython may not be available
  # see https://github.com/pypa/manylinux/issues/484
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module NumPy)
endif()

execute_process(
  COMMAND ${Python_EXECUTABLE} -c "import versioneer as v; print('.'.join(v.get_version().split('.')[:3]))"
  OUTPUT_VARIABLE PANDAS_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
project(pandas VERSION ${PANDAS_VERSION})

add_compile_definitions(NPY_NO_DEPRECATED_API=0)
add_subdirectory("pandas/_libs")
add_subdirectory("pandas/io/sas")
