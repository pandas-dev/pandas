# only add PUBLIC dependencies as well
#   https://cmake.org/cmake/help/latest/manual/cmake-packages.7.html#creating-a-package-configuration-file
include(CMakeFindDependencyMacro)

# Search in <PackageName>_ROOT:
#   https://cmake.org/cmake/help/v3.12/policy/CMP0074.html
if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()

# locate the installed FindABC.cmake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/Modules")

# this section stores which configuration options were set
set(HAVE_THREADS ON)
set(HAVE_IPP )
set(HAVE_ZLIB_NG TRUE)
set(DEACTIVATE_IPP ON)
set(DEACTIVATE_ZLIB OFF)
set(DEACTIVATE_ZSTD OFF)
set(PREFER_EXTERNAL_LZ4 OFF)
set(PREFER_EXTERNAL_ZLIB OFF)
set(PREFER_EXTERNAL_ZSTD OFF)

# find dependencies and their targets, which are used in our Blosc2Targets.cmake
#   additionally, the Blosc2_..._FOUND variables are used to support
#   find_package(Blosc2 ... COMPONENTS ... ...)
#   this enables downstream projects to express the need for specific features.
set(CMAKE_THREAD_PREFER_PTHREAD TRUE) # pre 3.1
set(THREADS_PREFER_PTHREAD_FLAG TRUE) # CMake 3.1+
if(HAVE_THREADS)
    find_dependency(Threads)
    set(Blosc2_THREADS_FOUND TRUE)
else()
    set(Blosc2_THREADS_FOUND FALSE)
endif()

if(NOT DEACTIVATE_IPP AND HAVE_IPP)
    find_dependency(IPP)
    set(Blosc2_IPP_FOUND FALSE)
else()
    set(Blosc2_IPP_FOUND TRUE)
endif()

if(PREFER_EXTERNAL_LZ4)
    find_dependency(LZ4)
endif()
set(Blosc2_LZ4_FOUND TRUE)

if(DEACTIVATE_ZLIB)
    set(Blosc2_ZLIB_FOUND FALSE)
elseif(NOT DEACTIVATE_ZLIB AND PREFER_EXTERNAL_ZLIB)
    if(HAVE_ZLIB_NG)
        find_dependency(ZLIB_NG)
    else()
        find_dependency(ZLIB)
    endif()
    set(Blosc2_ZLIB_FOUND TRUE)
endif()

if(DEACTIVATE_ZSTD)
    set(Blosc2_ZSTD_FOUND FALSE)
elseif(NOT PREFER_EXTERNAL_ZSTD AND PREFER_EXTERNAL_ZSTD)
    find_dependency(ZSTD)
    set(Blosc2_ZSTD_FOUND TRUE)
endif()

# define central Blosc2::blosc2_shared/static targets
include("${CMAKE_CURRENT_LIST_DIR}/Blosc2Targets.cmake")

# check if components are fulfilled and set Blosc2_<COMPONENT>_FOUND vars
#   Blosc2_FIND_COMPONENTS is a list set by find_package(... COMPONENTS ... ...)
#   likewise Blosc2_FIND_REQUIRED_... per component specified
foreach(comp ${Blosc2_FIND_COMPONENTS})
    if(NOT Blosc2_${comp}_FOUND)
        if(Blosc2_FIND_REQUIRED_${comp})
            set(Blosc2_FOUND FALSE)
        endif()
    endif()
endforeach()

