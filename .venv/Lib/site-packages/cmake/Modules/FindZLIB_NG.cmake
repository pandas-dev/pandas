find_path(ZLIB_NG_INCLUDE_DIR NAMES zlib-ng.h)

if(ZLIB_INCLUDE_DIRS)
  set(ZLIB_NG_LIBRARY_DIRS ${ZLIB_NG_INCLUDE_DIR})

  if("${ZLIB_NG_LIBRARY_DIRS}" MATCHES "/include$")
    # Strip off the trailing "/include" in the path.
    GET_FILENAME_COMPONENT(ZLIB_NG_LIBRARY_DIRS ${ZLIB_NG_LIBRARY_DIRS} PATH)
  endif("${ZLIB_NG_LIBRARY_DIRS}" MATCHES "/include$")

  if(EXISTS "${ZLIB_NG_LIBRARY_DIRS}/lib")
    set(ZLIB_NG_LIBRARY_DIRS ${ZLIB_NG_LIBRARY_DIRS}/lib)
  endif(EXISTS "${ZLIB_NG_LIBRARY_DIRS}/lib")
endif()

find_library(ZLIB_NG_LIBRARY NAMES z-ng libz-ng zlib-ng libz-ng.a)

set(ZLIB_NG_LIBRARIES ${ZLIB_NG_LIBRARY})
set(ZLIB_NG_INCLUDE_DIR ${ZLIB_NG_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ZLIB_NG DEFAULT_MSG ZLIB_NG_LIBRARY ZLIB_NG_INCLUDE_DIR)

if(ZLIB_NG_INCLUDE_DIR AND ZLIB_NG_LIBRARIES)
  set(ZLIB_NG_FOUND TRUE)
else(ZLIB_NG_INCLUDE_DIR AND ZLIB_NG_LIBRARIES)
  set(ZLIB_NG_FOUND FALSE)
endif(ZLIB_NG_INCLUDE_DIR AND ZLIB_NG_LIBRARIES)

# generate CMake target
if(ZLIB_NG_FOUND)
  message(STATUS "Found zlib-ng: ${ZLIB_NG_LIBRARIES}, ${ZLIB_NG_INCLUDE_DIR}")

  set(ZLIB_NG_INCLUDE_DIRS ${ZLIB_NG_INCLUDE_DIR})

  if(NOT TARGET ZLIB_NG::ZLIB_NG)
    add_library(ZLIB_NG::ZLIB_NG UNKNOWN IMPORTED)
    set_target_properties(ZLIB_NG::ZLIB_NG PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${ZLIB_NG_INCLUDE_DIRS}")

    set_property(TARGET ZLIB_NG::ZLIB_NG APPEND PROPERTY
      IMPORTED_LOCATION "${ZLIB_NG_LIBRARIES}")
  endif()
endif()

#[[
Copyright https://github.com/zlib-ng/minizip-ng, 2021

Condition of use and distribution are the same as zlib:

This software is provided 'as-is', without any express or implied
warranty.  In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgement in the product documentation would be
   appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
]]#
