# - Find GFLAGS (gflags.h, libgflags.a, libgflags.so, and libgflags.so.0)
# This module defines
#  GFLAGS_INCLUDE_DIR, directory containing headers
#  GFLAGS_LIBS, directory containing gflag libraries
#  GFLAGS_STATIC_LIB, path to libgflags.a
#  GFLAGS_SHARED_LIB, path to libgflags' shared library
#  GFLAGS_FOUND, whether gflags has been found

set(GFLAGS_SEARCH_HEADER_PATHS
  $ENV{GFLAGS_PREFIX}/include
)

set(GFLAGS_SEARCH_LIB_PATH
  $ENV{GFLAGS_PREFIX}/lib
)

find_path(GFLAGS_INCLUDE_DIR gflags/gflags.h PATHS
  ${GFLAGS_SEARCH_HEADER_PATHS}
  # make sure we don't accidentally pick up a different version
  NO_DEFAULT_PATH
)

find_library(GFLAGS_LIB_PATH NAMES gflags PATHS ${GFLAGS_SEARCH_LIB_PATH} NO_DEFAULT_PATH)

if (GFLAGS_INCLUDE_DIR AND GFLAGS_LIB_PATH)
  set(GFLAGS_FOUND TRUE)
  set(GFLAGS_LIB_NAME libgflags)
  set(GFLAGS_LIBS ${GFLAGS_SEARCH_LIB_PATH})
  set(GFLAGS_STATIC_LIB ${GFLAGS_SEARCH_LIB_PATH}/${GFLAGS_LIB_NAME}.a)
  set(GFLAGS_SHARED_LIB ${GFLAGS_LIBS}/${GFLAGS_LIB_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX})
else ()
  set(GFLAGS_FOUND FALSE)
endif ()

if (GFLAGS_FOUND)
  if (NOT GFlags_FIND_QUIETLY)
    message(STATUS "Found the GFlags library: ${GFLAGS_LIB_PATH}")
  endif ()
else ()
  if (NOT GFlags_FIND_QUIETLY)
    set(GFLAGS_ERR_MSG "Could not find the GFlags library. Looked for headers")
    set(GFLAGS_ERR_MSG "${GFLAGS_ERR_MSG} in ${GFLAGS_SEARCH_HEADER_PATHS}, and for libs")
    set(GFLAGS_ERR_MSG "${GFLAGS_ERR_MSG} in ${GFLAGS_SEARCH_LIB_PATH}")
    if (GFlags_FIND_REQUIRED)
      message(FATAL_ERROR "${GFLAGS_ERR_MSG}")
    else (GFlags_FIND_REQUIRED)
      message(STATUS "${GFLAGS_ERR_MSG}")
    endif (GFlags_FIND_REQUIRED)
  endif ()
endif ()

mark_as_advanced(
  GFLAGS_INCLUDE_DIR
  GFLAGS_LIBS
  GFLAGS_STATIC_LIB
  GFLAGS_SHARED_LIB
)
