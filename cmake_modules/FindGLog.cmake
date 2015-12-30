# - Find GLOG (logging.h, libglog.a, libglog.so, and libglog.so.0)
# This module defines
#  GLOG_INCLUDE_DIR, directory containing headers
#  GLOG_LIBS, directory containing glog libraries
#  GLOG_STATIC_LIB, path to libglog.a
#  GLOG_SHARED_LIB, path to libglog's shared library
#  GLOG_FOUND, whether glog has been found

set(GLOG_SEARCH_HEADER_PATHS $ENV{GLOG_PREFIX}/include)
set(GLOG_SEARCH_LIB_PATH $ENV{GLOG_PREFIX}/lib)

find_path(GLOG_INCLUDE_DIR glog/logging.h PATHS
  ${GLOG_SEARCH_HEADER_PATHS}
  # make sure we don't accidentally pick up a different version
  NO_DEFAULT_PATH
)

find_library(GLOG_LIB_PATH NAMES glog PATHS ${GLOG_SEARCH_LIB_PATH} NO_DEFAULT_PATH)

if (GLOG_INCLUDE_DIR AND GLOG_LIB_PATH)
  set(GLOG_FOUND TRUE)
  set(GLOG_LIBS ${GLOG_SEARCH_LIB_PATH})
  set(GLOG_LIB_NAME libglog)
  set(GLOG_STATIC_LIB ${GLOG_SEARCH_LIB_PATH}/${GLOG_LIB_NAME}.a)
  set(GLOG_SHARED_LIB ${GLOG_LIBS}/${GLOG_LIB_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX})
else ()
  set(GLOG_FOUND FALSE)
endif ()

if (GLOG_FOUND)
  if (NOT GLog_FIND_QUIETLY)
    message(STATUS "Found the GLog library: ${GLOG_LIB_PATH}")
  endif ()
else ()
  if (NOT GLog_FIND_QUIETLY)
    set(GLOG_ERR_MSG "Could not find the GLog library. Looked for headers")
    set(GLOG_ERR_MSG "${GLOG_ERR_MSG} in ${GLOG_SEARCH_HEADER_PATHS}, and for libs")
    set(GLOG_ERR_MSG "${GLOG_ERR_MSG} in ${GLOG_SEARCH_LIB_PATH}")
    if (GLog_FIND_REQUIRED)
      message(FATAL_ERROR "${GLOG_ERR_MSG}")
    else (GLog_FIND_REQUIRED)
      message(STATUS "${GLOG_ERR_MSG}")
    endif (GLog_FIND_REQUIRED)
  endif ()
endif ()

mark_as_advanced(
  GLOG_INCLUDE_DIR
  GLOG_LIBS
  GLOG_STATIC_LIB
  GLOG_SHARED_LIB
)
