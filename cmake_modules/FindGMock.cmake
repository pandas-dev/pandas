# Copyright (c) 2009-2010 Volvox Development Team
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# Author: Konstantin Lepa <konstantin.lepa@gmail.com>
#
# Find the Google Mock Framework, heavily cribbed from FindGTest.cmake.
# gmock ships a copy of gtest and bundles it in its libraries, so this also
# finds the gtest headers.
#
# This module defines
# GMOCK_INCLUDE_DIR, where to find gmock include files, etc.
# GTEST_INCLUDE_DIR, where to find gtest include files
# GMock_FOUND, If false, do not try to use gmock.
# GMOCK_STATIC_LIBRARY, Location of libgmock.a
# GMOCK_SHARED_LIBRARY, Location of libttest's shared library

# also defined, but not for general use are
# GMOCK_LIBRARY, where to find the GMock library.

set(GMOCK_SEARCH_PATH $ENV{GOOGLETEST_PREFIX})

set(GMOCK_H gmock/gmock.h)
set(GTEST_H gtest/gtest.h)

find_path(GMOCK_INCLUDE_DIR ${GMOCK_H}
  PATHS ${GMOCK_SEARCH_PATH}/include
        NO_DEFAULT_PATH
  DOC   "Path to the ${GMOCK_H} file"
)
find_path(GTEST_INCLUDE_DIR ${GTEST_H}
  PATHS ${GMOCK_SEARCH_PATH}/include
        NO_DEFAULT_PATH
  DOC   "Path to the ${GTEST_H} file"
)
find_library(GMOCK_LIBRARY
  NAMES gmock
  PATHS ${GMOCK_SEARCH_PATH}/lib/
        NO_DEFAULT_PATH
  DOC   "Google's framework for writing C++ tests (gmock)"
)

set(GMOCK_LIB_NAME libgmock)
if(GMOCK_INCLUDE_DIR AND GTEST_INCLUDE_DIR AND GMOCK_LIBRARY)
  set(GMOCK_STATIC_LIBRARY ${GMOCK_SEARCH_PATH}/lib/${GMOCK_LIB_NAME}.a)
  if(EXISTS "${GMOCK_STATIC_LIBRARY}")
    set(GMOCK_FOUND TRUE)
  endif()
else()
  set(GMOCK_FOUND FALSE)
endif()

if(GMOCK_FOUND)
  if(NOT GMock_FIND_QUIETLY)
    message(STATUS "Found the GMock library: ${GMOCK_STATIC_LIBRARY}")
  endif(NOT GMock_FIND_QUIETLY)
else(GMOCK_FOUND)
  if(NOT GMock_FIND_QUIETLY)
    if(GMock_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find the GMock library")
    else(GMock_FIND_REQUIRED)
      message(STATUS "Could not find the GMock library")
    endif(GMock_FIND_REQUIRED)
  endif(NOT GMock_FIND_QUIETLY)
endif(GMOCK_FOUND)

mark_as_advanced(
  GMOCK_INCLUDE_DIR
  GTEST_INCLUDE_DIR
  GMOCK_STATIC_LIBRARY
)
