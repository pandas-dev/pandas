# This file is a part of pandas. See LICENSE for details about reuse and
# copyright holders

# pandas_util

set(UTIL_SRCS
  bitarray.cc
)

set(UTIL_LIBS
  rt)


add_library(pandas_util STATIC
  ${UTIL_SRCS}
)
target_link_libraries(pandas_util ${UTIL_LIBS})
SET_TARGET_PROPERTIES(pandas_util PROPERTIES LINKER_LANGUAGE CXX)

# pandas_test_util

add_library(pandas_test_util)
target_link_libraries(pandas_test_util
  gflags
  glog
  gmock
  pandas_util)

SET_TARGET_PROPERTIES(pandas_test_util PROPERTIES LINKER_LANGUAGE CXX)

# pandas_test_main

add_library(pandas_test_main
  test_main.cc)

target_link_libraries(pandas_test_main
  gflags
  glog
  gmock
  pandas_util
  pandas_test_util
  dl
  rt)

ADD_PANDAS_TEST(bitarray-test)
