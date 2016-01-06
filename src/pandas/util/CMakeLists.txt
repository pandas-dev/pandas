# This file is a part of pandas. See LICENSE for details about reuse and
# copyright holders

# pandas_util

set(UTIL_SRCS
  bitarray.cc
  status.cc
)

set(UTIL_LIBS
)

add_library(pandas_util STATIC
  ${UTIL_SRCS}
)
target_link_libraries(pandas_util ${UTIL_LIBS})
SET_TARGET_PROPERTIES(pandas_util PROPERTIES LINKER_LANGUAGE CXX)

# pandas_test_util

add_library(pandas_test_util)
target_link_libraries(pandas_test_util
  pandas_util)

SET_TARGET_PROPERTIES(pandas_test_util PROPERTIES LINKER_LANGUAGE CXX)

# pandas_test_main
add_library(pandas_test_main
  test_main.cc
  ${THIRDPARTY_DIR}/gtest/gtest-all.cc)

target_link_libraries(pandas_test_main
  pandas_util
  pandas_test_util
  python
  dl)

if(APPLE)
  set_target_properties(pandas_test_main
		PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif()

ADD_PANDAS_TEST(bitarray-test)
