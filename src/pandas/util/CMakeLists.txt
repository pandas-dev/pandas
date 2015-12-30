# Copyright 2014 Cloudera, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# pandas_util

set(UTIL_SRCS
)

set(UTIL_LIBS
  glog
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
