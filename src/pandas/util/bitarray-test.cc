// Copyright 2015 Cloudera, Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <gtest/gtest.h>

#include "pandas/status.h"
#include "pandas/test-util.h"
#include "pandas/util/bitarray.h"

namespace pandas {

TEST(BitArrayTests, TestBasics) {
  BitArray arr;
  size_t length = 100;

  ASSERT_OK(arr.Init(length));

  ASSERT_EQ(0, arr.set_count());
  arr.Set(5);
  ASSERT_TRUE(arr.IsSet(5));
  ASSERT_EQ(1, arr.set_count());

  ASSERT_FALSE(arr.IsSet(0));

  arr.Unset(5);

  ASSERT_FALSE(arr.IsSet(5));
}

} // namespace pandas
