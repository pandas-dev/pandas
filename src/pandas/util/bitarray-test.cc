// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

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
