// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#include <cstdlib>
#include <vector>

#include <gtest/gtest.h>

#include "pandas/util.h"

namespace pandas {

TEST(UtilTests, TestNextPower2) {
  using util::next_power2;

  ASSERT_EQ(8, next_power2(6));
  ASSERT_EQ(8, next_power2(8));

  ASSERT_EQ(1, next_power2(1));
  ASSERT_EQ(256, next_power2(131));

  ASSERT_EQ(1024, next_power2(1000));

  ASSERT_EQ(4096, next_power2(4000));

  ASSERT_EQ(65536, next_power2(64000));

  ASSERT_EQ(1ULL << 32, next_power2((1ULL << 32) - 1));
  ASSERT_EQ(1ULL << 31, next_power2((1ULL << 31) - 1));
  ASSERT_EQ(1ULL << 63, next_power2((1ULL << 63) - 1));
}

} // namespace pandas
