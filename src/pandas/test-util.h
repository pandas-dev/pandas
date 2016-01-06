// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#ifndef PANDAS_TEST_UTIL_H_
#define PANDAS_TEST_UTIL_H_

#include <gtest/gtest.h>

#define ASSERT_RAISES(ENUM, expr)               \
  do {                                          \
    Status s = (expr);                          \
    ASSERT_TRUE(s.Is##ENUM());                  \
  } while (0)


#define ASSERT_OK(expr)                         \
  do {                                          \
    Status s = (expr);                          \
    ASSERT_TRUE(s.ok());                        \
  } while (0)


#define EXPECT_OK(expr)                         \
  do {                                          \
    Status s = (expr);                          \
    EXPECT_TRUE(s.ok());                        \
  } while (0)

#endif // PANDAS_TEST_UTIL_H_
