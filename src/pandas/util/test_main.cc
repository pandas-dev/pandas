// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gtest/gtest.h>

int main(int argc, char **argv) {
  google::InstallFailureSignalHandler();
  // InitGoogleTest() must precede ParseCommandLineFlags(), as the former
  // removes gtest-related flags from argv that would trip up the latter.
  ::testing::InitGoogleTest(&argc, argv);
  google::ParseCommandLineFlags(&argc, &argv, true);

  int ret = RUN_ALL_TESTS();

  return ret;
}
