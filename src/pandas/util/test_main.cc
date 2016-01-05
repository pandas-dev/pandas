// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#include <Python.h>

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gtest/gtest.h>

#include "pandas/do_import_numpy.h"

#include "pandas/numpy_interop.h"

int main(int argc, char **argv) {
  google::InstallFailureSignalHandler();
  // InitGoogleTest() must precede ParseCommandLineFlags(), as the former
  // removes gtest-related flags from argv that would trip up the latter.
  ::testing::InitGoogleTest(&argc, argv);
  google::ParseCommandLineFlags(&argc, &argv, true);

  Py_Initialize();
  pandas::import_numpy();

  int ret = RUN_ALL_TESTS();

  Py_Finalize();

  return ret;
}
