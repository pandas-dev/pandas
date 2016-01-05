// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#include <Python.h>

#include <gtest/gtest.h>

#include "pandas/do_import_numpy.h"

#include "pandas/numpy_interop.h"

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  Py_Initialize();
  pandas::import_numpy();

  int ret = RUN_ALL_TESTS();

  Py_Finalize();

  return ret;
}
