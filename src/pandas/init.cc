// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

// Strategy adopted from dynd-python
#define NUMPY_IMPORT_ARRAY
#include "pandas/numpy_interop.h"

#include "pandas/init.h"

void pandas::libpandas_init() {
  import_numpy();
}
