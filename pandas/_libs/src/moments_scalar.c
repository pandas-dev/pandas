/*
Copyright (c) 2026, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#include "pandas/moments.h"
#include <math.h>
#include <stddef.h>
#include <stdint.h>

Moments moments_reduce(size_t n, const double *values, bool skipna,
                       const uint8_t *mask, int max_moment) {
  Moments acc = {0};
  for (size_t i = 0; i < n; ++i) {
    double val = values[i];
    if (mask && mask[i])
      val = NAN;

    if (!isnan(val)) {
      moments_add_valuem(&acc, val, max_moment);
    } else if (!skipna) {
      return (Moments){(uint64_t)n, NAN, NAN, NAN, NAN};
    }
  }
  return acc;
}
