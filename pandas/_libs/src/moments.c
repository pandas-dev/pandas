/*
Copyright (c) 2026, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#include "pandas/moments.h"
#include <math.h>

Moments accumulate_moments_scalar(const double *values, int64_t n, int skipna,
                                  const uint8_t *mask, int max_moment) {
  Moments result = {0};

#ifdef _OPENMP
#  pragma omp parallel
#endif
  {
    Moments moments_local = {0};

#ifdef _OPENMP
#  pragma omp for nowait
#endif
    for (int64_t i = 0; i < n; i++) {
      double val = values[i];
      if (mask && mask[i]) {
        val = NAN;
      }
      if (skipna && isnan(val)) {
        continue;
      }
      moments_add_value(&moments_local, val, max_moment);
    }

#ifdef _OPENMP
#  pragma omp critical
#endif
    {
      result = moments_merge(result, moments_local, max_moment);
    }
  }

  return result;
}
