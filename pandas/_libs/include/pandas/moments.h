/*
Copyright (c) 2026, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  uint64_t n;
  double mean;
  double m2;
  double m3;
  double m4;
} Moments;

static inline void moments_add_valuem(Moments *moments, double val,
                                      int max_moment) {
  double delta = val - moments->mean;
  moments->n++;
  double n = (double)moments->n;
  double delta_n = delta / n;
  double term1 = delta * delta_n * (n - 1.0);

  if (max_moment >= 4) {
    moments->m4 +=
        delta_n *
        (-4.0 * moments->m3 +
         delta_n * (6.0 * moments->m2 + term1 * (n * n - 3.0 * n + 3.0)));
  }
  if (max_moment >= 3) {
    moments->m3 += delta_n * (term1 * (n - 2.0) - 3.0 * moments->m2);
  }
  moments->m2 += term1;
  moments->mean += delta_n;
}

void moments_merge(Moments *acc, const Moments *src, int max_moment);

/// Compute central moments until `max_moment` using `n` elements from `values`.
Moments moments_reduce(size_t n, const double *values, bool skipna,
                       const uint8_t *mask, int max_moment);
#ifdef __cplusplus
}
#endif
