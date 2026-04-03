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
// prevent name mangling
extern "C" {
#endif

// NOTE: It may be better to use the Moments struct
//       instead of using multiple double references
//       and a signed value to nobs.
typedef struct {
  int64_t n;
  double mean;
  double m2;
  double m3;
  double m4;
} Moments;

static inline void moments_add_value(double val, int64_t *nobs, double *mean,
                                     double *m2, double *m3, double *m4,
                                     int max_moment) {
  double delta = val - *mean;
  (*nobs)++;
  double n = (double)(*nobs);
  double delta_n = delta / n;
  double term1 = delta * delta_n * (n - 1.0);

  if (max_moment >= 4) {
    *m4 += delta_n * (-4.0 * *m3 +
                      delta_n * (6.0 * *m2 + term1 * (n * n - 3.0 * n + 3.0)));
  }
  if (max_moment >= 3) {
    *m3 += delta_n * (term1 * (n - 2.0) - 3.0 * *m2);
  }
  *m2 += term1;
  *mean += delta_n;
}

Moments moments_merge(Moments a, Moments b, int max_moment);

/// Compute central moments until `max_moment` using `n` elements from `values`.
void accumulate_moments_scalar(size_t n, const double *values, bool skipna,
                               const uint8_t *mask, int64_t *nobs, double *mean,
                               double *m2, double *m3, double *m4,
                               int max_moment);
#ifdef __cplusplus
}
#endif
