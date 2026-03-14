/*
Copyright (c) 2026, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#pragma once

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
// prevent name mangling
extern "C" {
#endif

typedef struct {
  int64_t n;
  double mean;
  double m2;
  double m3;
  double m4;
} Moments;

static inline void moments_add_value(Moments *moments, double val,
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

static inline Moments moments_merge(Moments a, Moments b, int max_moment) {

  if (a.n == 0) {
    return b;
  }
  if (b.n == 0) {
    return a;
  }

  Moments result;

  result.n = a.n + b.n;
  double n_a = (double)a.n;
  double n_b = (double)b.n;
  double delta = b.mean - a.mean;
  double delta_n = delta / (double)result.n;
  double term1 = delta * delta_n * n_a * n_b;

  if (max_moment >= 4) {
    result.m4 =
        a.m4 + b.m4 +
        delta_n * (4.0 * (n_a * b.m3 - n_b * a.m3) +
                   delta_n * (6.0 * (n_a * n_a * b.m2 + n_b * n_b * a.m2) +
                              term1 * (n_a * n_a - n_a * n_b + n_b * n_b)));
  }
  if (max_moment >= 3) {
    result.m3 =
        a.m3 + b.m3 +
        delta_n * (3.0 * (n_a * b.m2 - n_b * a.m2) + term1 * (n_a - n_b));
  }
  result.m2 = a.m2 + b.m2 + term1;
  result.mean = a.mean + delta_n * n_b;

  return result;
}

/// Compute central moments until `max_moment` using `n` elements from `values`.
/// The size is represented as signed integer for MSVC compatibility with
/// OpenMP.
Moments accumulate_moments_scalar(const double *values, int64_t n, int skipna,
                                  const uint8_t *mask, int max_moment);

static inline double calc_skew(Moments moments) {
  if (moments.n < 3) {
    return NAN;
  }
  double dnobs = (double)moments.n;
  double moments_ratio = moments.m3 / (moments.m2 * sqrt(moments.m2));
  double correction = (dnobs * sqrt(dnobs - 1.0)) / (dnobs - 2.0);
  return moments_ratio * correction;
}

static inline double calc_kurt(Moments moments) {
  if (moments.n < 4) {
    return NAN;
  }
  double dnobs = (double)moments.n;
  double moments_ratio = moments.m4 / (moments.m2 * moments.m2);
  double term1 = dnobs * (dnobs + 1.0) * moments_ratio;
  double term2 = 3.0 * (dnobs - 1.0);
  double inner = term1 - term2;
  double correction = (dnobs - 1.0) / ((dnobs - 2.0) * (dnobs - 3.0));
  return correction * inner;
}

#ifdef __cplusplus
}
#endif
