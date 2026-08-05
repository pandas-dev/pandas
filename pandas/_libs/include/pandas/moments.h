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
  double mean;
  double m2;
  double m3;
  double m4;
  size_t n;
} Moments;

/// Compute central moments until `max_moment` using `n` elements from `values`.
Moments moments_reduce(const double *values, size_t n, bool skipna,
                       const uint8_t *mask, int max_moment);
#ifdef __cplusplus
}
#endif
