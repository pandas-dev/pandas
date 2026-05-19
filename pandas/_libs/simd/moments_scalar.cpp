/*
Copyright (c) 2026, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#include "pandas/moments.h"
#include "pandas/simd/moments_simd.hpp"
#include <optional>
#include <span>

Moments moments_reduce(const double *values, size_t n, bool skipna,
                       const uint8_t *mask, int max_moment) {
  std::span<const double> values_span(values, n);
  std::optional<std::span<const uint8_t>> mask_span{};
  if (mask != nullptr)
    mask_span = std::span(mask, n);
  return pandas::moments::accumulate_moments_simd{}(
      xsimd::common{}, values_span, skipna, mask_span, max_moment);
}
