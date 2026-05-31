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

namespace pandas::moments {

template <>
Moments accumulate_moments_simd::operator()<xsimd::common>(
    xsimd::common /*arch*/, std::span<const double> values, bool skipna,
    std::optional<std::span<const uint8_t>> mask, int /*max_moment*/) noexcept {

  detail::FirstPassAcc total_acc{};

  for (std::size_t i = 0; i < values.size(); ++i) {
    const auto val = values[i];
    const bool isna = mask.has_value() ? (*mask)[i] != 0 : std::isnan(val);
    if (isna) {
      if (!skipna) {
        return {
            .mean = NAN, .m2 = NAN, .m3 = NAN, .m4 = NAN, .n = values.size()};
      }
      continue;
    }
    total_acc.sum += val;
    total_acc.count++;
  }

  if (total_acc.count == 0) {
    return {};
  }

  const double trial_mean =
      total_acc.sum / static_cast<double>(total_acc.count);
  detail::SecondPassAcc central_diffs{};

  for (std::size_t i = 0; i < values.size(); ++i) {
    const auto val = values[i];
    const bool isna = mask.has_value() ? (*mask)[i] != 0 : std::isnan(val);
    if (isna)
      continue;

    const double diff = val - trial_mean;
    const double diff2 = diff * diff;

    central_diffs.m1 += diff;
    central_diffs.m2 += diff2;
    central_diffs.m3 += diff2 * diff;
    central_diffs.m4 += diff2 * diff2;
  }

  return detail::compute_moments_with_correction(total_acc, central_diffs);
}

} // namespace pandas::moments

Moments moments_reduce(const double *values, size_t n, bool skipna,
                       const uint8_t *mask, int max_moment) {
  std::span<const double> values_span(values, n);
  std::optional<std::span<const uint8_t>> mask_span{};
  if (mask != nullptr)
    mask_span = std::span(mask, n);
  return pandas::moments::accumulate_moments_simd{}(
      xsimd::common{}, values_span, skipna, mask_span, max_moment);
}
