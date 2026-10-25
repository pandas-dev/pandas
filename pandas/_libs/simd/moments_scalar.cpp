/*
Copyright (c) 2026, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#include "pandas/moments.h"
#include "pandas/moments_common.hpp"
#include <cmath>
#include <optional>
#include <span>

using pandas::moments::detail::accumulate_central_diffs_scalar_direct;
using pandas::moments::detail::accumulate_mean_scalar_direct;
using pandas::moments::detail::CentralDiffs;
using pandas::moments::detail::compute_moments_with_correction;
using pandas::moments::detail::MeanAcc;

namespace {

std::optional<MeanAcc<double, std::size_t>>
accumulate_mean(const std::span<const double> values,
                const std::optional<std::span<const uint8_t>> mask,
                const bool skipna) {
  constexpr std::size_t leaf_size = 128;

  if (values.size() <= leaf_size) {
    return accumulate_mean_scalar_direct(values, mask, skipna);
  }

  const std::size_t mid = values.size() / 2;
  const std::optional<std::span<const uint8_t>> left_mask =
      mask.has_value() ? std::make_optional(mask->first(mid)) : std::nullopt;
  const std::optional<std::span<const uint8_t>> right_mask =
      mask.has_value() ? std::make_optional(mask->last(values.size() - mid))
                       : std::nullopt;

  const auto left = accumulate_mean(values.first(mid), left_mask, skipna);
  if (!left) {
    return std::nullopt;
  }

  const auto right =
      accumulate_mean(values.last(values.size() - mid), right_mask, skipna);
  if (!right) {
    return std::nullopt;
  }

  return MeanAcc<double, std::size_t>{.sum = left->sum + right->sum,
                                      .count = left->count + right->count};
}

CentralDiffs<double>
accumulate_central_diffs(const std::span<const double> values,
                         const std::optional<std::span<const uint8_t>> mask,
                         double mean, int max_moment) {
  constexpr std::size_t leaf_size = 128;
  if (values.size() <= leaf_size) {
    return accumulate_central_diffs_scalar_direct(values, mask, mean,
                                                  max_moment);
  }

  const std::size_t mid = values.size() / 2;

  const std::optional<std::span<const uint8_t>> left_mask =
      mask.has_value() ? std::make_optional(mask->first(mid)) : std::nullopt;
  const std::optional<std::span<const uint8_t>> right_mask =
      mask.has_value() ? std::make_optional(mask->last(values.size() - mid))
                       : std::nullopt;

  const auto left =
      accumulate_central_diffs(values.first(mid), left_mask, mean, max_moment);
  const auto right = accumulate_central_diffs(values.last(values.size() - mid),
                                              right_mask, mean, max_moment);

  return {.m1 = left.m1 + right.m1,
          .m2 = left.m2 + right.m2,
          .m3 = left.m3 + right.m3,
          .m4 = left.m4 + right.m4};
}

} // namespace

Moments moments_reduce(const double *values, size_t n, bool skipna,
                       const uint8_t *mask, int max_moment) {
  const std::span<const double> values_span(values, n);
  const auto mask_span =
      mask == nullptr
          ? std::nullopt
          : std::optional<std::span<const uint8_t>>(std::span(mask, n));

  const auto total_acc = accumulate_mean(values_span, mask_span, skipna);

  if (!total_acc.has_value()) {
    return {.mean = NAN, .m2 = NAN, .m3 = NAN, .m4 = NAN, .n = n};
  }

  if (total_acc->count == 0) {
    return {};
  }

  const auto count_double = static_cast<double>(total_acc->count);
  double trial_mean = total_acc->sum / count_double;

  const auto central_diffs =
      accumulate_central_diffs(values_span, mask_span, trial_mean, max_moment);

  return compute_moments_with_correction(*total_acc, central_diffs, max_moment);
}
