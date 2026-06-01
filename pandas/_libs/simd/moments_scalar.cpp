/*
Copyright (c) 2026, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#include "pandas/moments.h"
#include "pandas/simd/moments_simd.hpp"
#include <cmath>
#include <optional>
#include <span>

namespace pandas::moments {

namespace {

std::optional<detail::FirstPassAcc<double, std::size_t>>
accumulate_first_pass_scalar(const std::span<const double> values,
                             const std::optional<std::span<const uint8_t>> mask,
                             const bool skipna) {
  constexpr std::size_t leaf_size = 128;
  if (values.size() <= leaf_size) {
    double sum = 0.0;
    std::size_t count = 0;
    for (std::size_t i = 0; i < values.size(); ++i) {
      const double val = values[i];
      const bool isna = mask.has_value() ? (*mask)[i] != 0 : std::isnan(val);
      if (isna) {
        if (!skipna) {
          return std::nullopt;
        }
        continue;
      }
      sum += val;
      count++;
    }
    return detail::FirstPassAcc{.sum = sum, .count = count};
  }

  const std::size_t mid = values.size() / 2;
  const std::optional<std::span<const uint8_t>> left_mask =
      mask.has_value() ? std::make_optional((*mask).first(mid)) : std::nullopt;
  const std::optional<std::span<const uint8_t>> right_mask =
      mask.has_value() ? std::make_optional((*mask).last(values.size() - mid))
                       : std::nullopt;

  const auto left =
      accumulate_first_pass_scalar(values.first(mid), left_mask, skipna);
  if (!left) {
    return std::nullopt;
  }
  const auto right = accumulate_first_pass_scalar(
      values.last(values.size() - mid), right_mask, skipna);
  if (!right) {
    return std::nullopt;
  }

  return detail::FirstPassAcc{.sum = left->sum + right->sum,
                              .count = left->count + right->count};
}

detail::SecondPassAcc<double> accumulate_second_pass_scalar(
    const std::span<const double> values,
    const std::optional<std::span<const uint8_t>> mask, const double mean) {
  constexpr std::size_t leaf_size = 128;
  if (values.size() <= leaf_size) {
    detail::SecondPassAcc<double> acc{};
    for (std::size_t i = 0; i < values.size(); ++i) {
      const double val = values[i];
      const bool isna = mask.has_value() ? (*mask)[i] != 0 : std::isnan(val);
      if (isna) {
        continue;
      }
      const double diff = val - mean;
      const double diff2 = diff * diff;
      acc.m1 += diff;
      acc.m2 += diff2;
      acc.m3 += diff2 * diff;
      acc.m4 += diff2 * diff2;
    }
    return acc;
  }

  const std::size_t mid = values.size() / 2;
  const std::optional<std::span<const uint8_t>> left_mask =
      mask.has_value() ? std::make_optional((*mask).first(mid)) : std::nullopt;
  const std::optional<std::span<const uint8_t>> right_mask =
      mask.has_value() ? std::make_optional((*mask).last(values.size() - mid))
                       : std::nullopt;

  const auto left =
      accumulate_second_pass_scalar(values.first(mid), left_mask, mean);
  const auto right = accumulate_second_pass_scalar(
      values.last(values.size() - mid), right_mask, mean);

  return {.m1 = left.m1 + right.m1,
          .m2 = left.m2 + right.m2,
          .m3 = left.m3 + right.m3,
          .m4 = left.m4 + right.m4};
}

} // namespace

template <>
Moments accumulate_moments_simd::operator()<xsimd::common>(
    xsimd::common /*arch*/, const std::span<const double> values,
    const bool skipna, const std::optional<std::span<const uint8_t>> mask,
    int /*max_moment*/) noexcept {

  const auto total_acc_opt = accumulate_first_pass_scalar(values, mask, skipna);

  if (!total_acc_opt.has_value()) {
    return {.mean = NAN, .m2 = NAN, .m3 = NAN, .m4 = NAN, .n = values.size()};
  }

  const auto total_acc = *total_acc_opt;
  if (total_acc.count == 0) {
    return {};
  }

  const double trial_mean =
      total_acc.sum / static_cast<double>(total_acc.count);
  const detail::SecondPassAcc<double> central_diffs =
      accumulate_second_pass_scalar(values, mask, trial_mean);

  return detail::compute_moments_with_correction(total_acc, central_diffs);
}

} // namespace pandas::moments

Moments moments_reduce(const double *values, size_t n, bool skipna,
                       const uint8_t *mask, int max_moment) {
  const std::span<const double> values_span(values, n);
  std::optional<std::span<const uint8_t>> mask_span{};
  if (mask != nullptr) {
    mask_span = std::span(mask, n);
  }
  return pandas::moments::accumulate_moments_simd{}(
      xsimd::common{}, values_span, skipna, mask_span, max_moment);
}
