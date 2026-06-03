/*
Copyright (c) 2026, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#pragma once

#include "pandas/moments.h"
#include "xsimd/xsimd.hpp"
#include <array>
#include <cassert>
#include <cmath>
#include <optional>
#include <span>

namespace pandas::moments {

namespace detail {

template <typename T, typename U> struct FirstPassAcc {
  T sum{};
  U count{};
};

template <typename T> struct SecondPassAcc {
  T m1{}, m2{}, m3{}, m4{};
};

/// Compute Moments struct using two-pass
/// with correction presented in Pébay (2016).
///
/// References:
/// Pébay, P., Terriberry, T. B., Kolla, H., & Bennett, J. (2016).
/// Numerically stable, scalable formulas for parallel
/// and online computation of higher-order multivariate central moments with
/// arbitrary weights. Comput. Stat., 31(4), 1305–1325.
/// https://doi.org/10.1007/s00180-015-0637-z
inline Moments
compute_moments_with_correction(FirstPassAcc<double, std::size_t> total_acc,
                                SecondPassAcc<double> central_diffs) {
  const auto count_double = static_cast<double>(total_acc.count);
  const double trial_mean = total_acc.sum / count_double;
  const double mean = trial_mean + (central_diffs.m1 / count_double);

  const double correction_term = central_diffs.m1 / count_double;
  const double correction_term2 = correction_term * correction_term;
  const double term1 = central_diffs.m1 * correction_term;
  const double m2 = central_diffs.m2 - term1;
  const double m3 = central_diffs.m3 -
                    (3.0 * central_diffs.m2 * correction_term) +
                    (2.0 * correction_term * term1);
  const double m4 = central_diffs.m4 -
                    (4.0 * central_diffs.m3 * correction_term) +
                    (6.0 * central_diffs.m2 * correction_term2) -
                    (3.0 * term1 * correction_term2);
  return {.mean = mean, .m2 = m2, .m3 = m3, .m4 = m4, .n = total_acc.count};
}

template <class Arch>
std::optional<
    FirstPassAcc<xsimd::batch<double, Arch>, xsimd::batch<std::uint64_t, Arch>>>
accumulate_first_pass_internal(std::span<const double> values, bool skipna) {
  using batch_type = xsimd::batch<double, Arch>;
  constexpr std::size_t step = batch_type::size;
  constexpr std::size_t leaf_size = 128 * step;

  assert(values.size() % step == 0);

  if (values.size() <= leaf_size) {
    xsimd::batch<std::uint64_t, Arch> vec_count{};
    batch_type vec_sum{};
    for (std::size_t i = 0; i < values.size(); i += step) {
      const auto val = xsimd::load_unaligned<Arch>(&values[i]);
      const xsimd::batch_bool<double, Arch> isna = xsimd::isnan(val);

      if (!skipna && xsimd::any(isna)) {
        return std::nullopt;
      }

      vec_count = xsimd::incr_if(vec_count,
                                 xsimd::batch_bool_cast<uint64_t, double, Arch>(
                                     xsimd::bitwise_not(isna)));
      vec_sum += xsimd::select(isna, batch_type(0), val);
    }
    return FirstPassAcc{.sum = vec_sum, .count = vec_count};
  }

  const std::size_t mid = (values.size() / 2 / step) * step;

  const auto left =
      accumulate_first_pass_internal<Arch>(values.first(mid), skipna);
  if (!left) {
    return std::nullopt;
  }

  const auto right = accumulate_first_pass_internal<Arch>(
      values.last(values.size() - mid), skipna);
  if (!right) {
    return std::nullopt;
  }

  return FirstPassAcc{.sum = left->sum + right->sum,
                      .count = left->count + right->count};
}

template <class Arch>
std::optional<
    FirstPassAcc<xsimd::batch<double, Arch>, xsimd::batch<std::uint64_t, Arch>>>
accumulate_first_pass_internal(std::span<const double> values,
                               std::span<const uint8_t> mask, bool skipna) {
  using value_batch_type = xsimd::batch<double, Arch>;
  using mask_batch_type = xsimd::batch<uint8_t, Arch>;
  constexpr std::size_t mask_step = mask_batch_type::size;
  constexpr std::size_t value_step = value_batch_type::size;

  // Use multiple accumulators to increase throughput.
  constexpr std::size_t naccumulators = 4;
  constexpr std::size_t leaf_size = 128 * naccumulators * value_step;
  static_assert(leaf_size % mask_step == 0);

  assert(values.size() == mask.size());
  assert(values.size() % mask_step == 0);

  if (values.size() <= leaf_size) {
    std::array<value_batch_type, naccumulators> sum{};
    xsimd::batch<std::uint64_t, Arch> count{};

    for (std::size_t i = 0; i < values.size(); i += mask_step) {
      const auto vector_mask = xsimd::load_unaligned<Arch>(&mask[i]);
      const typename mask_batch_type::batch_bool_type isna =
          vector_mask != mask_batch_type(0);

      if (!skipna && xsimd::any(isna)) {
        return std::nullopt;
      }

      const auto isna_bitmask = isna.mask();
      if (isna_bitmask == 0) {
        for (std::size_t j = 0; j < mask_step;
             j += naccumulators * value_step) {
          for (std::size_t k = 0; k < naccumulators; k++) {
            sum[k] +=
                xsimd::load_unaligned<Arch>(&values[i + j + (k * value_step)]);
          }
        }

        count += xsimd::batch<std::uint64_t, Arch>(8);
        continue;
      }

      constexpr uint64_t batch_mask = (1ULL << value_step) - 1;
      for (std::size_t j = 0; j < mask_step; j += naccumulators * value_step) {
        for (std::size_t k = 0; k < naccumulators; k++) {
          std::size_t idx = j + (k * value_step);
          const std::uint64_t batch_isna_bitmask =
              (isna_bitmask >> idx) & batch_mask;

          const auto isna_batch =
              xsimd::batch_bool<uint64_t, Arch>::from_mask(batch_isna_bitmask);
          const auto val = xsimd::load_unaligned<Arch>(&values[i + idx]);

          count = xsimd::incr_if(count, !isna_batch);
          sum[k] += xsimd::select(
              xsimd::batch_bool_cast<double, uint64_t, Arch>(isna_batch),
              value_batch_type(0), val);
        }
      }
    }

    for (std::size_t k = naccumulators / 2; k > 0; k /= 2) {
      for (std::size_t i = 0; i < k; i++) {
        sum[i] += sum[i + k];
      }
    }

    return FirstPassAcc{.sum = sum[0], .count = count};
  }

  const std::size_t mid = (values.size() / 2 / mask_step) * mask_step;
  const auto left = accumulate_first_pass_internal<Arch>(
      values.first(mid), mask.first(mid), skipna);
  if (!left)
    return std::nullopt;
  const auto right = accumulate_first_pass_internal<Arch>(
      values.last(values.size() - mid), mask.last(mask.size() - mid), skipna);
  if (!right)
    return std::nullopt;

  return FirstPassAcc{left->sum + right->sum, left->count + right->count};
}

// TODO: use std::expected when we start to use c++23
template <class Arch>
std::optional<FirstPassAcc<double, std::size_t>>
accumulate_first_pass(std::span<const double> values, bool skipna) {
  using batch_type = xsimd::batch<double, Arch>;
  constexpr std::size_t step = batch_type::size;
  const std::size_t vec_size = (values.size() / step) * step;

  const auto vec_acc =
      accumulate_first_pass_internal<Arch>(values.first(vec_size), skipna);
  if (!vec_acc) {
    return std::nullopt;
  }

  double tail_sum = 0.0;
  std::size_t tail_count = 0;
  for (std::size_t i = vec_size; i < values.size(); i++) {
    const auto val = values[i];
    const auto isna = std::isnan(val);

    if (isna) {
      if (!skipna) {
        return std::nullopt;
      }
      continue;
    }

    tail_count++;
    tail_sum += val;
  }

  return FirstPassAcc<double, std::size_t>{
      .sum = xsimd::reduce_add(vec_acc->sum) + tail_sum,
      .count = xsimd::reduce_add(vec_acc->count) + tail_count};
}

template <class Arch>
std::optional<FirstPassAcc<double, std::size_t>>
accumulate_first_pass(std::span<const double> values,
                      std::span<const uint8_t> mask, bool skipna) {
  assert(values.size() == mask.size());
  using mask_batch_type = xsimd::batch<uint8_t, Arch>;

  constexpr std::size_t mask_step = mask_batch_type::size;

  const std::size_t vec_size = (values.size() / mask_step) * mask_step;

  const auto vec_acc = accumulate_first_pass_internal<Arch>(
      values.first(vec_size), mask.first(vec_size), skipna);
  if (!vec_acc) {
    return std::nullopt;
  }

  double tail_sum = 0.0;
  std::size_t tail_count = 0;
  for (std::size_t i = vec_size; i < values.size(); i++) {
    const auto val = values[i];
    const auto isna = mask[i] != 0;

    if (isna) {
      if (!skipna) {
        return std::nullopt;
      }
      continue;
    }

    tail_count++;
    tail_sum += val;
  }

  return FirstPassAcc<double, std::size_t>{
      .sum = xsimd::reduce_add(vec_acc->sum) + tail_sum,
      .count = xsimd::reduce_add(vec_acc->count) + tail_count};
}

template <class Arch>
SecondPassAcc<xsimd::batch<double, Arch>>
accumulate_second_pass_internal(std::span<const double> values,
                                const xsimd::batch<double, Arch> &mean_vector) {
  using batch_type = xsimd::batch<double, Arch>;
  constexpr std::size_t step = batch_type::size;
  constexpr std::size_t leaf_size = 128 * step;

  if (values.size() <= leaf_size) {
    SecondPassAcc<batch_type> acc{};
    for (std::size_t i = 0; i < values.size(); i += step) {
      const auto val = xsimd::load_unaligned<Arch>(&values[i]);
      const auto isna = xsimd::isnan(val);

      const auto diff = xsimd::select(isna, batch_type(0), val - mean_vector);
      const auto diff2 = diff * diff;

      acc.m1 += diff;
      acc.m2 += diff2;
      acc.m3 += diff2 * diff;
      acc.m4 += diff2 * diff2;
    }
    return acc;
  }

  const std::size_t mid = (values.size() / 2 / step) * step;
  const auto left =
      accumulate_second_pass_internal<Arch>(values.first(mid), mean_vector);
  const auto right = accumulate_second_pass_internal<Arch>(
      values.last(values.size() - mid), mean_vector);

  return SecondPassAcc<batch_type>{
      .m1 = left.m1 + right.m1,
      .m2 = left.m2 + right.m2,
      .m3 = left.m3 + right.m3,
      .m4 = left.m4 + right.m4,
  };
}

template <class Arch>
SecondPassAcc<xsimd::batch<double, Arch>>
accumulate_second_pass_masked_internal(
    std::span<const double> values, std::span<const uint8_t> mask,
    const xsimd::batch<double, Arch> &mean_vector) {
  using value_batch_type = xsimd::batch<double, Arch>;
  using mask_batch_type = xsimd::batch<uint8_t, Arch>;
  constexpr std::size_t mask_step = mask_batch_type::size;
  constexpr std::size_t value_step = value_batch_type::size;
  constexpr std::size_t leaf_size = 128 * mask_step;

  if (values.size() <= leaf_size) {
    SecondPassAcc<value_batch_type> acc{};
    const uint64_t batch_mask = (1ULL << value_step) - 1;

    for (std::size_t i = 0; i < values.size(); i += mask_step) {
      const auto vector_mask = xsimd::load_unaligned<Arch>(&mask[i]);
      const typename mask_batch_type::batch_bool_type isna =
          vector_mask != mask_batch_type(0);
      const auto isna_bitmask = isna.mask();

      for (std::size_t k = 0; k < mask_step; k += value_step) {
        const std::uint64_t batch_isna_bitmask =
            (isna_bitmask >> k) & batch_mask;
        const auto isna_batch =
            xsimd::batch_bool<double, Arch>::from_mask(batch_isna_bitmask);
        const auto val = xsimd::load_unaligned<Arch>(&values[i + k]);

        const auto diff =
            xsimd::select(isna_batch, value_batch_type(0), val - mean_vector);
        const auto diff2 = diff * diff;

        acc.m1 += diff;
        acc.m2 += diff2;
        acc.m3 += diff2 * diff;
        acc.m4 += diff2 * diff2;
      }
    }
    return acc;
  }

  const std::size_t mid = (values.size() / 2 / mask_step) * mask_step;
  const auto left = accumulate_second_pass_masked_internal<Arch>(
      values.first(mid), mask.first(mid), mean_vector);
  const auto right = accumulate_second_pass_masked_internal<Arch>(
      values.last(values.size() - mid), mask.last(values.size() - mid),
      mean_vector);

  return {
      left.m1 + right.m1,
      left.m2 + right.m2,
      left.m3 + right.m3,
      left.m4 + right.m4,
  };
}

template <class Arch>
SecondPassAcc<double> accumulate_second_pass(std::span<const double> values,
                                             double mean) {
  using batch_type = xsimd::batch<double, Arch>;
  constexpr std::size_t step = batch_type::size;
  const std::size_t vec_size = (values.size() / step) * step;

  batch_type mean_vector = xsimd::broadcast<double, Arch>(mean);

  const auto vec_acc = accumulate_second_pass_internal<Arch>(
      values.first(vec_size), mean_vector);

  double tail_m1{};
  double tail_m2{};
  double tail_m3{};
  double tail_m4{};
  for (std::size_t i = vec_size; i < values.size(); i++) {
    const auto val = values[i];
    const auto isna = std::isnan(val);

    if (isna) {
      // skipna is guaranteed on the second pass
      continue;
    }

    const auto diff = val - mean;
    const auto diff2 = diff * diff;

    tail_m1 += diff;
    tail_m2 += diff2;
    tail_m3 += diff2 * diff;
    tail_m4 += diff2 * diff2;
  }

  return {.m1 = xsimd::reduce_add(vec_acc.m1) + tail_m1,
          .m2 = xsimd::reduce_add(vec_acc.m2) + tail_m2,
          .m3 = xsimd::reduce_add(vec_acc.m3) + tail_m3,
          .m4 = xsimd::reduce_add(vec_acc.m4) + tail_m4};
}

template <class Arch>
SecondPassAcc<double> accumulate_second_pass(std::span<const double> values,
                                             std::span<const uint8_t> mask,
                                             double mean) {
  assert(values.size() == mask.size());
  using value_batch_type = xsimd::batch<double, Arch>;
  using mask_batch_type = xsimd::batch<uint8_t, Arch>;
  constexpr std::size_t mask_step = mask_batch_type::size;

  const std::size_t vec_size = (values.size() / mask_step) * mask_step;

  value_batch_type mean_vector = xsimd::broadcast<double, Arch>(mean);

  const auto vec_acc = accumulate_second_pass_masked_internal<Arch>(
      values.first(vec_size), mask.first(vec_size), mean_vector);

  double tail_m1{};
  double tail_m2{};
  double tail_m3{};
  double tail_m4{};
  for (std::size_t i = vec_size; i < values.size(); i++) {
    const auto val = values[i];
    const auto isna = mask[i] != 0;

    if (isna) {
      continue;
    }

    const auto diff = val - mean;
    const auto diff2 = diff * diff;

    tail_m1 += diff;
    tail_m2 += diff2;
    tail_m3 += diff2 * diff;
    tail_m4 += diff2 * diff2;
  }

  return {.m1 = xsimd::reduce_add(vec_acc.m1) + tail_m1,
          .m2 = xsimd::reduce_add(vec_acc.m2) + tail_m2,
          .m3 = xsimd::reduce_add(vec_acc.m3) + tail_m3,
          .m4 = xsimd::reduce_add(vec_acc.m4) + tail_m4};
}

template <class Arch>
Moments accumulate_moments_simd_impl(std::span<const double> values,
                                     bool skipna) {
  const auto total_acc_opt = accumulate_first_pass<Arch>(values, skipna);
  if (!total_acc_opt.has_value()) {
    return {.mean = NAN, .m2 = NAN, .m3 = NAN, .m4 = NAN, .n = values.size()};
  }

  const auto total_acc = *total_acc_opt;
  const auto count_double = static_cast<double>(total_acc.count);
  double trial_mean = total_acc.sum / count_double;

  const SecondPassAcc<double> central_diffs =
      accumulate_second_pass<Arch>(values, trial_mean);

  return compute_moments_with_correction(total_acc, central_diffs);
}

template <class Arch>
Moments accumulate_moments_simd_impl(std::span<const double> values,
                                     std::span<const uint8_t> mask,
                                     bool skipna) {
  const auto total_acc = accumulate_first_pass<Arch>(values, mask, skipna);
  if (!total_acc.has_value()) {
    return {.mean = NAN, .m2 = NAN, .m3 = NAN, .m4 = NAN, .n = values.size()};
  }

  const auto count_double = static_cast<double>(total_acc->count);
  double trial_mean = total_acc->sum / count_double;

  const SecondPassAcc<double> central_diffs =
      accumulate_second_pass<Arch>(values, mask, trial_mean);

  return compute_moments_with_correction(*total_acc, central_diffs);
}

} // namespace detail

struct accumulate_moments_simd {
  template <class Arch>
  Moments operator()(Arch /*arch*/, std::span<const double> values, bool skipna,
                     std::optional<std::span<const uint8_t>> mask,
                     int max_moment) noexcept;
};

template <class Arch>
Moments accumulate_moments_simd::operator()(
    Arch /*arch*/, std::span<const double> values, bool skipna,
    std::optional<std::span<const uint8_t>> mask, int /*max_moment*/) noexcept {
  if (mask.has_value()) {
    return detail::accumulate_moments_simd_impl<Arch>(values, *mask, skipna);
  }
  return detail::accumulate_moments_simd_impl<Arch>(values, skipna);
}

} // namespace pandas::moments
