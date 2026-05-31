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

struct FirstPassAcc {
  double sum{};
  std::size_t count{};
};

struct SecondPassAcc {
  double m1, m2, m3, m4;
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
inline Moments compute_moments_with_correction(FirstPassAcc total_acc,
                                               SecondPassAcc central_diffs) {
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

// TODO: use std::expected when we start to use c++23
template <class Arch>
std::optional<FirstPassAcc>
accumulate_first_pass(std::span<const double> values, bool skipna) {
  using batch_type = xsimd::batch<double, Arch>;
  constexpr std::size_t step = batch_type::size;
  const std::size_t tail_size = values.size() % step;
  const std::size_t vec_size = values.size() - tail_size;

  xsimd::batch<std::uint64_t, Arch> vec_count{};
  batch_type vec_sum{};
  for (std::size_t i = 0; i < vec_size; i += step) {
    const auto val = xsimd::load_unaligned(&values[i]);
    const xsimd::batch_bool<double, Arch> isna = xsimd::isnan(val);

    if (!skipna && xsimd::any(isna)) {
      return std::nullopt;
    }

    vec_count = xsimd::incr_if(vec_count,
                               xsimd::batch_bool_cast<uint64_t, double, Arch>(
                                   xsimd::bitwise_not(isna)));
    vec_sum += xsimd::select(isna, batch_type(0), val);
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

  const std::size_t count = xsimd::reduce_add(vec_count) + tail_count;
  const double sum = xsimd::reduce_add(vec_sum) + tail_sum;
  return FirstPassAcc{.sum = sum, .count = count};
}

template <class Arch>
std::optional<FirstPassAcc>
accumulate_first_pass(std::span<const double> values,
                      std::span<const uint8_t> mask, bool skipna) {
  assert(values.size() == mask.size());
  using value_batch_type = xsimd::batch<double, Arch>;
  using mask_batch_type = xsimd::batch<uint8_t, Arch>;

  constexpr std::size_t mask_step = mask_batch_type::size;
  constexpr std::size_t value_step = value_batch_type::size;
  constexpr std::size_t num_accumulators = mask_step / value_step;

  const std::size_t tail_size = values.size() % mask_step;
  const std::size_t vec_size = values.size() - tail_size;

  std::array<xsimd::batch<std::uint64_t, Arch>, num_accumulators> vec_counts{};
  std::array<value_batch_type, num_accumulators> vec_sums{};

  for (std::size_t i = 0; i < vec_size; i += mask_step) {
    const auto vector_mask = xsimd::load_unaligned(&mask[i]);
    const typename mask_batch_type::batch_bool_type isna =
        vector_mask != mask_batch_type(0);

    if (!skipna && xsimd::any(isna)) {
      return std::nullopt;
    }

    const auto isna_bitmask = isna.mask();

    if (isna_bitmask == 0) {
      for (std::size_t k = 0; k < num_accumulators; ++k) {
        vec_sums[k] += xsimd::load_unaligned(&values[i + (k * value_step)]);
        vec_counts[k] = xsimd::incr(vec_counts[k]);
      }
      continue;
    }

    constexpr uint64_t all_na_mask = (1ULL << mask_step) - 1;
    if (isna_bitmask == all_na_mask) {
      continue;
    }

    constexpr uint64_t batch_mask = (1ULL << value_step) - 1;

    static_assert(mask_step % value_step == 0);

    for (std::size_t k = 0; k < num_accumulators; ++k) {
      const std::size_t value_offset = k * value_step;
      const std::uint64_t batch_isna_bitmask =
          (isna_bitmask >> value_offset) & batch_mask;

      const auto isna_batch =
          xsimd::batch_bool<uint64_t, Arch>::from_mask(batch_isna_bitmask);
      const auto val = xsimd::load_unaligned(&values[i + value_offset]);

      vec_counts[k] = xsimd::incr_if(vec_counts[k], !isna_batch);
      vec_sums[k] += xsimd::select(
          xsimd::batch_bool_cast<double, uint64_t, Arch>(isna_batch),
          value_batch_type(0), val);
    }
  }

  // Reduce accumulators to the first one
  for (std::size_t k = num_accumulators / 2; k > 0; k /= 2) {
    for (std::size_t i = 0; i < k; ++i) {
      vec_counts[i] += vec_counts[i + k];
      vec_sums[i] += vec_sums[i + k];
    }
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

  const std::size_t count = xsimd::reduce_add(vec_counts[0]) + tail_count;
  const double sum = xsimd::reduce_add(vec_sums[0]) + tail_sum;
  return FirstPassAcc{.sum = sum, .count = count};
}

template <class Arch>
SecondPassAcc accumulate_second_pass(std::span<const double> values,
                                     double mean) {
  using batch_type = xsimd::batch<double, Arch>;
  constexpr std::size_t step = batch_type::size;
  const std::size_t tail_size = values.size() % step;
  const std::size_t vec_size = values.size() - tail_size;

  batch_type mean_vector = xsimd::broadcast(mean);

  batch_type m1_vector{};
  batch_type m2_vector{};
  batch_type m3_vector{};
  batch_type m4_vector{};
  for (std::size_t i = 0; i < vec_size; i += step) {
    const auto val = xsimd::load_unaligned(&values[i]);
    const auto isna = xsimd::isnan(val);

    const auto diff = xsimd::select(isna, batch_type(0), val - mean_vector);
    const auto diff2 = diff * diff;

    m1_vector += diff;
    m2_vector += diff2;
    m3_vector += diff2 * diff;
    m4_vector += diff2 * diff2;
  }

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

  return {xsimd::reduce_add(m1_vector) + tail_m1,
          xsimd::reduce_add(m2_vector) + tail_m2,
          xsimd::reduce_add(m3_vector) + tail_m3,
          xsimd::reduce_add(m4_vector) + tail_m4};
}

template <class Arch>
SecondPassAcc accumulate_second_pass(std::span<const double> values,
                                     std::span<const uint8_t> mask,
                                     double mean) {
  assert(values.size() == mask.size());
  using value_batch_type = xsimd::batch<double, Arch>;
  using mask_batch_type = xsimd::batch<uint8_t, Arch>;
  constexpr std::size_t mask_step = mask_batch_type::size;
  constexpr std::size_t value_step = value_batch_type::size;

  const std::size_t tail_size = values.size() % mask_step;
  const std::size_t vec_size = values.size() - tail_size;

  value_batch_type mean_vector = xsimd::broadcast(mean);

  value_batch_type m1_vector{};
  value_batch_type m2_vector{};
  value_batch_type m3_vector{};
  value_batch_type m4_vector{};
  for (std::size_t i = 0; i < vec_size; i += mask_step) {
    const auto vector_mask = xsimd::load_unaligned(&mask[i]);
    const typename mask_batch_type::batch_bool_type isna =
        vector_mask != mask_batch_type(0);

    const auto isna_bitmask = isna.mask();
    const uint64_t batch_mask = (1ULL << value_step) - 1;

    static_assert(mask_step % value_step == 0);
    for (std::size_t j = 0; j < mask_step; j += value_step) {
      const std::uint64_t batch_isna_bitmask = (isna_bitmask >> j) & batch_mask;
      const auto isna_batch =
          xsimd::batch_bool<double, Arch>::from_mask(batch_isna_bitmask);
      const auto val = xsimd::load_unaligned(&values[i + j]);

      const auto diff =
          xsimd::select(isna_batch, value_batch_type(0), val - mean_vector);
      const auto diff2 = diff * diff;

      m1_vector += diff;
      m2_vector += diff2;
      m3_vector += diff2 * diff;
      m4_vector += diff2 * diff2;
    }
  }

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

  return {xsimd::reduce_add(m1_vector) + tail_m1,
          xsimd::reduce_add(m2_vector) + tail_m2,
          xsimd::reduce_add(m3_vector) + tail_m3,
          xsimd::reduce_add(m4_vector) + tail_m4};
}

template <class Arch>
Moments accumulate_moments_simd_impl(std::span<const double> values,
                                     bool skipna) {
  const std::optional<FirstPassAcc> total_acc_opt =
      accumulate_first_pass<Arch>(values, skipna);
  if (!total_acc_opt.has_value()) {
    return {.mean = NAN, .m2 = NAN, .m3 = NAN, .m4 = NAN, .n = values.size()};
  }

  const auto total_acc = *total_acc_opt;
  const auto count_double = static_cast<double>(total_acc.count);
  double trial_mean = total_acc.sum / count_double;

  const SecondPassAcc central_diffs =
      accumulate_second_pass<Arch>(values, trial_mean);

  return compute_moments_with_correction(total_acc, central_diffs);
}

template <class Arch>
Moments accumulate_moments_simd_impl(std::span<const double> values,
                                     std::span<const uint8_t> mask,
                                     bool skipna) {
  const std::optional<FirstPassAcc> total_acc_opt =
      accumulate_first_pass<Arch>(values, mask, skipna);
  if (!total_acc_opt.has_value()) {
    return {.mean = NAN, .m2 = NAN, .m3 = NAN, .m4 = NAN, .n = values.size()};
  }

  const auto total_acc = *total_acc_opt;
  const auto count_double = static_cast<double>(total_acc.count);
  double trial_mean = total_acc.sum / count_double;

  const SecondPassAcc central_diffs =
      accumulate_second_pass<Arch>(values, mask, trial_mean);

  return compute_moments_with_correction(total_acc, central_diffs);
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
