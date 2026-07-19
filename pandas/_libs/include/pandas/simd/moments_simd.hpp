/*
Copyright (c) 2026, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#pragma once

#include "pandas/moments.h"
#include "xsimd/xsimd.hpp" // IWYU pragma: keep
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <optional>
#include <span>

namespace pandas::moments {

namespace detail {

template <typename T, typename U> struct MeanAcc {
  T sum{};
  U count{};
};

template <typename T> struct CentralDiffs {
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
static inline Moments
compute_moments_with_correction(MeanAcc<double, std::size_t> total_acc,
                                CentralDiffs<double> central_diffs,
                                int max_moment) {
  const auto count_double = static_cast<double>(total_acc.count);
  const double trial_mean = total_acc.sum / count_double;
  const double correction_term = central_diffs.m1 / count_double;
  const double correction_term2 = correction_term * correction_term;
  const double term1 = central_diffs.m1 * correction_term;

  Moments result{};

  result.n = total_acc.count;
  result.mean = trial_mean + correction_term;
  result.m2 = central_diffs.m2 - term1;
  result.m2 = std::max(result.m2, 0.0);

  if (max_moment >= 3) {
    result.m3 = central_diffs.m3 - (3.0 * central_diffs.m2 * correction_term) +
                (2.0 * correction_term * term1);
  }
  if (max_moment >= 4) {
    result.m4 = central_diffs.m4 - (4.0 * central_diffs.m3 * correction_term) +
                (6.0 * central_diffs.m2 * correction_term2) -
                (3.0 * term1 * correction_term2);
    result.m4 = std::max(result.m4, 0.0);
  }

  return result;
}

static inline std::optional<detail::MeanAcc<double, std::size_t>>
accumulate_mean_scalar_direct(
    const std::span<const double> values,
    const std::optional<std::span<const uint8_t>> mask, const bool skipna) {
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
  // Explicit template arguments: CTAD for aggregates (P1816) is unavailable
  // in Apple clang < 17 (Xcode < 16).
  return MeanAcc<double, std::size_t>{.sum = sum, .count = count};
}

template <class Arch>
std::optional<
    MeanAcc<xsimd::batch<double, Arch>, xsimd::batch<std::uint64_t, Arch>>>
accumulate_mean_direct(std::span<const double> values, bool skipna) {
  using batch_type = xsimd::batch<double, Arch>;
  constexpr std::size_t step = batch_type::size;

  xsimd::batch<std::uint64_t, Arch> vec_count{};
  batch_type vec_sum{};

  for (std::size_t i = 0; i < values.size(); i += step) {
    const auto val = xsimd::load_unaligned<Arch>(&values[i]);
    const auto isna = xsimd::isnan(val);

    if (!skipna && xsimd::any(isna)) {
      return std::nullopt;
    }

    const auto is_not_na_si = xsimd::batch_bool_cast<uint64_t>(!isna);

    vec_count = xsimd::incr_if(vec_count, is_not_na_si);
    vec_sum += xsimd::select(isna, batch_type(0), val);
  }
  return MeanAcc<batch_type, xsimd::batch<std::uint64_t, Arch>>{vec_sum,
                                                                vec_count};
}

template <class Arch, std::size_t NAcc>
std::optional<
    MeanAcc<xsimd::batch<double, Arch>, xsimd::batch<std::uint64_t, Arch>>>
accumulate_mean_direct(std::span<const double> values,
                       std::span<const uint8_t> mask, bool skipna) {
  using value_batch_type = xsimd::batch<double, Arch>;
  using mask_batch_type = xsimd::batch<uint8_t, Arch>;
  constexpr std::size_t mask_step = mask_batch_type::size;
  constexpr std::size_t value_step = value_batch_type::size;

  static_assert(mask_step % (NAcc * value_step) == 0);

  std::array<value_batch_type, NAcc> sum{};
  xsimd::batch<std::uint64_t, Arch> count{};

  for (std::size_t i = 0; i < values.size(); i += mask_step) {
    const auto vector_mask = xsimd::load_unaligned<Arch>(&mask[i]);
    const auto isna = vector_mask != mask_batch_type(0);

    if (!skipna && xsimd::any(isna)) {
      return std::nullopt;
    }

    const std::uint64_t isna_bitmask = isna.mask();
    if (isna_bitmask == 0) {
      for (std::size_t j = 0; j < mask_step; j += NAcc * value_step) {
        for (std::size_t k = 0; k < NAcc; k++) {
          sum[k] +=
              xsimd::load_unaligned<Arch>(&values[i + j + (k * value_step)]);
        }
      }

      count += xsimd::batch<std::uint64_t, Arch>(mask_step / value_step);
      continue;
    }

    constexpr uint64_t batch_mask = (1ULL << value_step) - 1;
    for (std::size_t j = 0; j < mask_step; j += NAcc * value_step) {
      for (std::size_t k = 0; k < NAcc; k++) {
        std::size_t idx = j + (k * value_step);
        const std::uint64_t batch_isna_bitmask =
            (isna_bitmask >> idx) & batch_mask;

        const auto isna_batch =
            xsimd::batch_bool<uint64_t, Arch>::from_mask(batch_isna_bitmask);
        const auto isna_batch_pd = xsimd::batch_bool_cast<double>(isna_batch);
        const auto val = xsimd::load_unaligned<Arch>(&values[i + idx]);

        count = xsimd::incr_if(count, !isna_batch);
        sum[k] += xsimd::select(isna_batch_pd, value_batch_type(0), val);
      }
    }
  }

  for (std::size_t k = NAcc / 2; k > 0; k /= 2) {
    for (std::size_t i = 0; i < k; i++) {
      sum[i] += sum[i + k];
    }
  }

  return MeanAcc<value_batch_type, xsimd::batch<std::uint64_t, Arch>>{
      .sum = sum[0], .count = count};
}

static inline CentralDiffs<double> accumulate_central_diffs_scalar_direct(
    std::span<const double> values,
    std::optional<std::span<const uint8_t>> mask, double mean, int max_moment) {
  assert(!mask.has_value() || mask->size() == values.size());

  CentralDiffs<double> acc{};
  for (std::size_t i = 0; i < values.size(); i++) {
    const auto isna = mask ? ((*mask)[i] != 0) : std::isnan(values[i]);

    if (isna) {
      continue;
    }

    const auto diff = values[i] - mean;
    const auto diff2 = diff * diff;

    acc.m1 += diff;
    acc.m2 += diff2;
    if (max_moment >= 3) {
      acc.m3 += diff2 * diff;
    }
    if (max_moment >= 4) {
      acc.m4 += diff2 * diff2;
    }
  }

  return acc;
}

template <class Arch>
CentralDiffs<xsimd::batch<double, Arch>>
accumulate_central_diffs_direct(std::span<const double> values,
                                const xsimd::batch<double, Arch> &mean_vector,
                                int max_moment) {
  using batch_type = xsimd::batch<double, Arch>;
  constexpr std::size_t step = batch_type::size;

  CentralDiffs<batch_type> acc{};
  for (std::size_t i = 0; i < values.size(); i += step) {
    const auto val = xsimd::load_unaligned<Arch>(&values[i]);
    const auto isna = xsimd::isnan(val);

    const auto diff = xsimd::select(isna, batch_type(0), val - mean_vector);
    const auto diff2 = diff * diff;

    acc.m1 += diff;
    acc.m2 += diff2;
    if (max_moment >= 3) {
      acc.m3 += diff2 * diff;
    }
    if (max_moment >= 4) {
      acc.m4 += diff2 * diff2;
    }
  }
  return acc;
}

template <class Arch>
CentralDiffs<xsimd::batch<double, Arch>> accumulate_central_diffs_direct(
    std::span<const double> values, std::span<const uint8_t> mask,
    const xsimd::batch<double, Arch> &mean_vector, int max_moment) {
  using value_batch_type = xsimd::batch<double, Arch>;
  using mask_batch_type = xsimd::batch<uint8_t, Arch>;
  constexpr std::size_t mask_step = mask_batch_type::size;
  constexpr std::size_t value_step = value_batch_type::size;

  CentralDiffs<value_batch_type> acc{};
  constexpr uint64_t batch_mask = (1ULL << value_step) - 1;

  for (std::size_t i = 0; i < values.size(); i += mask_step) {
    const auto vector_mask = xsimd::load_unaligned<Arch>(&mask[i]);
    const auto isna = vector_mask != mask_batch_type(0);
    const std::uint64_t isna_bitmask = isna.mask();

    for (std::size_t k = 0; k < mask_step; k += value_step) {
      const std::uint64_t batch_isna_bitmask = (isna_bitmask >> k) & batch_mask;
      const auto isna_batch =
          xsimd::batch_bool<double, Arch>::from_mask(batch_isna_bitmask);
      const auto val = xsimd::load_unaligned<Arch>(&values[i + k]);

      const auto diff =
          xsimd::select(isna_batch, value_batch_type(0), val - mean_vector);
      const auto diff2 = diff * diff;

      acc.m1 += diff;
      acc.m2 += diff2;
      if (max_moment >= 3) {
        acc.m3 += diff2 * diff;
      }
      if (max_moment >= 4) {
        acc.m4 += diff2 * diff2;
      }
    }
  }
  return acc;
}

template <class Arch>
std::optional<
    MeanAcc<xsimd::batch<double, Arch>, xsimd::batch<std::uint64_t, Arch>>>
accumulate_mean_pairwise(std::span<const double> values,
                         std::optional<std::span<const uint8_t>> mask,
                         bool skipna) {
  assert(!mask || mask->size() == values.size());
  using value_batch_type = xsimd::batch<double, Arch>;
  using mask_batch_type = xsimd::batch<uint8_t, Arch>;

  constexpr std::size_t naccumulators = 4;
  const std::size_t step =
      mask ? mask_batch_type::size : value_batch_type::size;

  std::size_t leaf_size;
  if (mask) {
    leaf_size = 128 * naccumulators * value_batch_type::size;
  } else {
    leaf_size = 128 * value_batch_type::size;
  }

  if (values.size() <= leaf_size) {
    if (mask) {
      return accumulate_mean_direct<Arch, naccumulators>(values, *mask, skipna);
    }
    return accumulate_mean_direct<Arch>(values, skipna);
  }

  const std::size_t mid = (values.size() / 2 / step) * step;

  const auto left_mask =
      mask ? std::make_optional(mask->first(mid)) : std::nullopt;
  const auto left =
      accumulate_mean_pairwise<Arch>(values.first(mid), left_mask, skipna);
  if (!left) {
    return std::nullopt;
  }

  const auto right_mask =
      mask ? std::make_optional(mask->last(mask->size() - mid)) : std::nullopt;
  const auto right = accumulate_mean_pairwise<Arch>(
      values.last(values.size() - mid), right_mask, skipna);
  if (!right) {
    return std::nullopt;
  }

  return MeanAcc<value_batch_type, xsimd::batch<std::uint64_t, Arch>>{
      .sum = left->sum + right->sum, .count = left->count + right->count};
}

template <class Arch>
CentralDiffs<xsimd::batch<double, Arch>>
accumulate_central_diffs_pairwise(std::span<const double> values,
                                  std::optional<std::span<const uint8_t>> mask,
                                  const xsimd::batch<double, Arch> &mean_vector,
                                  int max_moment) {
  assert(!mask || mask->size() == values.size());

  using value_batch_type = xsimd::batch<double, Arch>;
  using mask_batch_type = xsimd::batch<uint8_t, Arch>;

  const std::size_t step =
      mask ? mask_batch_type::size : value_batch_type::size;
  const std::size_t leaf_size = 128 * step;

  if (values.size() <= leaf_size) {
    if (mask) {
      return accumulate_central_diffs_direct<Arch>(values, *mask, mean_vector,
                                                   max_moment);
    }
    return accumulate_central_diffs_direct<Arch>(values, mean_vector,
                                                 max_moment);
  }

  const std::size_t mid = (values.size() / 2 / step) * step;

  const auto left_mask =
      mask ? std::make_optional(mask->first(mid)) : std::nullopt;
  const auto left = accumulate_central_diffs_pairwise<Arch>(
      values.first(mid), left_mask, mean_vector, max_moment);

  const auto right_mask =
      mask ? std::make_optional(mask->last(mask->size() - mid)) : std::nullopt;
  const auto right = accumulate_central_diffs_pairwise<Arch>(
      values.last(values.size() - mid), right_mask, mean_vector, max_moment);

  return {
      .m1 = left.m1 + right.m1,
      .m2 = left.m2 + right.m2,
      .m3 = max_moment >= 3 ? left.m3 + right.m3 : value_batch_type{},
      .m4 = max_moment >= 4 ? left.m4 + right.m4 : value_batch_type{},
  };
}

// TODO: use std::expected when we start to use c++23
template <class Arch>
std::optional<MeanAcc<double, std::size_t>>
accumulate_mean(std::span<const double> values,
                std::optional<std::span<const uint8_t>> mask, bool skipna) {
  using value_batch_type = xsimd::batch<double, Arch>;
  using mask_batch_type = xsimd::batch<uint8_t, Arch>;

  const std::size_t step =
      mask ? mask_batch_type::size : value_batch_type::size;
  const std::size_t vec_size = (values.size() / step) * step;
  const auto vec_mask =
      mask ? std::make_optional(mask->first(vec_size)) : std::nullopt;
  const auto vec_acc =
      accumulate_mean_pairwise<Arch>(values.first(vec_size), vec_mask, skipna);
  if (!vec_acc) {
    return std::nullopt;
  }

  const std::optional<std::span<const uint8_t>> tail_mask =
      mask.has_value()
          ? std::make_optional(mask->last(values.size() - vec_size))
          : std::nullopt;

  const auto tail_total = accumulate_mean_scalar_direct(
      values.last(values.size() - vec_size), tail_mask, skipna);

  if (!tail_total) {
    return std::nullopt;
  }

  return MeanAcc<double, std::size_t>{
      .sum = xsimd::reduce_add(vec_acc->sum) + tail_total->sum,
      .count = xsimd::reduce_add(vec_acc->count) + tail_total->count};
}

template <class Arch>
CentralDiffs<double>
accumulate_central_diffs(std::span<const double> values,
                         std::optional<std::span<const uint8_t>> mask,
                         double mean, int max_moment) {
  using value_batch_type = xsimd::batch<double, Arch>;
  using mask_batch_type = xsimd::batch<uint8_t, Arch>;

  const auto mean_vector = xsimd::broadcast<double, Arch>(mean);

  const std::size_t step =
      mask ? mask_batch_type::size : value_batch_type::size;
  const std::size_t vec_size = (values.size() / step) * step;
  const auto vec_mask =
      mask ? std::make_optional(mask->first(vec_size)) : std::nullopt;
  const auto vec_acc = accumulate_central_diffs_pairwise<Arch>(
      values.first(vec_size), vec_mask, mean_vector, max_moment);

  const std::size_t tail_size = values.size() - vec_size;
  const auto tail_mask = mask.has_value()
                             ? std::make_optional(mask->last(tail_size))
                             : std::nullopt;

  const auto tail_acc = accumulate_central_diffs_scalar_direct(
      values.last(tail_size), tail_mask, mean, max_moment);

  return {
      .m1 = xsimd::reduce_add(vec_acc.m1) + tail_acc.m1,
      .m2 = xsimd::reduce_add(vec_acc.m2) + tail_acc.m2,
      .m3 = max_moment >= 3 ? xsimd::reduce_add(vec_acc.m3) + tail_acc.m3 : 0.0,
      .m4 = max_moment >= 4 ? xsimd::reduce_add(vec_acc.m4) + tail_acc.m4 : 0.0,
  };
}

template <class Arch>
Moments
accumulate_moments_simd_impl(std::span<const double> values,
                             std::optional<std::span<const uint8_t>> mask,
                             bool skipna, int max_moment) {
  const auto total_acc = accumulate_mean<Arch>(values, mask, skipna);

  if (!total_acc.has_value()) {
    return {.mean = NAN, .m2 = NAN, .m3 = NAN, .m4 = NAN, .n = values.size()};
  }

  if (total_acc->count == 0) {
    return {};
  }

  const auto count_double = static_cast<double>(total_acc->count);
  double trial_mean = total_acc->sum / count_double;

  const auto central_diffs =
      accumulate_central_diffs<Arch>(values, mask, trial_mean, max_moment);

  return compute_moments_with_correction(*total_acc, central_diffs, max_moment);
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
    std::optional<std::span<const uint8_t>> mask, int max_moment) noexcept {
  return detail::accumulate_moments_simd_impl<Arch>(values, mask, skipna,
                                                    max_moment);
}

} // namespace pandas::moments
