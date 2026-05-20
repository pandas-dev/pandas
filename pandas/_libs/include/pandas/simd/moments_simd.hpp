/*
Copyright (c) 2026, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#pragma once

#include "pandas/moments.h"
#include "pandas_simd_config.h"
#include "xsimd/xsimd.hpp"
#include <cassert>
#include <cmath>
#include <optional>
#include <span>

namespace pandas::moments {

namespace detail {

template <class A> struct MomentsBatch {
  xsimd::batch<double, A> mean, m2, m3, m4, nobs;
};

static inline void moments_add_value(Moments &moments, double val,
                                     int max_moment) {
  const auto delta = val - moments.mean;
  moments.n++;
  const auto n = static_cast<double>(moments.n);
  const auto delta_n = delta / n;
  const auto term1 = delta * delta_n * (n - 1.0);

  if (max_moment >= 4) {
    const auto m3_term = -4.0 * moments.m3;
    const auto m2_term = 6.0 * moments.m2;
    const auto m0_term = (n * (n - 3.0)) + 3.0;
    moments.m4 +=
        delta_n * (m3_term + (delta_n * (m2_term + (term1 * m0_term))));
  }
  if (max_moment >= 3) {
    const auto m2_term = -3.0 * moments.m2;
    const auto m0_term = n - 2.0;
    moments.m3 += delta_n * (m2_term + (term1 * m0_term));
  }
  moments.m2 += term1;
  moments.mean += delta_n;
}

template <class Arch>
inline void
update_moments_batch(MomentsBatch<Arch> &acc, xsimd::batch<double, Arch> val,
                     xsimd::batch_bool<double, Arch> nan_mask, int max_moment) {
  using batch_type = xsimd::batch<double, Arch>;
  const batch_type zero(0.0);
  const batch_type one(1.0);
  const batch_type three(3.0);

  const auto nobs_increment = xsimd::select(nan_mask, zero, one);
  acc.nobs += nobs_increment;

  const auto n_nonzero = xsimd::max(acc.nobs, one);
  const auto delta = xsimd::select(nan_mask, zero, val - acc.mean);
  const auto delta_n = delta / n_nonzero;
  const auto delta_n2 = delta * delta_n;
  const auto term1 = delta_n2 * (acc.nobs - one);

  if (max_moment >= 4) {
    const auto m3_term = batch_type(-4.0) * acc.m3;
    const auto m2_term = batch_type(6.0) * acc.m2;
    const auto m0_term = (acc.nobs * (acc.nobs - three)) + three;
    acc.m4 += delta_n * (m3_term + (delta_n * (m2_term + (term1 * m0_term))));
  }

  if (max_moment >= 3) {
    const auto m2_term = three * acc.m2;
    const auto m0_term = acc.nobs - batch_type(2.0);
    acc.m3 += delta_n * ((term1 * m0_term) - m2_term);
  }

  acc.m2 += term1;
  acc.mean += delta_n;
}

/// Merge results from moments accumulators.
/// It uses the formula for merging central moments:
/// $M_{p; N} = \sum_{k=1}^l \sum_{j=0}^p
/// \binom{p}{j} M_{p-j; k} (-\frac{\delta_k}{n})^j$
/// where $\delta_k = \sum_{j=1}^k n_j * (\bar{x}_j - \bar{x}_k)$.
template <class A>
Moments merge_batches(const MomentsBatch<A> &acc, int max_moment) {
  using batch_type = xsimd::batch<double, A>;
  constexpr std::size_t step = batch_type::size;
  Moments result{};

  const auto total_n = xsimd::reduce_add(acc.nobs);
  assert(total_n >= 0);
  result.n = static_cast<std::size_t>(total_n);

  if (result.n == 0) {
    return result;
  }

  auto mean_other = xsimd::rotate_left<1>(acc.mean);
  auto nobs_other = xsimd::rotate_left<1>(acc.nobs);
  batch_type delta(0.0);
  for (std::size_t i = 0; i + 1 < step; ++i) {
    delta += nobs_other * (mean_other - acc.mean);
    mean_other = xsimd::rotate_left<1>(mean_other);
    nobs_other = xsimd::rotate_left<1>(nobs_other);
  }

  const batch_type total_n_v(total_n);
  const auto delta_n = delta / total_n_v;
  const auto delta2_n = delta_n * delta_n;

  if (max_moment >= 4) {
    const auto m3_term = batch_type(-4.0) * acc.m3;
    const auto m2_term = batch_type(6.0) * acc.m2;
    const auto m4_acc =
        acc.m4 +
        (delta_n * (m3_term + (delta_n * (m2_term + (delta2_n * acc.nobs)))));
    result.m4 = xsimd::reduce_add(m4_acc);
  }

  if (max_moment >= 3) {
    const auto m2_term = batch_type(-3.0) * acc.m2;
    const auto m3_acc = acc.m3 + (delta_n * (m2_term - (acc.nobs * delta2_n)));
    result.m3 = xsimd::reduce_add(m3_acc);
  }

  const auto m2_acc = acc.m2 + (acc.nobs * delta2_n);
  result.m2 = xsimd::reduce_add(m2_acc);

  const auto mean_v = acc.mean + delta_n;
  result.mean = mean_v.first();

  [[maybe_unused]] constexpr double rtol = 1e-12;
  [[maybe_unused]] constexpr double atol = 1e-8;
  assert((xsimd::all(xsimd::isnan(mean_v)) ||
          xsimd::all(xsimd::abs(mean_v - result.mean) <
                     ((rtol * xsimd::abs(result.mean)) + atol))) &&
         "mean lanes aren't homogeneous after merge");

  return result;
}

/// https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Higher-order_statistics
static inline void moments_merge(Moments &acc, const Moments &src,
                                 int max_moment) {
  if (acc.n == 0) {
    acc = src;
    return;
  }
  if (src.n == 0) {
    return;
  }

  const auto n_a = static_cast<double>(acc.n);
  const auto n_b = static_cast<double>(src.n);
  acc.n += src.n;

  const auto delta = src.mean - acc.mean;
  const auto delta_n = delta / static_cast<double>(acc.n);
  const auto term1 = delta * delta_n * n_a * n_b;

  if (max_moment >= 4) {
    const auto m3_term = 4.0 * ((n_a * src.m3) - (n_b * acc.m3));
    const auto m2_term = 6.0 * ((n_a * n_a * src.m2) + (n_b * n_b * acc.m2));
    const auto m0_term = (n_a * n_a) - (n_a * n_b) + (n_b * n_b);
    acc.m4 += src.m4 +
              (delta_n * (m3_term + (delta_n * (m2_term + (term1 * m0_term)))));
  }

  if (max_moment >= 3) {
    const auto m2_term = 3.0 * ((n_a * src.m2) - (n_b * acc.m2));
    const auto m0_term = n_a - n_b;
    acc.m3 += src.m3 + (delta_n * (m2_term + (term1 * m0_term)));
  }

  acc.m2 += src.m2 + term1;
  acc.mean += delta_n * n_b;
}

template <class A>
static inline void set_moments_nan(MomentsBatch<A> &acc, std::size_t total_n) {
  using batch_type = xsimd::batch<double, A>;
  acc.nobs = batch_type(static_cast<double>(total_n));
  acc.mean = batch_type(NAN);
  acc.m2 = batch_type(NAN);
  acc.m3 = batch_type(NAN);
  acc.m4 = batch_type(NAN);
}

template <class Arch>
void accumulate_moments_simd_impl(MomentsBatch<Arch> &acc, int max_moment,
                                  std::span<const double> values, bool skipna) {
  using batch_type = xsimd::batch<double, Arch>;
  constexpr std::size_t step = batch_type::size;
  assert(values.size() % step == 0);

  for (std::size_t i = 0; i < values.size(); i += step) {
    auto val = xsimd::load_unaligned<Arch>(&values[i]);
    auto nan_mask = xsimd::isnan(val);

    if (!skipna && xsimd::any(nan_mask)) [[unlikely]] {
      const std::size_t nobs_per_lane = values.size() / step;
      set_moments_nan(acc, nobs_per_lane);
      return;
    }

    detail::update_moments_batch(acc, val, nan_mask, max_moment);
  }
}

template <class A>
void accumulate_moments_simd_masked_impl(detail::MomentsBatch<A> &acc,
                                         int max_moment,
                                         std::span<const double> values,
                                         std::span<const uint8_t> mask,
                                         bool skipna) {
  using mask_batch_type = xsimd::batch<uint8_t, A>;
  using value_batch_type = xsimd::batch<double, A>;
  constexpr std::size_t mask_step = mask_batch_type::size;
  constexpr std::size_t val_step = value_batch_type::size;

  assert(values.size() == mask.size());
  assert(mask.size() % mask_step == 0);

  std::size_t left = 0;
  for (std::size_t right = 0; right < mask.size(); right += mask_step) {
    const mask_batch_type mask_batch = xsimd::load_unaligned<A>(&mask[right]);
    const auto is_masked = mask_batch != mask_batch_type(0U);

    if (!xsimd::any(is_masked)) {
      continue;
    }

    if (!skipna) {
      const std::size_t nobs_per_lane = values.size() / val_step;
      set_moments_nan(acc, nobs_per_lane);
      return;
    }

    // NaN values aren't skipped when there is a mask
    accumulate_moments_simd_impl<A>(acc, max_moment,
                                    values.subspan(left, right - left),
                                    /*skipna=*/false);

    const std::uint64_t is_masked_bitmask = is_masked.mask();
    static_assert(mask_step % val_step == 0);

    const std::uint64_t lane_values_mask = ((1 << val_step) - 1);
    for (std::size_t i = 0; i < mask_step; i += val_step) {
      const std::uint64_t lane_mask_bits =
          (is_masked_bitmask >> i) & lane_values_mask;
      const auto isna_pd =
          xsimd::batch_bool<double, A>::from_mask(lane_mask_bits);

      const value_batch_type val = xsimd::load_unaligned<A>(&values[right + i]);

      detail::update_moments_batch(acc, val, isna_pd, max_moment);
    }

    left = right + mask_step;
  }

  accumulate_moments_simd_impl<A>(acc, max_moment,
                                  values.last(values.size() - left), false);
}
} // namespace detail

struct accumulate_moments_simd {
  template <class Arch>
  Moments operator()(Arch, std::span<const double> values, bool skipna,
                     std::optional<std::span<const uint8_t>> mask,
                     int max_moment) noexcept;
};

template <>
inline Moments accumulate_moments_simd::operator()<xsimd::common>(
    xsimd::common, std::span<const double> values, bool skipna,
    std::optional<std::span<const uint8_t>> mask, int max_moment) noexcept {
  Moments acc{};
  for (std::size_t i = 0; i < values.size(); i++) {
    const auto val = values[i];
    const auto isna_entry =
        mask.has_value() ? (*mask)[i] != 0 : std::isnan(val);

    if (skipna && isna_entry) {
      continue;
    }
    if (isna_entry) [[unlikely]] {
      acc = {.mean = NAN, .m2 = NAN, .m3 = NAN, .m4 = NAN, .n = values.size()};
      break;
    }
    detail::moments_add_value(acc, val, max_moment);
  }
  return acc;
}

template <class Arch>
Moments accumulate_moments_simd::operator()(
    Arch, std::span<const double> values, bool skipna,
    std::optional<std::span<const uint8_t>> mask, int max_moment) noexcept {
  using values_batch_type = xsimd::batch<double, Arch>;
  detail::MomentsBatch<Arch> acc_simd{};

  std::size_t vec_size;
  std::size_t tail_size;

  if (mask.has_value()) {
    using mask_batch_type = xsimd::batch<uint8_t, Arch>;
    constexpr std::size_t batch_size = mask_batch_type::size;
    tail_size = values.size() % batch_size;
    vec_size = values.size() - tail_size;

    detail::accumulate_moments_simd_masked_impl<Arch>(
        acc_simd, max_moment, values.first(vec_size), mask->first(vec_size),
        skipna);
  } else {
    constexpr std::size_t batch_size = values_batch_type::size;
    tail_size = values.size() % batch_size;
    vec_size = values.size() - tail_size;

    detail::accumulate_moments_simd_impl<Arch>(acc_simd, max_moment,
                                               values.first(vec_size), skipna);
  }

  Moments moments_acc = detail::merge_batches(acc_simd, max_moment);

  auto values_tail = values.last(tail_size);
  std::optional<std::span<const uint8_t>> mask_tail{};
  if (mask.has_value()) {
    mask_tail = mask->last(tail_size);
  }

  Moments tail = accumulate_moments_simd{}(xsimd::common{}, values_tail, skipna,
                                           mask_tail, max_moment);
  detail::moments_merge(moments_acc, tail, max_moment);

  return moments_acc;
}

} // namespace pandas::moments
