/*
Copyright (c) 2026, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#pragma once

#include "pandas/moments.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <optional>
#include <span>

/// Common functions used for simd and scalar versions.

namespace pandas::moments::detail {

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

static inline std::optional<MeanAcc<double, std::size_t>>
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

} // namespace pandas::moments::detail
