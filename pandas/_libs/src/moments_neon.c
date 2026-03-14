/*
Copyright (c) 2026, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#include "pandas/moments.h"
#include <arm_neon.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>

Moments moments_reduce(size_t n, const double *values, bool skipna,
                       const uint8_t *mask, int max_moment) {
  float64x2_t mean = vdupq_n_f64(0.0);
  float64x2_t m2 = vdupq_n_f64(0.0);
  float64x2_t m3 = vdupq_n_f64(0.0);
  float64x2_t m4 = vdupq_n_f64(0.0);
  float64x2_t nobs = vdupq_n_f64(0.0);

  float64x2_t one = vdupq_n_f64(1.0);
  float64x2_t two = vdupq_n_f64(2.0);
  float64x2_t three = vdupq_n_f64(3.0);
  float64x2_t six = vdupq_n_f64(6.0);
  float64x2_t minus_four = vdupq_n_f64(-4.0);

  size_t vecsize = n - (n % 2);
  for (size_t i = 0; i < vecsize; i += 2) {
    float64x2_t val = vld1q_f64(values + i);

    if (mask) {
      double m_hi = mask[i + 1] == 0 ? 0.0 : NAN;
      double m_lo = mask[i] == 0 ? 0.0 : NAN;
      float64x2_t m_nan_mask = {m_lo, m_hi};
      // NaN will propagate to val
      val = vaddq_f64(val, m_nan_mask);
    }

    // vceqq_f64 returns all ones if val == val (not NaN)
    uint64x2_t is_valid_mask = vceqq_f64(val, val);

    if (!skipna) {
      if (!vgetq_lane_u64(is_valid_mask, 0) ||
          !vgetq_lane_u64(is_valid_mask, 1)) {
        return (Moments){(uint64_t)n, NAN, NAN, NAN, NAN};
      }
    }

    float64x2_t n_inc = vreinterpretq_f64_u64(
        vandq_u64(is_valid_mask, vreinterpretq_u64_f64(one)));
    nobs = vaddq_f64(nobs, n_inc);

    // avoid division by zero
    float64x2_t nobs_nonzero = vmaxq_f64(nobs, one);

    float64x2_t delta = vsubq_f64(val, mean);
    // Zero out delta for skipped elements so they don't affect moments
    delta = vreinterpretq_f64_u64(
        vandq_u64(is_valid_mask, vreinterpretq_u64_f64(delta)));

    float64x2_t delta_n = vdivq_f64(delta, nobs_nonzero);
    float64x2_t term1 =
        vmulq_f64(delta, vmulq_f64(delta_n, vsubq_f64(nobs, one)));

    if (max_moment >= 4) {
      float64x2_t n2 = vmulq_f64(nobs, nobs);
      float64x2_t n_term =
          vaddq_f64(vsubq_f64(n2, vmulq_f64(three, nobs)), three);
      float64x2_t m4_update = vmulq_f64(
          delta_n,
          vaddq_f64(vmulq_f64(minus_four, m3),
                    vmulq_f64(delta_n, vaddq_f64(vmulq_f64(six, m2),
                                                 vmulq_f64(term1, n_term)))));
      m4 = vaddq_f64(m4, m4_update);
    }

    if (max_moment >= 3) {
      float64x2_t m3_update =
          vmulq_f64(delta_n, vsubq_f64(vmulq_f64(term1, vsubq_f64(nobs, two)),
                                       vmulq_f64(three, m2)));
      m3 = vaddq_f64(m3, m3_update);
    }

    m2 = vaddq_f64(m2, term1);
    mean = vaddq_f64(mean, delta_n);
  }

  double n_arr[2], mean_arr[2], m2_arr[2], m3_arr[2], m4_arr[2];
  vst1q_f64(n_arr, nobs);
  vst1q_f64(mean_arr, mean);
  vst1q_f64(m2_arr, m2);
  vst1q_f64(m3_arr, m3);
  vst1q_f64(m4_arr, m4);

  Moments moments_arr[2];
  for (int j = 0; j < 2; j++) {
    moments_arr[j] = (Moments){(uint64_t)n_arr[j], mean_arr[j], m2_arr[j],
                               m3_arr[j], m4_arr[j]};
  }

  // Handle last element
  if (n & 1ULL) {
    double val = values[vecsize];
    if (mask && mask[vecsize])
      val = NAN;

    if (!isnan(val)) {
      moments_add_valuem(&moments_arr[1], val, max_moment);
    } else if (!skipna) {
      return (Moments){(uint64_t)n, NAN, NAN, NAN, NAN};
    }
  }

  moments_merge(&moments_arr[0], &moments_arr[1], max_moment);
  return moments_arr[0];
}
