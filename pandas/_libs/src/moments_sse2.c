/*
Copyright (c) 2026, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#include "pandas/moments.h"
#include <emmintrin.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>

static __m128d get_nan_mask(const uint8_t *mask) {
  double m_hi = mask[1] == 0 ? 0.0 : NAN;
  double m_lo = mask[0] == 0 ? 0.0 : NAN;
  return _mm_set_pd(m_hi, m_lo);
}

static void update_accumulators(int max_moment, __m128d val, __m128d skip_mask,
                                __m128d *nobs, __m128d *mean, __m128d *m2,
                                __m128d *m3, __m128d *m4) {
  const __m128d one = _mm_set1_pd(1.0);
  const __m128d two = _mm_set1_pd(2.0);
  const __m128d three = _mm_set1_pd(3.0);
  const __m128d six = _mm_set1_pd(6.0);
  const __m128d minus_four = _mm_set1_pd(-4.0);

  const __m128d n_inc = _mm_andnot_pd(skip_mask, one);
  *nobs = _mm_add_pd(*nobs, n_inc);

  // avoid division by zero
  const __m128d nobs_nonzero = _mm_max_pd(*nobs, one);

  __m128d delta = _mm_sub_pd(val, *mean);
  // Zero out delta for skipped elements so they don't affect moments
  delta = _mm_andnot_pd(skip_mask, delta);

  const __m128d delta_n = _mm_div_pd(delta, nobs_nonzero);
  const __m128d term1 =
      _mm_mul_pd(delta, _mm_mul_pd(delta_n, _mm_sub_pd(*nobs, one)));

  if (max_moment >= 4) {
    const __m128d n2 = _mm_mul_pd(*nobs, *nobs);
    const __m128d n_term =
        _mm_add_pd(_mm_sub_pd(n2, _mm_mul_pd(three, *nobs)), three);
    const __m128d m4_update = _mm_mul_pd(
        delta_n,
        _mm_add_pd(_mm_mul_pd(minus_four, *m3),
                   _mm_mul_pd(delta_n, _mm_add_pd(_mm_mul_pd(six, *m2),
                                                  _mm_mul_pd(term1, n_term)))));
    *m4 = _mm_add_pd(*m4, m4_update);
  }

  if (max_moment >= 3) {
    const __m128d m3_update = _mm_mul_pd(
        delta_n, _mm_sub_pd(_mm_mul_pd(term1, _mm_sub_pd(*nobs, two)),
                            _mm_mul_pd(three, *m2)));
    *m3 = _mm_add_pd(*m3, m3_update);
  }

  *m2 = _mm_add_pd(*m2, term1);
  *mean = _mm_add_pd(*mean, delta_n);
}

static void store_accumulators_in_struct(Moments dst[2], __m128d nobs,
                                         __m128d mean, __m128d m2, __m128d m3,
                                         __m128d m4) {
  double n_arr[2], mean_arr[2], m2_arr[2], m3_arr[2], m4_arr[2];
  _mm_storeu_pd(n_arr, nobs);
  _mm_storeu_pd(mean_arr, mean);
  _mm_storeu_pd(m2_arr, m2);
  _mm_storeu_pd(m3_arr, m3);
  _mm_storeu_pd(m4_arr, m4);

  for (int j = 0; j < 2; j++) {
    dst[j].n = (uint64_t)n_arr[j];
    dst[j].mean = mean_arr[j];
    dst[j].m2 = m2_arr[j];
    dst[j].m3 = m3_arr[j];
    dst[j].m4 = m4_arr[j];
  }
}

Moments moments_reduce(size_t n, const double *values, bool skipna,
                       const uint8_t *mask, int max_moment) {
  __m128d nobs = _mm_setzero_pd();
  __m128d mean = _mm_setzero_pd();
  __m128d m2 = _mm_setzero_pd();
  __m128d m3 = _mm_setzero_pd();
  __m128d m4 = _mm_setzero_pd();

  size_t vecsize = n - (n % 2);
  for (size_t i = 0; i < vecsize; i += 2) {
    __m128d val = _mm_loadu_pd(values + i);

    if (mask) {
      // NaN will propagate to val
      const __m128d m_nan_mask = get_nan_mask(mask + i);
      val = _mm_add_pd(val, m_nan_mask);
    }

    __m128d is_nan = _mm_cmpneq_pd(val, val);

    if (!skipna) {
      const bool contains_nan = _mm_movemask_pd(is_nan) != 0;
      if (contains_nan) {
        return (Moments){(uint64_t)n, NAN, NAN, NAN, NAN};
      }
    }

    const __m128d skip_mask = is_nan;
    update_accumulators(max_moment, val, skip_mask, &nobs, &mean, &m2, &m3,
                        &m4);
  }

  Moments moments_arr[2];
  store_accumulators_in_struct(moments_arr, nobs, mean, m2, m3, m4);

  // Handle last element
  if (n % 2 == 1) {
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
