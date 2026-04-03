/*
Copyright (c) 2026, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#include "pandas/moments.h"
#include <math.h>
#include <stdint.h>

static inline void moments_add_valuem(Moments *moments, double val,
                                      int max_moment) {
  double delta = val - moments->mean;
  moments->n++;
  double n = (double)moments->n;
  double delta_n = delta / n;
  double term1 = delta * delta_n * (n - 1.0);

  if (max_moment >= 4) {
    moments->m4 +=
        delta_n *
        (-4.0 * moments->m3 +
         delta_n * (6.0 * moments->m2 + term1 * (n * n - 3.0 * n + 3.0)));
  }
  if (max_moment >= 3) {
    moments->m3 += delta_n * (term1 * (n - 2.0) - 3.0 * moments->m2);
  }
  moments->m2 += term1;
  moments->mean += delta_n;
}

Moments moments_merge(Moments a, Moments b, int max_moment) {

  if (a.n == 0) {
    return b;
  }
  if (b.n == 0) {
    return a;
  }

  Moments result;

  result.n = a.n + b.n;
  double n_a = (double)a.n;
  double n_b = (double)b.n;
  double delta = b.mean - a.mean;
  double delta_n = delta / (double)result.n;
  double term1 = delta * delta_n * n_a * n_b;

  if (max_moment >= 4) {
    result.m4 =
        a.m4 + b.m4 +
        delta_n * (4.0 * (n_a * b.m3 - n_b * a.m3) +
                   delta_n * (6.0 * (n_a * n_a * b.m2 + n_b * n_b * a.m2) +
                              term1 * (n_a * n_a - n_a * n_b + n_b * n_b)));
  }
  if (max_moment >= 3) {
    result.m3 =
        a.m3 + b.m3 +
        delta_n * (3.0 * (n_a * b.m2 - n_b * a.m2) + term1 * (n_a - n_b));
  }
  result.m2 = a.m2 + b.m2 + term1;
  result.mean = a.mean + delta_n * n_b;

  return result;
}

#ifndef __has_attribute
#  define __has_attribute(x) 0
#endif // __has_attribute

#if defined(__x86_64__) && __has_attribute(target_clones) && defined(__GLIBC__)
#  define PANDAS_SIMD_TARGETS __attribute__((target_clones("avx2", "default")))
#else
#  define PANDAS_SIMD_TARGETS
#endif // x86_64 + glibc + target_clones

/* --- SIMD Implementation --- */
#if __has_attribute(ext_vector_type)
typedef double v4d __attribute__((ext_vector_type(4), aligned(1)));
typedef long long v4si __attribute__((ext_vector_type(4), aligned(1)));
typedef uint8_t v4u8 __attribute__((ext_vector_type(4), aligned(1)));
#  define PANDAS_HAS_SIMD 1
#elif __has_attribute(vector_size)
typedef double v4d __attribute__((vector_size(4 * sizeof(double)), aligned(1)));
typedef long long v4si
    __attribute__((vector_size(4 * sizeof(long long)), aligned(1)));
typedef uint8_t v4u8
    __attribute__((vector_size(4 * sizeof(uint8_t)), aligned(1)));
#  define PANDAS_HAS_SIMD 1
#endif // __has_attribute(ext_vector_type)

#ifdef PANDAS_HAS_SIMD
/* Vector select: returns (mask ? a : b) where mask is all-ones or all-zeros */
#  define v_select(mask, a, b)                                                 \
    ((v4d)(((v4si)(mask) & (v4si)(a)) | (~(v4si)(mask) & (v4si)(b))))

PANDAS_SIMD_TARGETS
Moments accumulate_moments_simd(size_t n, const double *values, int skipna,
                                const uint8_t *mask, int max_moment) {
  v4d v_mean = {0.0, 0.0, 0.0, 0.0};
  v4d v_m2 = {0.0, 0.0, 0.0, 0.0};
  v4d v_m3 = {0.0, 0.0, 0.0, 0.0};
  v4d v_m4 = {0.0, 0.0, 0.0, 0.0};
  // Using double for v_n to avoid epi64 -> pd conversions.
  // The results are exact up to n < 2^53.
  v4d v_n = {0.0, 0.0, 0.0, 0.0};

  v4d v_one = {1.0, 1.0, 1.0, 1.0};
  v4d v_zero = {0.0, 0.0, 0.0, 0.0};
  v4si v_zerosi = {0, 0, 0, 0};
  v4d v_nan = {NAN, NAN, NAN, NAN};

  size_t i = 0;
  for (; i + 3 < n; i += 4) {
    v4d v_val = *(v4d *)(values + i);

    if (mask) {
      v4u8 mask_vec = *(v4u8 *)(mask + i);
      v4si mask_si = __builtin_convertvector(mask_vec, v4si);
      // replace val with NAN where mask is 1.
      v_val = v_select(mask_si != v_zerosi, v_nan, v_val);
    }

    // skip_mask is 1.0 if we should skip, 0.0 otherwise
    v4d v_is_nan = v_select(v_val == v_val, v_zero, v_one);
    v4d v_skip_mask = skipna ? v_is_nan : v_zero;

    // Increment 1 where we do NOT skip
    v4d v_n_increment = v_one - v_skip_mask;
    v_n += v_n_increment;
    v4d v_n_nonzero = v_select(v_n < v_one, v_one, v_n);

    v4d v_delta = v_val - v_mean;

    // replace delta with zero when skipping to don't update moments
    v_delta = v_select(v_skip_mask == v_zero, v_delta, v_zero);
    v4d v_delta_n = v_delta / v_n_nonzero;

    v4d v_delta_n2 = v_delta * v_delta_n;
    v4d v_term1 = v_delta_n2 * (v_n - v_one);

    if (max_moment >= 4) {
      v4d v_n2 = v_n * v_n;
      // n * n - 3.0 * n + 3.0
      v4d v_n_term = (v_n2 - (3.0 * v_n)) + 3.0;
      v_m4 += v_delta_n * ((-4.0 * v_m3) +
                           (v_delta_n * ((6.0 * v_m2) + (v_term1 * v_n_term))));
    }

    if (max_moment >= 3) {
      v_m3 += v_delta_n * ((v_term1 * (v_n - 2.0)) - (3.0 * v_m2));
    }

    v_m2 += v_term1;
    v_mean += v_delta_n;
  }

  double n_arrd[4], mean_arr[4], m2_arr[4], m3_arr[4], m4_arr[4];
  *(v4d *)n_arrd = v_n;
  *(v4d *)mean_arr = v_mean;
  *(v4d *)m2_arr = v_m2;
  *(v4d *)m3_arr = v_m3;
  *(v4d *)m4_arr = v_m4;

  int64_t n_arr[4];
  for (int j = 0; j < 4; j++) {
    n_arr[j] = (int64_t)n_arrd[j];
  }

  // Distribute remaining values across chunks
  for (; i < n; i++) {
    double val = values[i];
    if (mask && mask[i])
      val = NAN;
    if (skipna && isnan(val))
      continue;
    moments_add_value(val, &n_arr[i % 4], &mean_arr[i % 4], &m2_arr[i % 4],
                      &m3_arr[i % 4], &m4_arr[i % 4], max_moment);
  }
  Moments moments_arr[4];
  for (int j = 0; j < 4; j++) {
    moments_arr[j] = (Moments){.n = (uint64_t)n_arr[j],
                               .mean = mean_arr[j],
                               .m2 = m2_arr[j],
                               .m3 = m3_arr[j],
                               .m4 = m4_arr[j]};
  }

  // pairwise merge for numerical stability
  Moments m_01 = moments_merge(moments_arr[0], moments_arr[1], max_moment);
  Moments m_23 = moments_merge(moments_arr[2], moments_arr[3], max_moment);

  return moments_merge(m_01, m_23, max_moment);
}
#  undef v_select
#endif // PANDAS_HAS_SIMD

/* --- Scalar Fallback Implementation --- */

Moments accumulate_moments_scalar_block(size_t n, const double *values,
                                        int skipna, const uint8_t *mask,
                                        int max_moment) {
  Moments moments = {0};
  for (size_t i = 0; i < n; i++) {
    double val = values[i];
    if (mask && mask[i])
      val = NAN;
    if (skipna && isnan(val))
      continue;
    moments_add_valuem(&moments, val, max_moment);
  }
  return moments;
}

/* --- Accumulation Dispatch (Choose SIMD or Scalar) --- */

Moments accumulate_moments_dispatch(size_t n, const double *values, int skipna,
                                    const uint8_t *mask, int max_moment) {
#if defined(PANDAS_HAS_SIMD)
  return accumulate_moments_simd(n, values, skipna, mask, max_moment);
#endif // defined(PANDAS_HAS_SIMD)
  return accumulate_moments_scalar_block(n, values, skipna, mask, max_moment);
}

/* --- Moments 1D Accumulator Implementation --- */

void accumulate_moments_scalar(size_t n, const double *values, bool skipna,
                               const uint8_t *mask, int64_t *nobs, double *mean,
                               double *m2, double *m3, double *m4,
                               int max_moment) {
  // PERF: It's possible to parallelize moment reductions
  // and call `moments_merge` to join the results.
  Moments acc =
      accumulate_moments_dispatch(n, values, skipna, mask, max_moment);

  Moments tmp = (Moments){
      .n = (uint64_t)*nobs, .mean = *mean, .m2 = *m2, .m3 = 0.0, .m4 = 0.0};
  if (max_moment >= 4) {
    tmp.m4 = *m4;
  }
  if (max_moment >= 3) {
    tmp.m3 = *m3;
  }
  Moments result = moments_merge(tmp, acc, max_moment);

  if (max_moment >= 4) {
    *m4 = result.m4;
  }
  if (max_moment >= 3) {
    *m3 = result.m3;
  }
  *m2 = result.m2;
  *mean = result.mean;
  *nobs = (int64_t)result.n;
}
