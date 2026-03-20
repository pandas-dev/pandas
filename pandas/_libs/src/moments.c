/*
Copyright (c) 2026, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#include "pandas/moments.h"
#include <math.h>
#include <stdint.h>

#ifdef _OPENMP
#  include <omp.h>
#else
static inline int omp_get_thread_num() { return 0; }
static inline int omp_get_num_threads() { return 1; }
#endif

#ifndef __has_attribute
#  define __has_attribute(x) 0
#endif

#if defined(__x86_64__) && __has_attribute(target_clones) && defined(__GLIBC__)
#  define PANDAS_SIMD_TARGETS __attribute__((target_clones("avx2", "default")))
#else
#  define PANDAS_SIMD_TARGETS
#endif

/* --- SIMD Implementation --- */
#if defined(__clang__) && __has_attribute(ext_vector_type)
typedef double v4d __attribute__((ext_vector_type(4), aligned(1)));
typedef long long v4si __attribute__((ext_vector_type(4), aligned(1)));
#  define PANDAS_HAS_SIMD 1
#elif defined(__GNUC__) && __has_attribute(vector_size)
typedef double v4d __attribute__((vector_size(32), aligned(1)));
typedef long long v4si __attribute__((vector_size(32), aligned(1)));
#  define PANDAS_HAS_SIMD 1
#endif

#ifdef PANDAS_HAS_SIMD
/* Vector select: returns (mask ? a : b) where mask is all-ones or all-zeros */
#  if defined(__clang__)
#    define v_select(mask, a, b) ((mask) ? (a) : (b))
#  else
// gcc doesn't have ternary operators when compiling C files.
#    define v_select(mask, a, b)                                               \
      ((v4d)(((v4si)(mask) & (v4si)(a)) | (~(v4si)(mask) & (v4si)(b))))
#  endif

PANDAS_SIMD_TARGETS Moments accumulate_moments_simd(const double *values,
                                                    int64_t n, int skipna,
                                                    const uint8_t *mask,
                                                    int max_moment) {
  v4d v_mean = {0.0, 0.0, 0.0, 0.0};
  v4d v_m2 = {0.0, 0.0, 0.0, 0.0};
  v4d v_m3 = {0.0, 0.0, 0.0, 0.0};
  v4d v_m4 = {0.0, 0.0, 0.0, 0.0};
  // Using double for v_n to avoid epi64 -> pd conversions.
  // The results are exact up to n < 2^53.
  v4d v_n = {0.0, 0.0, 0.0, 0.0};

  v4d v_one = {1.0, 1.0, 1.0, 1.0};
  v4d v_zero = {0.0, 0.0, 0.0, 0.0};
  v4d v_nan = {NAN, NAN, NAN, NAN};

  int64_t i = 0;
  for (; i < n - 4; i += 4) {
    v4d v_val = *(v4d *)(values + i);

    if (mask) {
      v4d v_m_val = {mask[i + 0] ? 1.0 : 0.0, mask[i + 1] ? 1.0 : 0.0,
                     mask[i + 2] ? 1.0 : 0.0, mask[i + 3] ? 1.0 : 0.0};
      // replace val with NAN where mask is 1.
      v_val = v_select(v_m_val == v_one, v_nan, v_val);
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

  double n_arr[4], mean_arr[4], m2_arr[4], m3_arr[4], m4_arr[4];
  *(v4d *)n_arr = v_n;
  *(v4d *)mean_arr = v_mean;
  *(v4d *)m2_arr = v_m2;
  *(v4d *)m3_arr = v_m3;
  *(v4d *)m4_arr = v_m4;

  Moments moments_arr[4];
  for (int j = 0; j < 4; j++) {
    moments_arr[j] = (Moments){(int64_t)round(n_arr[j]), mean_arr[j], m2_arr[j],
                               m3_arr[j], m4_arr[j]};
  }

  // Distribute remaining values across chunks
  for (; i < n; i++) {
    double val = values[i];
    if (mask && mask[i])
      val = NAN;
    if (skipna && isnan(val))
      continue;
    moments_add_value(&moments_arr[i % 4], val, max_moment);
  }

  // pairwise merge for numerical stability
  Moments m_01 = moments_merge(moments_arr[0], moments_arr[1], max_moment);
  Moments m_23 = moments_merge(moments_arr[2], moments_arr[3], max_moment);

  return moments_merge(m_01, m_23, max_moment);
}
#endif

/* --- Scalar Fallback Implementation --- */

static inline Moments accumulate_moments_scalar_block(const double *values,
                                                      int64_t n, int skipna,
                                                      const uint8_t *mask,
                                                      int max_moment) {
  Moments moments = {0};
  for (int64_t i = 0; i < n; i++) {
    double val = values[i];
    if (mask && mask[i])
      val = NAN;
    if (skipna && isnan(val))
      continue;
    moments_add_value(&moments, val, max_moment);
  }
  return moments;
}

/* --- Accumulation Dispatch (Choose SIMD or Scalar) --- */

static inline Moments accumulate_moments_dispatch(const double *values,
                                                  int64_t n, int skipna,
                                                  const uint8_t *mask,
                                                  int max_moment) {
#if defined(PANDAS_HAS_SIMD)
  return accumulate_moments_simd(values, n, skipna, mask, max_moment);
#endif
  return accumulate_moments_scalar_block(values, n, skipna, mask, max_moment);
}

/* --- Public API (Orchestrates OpenMP Parallelism) --- */

Moments accumulate_moments_scalar(const double *values, int64_t n, int skipna,
                                  const uint8_t *mask, int max_moment) {
  Moments result = {0};

#ifdef _OPENMP
  // parallel threshold chosen based on `perf report` where libgomp wasn't the
  // overhead
#  pragma omp parallel if (n > 10000)
#endif
  {
    int tid = omp_get_thread_num();
    int num_threads = omp_get_num_threads();

    int64_t start = (tid * n) / num_threads;
    int64_t end = tid == num_threads - 1 ? n : ((tid + 1) * n) / num_threads;

    Moments moments_local =
        accumulate_moments_dispatch(values + start, end - start, skipna,
                                    mask ? mask + start : NULL, max_moment);

#ifdef _OPENMP
#  pragma omp critical
#endif
    {
      result = moments_merge(result, moments_local, max_moment);
    }
  }

  return result;
}
