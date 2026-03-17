/*
Copyright (c) 2026, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#include "pandas/moments.h"
#include <math.h>
#include <stdint.h>
#include <string.h>

#ifdef _OPENMP
#  include <omp.h>
#else
static inline int omp_get_thread_num() { return 0; }
static inline int omp_get_num_threads() { return 1; }
#endif

/* --- AVX2 SIMD Implementation --- */

#ifdef __x86_64__
#  if defined(__clang__)
#    pragma clang attribute push(__attribute__((target("avx2"))),              \
                                 apply_to = function)
#    define PANDAS_HAS_AVX2 1
#  elif defined(__GNUC__)
#    pragma GCC push_options
#    pragma GCC target("avx2")
#    define PANDAS_HAS_AVX2 1
#  endif // defined(__GNUC__)
#endif   // ifdef __x86_64__

#ifdef PANDAS_HAS_AVX2
#  include <immintrin.h>

// Vectorization of [moments_add_value]
static inline Moments accumulate_moments_avx2(const double *values, int64_t n,
                                              int skipna, const uint8_t *mask,
                                              int max_moment) {
  __m256d v_mean = _mm256_setzero_pd();
  __m256d v_m2 = _mm256_setzero_pd();
  __m256d v_m3 = _mm256_setzero_pd();
  __m256d v_m4 = _mm256_setzero_pd();
  // Using double for v_n to avoid epi64 -> pd conversions.
  // The results are exact up to n < 2^53.
  __m256d v_n = _mm256_setzero_pd();

  __m256d v_one = _mm256_set1_pd(1.0);
  __m256d v_zero = _mm256_setzero_pd();
  __m256d v_nan = _mm256_set1_pd(NAN);

  int64_t i = 0;
  for (; i < n - 4; i += 4) {
    __m256d v_val = _mm256_loadu_pd(values + i);

    if (mask) {
      uint32_t m32;
      memcpy(&m32, mask + i, sizeof(uint32_t));
      __m128i v_mask_bytes = _mm_cvtsi32_si128((int)m32);
      __m256i v_mask_epi64 = _mm256_cvtepu8_epi64(v_mask_bytes);
      __m256i v_mask_is_zero =
          _mm256_cmpeq_epi64(v_mask_epi64, _mm256_setzero_si256());
      // replace val with NAN where mask is 1.
      v_val =
          _mm256_blendv_pd(v_nan, v_val, _mm256_castsi256_pd(v_mask_is_zero));
    }

    __m256d v_skip_mask =
        skipna ? _mm256_cmp_pd(v_val, v_val, _CMP_UNORD_Q) : v_zero;
    // Increment 1 where we do not skip
    __m256d v_n_increment = _mm256_andnot_pd(v_skip_mask, v_one);
    v_n = _mm256_add_pd(v_n, v_n_increment);
    __m256d v_n_nonzero = _mm256_max_pd(v_n, v_one);

    __m256d v_delta = _mm256_sub_pd(v_val, v_mean);

    // replace delta with zero when skipping to don't update moments
    v_delta = _mm256_blendv_pd(v_delta, v_zero, v_skip_mask);
    __m256d v_delta_n = _mm256_div_pd(v_delta, v_n_nonzero);

    __m256d v_delta_n2 = _mm256_mul_pd(v_delta, v_delta_n);
    __m256d v_term1 = _mm256_mul_pd(v_delta_n2, _mm256_sub_pd(v_n, v_one));

    if (max_moment >= 4) {
      __m256d v_n2 = _mm256_mul_pd(v_n, v_n);
      // n * n - 3.0 * n + 3.0
      __m256d v_n_term = _mm256_add_pd(
          _mm256_sub_pd(v_n2, _mm256_mul_pd(_mm256_set1_pd(3.0), v_n)),
          _mm256_set1_pd(3.0));
      __m256d v_m4_update = _mm256_mul_pd(
          v_delta_n,
          _mm256_add_pd(
              _mm256_mul_pd(_mm256_set1_pd(-4.0), v_m3),
              _mm256_mul_pd(
                  v_delta_n,
                  _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(6.0), v_m2),
                                _mm256_mul_pd(v_term1, v_n_term)))));
      v_m4 = _mm256_add_pd(v_m4, v_m4_update);
    }

    if (max_moment >= 3) {
      __m256d v_m3_term = _mm256_mul_pd(
          v_delta_n,
          _mm256_sub_pd(
              _mm256_mul_pd(v_term1, _mm256_sub_pd(v_n, _mm256_set1_pd(2.0))),
              _mm256_mul_pd(_mm256_set1_pd(3.0), v_m2)));
      v_m3 = _mm256_add_pd(v_m3, v_m3_term);
    }

    v_m2 = _mm256_add_pd(v_m2, v_term1);
    v_mean = _mm256_add_pd(v_mean, v_delta_n);
  }

  double n_arr[4], mean_arr[4], m2_arr[4], m3_arr[4], m4_arr[4];
  _mm256_storeu_pd(n_arr, v_n);
  _mm256_storeu_pd(mean_arr, v_mean);
  _mm256_storeu_pd(m2_arr, v_m2);
  _mm256_storeu_pd(m3_arr, v_m3);
  _mm256_storeu_pd(m4_arr, v_m4);

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

#ifdef __x86_64__
#  if defined(__clang__)
#    pragma clang attribute pop
#  elif defined(__GNUC__)
#    pragma GCC pop_options
#  endif
#endif // ifdef __x86_64__

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

/* --- Accumulation Dispatch (Choose AVX2 or Scalar) --- */

static inline Moments accumulate_moments_dispatch(const double *values,
                                                  int64_t n, int skipna,
                                                  const uint8_t *mask,
                                                  int max_moment) {
#if defined(PANDAS_HAS_AVX2) && (defined(__GNUC__) || defined(__clang__))
  if (__builtin_cpu_supports("avx2")) {
    return accumulate_moments_avx2(values, n, skipna, mask, max_moment);
  }
#endif
  return accumulate_moments_scalar_block(values, n, skipna, mask, max_moment);
}
#undef PANDAS_HAS_AVX2

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
