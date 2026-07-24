/*
 * Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
 * All rights reserved.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * This header declares the small backend-selection layer for matmul block
 * acceleration. The portable naive kernel remains in blosc2_ext.pyx, while
 * platform BLAS backends such as Accelerate and runtime-discovered CBLAS
 * providers are exposed through this module.
 */

#ifndef PYTHON_BLOSC2_MATMUL_KERNELS_H
#define PYTHON_BLOSC2_MATMUL_KERNELS_H

#ifdef __cplusplus
extern "C" {
#endif

typedef enum b2_matmul_backend {
    B2_MATMUL_BACKEND_AUTO = 0,
    B2_MATMUL_BACKEND_NAIVE = 1,
    B2_MATMUL_BACKEND_ACCELERATE = 2,
    B2_MATMUL_BACKEND_CBLAS = 3,
} b2_matmul_backend;

int b2_has_accelerate(void);
int b2_has_cblas(void);
void b2_clear_cblas_candidates(void);
int b2_add_cblas_candidate(const char *path);
int b2_init_cblas(void);
void b2_set_matmul_backend(int backend);
int b2_get_matmul_backend(void);
int b2_get_selected_matmul_backend(void);
const char *b2_get_matmul_backend_name(void);
const char *b2_get_selected_matmul_backend_name(void);
const char *b2_get_loaded_cblas_path(void);

int b2_gemm_accelerate_f32(const float *A, const float *B, float *C, int M, int K, int N);
int b2_gemm_accelerate_f64(const double *A, const double *B, double *C, int M, int K, int N);
int b2_gemm_cblas_f32(const float *A, const float *B, float *C, int M, int K, int N);
int b2_gemm_cblas_f64(const double *A, const double *B, double *C, int M, int K, int N);

#ifdef __cplusplus
}
#endif

#endif
