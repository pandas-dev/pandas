/*
 * Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
 * All rights reserved.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

/* This module hosts the backend-selection and accelerator-specific GEMM
 * wrappers used by the matmul fast path. The portable naive kernel stays in
 * blosc2_ext.pyx, and external BLAS-style integrations live here.
 */

#include "matmul_kernels.h"
#include "blosc2.h"
#include <stdlib.h>
#include <string.h>

#if defined(_WIN32)
#include <windows.h>
#else
#include <dlfcn.h>
#endif

static int g_b2_matmul_backend = B2_MATMUL_BACKEND_AUTO;

typedef enum b2_cblas_order {
    B2_CBLAS_ROW_MAJOR = 101,
} b2_cblas_order;

typedef enum b2_cblas_transpose {
    B2_CBLAS_NO_TRANS = 111,
} b2_cblas_transpose;

typedef void (*b2_cblas_sgemm_fn)(
    int Order, int TransA, int TransB, int M, int N, int K,
    float alpha, const float *A, int lda, const float *B, int ldb,
    float beta, float *C, int ldc
);
typedef void (*b2_cblas_dgemm_fn)(
    int Order, int TransA, int TransB, int M, int N, int K,
    double alpha, const double *A, int lda, const double *B, int ldb,
    double beta, double *C, int ldc
);

#define B2_CBLAS_MAX_CANDIDATES 16

static char *g_b2_cblas_candidates[B2_CBLAS_MAX_CANDIDATES];
static int g_b2_cblas_ncandidates = 0;
static int g_b2_cblas_initialized = 0;
static int g_b2_cblas_available = 0;
static const char *g_b2_cblas_path = NULL;

#if defined(_WIN32)
static HMODULE g_b2_cblas_handle = NULL;
#else
static void *g_b2_cblas_handle = NULL;
#endif

static b2_cblas_sgemm_fn g_b2_cblas_sgemm = NULL;
static b2_cblas_dgemm_fn g_b2_cblas_dgemm = NULL;

static void b2_reset_cblas_state(void) {
    g_b2_cblas_initialized = 0;
    g_b2_cblas_available = 0;
    g_b2_cblas_path = NULL;
    g_b2_cblas_sgemm = NULL;
    g_b2_cblas_dgemm = NULL;
    if (g_b2_cblas_handle != NULL) {
#if defined(_WIN32)
        FreeLibrary(g_b2_cblas_handle);
#else
        dlclose(g_b2_cblas_handle);
#endif
        g_b2_cblas_handle = NULL;
    }
}

void b2_clear_cblas_candidates(void) {
    int i;
    for (i = 0; i < g_b2_cblas_ncandidates; ++i) {
        free(g_b2_cblas_candidates[i]);
        g_b2_cblas_candidates[i] = NULL;
    }
    g_b2_cblas_ncandidates = 0;
    b2_reset_cblas_state();
}

int b2_add_cblas_candidate(const char *path) {
    size_t len;
    char *copy;
    if (path == NULL || path[0] == '\0') {
        return -1;
    }
    if (g_b2_cblas_ncandidates >= B2_CBLAS_MAX_CANDIDATES) {
        return -1;
    }
    len = strlen(path);
    copy = (char *)malloc(len + 1);
    if (copy == NULL) {
        return -1;
    }
    memcpy(copy, path, len + 1);
    g_b2_cblas_candidates[g_b2_cblas_ncandidates++] = copy;
    BLOSC_TRACE_INFO("matmul cblas candidate: %s", copy);
    b2_reset_cblas_state();
    return 0;
}

int b2_init_cblas(void) {
    int i;
    if (g_b2_cblas_initialized) {
        return g_b2_cblas_available;
    }
    g_b2_cblas_initialized = 1;
    for (i = 0; i < g_b2_cblas_ncandidates; ++i) {
#if defined(_WIN32)
        HMODULE handle = LoadLibraryA(g_b2_cblas_candidates[i]);
        if (handle == NULL) {
            BLOSC_TRACE_INFO("matmul cblas reject: failed to load %s", g_b2_cblas_candidates[i]);
            continue;
        }
        g_b2_cblas_sgemm = (b2_cblas_sgemm_fn)GetProcAddress(handle, "cblas_sgemm");
        g_b2_cblas_dgemm = (b2_cblas_dgemm_fn)GetProcAddress(handle, "cblas_dgemm");
        if (g_b2_cblas_sgemm != NULL && g_b2_cblas_dgemm != NULL) {
            g_b2_cblas_handle = handle;
            g_b2_cblas_path = g_b2_cblas_candidates[i];
            g_b2_cblas_available = 1;
            BLOSC_TRACE_INFO("matmul cblas selected: %s", g_b2_cblas_path);
            return 1;
        }
        BLOSC_TRACE_INFO("matmul cblas reject: missing cblas_sgemm/cblas_dgemm in %s", g_b2_cblas_candidates[i]);
        FreeLibrary(handle);
#else
        void *handle = dlopen(g_b2_cblas_candidates[i], RTLD_NOW | RTLD_LOCAL);
        if (handle == NULL) {
            BLOSC_TRACE_INFO("matmul cblas reject: failed to load %s", g_b2_cblas_candidates[i]);
            continue;
        }
        g_b2_cblas_sgemm = (b2_cblas_sgemm_fn)dlsym(handle, "cblas_sgemm");
        g_b2_cblas_dgemm = (b2_cblas_dgemm_fn)dlsym(handle, "cblas_dgemm");
        if (g_b2_cblas_sgemm != NULL && g_b2_cblas_dgemm != NULL) {
            g_b2_cblas_handle = handle;
            g_b2_cblas_path = g_b2_cblas_candidates[i];
            g_b2_cblas_available = 1;
            BLOSC_TRACE_INFO("matmul cblas selected: %s", g_b2_cblas_path);
            return 1;
        }
        BLOSC_TRACE_INFO("matmul cblas reject: missing cblas_sgemm/cblas_dgemm in %s", g_b2_cblas_candidates[i]);
        dlclose(handle);
#endif
        g_b2_cblas_sgemm = NULL;
        g_b2_cblas_dgemm = NULL;
    }
    BLOSC_TRACE_INFO("matmul cblas unavailable after probing %d candidate(s)", g_b2_cblas_ncandidates);
    return 0;
}

int b2_has_accelerate(void) {
#if defined(__APPLE__)
    return 1;
#else
    return 0;
#endif
}

int b2_has_cblas(void) {
#if defined(__APPLE__) || defined(EMSCRIPTEN)
    return 0;
#else
    return b2_init_cblas();
#endif
}

void b2_set_matmul_backend(int backend) {
    switch (backend) {
        case B2_MATMUL_BACKEND_AUTO:
        case B2_MATMUL_BACKEND_NAIVE:
        case B2_MATMUL_BACKEND_ACCELERATE:
        case B2_MATMUL_BACKEND_CBLAS:
            g_b2_matmul_backend = backend;
            break;
        default:
            g_b2_matmul_backend = B2_MATMUL_BACKEND_AUTO;
            break;
    }
}

int b2_get_matmul_backend(void) {
    return g_b2_matmul_backend;
}

int b2_get_selected_matmul_backend(void) {
    if (g_b2_matmul_backend == B2_MATMUL_BACKEND_ACCELERATE && !b2_has_accelerate()) {
        BLOSC_TRACE_INFO("matmul backend fallback: accelerate unavailable, using naive");
        return B2_MATMUL_BACKEND_NAIVE;
    }
    if (g_b2_matmul_backend == B2_MATMUL_BACKEND_CBLAS && !b2_has_cblas()) {
        BLOSC_TRACE_INFO("matmul backend fallback: cblas unavailable, using naive");
        return B2_MATMUL_BACKEND_NAIVE;
    }
    if (g_b2_matmul_backend == B2_MATMUL_BACKEND_AUTO) {
        if (b2_has_accelerate()) {
            return B2_MATMUL_BACKEND_ACCELERATE;
        }
        if (b2_has_cblas()) {
            return B2_MATMUL_BACKEND_CBLAS;
        }
        BLOSC_TRACE_INFO("matmul backend fallback: auto found no accelerate/cblas backend, using naive");
        return B2_MATMUL_BACKEND_NAIVE;
    }
    return g_b2_matmul_backend;
}

const char *b2_get_loaded_cblas_path(void) {
    if (!b2_has_cblas()) {
        return NULL;
    }
    return g_b2_cblas_path;
}

const char *b2_get_matmul_backend_name(void) {
    switch (g_b2_matmul_backend) {
        case B2_MATMUL_BACKEND_NAIVE:
            return "naive";
        case B2_MATMUL_BACKEND_ACCELERATE:
            return "accelerate";
        case B2_MATMUL_BACKEND_CBLAS:
            return "cblas";
        case B2_MATMUL_BACKEND_AUTO:
        default:
            return "auto";
    }
}

const char *b2_get_selected_matmul_backend_name(void) {
    switch (b2_get_selected_matmul_backend()) {
        case B2_MATMUL_BACKEND_ACCELERATE:
            return "accelerate";
        case B2_MATMUL_BACKEND_CBLAS:
            return "cblas";
        case B2_MATMUL_BACKEND_NAIVE:
        default:
            return "naive";
    }
}

#if defined(__APPLE__)
#include <Accelerate/Accelerate.h>

int b2_gemm_accelerate_f32(const float *A, const float *B, float *C, int M, int K, int N) {
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0f, A, K, B, N, 1.0f, C, N);
    return 0;
}

int b2_gemm_accelerate_f64(const double *A, const double *B, double *C, int M, int K, int N) {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0, A, K, B, N, 1.0, C, N);
    return 0;
}
#else
int b2_gemm_accelerate_f32(const float *A, const float *B, float *C, int M, int K, int N) {
    (void)A;
    (void)B;
    (void)C;
    (void)M;
    (void)K;
    (void)N;
    return -1;
}

int b2_gemm_accelerate_f64(const double *A, const double *B, double *C, int M, int K, int N) {
    (void)A;
    (void)B;
    (void)C;
    (void)M;
    (void)K;
    (void)N;
    return -1;
}
#endif

int b2_gemm_cblas_f32(const float *A, const float *B, float *C, int M, int K, int N) {
    if (!b2_has_cblas()) {
        return -1;
    }
    g_b2_cblas_sgemm(
        B2_CBLAS_ROW_MAJOR, B2_CBLAS_NO_TRANS, B2_CBLAS_NO_TRANS,
        M, N, K, 1.0f, A, K, B, N, 1.0f, C, N
    );
    return 0;
}

int b2_gemm_cblas_f64(const double *A, const double *B, double *C, int M, int K, int N) {
    if (!b2_has_cblas()) {
        return -1;
    }
    g_b2_cblas_dgemm(
        B2_CBLAS_ROW_MAJOR, B2_CBLAS_NO_TRANS, B2_CBLAS_NO_TRANS,
        M, N, K, 1.0, A, K, B, N, 1.0, C, N
    );
    return 0;
}
