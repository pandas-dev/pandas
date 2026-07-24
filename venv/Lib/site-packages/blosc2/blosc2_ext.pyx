#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################

#cython: language_level=3

import glob
import os
import dataclasses
import ast
import atexit
import pathlib
import sys
import time
import warnings

import _ctypes

import cython
from cpython cimport (
    Py_buffer,
    PyBUF_SIMPLE,
    PyBuffer_Release,
    PyBytes_FromStringAndSize,
    PyObject_GetBuffer,
)
from cpython.ref cimport Py_INCREF, Py_DECREF
from cpython.pycapsule cimport PyCapsule_GetPointer, PyCapsule_New
from cython.operator cimport dereference
from libc.stdint cimport uintptr_t
from libc.stdlib cimport free, malloc, realloc, calloc
from libc.stdlib cimport abs as c_abs
from libc.string cimport memcpy, memset, strcpy, strdup, strlen
from libcpp cimport bool as c_bool

from enum import Enum

import numpy as np
from msgpack import packb, unpackb

import blosc2

cimport numpy as np

np.import_array()

cdef extern from "<stdint.h>":
    ctypedef   signed char  int8_t
    ctypedef   signed short int16_t
    ctypedef   signed int   int32_t
    ctypedef   signed long  int64_t
    ctypedef unsigned char  uint8_t
    ctypedef unsigned short uint16_t
    ctypedef unsigned int   uint32_t
    ctypedef unsigned long long uint64_t

ctypedef fused T:
    float
    double
    int32_t
    int64_t


cdef extern from "<stdio.h>":
    int printf(const char *format, ...) nogil

cdef extern from "matmul_kernels.h":
    ctypedef enum b2_matmul_backend:
        B2_MATMUL_BACKEND_AUTO
        B2_MATMUL_BACKEND_NAIVE
        B2_MATMUL_BACKEND_ACCELERATE
        B2_MATMUL_BACKEND_CBLAS

    int b2_has_accelerate() nogil
    int b2_has_cblas() nogil
    void b2_clear_cblas_candidates()
    int b2_add_cblas_candidate(const char *path)
    int b2_init_cblas()
    void b2_set_matmul_backend(int backend) nogil
    int b2_get_matmul_backend() nogil
    int b2_get_selected_matmul_backend() nogil
    const char *b2_get_matmul_backend_name() nogil
    const char *b2_get_selected_matmul_backend_name() nogil
    const char *b2_get_loaded_cblas_path() nogil
    int b2_gemm_accelerate_f32(const float *a, const float *b, float *c, int m, int k, int n) nogil
    int b2_gemm_accelerate_f64(const double *a, const double *b, double *c, int m, int k, int n) nogil
    int b2_gemm_cblas_f32(const float *a, const float *b, float *c, int m, int k, int n) nogil
    int b2_gemm_cblas_f64(const double *a, const double *b, double *c, int m, int k, int n) nogil


def _discover_matmul_cblas_candidates():
    if sys.platform == "darwin":
        return []

    prefix = pathlib.Path(sys.prefix)
    if sys.platform.startswith("win"):
        libdirs = [prefix / "Library" / "bin", prefix / "Library" / "lib", prefix / "DLLs"]
        patterns = [
            "mkl_rt.dll",
            "libopenblas*.dll",
            "openblas*.dll",
            "cblas.dll",
            "blas.dll",
        ]
    else:
        libdirs = [prefix / "lib", prefix / "lib64"]
        patterns = [
            "libcblas.so",
            "libcblas.so.*",
            "libopenblas.so",
            "libopenblas.so.*",
            "libflexiblas.so",
            "libflexiblas.so.*",
            "libblis.so",
            "libblis.so.*",
            "libmkl_rt.so",
            "libmkl_rt.so.*",
            "libblas.so",
            "libblas.so.*",
        ]

    try:
        config = np.show_config(mode="dicts")
        blas_cfg = config.get("Build Dependencies", {}).get("blas", {})
        libdir = blas_cfg.get("lib directory")
        if libdir:
            libdirs.insert(0, pathlib.Path(libdir))
    except Exception:
        pass

    candidates = []
    seen = set()
    for libdir in libdirs:
        if not libdir.exists():
            continue
        for pattern in patterns:
            for match in sorted(glob.glob(str(libdir / pattern))):
                resolved = pathlib.Path(match).resolve()
                path = str(resolved)
                if path not in seen and resolved.exists():
                    seen.add(path)
                    candidates.append(path)
    return candidates


def _configure_matmul_cblas_backend():
    cdef bytes path_bytes

    b2_clear_cblas_candidates()
    for path in _discover_matmul_cblas_candidates():
        path_bytes = os.fsencode(path)
        b2_add_cblas_candidate(path_bytes)
    b2_init_cblas()


_configure_matmul_cblas_backend()

cdef extern from "blosc2.h":

    ctypedef enum:
        BLOSC2_MAX_FILTERS
        BLOSC2_DEFINED_FILTERS_START
        BLOSC2_DEFINED_FILTERS_STOP
        BLOSC2_GLOBAL_REGISTERED_FILTERS_START
        BLOSC2_GLOBAL_REGISTERED_FILTERS_STOP
        BLOSC2_GLOBAL_REGISTERED_FILTERS
        BLOSC2_USER_REGISTERED_FILTERS_START
        BLOSC2_USER_REGISTERED_FILTERS_STOP
        BLOSC2_MAX_UDFILTERS
        BLOSC2_MAX_METALAYERS
        BLOSC2_MAX_VLMETALAYERS
        BLOSC2_PREFILTER_INPUTS_MAX
        BLOSC_MAX_CODECS
        BLOSC_MIN_HEADER_LENGTH
        BLOSC_EXTENDED_HEADER_LENGTH
        BLOSC2_MAX_OVERHEAD
        BLOSC2_MAX_BUFFERSIZE
        BLOSC2_MAXBLOCKSIZE
        BLOSC2_MAXTYPESIZE
        BLOSC_MAX_TYPESIZE
        BLOSC_MIN_BUFFERSIZE

    ctypedef enum:
        BLOSC2_SPECIAL_ZERO
        BLOSC2_SPECIAL_NAN
        BLOSC2_SPECIAL_UNINIT

    ctypedef enum:
        BLOSC2_VERSION_STRING
        BLOSC2_VERSION_REVISION
        BLOSC2_VERSION_DATE

    ctypedef enum:
        BLOSC2_ERROR_SUCCESS
        BLOSC2_ERROR_FAILURE
        BLOSC2_ERROR_STREAM
        BLOSC2_ERROR_DATA
        BLOSC2_ERROR_MEMORY_ALLOC
        BLOSC2_ERROR_READ_BUFFER
        BLOSC2_ERROR_WRITE_BUFFER
        BLOSC2_ERROR_CODEC_SUPPORT
        BLOSC2_ERROR_CODEC_PARAM
        BLOSC2_ERROR_CODEC_DICT
        BLOSC2_ERROR_VERSION_SUPPORT
        BLOSC2_ERROR_INVALID_HEADER
        BLOSC2_ERROR_INVALID_PARAM
        BLOSC2_ERROR_FILE_READ
        BLOSC2_ERROR_FILE_WRITE
        BLOSC2_ERROR_FILE_OPEN
        BLOSC2_ERROR_NOT_FOUND
        BLOSC2_ERROR_RUN_LENGTH
        BLOSC2_ERROR_FILTER_PIPELINE
        BLOSC2_ERROR_CHUNK_INSERT
        BLOSC2_ERROR_CHUNK_APPEND
        BLOSC2_ERROR_CHUNK_UPDATE
        BLOSC2_ERROR_2GB_LIMIT
        BLOSC2_ERROR_SCHUNK_COPY
        BLOSC2_ERROR_FRAME_TYPE
        BLOSC2_ERROR_FILE_TRUNCATE
        BLOSC2_ERROR_THREAD_CREATE
        BLOSC2_ERROR_POSTFILTER
        BLOSC2_ERROR_FRAME_SPECIAL
        BLOSC2_ERROR_SCHUNK_SPECIAL
        BLOSC2_ERROR_PLUGIN_IO
        BLOSC2_ERROR_FILE_REMOVE

    ctypedef enum:
        BLOSC2_DEFINED_CODECS_START
        BLOSC2_DEFINED_CODECS_STOP
        BLOSC2_GLOBAL_REGISTERED_CODECS_START
        BLOSC2_GLOBAL_REGISTERED_CODECS_STOP
        BLOSC2_GLOBAL_REGISTERED_CODECS
        BLOSC2_USER_REGISTERED_CODECS_START
        BLOSC2_USER_REGISTERED_CODECS_STOP

    ctypedef enum:
        BLOSC2_IO_FILESYSTEM
        BLOSC2_IO_FILESYSTEM_MMAP
        BLOSC_IO_LAST_BLOSC_DEFINED
        BLOSC_IO_LAST_REGISTERED

    cdef int INT_MAX

    void blosc2_init()
    void blosc2_destroy()

    int blosc1_compress(int clevel, int doshuffle, size_t typesize,
                       size_t nbytes, const void* src, void* dest,
                       size_t destsize)

    int blosc1_decompress(const void* src, void* dest, size_t destsize)

    int blosc1_getitem(const void* src, int start, int nitems, void* dest)

    int blosc2_getitem(const void* src, int32_t srcsize, int start, int nitems,
                       void* dest, int32_t destsize)

    ctypedef void(*blosc2_threads_callback)(void *callback_data, void (*dojob)(void *), int numjobs,
                                            size_t jobdata_elsize, void *jobdata)

    void blosc2_set_threads_callback(blosc2_threads_callback callback, void *callback_data)

    int16_t blosc2_set_nthreads(int16_t nthreads)

    const char* blosc1_get_compressor()

    int blosc1_set_compressor(const char* compname)

    void blosc2_set_delta(int dodelta)

    int blosc2_compcode_to_compname(int compcode, const char** compname)

    int blosc2_compname_to_compcode(const char* compname)

    const char* blosc2_list_compressors()

    int blosc2_get_complib_info(const char* compname, char** complib,
                               char** version)

    int blosc2_free_resources()

    int blosc2_cbuffer_sizes(const void* cbuffer, int32_t* nbytes,
                             int32_t* cbytes, int32_t* blocksize) nogil

    int blosc1_cbuffer_validate(const void* cbuffer, size_t cbytes, size_t* nbytes)

    void blosc1_cbuffer_metainfo(const void* cbuffer, size_t* typesize, int* flags)

    void blosc1_cbuffer_versions(const void* cbuffer, int* version, int* versionlz)

    const char* blosc2_cbuffer_complib(const void* cbuffer)


    ctypedef struct blosc2_context:
        pass

    ctypedef struct blosc2_prefilter_params:
        void* user_data
        const uint8_t* input
        uint8_t* output
        int32_t output_size
        int32_t output_typesize
        int32_t output_offset
        int64_t nchunk
        int32_t nblock
        int32_t tid
        uint8_t* ttmp
        size_t ttmp_nbytes
        blosc2_context* ctx
        c_bool output_is_disposable

    ctypedef struct blosc2_postfilter_params:
        void *user_data
        const uint8_t *input
        uint8_t *output
        int32_t size
        int32_t typesize
        int32_t offset
        int64_t nchunk
        int32_t nblock
        int32_t tid
        uint8_t *ttmp
        size_t ttmp_nbytes
        blosc2_context *ctx

    ctypedef int(*blosc2_prefilter_fn)(blosc2_prefilter_params* params)

    ctypedef int(*blosc2_postfilter_fn)(blosc2_postfilter_params *params)

    ctypedef struct blosc2_cparams:
        uint8_t compcode
        uint8_t compcode_meta
        uint8_t clevel
        int use_dict
        int32_t typesize
        int16_t nthreads
        int32_t blocksize
        int32_t splitmode
        void *schunk
        uint8_t filters[BLOSC2_MAX_FILTERS]
        uint8_t filters_meta[BLOSC2_MAX_FILTERS]
        blosc2_prefilter_fn prefilter
        blosc2_prefilter_params* preparams
        int tuner_id
        void* tuner_params
        c_bool instr_codec
        void* codec_params
        void* filter_params[BLOSC2_MAX_FILTERS]

    cdef const blosc2_cparams BLOSC2_CPARAMS_DEFAULTS

    ctypedef struct blosc2_dparams:
        int16_t nthreads
        void* schunk
        blosc2_postfilter_fn postfilter
        blosc2_postfilter_params *postparams
        int32_t typesize

    cdef const blosc2_dparams BLOSC2_DPARAMS_DEFAULTS

    blosc2_context* blosc2_create_cctx(blosc2_cparams cparams) nogil

    blosc2_context* blosc2_create_dctx(blosc2_dparams dparams) nogil

    void blosc2_free_ctx(blosc2_context * context) nogil

    int blosc2_set_maskout(blosc2_context *ctx, c_bool *maskout, int nblocks)


    int blosc2_compress(int clevel, int doshuffle, int32_t typesize,
                        const void * src, int32_t srcsize, void * dest,
                        int32_t destsize) nogil

    int blosc2_decompress(const void * src, int32_t srcsize,
                          void * dest, int32_t destsize)

    int blosc2_compress_ctx(
            blosc2_context * context, const void * src, int32_t srcsize, void * dest,
            int32_t destsize) nogil

    int blosc2_vlcompress_ctx(
            blosc2_context * context, const void * const * srcs, const int32_t * srcsizes,
            int32_t nblocks, void * dest, int32_t destsize) nogil

    int blosc2_decompress_ctx(blosc2_context * context, const void * src,
                              int32_t srcsize, void * dest, int32_t destsize) nogil

    int blosc2_vldecompress_ctx(blosc2_context* context, const void* src,
                                int32_t srcsize, void** dests,
                                int32_t* destsizes, int32_t maxblocks)

    int blosc2_getitem_ctx(blosc2_context* context, const void* src,
                           int32_t srcsize, int start, int nitems, void* dest,
                           int32_t destsize) nogil



    ctypedef struct blosc2_storage:
        c_bool contiguous
        char* urlpath
        blosc2_cparams* cparams
        blosc2_dparams* dparams
        blosc2_io *io

    cdef const blosc2_storage BLOSC2_STORAGE_DEFAULTS

    ctypedef struct blosc2_frame:
        pass

    ctypedef struct blosc2_metalayer:
        char* name
        uint8_t* content
        int32_t content_len


    ctypedef struct blosc2_tuner:
        void(*init)(void *config, blosc2_context*cctx, blosc2_context*dctx)
        void (*next_blocksize)(blosc2_context *context)
        void(*next_cparams)(blosc2_context *context)
        void(*update)(blosc2_context *context, double ctime)
        void (*free)(blosc2_context *context)
        int id
        char *name

    ctypedef struct blosc2_io:
        uint8_t id
        const char *name
        void* params

    ctypedef struct blosc2_stdio_mmap:
        const char* mode
        int64_t initial_mapping_size
        c_bool needs_free

    cdef const blosc2_stdio_mmap BLOSC2_STDIO_MMAP_DEFAULTS

    ctypedef struct blosc2_stdio_params:
        c_bool locking

    ctypedef struct blosc2_schunk:
        uint8_t version
        uint8_t compcode
        uint8_t compcode_meta
        uint8_t clevel
        uint8_t splitmode
        int32_t typesize
        int32_t blocksize
        int32_t chunksize
        uint8_t filters[BLOSC2_MAX_FILTERS]
        uint8_t filters_meta[BLOSC2_MAX_FILTERS]
        int64_t nchunks
        int64_t current_nchunk
        int64_t nbytes
        int64_t cbytes
        uint8_t** data
        size_t data_len
        blosc2_storage* storage
        blosc2_frame* frame
        blosc2_context* cctx
        blosc2_context* dctx
        blosc2_metalayer *metalayers[BLOSC2_MAX_METALAYERS]
        uint16_t nmetalayers
        blosc2_metalayer *vlmetalayers[BLOSC2_MAX_VLMETALAYERS]
        int16_t nvlmetalayers
        int tuner_id
        void *tuner_params
        int8_t ndim
        int64_t *blockshape
        int64_t change_tick

    blosc2_schunk *blosc2_schunk_new(blosc2_storage *storage)
    blosc2_schunk *blosc2_schunk_copy(blosc2_schunk *schunk, blosc2_storage *storage)
    blosc2_schunk *blosc2_schunk_from_buffer(uint8_t *cframe, int64_t len, c_bool copy)
    blosc2_schunk *blosc2_schunk_open_offset(const char* urlpath, int64_t offset)
    blosc2_schunk* blosc2_schunk_open_offset_udio(const char* urlpath, int64_t offset, const blosc2_io *udio)

    int64_t blosc2_schunk_to_buffer(blosc2_schunk* schunk, uint8_t** cframe, c_bool* needs_free) nogil
    void blosc2_schunk_avoid_cframe_free(blosc2_schunk *schunk, c_bool avoid_cframe_free)
    int blosc2_schunk_lock(blosc2_schunk *schunk) nogil
    int blosc2_schunk_unlock(blosc2_schunk *schunk) nogil
    int blosc2_schunk_refresh(blosc2_schunk *schunk) nogil
    int64_t blosc2_schunk_to_file(blosc2_schunk* schunk, const char* urlpath)
    int64_t blosc2_schunk_free(blosc2_schunk *schunk) nogil
    int64_t blosc2_schunk_append_chunk(blosc2_schunk *schunk, uint8_t *chunk, c_bool copy)
    int64_t blosc2_schunk_update_chunk(blosc2_schunk *schunk, int64_t nchunk, uint8_t *chunk, c_bool copy)
    int64_t blosc2_schunk_insert_chunk(blosc2_schunk *schunk, int64_t nchunk, uint8_t *chunk, c_bool copy)
    int64_t blosc2_schunk_delete_chunk(blosc2_schunk *schunk, int64_t nchunk)
    int64_t blosc2_schunk_fill_special(blosc2_schunk *schunk, int64_t nitems, int special_value,
                                       int32_t chunksize);

    int64_t blosc2_schunk_append_buffer(blosc2_schunk *schunk, void *src, int32_t nbytes)
    int blosc2_schunk_decompress_chunk(blosc2_schunk *schunk, int64_t nchunk, void *dest, int32_t nbytes)

    int blosc2_schunk_get_chunk(blosc2_schunk *schunk, int64_t nchunk, uint8_t ** chunk,
                                c_bool *needs_free) nogil
    int blosc2_schunk_get_lazychunk(blosc2_schunk *schunk, int64_t nchunk, uint8_t ** chunk,
                                    c_bool *needs_free) nogil
    int blosc2_schunk_get_vlblock(blosc2_schunk *schunk, int64_t nchunk, int32_t nblock,
                                  uint8_t **dest, int32_t *destsize)
    int blosc2_schunk_get_slice_buffer(blosc2_schunk *schunk, int64_t start, int64_t stop, void *buffer)
    int blosc2_schunk_set_slice_buffer(blosc2_schunk *schunk, int64_t start, int64_t stop, void *buffer)
    int blosc2_schunk_get_cparams(blosc2_schunk *schunk, blosc2_cparams** cparams)
    int blosc2_schunk_get_dparams(blosc2_schunk *schunk, blosc2_dparams** dparams)
    int blosc2_schunk_reorder_offsets(blosc2_schunk *schunk, int64_t *offsets_order)
    int64_t blosc2_schunk_frame_len(blosc2_schunk* schunk)

    int blosc2_chunk_repeatval(blosc2_cparams cparams, const int32_t nbytes,
                               void *dest, int32_t destsize, const void *repeatval)

    int blosc2_meta_exists(blosc2_schunk *schunk, const char *name)
    int blosc2_meta_add(blosc2_schunk *schunk, const char *name, uint8_t *content,
                        int32_t content_len)
    int blosc2_meta_update(blosc2_schunk *schunk, const char *name, uint8_t *content,
                           int32_t content_len)
    int blosc2_meta_get(blosc2_schunk *schunk, const char *name, uint8_t **content,
                        int32_t *content_len)
    int blosc2_vlmeta_exists(blosc2_schunk *schunk, const char *name)
    int blosc2_vlmeta_add(blosc2_schunk *schunk, const char *name,
                          uint8_t *content, int32_t content_len, blosc2_cparams *cparams)
    int blosc2_vlmeta_update(blosc2_schunk *schunk, const char *name,
                             uint8_t *content, int32_t content_len, blosc2_cparams *cparams)
    int blosc2_vlmeta_get(blosc2_schunk *schunk, const char *name,
                          uint8_t **content, int32_t *content_len)
    int blosc2_vlmeta_delete(blosc2_schunk *schunk, const char *name)
    int blosc2_vlmeta_get_names(blosc2_schunk *schunk, char **names)
    int blosc2_vldecompress_block_ctx(blosc2_context* context, const void* src,
                                      int32_t srcsize, int32_t nblock, uint8_t** dest,
                                      int32_t* destsize)


    int blosc1_get_blocksize()
    void blosc1_set_blocksize(size_t blocksize)
    void blosc1_set_schunk(blosc2_schunk *schunk)

    int blosc2_remove_dir(const char *path)
    int blosc2_remove_urlpath(const char *path)

    ctypedef int(*blosc2_codec_encoder_cb)(const uint8_t *input, int32_t input_len, uint8_t *output, int32_t output_len,
                                          uint8_t meta, blosc2_cparams *cparams, const void *chunk)
    ctypedef int(*blosc2_codec_decoder_cb)(const uint8_t *input, int32_t input_len, uint8_t *output, int32_t output_len,
                                          uint8_t meta, blosc2_dparams *dparams, const void *chunk)

    ctypedef struct blosc2_codec:
        uint8_t compcode
        char* compname
        uint8_t complib
        uint8_t version
        blosc2_codec_encoder_cb encoder
        blosc2_codec_decoder_cb decoder

    int blosc2_register_codec(blosc2_codec *codec)

    ctypedef int(*blosc2_filter_forward_cb)(const uint8_t *, uint8_t *, int32_t, uint8_t, blosc2_cparams *, uint8_t)
    ctypedef int(*blosc2_filter_backward_cb)(const uint8_t *, uint8_t *, int32_t, uint8_t, blosc2_dparams *, uint8_t)

    ctypedef struct blosc2_filter:
        uint8_t id
        char* name
        blosc2_filter_forward_cb forward
        blosc2_filter_backward_cb backward

    int blosc2_register_filter(blosc2_filter *filter)

    int blosc2_get_slice_nchunks(blosc2_schunk * schunk, int64_t *start, int64_t *stop, int64_t ** chunks_idx)


cdef extern from "b2nd.h":
    ctypedef enum:
        B2ND_MAX_DIM
        B2ND_MAX_METALAYERS
        B2ND_DEFAULT_DTYPE_FORMAT

    cdef struct chunk_cache_s:
        uint8_t *data
        int64_t nchunk

    ctypedef struct b2nd_array_t:
        blosc2_schunk* sc
        int64_t shape[B2ND_MAX_DIM]
        int32_t chunkshape[B2ND_MAX_DIM]
        int64_t extshape[B2ND_MAX_DIM]
        int32_t blockshape[B2ND_MAX_DIM]
        int64_t extchunkshape[B2ND_MAX_DIM]
        int64_t nitems
        int32_t chunknitems
        int64_t extnitems
        int32_t blocknitems
        int64_t extchunknitems
        int8_t ndim
        chunk_cache_s chunk_cache
        int64_t item_array_strides[B2ND_MAX_DIM]
        int64_t item_chunk_strides[B2ND_MAX_DIM]
        int64_t item_extchunk_strides[B2ND_MAX_DIM]
        int64_t item_block_strides[B2ND_MAX_DIM]
        int64_t block_chunk_strides[B2ND_MAX_DIM]
        int64_t chunk_array_strides[B2ND_MAX_DIM]
        char *dtype
        int8_t dtype_format

    ctypedef struct b2nd_context_t:
        pass
    b2nd_context_t *b2nd_create_ctx(blosc2_storage *b2_storage, int8_t ndim, int64_t *shape,
                                    int32_t *chunkshape, int32_t *blockshape, char *dtype,
                                    int8_t dtype_format, blosc2_metalayer *metalayers, int32_t nmetalayers)
    int b2nd_free_ctx(b2nd_context_t *ctx)

    int b2nd_uninit(b2nd_context_t *ctx, b2nd_array_t ** array)

    int b2nd_nans(b2nd_context_t * ctx, b2nd_array_t ** array)

    int b2nd_empty(b2nd_context_t *ctx, b2nd_array_t **array)
    int b2nd_zeros(b2nd_context_t *ctx, b2nd_array_t **array)
    int b2nd_full(b2nd_context_t *ctx, b2nd_array_t ** array, void *fill_value)

    int b2nd_free(b2nd_array_t *array)
    int b2nd_get_slice_cbuffer(b2nd_array_t *array,
                               int64_t *start, int64_t *stop,
                               void *buffer, int64_t *buffershape, int64_t buffersize)
    int b2nd_set_slice_cbuffer(void *buffer, int64_t *buffershape, int64_t buffersize,
                               int64_t *start, int64_t *stop, b2nd_array_t *array)
    int b2nd_get_slice(b2nd_context_t *ctx, b2nd_array_t **array, b2nd_array_t *src, const int64_t *start,
                       const int64_t *stop)
    int b2nd_from_cbuffer(b2nd_context_t *ctx, b2nd_array_t **array, void *buffer, int64_t buffersize)
    int b2nd_to_cbuffer(b2nd_array_t *array, void *buffer, int64_t buffersize)
    int b2nd_get_sparse_cbuffer(b2nd_array_t *array, int64_t ncoords, const int64_t *coords,
                                void *buffer, int64_t buffersize)
    int b2nd_from_cframe(uint8_t *cframe, int64_t cframe_len, c_bool copy, b2nd_array_t ** array);
    int b2nd_to_cframe(const b2nd_array_t *array, uint8_t ** cframe, int64_t *cframe_len,
                       c_bool *needs_free);

    int b2nd_squeeze(b2nd_array_t *array, b2nd_array_t **view)
    int b2nd_squeeze_index(b2nd_array_t *array, b2nd_array_t **view, const c_bool *index)
    int b2nd_resize(b2nd_array_t *array, const int64_t *new_shape, const int64_t *start)
    int b2nd_refresh(b2nd_array_t *array)
    int b2nd_copy(b2nd_context_t *ctx, b2nd_array_t *src, b2nd_array_t **array)
    int b2nd_concatenate(b2nd_context_t *ctx, b2nd_array_t *src1, b2nd_array_t *src2,
                         int8_t axis, c_bool copy, b2nd_array_t **array)
    int b2nd_expand_dims(const b2nd_array_t *array, b2nd_array_t ** view, const c_bool *axis, const uint8_t final_dims)
    int b2nd_get_orthogonal_selection(const b2nd_array_t *array, int64_t ** selection,
                                      int64_t *selection_size, void *buffer,
                                      int64_t *buffershape, int64_t buffersize)
    int b2nd_set_orthogonal_selection(const b2nd_array_t *array, int64_t ** selection,
                                      int64_t *selection_size, void *buffer,
                                      int64_t *buffershape, int64_t buffersize)
    int b2nd_from_schunk(blosc2_schunk *schunk, b2nd_array_t **array)

    void blosc2_unidim_to_multidim(uint8_t ndim, int64_t *shape, int64_t i, int64_t *index) nogil
    int b2nd_copy_buffer2(int8_t ndim,
                          int32_t itemsize,
                          const void *src, const int64_t *src_pad_shape,
                          const int64_t *src_start, const int64_t *src_stop,
                          void *dst, const int64_t *dst_pad_shape,
                          const int64_t *dst_start)


# miniexpr C API declarations
cdef extern from "miniexpr.h":
    ctypedef enum me_dtype:
        ME_AUTO,
        ME_BOOL
        ME_INT8
        ME_INT16
        ME_INT32
        ME_INT64
        ME_UINT8
        ME_UINT16
        ME_UINT32
        ME_UINT64
        ME_FLOAT32
        ME_FLOAT64
        ME_COMPLEX64
        ME_COMPLEX128
        ME_STRING

    # typedef struct me_variable
    ctypedef struct me_variable:
        const char *name
        me_dtype dtype
        const void *address
        int type
        void *context
        size_t itemsize

    ctypedef struct me_expr:
        int type
        double value
        const double *bound
        const void *function
        void *output
        int nitems
        me_dtype dtype
        me_dtype input_dtype
        void *bytecode
        int ncode
        void *parameters[1]

    int me_compile_nd_jit(const char *expression, const me_variable *variables,
                          int var_count, me_dtype dtype, int ndims,
                          const int64_t *shape, const int32_t *chunkshape,
                          const int32_t *blockshape, int jit_mode,
                          int *error, me_expr **out)

    ctypedef enum me_compile_status:
        ME_COMPILE_SUCCESS
        ME_COMPILE_ERR_OOM
        ME_COMPILE_ERR_PARSE
        ME_COMPILE_ERR_INVALID_ARG
        ME_COMPILE_ERR_COMPLEX_UNSUPPORTED
        ME_COMPILE_ERR_REDUCTION_INVALID
        ME_COMPILE_ERR_VAR_MIXED
        ME_COMPILE_ERR_VAR_UNSPECIFIED
        ME_COMPILE_ERR_INVALID_ARG_TYPE
        ME_COMPILE_ERR_MIXED_TYPE_NESTED

    ctypedef enum me_simd_ulp_mode:
        ME_SIMD_ULP_DEFAULT
        ME_SIMD_ULP_1
        ME_SIMD_ULP_3_5

    ctypedef enum me_jit_mode:
        ME_JIT_DEFAULT
        ME_JIT_ON
        ME_JIT_OFF

    ctypedef struct me_eval_params:
        c_bool disable_simd
        me_simd_ulp_mode simd_ulp_mode
        me_jit_mode jit_mode

    int me_eval(const me_expr *expr, const void **vars_block,
                int n_vars, void *output_block, int chunk_nitems,
                const me_eval_params *params) nogil

    int me_eval_nd(const me_expr *expr, const void **vars_block,
                   int n_vars, void *output_block, int block_nitems,
                   int64_t nchunk, int64_t nblock, const me_eval_params *params) nogil

    int me_nd_valid_nitems(const me_expr *expr, int64_t nchunk, int64_t nblock, int64_t *valid_nitems) nogil

    void me_print(const me_expr *n) nogil
    void me_free(me_expr *n) nogil

    bint me_expr_has_jit_kernel(const me_expr *expr) nogil

    ctypedef int (*me_wasm_jit_instantiate_helper)(
        const unsigned char *wasm_bytes,
        int wasm_len,
        int bridge_lookup_fn_idx
    )
    ctypedef void (*me_wasm_jit_free_helper)(int fn_idx)
    void me_register_wasm_jit_helpers(me_wasm_jit_instantiate_helper instantiate_helper,
                                      me_wasm_jit_free_helper free_helper)


cdef extern from "miniexpr_numpy.h":
    me_dtype me_dtype_from_numpy(int numpy_type_num)

cdef extern from "pythread.h":
    ctypedef void* PyThread_type_lock
    PyThread_type_lock PyThread_allocate_lock() nogil
    int PyThread_acquire_lock(PyThread_type_lock lock, int waitflag) nogil
    void PyThread_release_lock(PyThread_type_lock lock) nogil
    void PyThread_free_lock(PyThread_type_lock lock) nogil


ctypedef struct user_filters_udata:
    char* py_func
    int input_cdtype
    int output_cdtype
    int32_t chunkshape

ctypedef struct filler_udata:
    char* py_func
    uintptr_t inputs_id
    int output_cdtype
    int32_t chunkshape

ctypedef struct udf_udata:
    char* py_func
    uintptr_t inputs_id
    int output_cdtype
    b2nd_array_t *array
    int64_t chunks_in_array[B2ND_MAX_DIM]
    int64_t blocks_in_chunk[B2ND_MAX_DIM]

ctypedef enum:
    ME_CACHE_EMPTY
    ME_CACHE_LOADING
    ME_CACHE_READY
    ME_CACHE_ERROR

ctypedef struct me_input_cache_s:
    uint8_t* data
    int64_t nchunk
    int state
    PyThread_type_lock state_lock
    PyThread_type_lock ready_lock

ctypedef struct me_udata:
    b2nd_array_t** inputs
    me_input_cache_s* input_chunk_caches
    int ninputs
    me_eval_params* eval_params
    b2nd_array_t* array
    void* aux_reduc_ptr
    int64_t chunks_in_array[B2ND_MAX_DIM]
    int64_t blocks_in_chunk[B2ND_MAX_DIM]
    me_expr* miniexpr_handle
    # Optional candidate-block bitmap (1-D arrays only).  When non-NULL, blocks
    # whose byte is 0 are skipped: the prefilter writes a zero (false) output
    # and returns without decompressing inputs or running miniexpr.  Indexed by
    # the global block number ``nchunk * blocks_in_chunk[0] + nblock``.  NULL
    # means "evaluate every block" (the default, fully inert).
    const uint8_t* candidate_blocks
    int64_t candidate_blocks_len

ctypedef struct mm_udata:
    b2nd_array_t** inputs
    b2nd_array_t* array
    int64_t chunks_strides[3][B2ND_MAX_DIM]
    int64_t blocks_strides[3][B2ND_MAX_DIM]
    int64_t el_strides[3][B2ND_MAX_DIM]

MAX_TYPESIZE = BLOSC2_MAXTYPESIZE
MAX_BUFFERSIZE = BLOSC2_MAX_BUFFERSIZE
MAX_BLOCKSIZE = BLOSC2_MAXBLOCKSIZE
MAX_OVERHEAD = BLOSC2_MAX_OVERHEAD
MAX_DIM = B2ND_MAX_DIM
VERSION_STRING = (<char*>BLOSC2_VERSION_STRING).decode("utf-8")
VERSION_DATE = (<char*>BLOSC2_VERSION_DATE).decode("utf-8")
MIN_HEADER_LENGTH = BLOSC_MIN_HEADER_LENGTH
EXTENDED_HEADER_LENGTH = BLOSC_EXTENDED_HEADER_LENGTH
DEFINED_CODECS_STOP = BLOSC2_DEFINED_CODECS_STOP
GLOBAL_REGISTERED_CODECS_STOP = BLOSC2_GLOBAL_REGISTERED_CODECS_STOP
USER_REGISTERED_CODECS_STOP = BLOSC2_USER_REGISTERED_CODECS_STOP
DEFAULT_DTYPE_FORMAT = B2ND_DEFAULT_DTYPE_FORMAT

cdef _check_comp_length(comp_name, comp_len):
    if comp_len < BLOSC_MIN_HEADER_LENGTH:
        raise ValueError(f"{comp_name} cannot be less than {BLOSC_MIN_HEADER_LENGTH} bytes")


blosc2_init()

@atexit.register
def destroy():
    blosc2_destroy()


def register_wasm_jit_helpers(uintptr_t instantiate_ptr, uintptr_t free_ptr):
    cdef me_wasm_jit_instantiate_helper instantiate_helper = (
        <me_wasm_jit_instantiate_helper>instantiate_ptr
    )
    cdef me_wasm_jit_free_helper free_helper = <me_wasm_jit_free_helper>free_ptr
    me_register_wasm_jit_helpers(instantiate_helper, free_helper)


cdef inline me_dtype _me_dtype_from_numpy_dtype(dtype_obj):
    dtype = np.dtype(dtype_obj)
    cdef int itemsize = <int>dtype.itemsize
    kind = dtype.kind
    if kind == "b":
        return ME_BOOL
    if kind == "i":
        if itemsize == 1:
            return ME_INT8
        if itemsize == 2:
            return ME_INT16
        if itemsize == 4:
            return ME_INT32
        if itemsize == 8:
            return ME_INT64
    elif kind == "u":
        if itemsize == 1:
            return ME_UINT8
        if itemsize == 2:
            return ME_UINT16
        if itemsize == 4:
            return ME_UINT32
        if itemsize == 8:
            return ME_UINT64
    elif kind == "f":
        if itemsize == 4:
            return ME_FLOAT32
        if itemsize == 8:
            return ME_FLOAT64
    elif kind == "c":
        if itemsize == 8:
            return ME_COMPLEX64
        if itemsize == 16:
            return ME_COMPLEX128
    elif kind == "U":
        # miniexpr string variables use fixed-size UCS4 (numpy unicode) storage.
        if itemsize <= 0 or itemsize % 4 != 0:
            raise TypeError(
                f"miniexpr string operands require unicode dtype with UCS4 itemsize; got '{dtype}'"
            )
        return ME_STRING
    return <me_dtype>-1


cdef inline str _me_compile_status_name(int rc):
    if rc == ME_COMPILE_SUCCESS:
        return "ME_COMPILE_SUCCESS"
    if rc == ME_COMPILE_ERR_OOM:
        return "ME_COMPILE_ERR_OOM"
    if rc == ME_COMPILE_ERR_PARSE:
        return "ME_COMPILE_ERR_PARSE"
    if rc == ME_COMPILE_ERR_INVALID_ARG:
        return "ME_COMPILE_ERR_INVALID_ARG"
    if rc == ME_COMPILE_ERR_COMPLEX_UNSUPPORTED:
        return "ME_COMPILE_ERR_COMPLEX_UNSUPPORTED"
    if rc == ME_COMPILE_ERR_REDUCTION_INVALID:
        return "ME_COMPILE_ERR_REDUCTION_INVALID"
    if rc == ME_COMPILE_ERR_VAR_MIXED:
        return "ME_COMPILE_ERR_VAR_MIXED"
    if rc == ME_COMPILE_ERR_VAR_UNSPECIFIED:
        return "ME_COMPILE_ERR_VAR_UNSPECIFIED"
    if rc == ME_COMPILE_ERR_INVALID_ARG_TYPE:
        return "ME_COMPILE_ERR_INVALID_ARG_TYPE"
    if rc == ME_COMPILE_ERR_MIXED_TYPE_NESTED:
        return "ME_COMPILE_ERR_MIXED_TYPE_NESTED"
    return "ME_COMPILE_ERR_UNKNOWN"


cdef inline str _me_compile_error_details(int rc, int error):
    cdef str details = f"{_me_compile_status_name(rc)} ({rc})"
    if rc == ME_COMPILE_ERR_PARSE and error > 0:
        details += f", parse_error_pos={error}"
    elif error != 0:
        details += f", error_pos={error}"
    return details


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def nearest_divisor(int64_t a, int64_t b, bint strict=False):
    """Find the divisor of `a` that is closest to `b`.

    Parameters
    ----------
    a : int
        The number for which to find divisors.
    b : int
        The reference value to compare divisors against.
    strict : bool, optional
        If True, always use the downward search algorithm.

    Returns
    -------
    int
        The divisor of `a` that is closest to `b`.

    Notes
    -----
    This is a *much* faster version than its Python counterpart.
    """
    cdef:
        int64_t i, closest, min_diff, diff
        bint found

    if a > 100_000 or strict:
        # For large numbers or when strict=True, search downwards from b
        i = b
        while i > 0:
            if a % i == 0:
                return i
            i -= 1
        return 1  # Fallback to 1, which is always a divisor

    # For smaller numbers, find the closest divisor
    closest = 1
    min_diff = a  # Initialize to a large value
    found = False

    # Search for divisors up to sqrt(a)
    i = 1
    while i * i <= a:
        if a % i == 0:
            # Check i as a divisor
            diff = c_abs(i - b)
            if diff < min_diff:
                min_diff = diff
                closest = i
                found = True

            # Check a/i as a divisor
            diff = c_abs(a // i - b)
            if diff < min_diff:
                min_diff = diff
                closest = a // i
                found = True
        i += 1

    return closest if found else 1


def cbuffer_sizes(src):
    cdef const uint8_t[:] typed_view_src
    mem_view_src = memoryview(src)
    typed_view_src = mem_view_src.cast('B')
    _check_comp_length('src', typed_view_src.nbytes)
    cdef int32_t nbytes
    cdef int32_t cbytes
    cdef int32_t blocksize
    blosc2_cbuffer_sizes(<void*>&typed_view_src[0], &nbytes, &cbytes, &blocksize)
    return nbytes, cbytes, blocksize


cpdef compress(src, int32_t typesize=8, int clevel=9, filter=blosc2.Filter.SHUFFLE, codec=blosc2.Codec.BLOSCLZ):
    set_compressor(codec)
    cdef int32_t len_src = <int32_t> len(src)
    cdef Py_buffer buf
    PyObject_GetBuffer(src, &buf, PyBUF_SIMPLE)
    dest = bytes(buf.len + BLOSC2_MAX_OVERHEAD)
    cdef int32_t len_dest =  <int32_t> len(dest)
    cdef int size
    cdef int filter_ = filter.value if isinstance(filter, Enum) else 0
    if RELEASEGIL:
        _dest = <void*> <char *> dest
        with nogil:
            size = blosc2_compress(clevel, filter_, <int32_t> typesize, buf.buf, <int32_t> buf.len, _dest, len_dest)
    else:
        size = blosc2_compress(clevel, filter_, <int32_t> typesize, buf.buf, <int32_t> buf.len, <void*> <char *> dest, len_dest)
    PyBuffer_Release(&buf)
    if size > 0:
        return dest[:size]
    else:
        raise ValueError("Cannot compress")


def decompress(src, dst=None, as_bytearray=False):
    cdef int32_t nbytes
    cdef int32_t cbytes
    cdef int32_t blocksize
    cdef const uint8_t[:] typed_view_src

    mem_view_src = memoryview(src)
    typed_view_src = mem_view_src.cast('B')
    _check_comp_length('src', len(typed_view_src))
    blosc2_cbuffer_sizes(<void*>&typed_view_src[0], &nbytes, &cbytes, &blocksize)
    cdef Py_buffer buf
    if dst is not None:
        PyObject_GetBuffer(dst, &buf, PyBUF_SIMPLE)
        if buf.len == 0:
            raise ValueError("The dst length must be greater than 0")
        size = blosc1_decompress(<void*>&typed_view_src[0], buf.buf, buf.len)
        PyBuffer_Release(&buf)
    else:
        dst = PyBytes_FromStringAndSize(NULL, nbytes)
        if dst is None:
            raise RuntimeError("Could not get a bytes object")
        size = blosc1_decompress(<void*>&typed_view_src[0], <void*> <char *> dst, len(dst))
        if as_bytearray:
            dst = bytearray(dst)
        if size >= 0:
            return dst
    if size < 0:
        raise RuntimeError("Cannot decompress")


def set_compressor(codec):
    codec = codec.name.lower().encode("utf-8")
    size = blosc1_set_compressor(codec)
    if size == -1:
        raise ValueError("The code is not available")
    else:
        return size

def free_resources():
    rc = blosc2_free_resources()
    if rc < 0:
        raise ValueError("Could not free the resources")

def set_nthreads(nthreads):
    if nthreads > INT_MAX:
        raise ValueError("nthreads must be less or equal than 2^31 - 1.")
    rc = blosc2_set_nthreads(nthreads)
    if rc < 0:
        raise ValueError("nthreads must be a positive integer.")
    else:
        return rc

def set_blocksize(size_t blocksize=0):
    blosc1_set_blocksize(blocksize)

def clib_info(codec):
    cdef char* clib
    cdef char* version
    codec = codec.name.lower().encode("utf-8")
    rc = blosc2_get_complib_info(codec, &clib, &version)
    if rc >= 0:
        return clib, version
    else:
        raise ValueError("The compression library is not supported.")

def get_clib(bytesobj):
    rc = blosc2_cbuffer_complib(<void *> <char*> bytesobj)
    if rc == NULL:
        raise ValueError("Cannot get the info for the compressor")
    else:
        return rc

def get_compressor():
    return blosc1_get_compressor()


cdef c_bool RELEASEGIL = False

def set_releasegil(c_bool gilstate):
    global RELEASEGIL
    oldstate = RELEASEGIL
    RELEASEGIL = gilstate
    return oldstate

def get_blocksize():
    return blosc1_get_blocksize()

cdef _check_cparams(blosc2_cparams *cparams):
    if cparams.nthreads > 1:
        if BLOSC2_USER_REGISTERED_CODECS_START <= cparams.compcode <= BLOSC2_USER_REGISTERED_CODECS_STOP\
                and cparams.compcode in blosc2.ucodecs_registry.keys():
            raise ValueError("Cannot use multi-threading with user defined Python codecs")

        ufilters = [BLOSC2_USER_REGISTERED_FILTERS_START <= filter <= BLOSC2_USER_REGISTERED_FILTERS_STOP
                  for filter in cparams.filters]
        for i in range(len(ufilters)):
            if ufilters[i] and cparams.filters[i] in blosc2.ufilters_registry.keys():
                raise ValueError("Cannot use multi-threading with user defined Python filters")

        if cparams.prefilter != NULL and (cparams.prefilter != <blosc2_prefilter_fn>miniexpr_prefilter and cparams.prefilter != <blosc2_prefilter_fn>matmul_prefilter):
            # Note: miniexpr_prefilter uses miniexpr C API which is thread-friendly,
            raise ValueError("`nthreads` must be 1 when a prefilter is set")

cdef _check_dparams(blosc2_dparams* dparams, blosc2_cparams* cparams=NULL):
    if cparams == NULL:
        return
    if dparams.nthreads > 1:
        if BLOSC2_USER_REGISTERED_CODECS_START <= cparams.compcode <= BLOSC2_USER_REGISTERED_CODECS_STOP\
                and cparams.compcode in blosc2.ucodecs_registry.keys():
            raise ValueError("Cannot use multi-threading with user defined Python codecs")

        ufilters = [BLOSC2_USER_REGISTERED_FILTERS_START <= filter <= BLOSC2_USER_REGISTERED_FILTERS_STOP
                    for filter in cparams.filters]
        for i in range(len(ufilters)):
            if ufilters[i] and cparams.filters[i] in blosc2.ufilters_registry.keys():
                raise ValueError("Cannot use multi-threading with user defined Python filters")

        if dparams.postfilter != NULL:
            raise ValueError("`nthreads` must be 1 when a postfilter is set")


cdef create_cparams_from_kwargs(blosc2_cparams *cparams, kwargs):
    if "compcode" in kwargs:
        raise NameError("`compcode` has been renamed to `codec`. Please go update your code.")
    if "shuffle" in kwargs:
        raise NameError("`shuffle` has been substituted by `filters`. Please go update your code.")
    codec = kwargs.get('codec', blosc2.cparams_dflts['codec'])
    cparams.compcode = codec if not isinstance(codec, blosc2.Codec) else codec.value
    cparams.compcode_meta = kwargs.get('codec_meta', blosc2.cparams_dflts['codec_meta'])
    cparams.clevel = kwargs.get('clevel', blosc2.cparams_dflts['clevel'])
    cparams.use_dict = kwargs.get('use_dict', blosc2.cparams_dflts['use_dict'])
    cparams.typesize = typesize = kwargs.get('typesize', blosc2.cparams_dflts['typesize'])
    cparams.nthreads = kwargs.get('nthreads', 1 if blosc2.IS_WASM else blosc2.nthreads)
    if blosc2.IS_WASM:
        cparams.nthreads = 1
    cparams.blocksize = kwargs.get('blocksize', blosc2.cparams_dflts['blocksize'])
    splitmode = kwargs.get('splitmode', blosc2.cparams_dflts['splitmode'])
    cparams.splitmode = splitmode.value
    # TODO: support the commented ones in the future
    #schunk_c = kwargs.get('schunk', blosc2.cparams_dflts['schunk'])
    #cparams.schunk = <void *> schunk_c
    cparams.schunk = NULL
    for i in range(BLOSC2_MAX_FILTERS):
        cparams.filters[i] = 0
        cparams.filters_meta[i] = 0

    filters = kwargs.get('filters', blosc2.cparams_dflts['filters'])
    if len(filters) > BLOSC2_MAX_FILTERS:
        raise ValueError(f"filters list cannot exceed {BLOSC2_MAX_FILTERS}")
    for i, filter in enumerate(filters):
        cparams.filters[i] = filter.value if isinstance(filter, Enum) else filter
        # Bytedelta does not work on typesize 1
        if cparams.filters[i] == blosc2.Filter.BYTEDELTA.value and typesize == 1:
            cparams.filters[i] = 0

    if "filters_meta" not in kwargs:
        # If not specified, we can still assign a 0 list to it
        filters_meta = [0] * len(filters)
    else:
        filters_meta = kwargs['filters_meta']
        if len(filters) != len(filters_meta):
            raise ValueError("filters and filters_meta lists must have same length")
    cdef int8_t meta_value
    for i, meta in enumerate(filters_meta):
        # We still may want to encode negative values
        meta_value = <int8_t>meta if meta < 0 else meta
        if meta_value == 0 and cparams.filters[i] == blosc2.Filter.BYTEDELTA.value:
            # bytedelta typesize cannot be zero when using compress2
            cparams.filters_meta[i] = <uint8_t>typesize
        else:
            cparams.filters_meta[i] = <uint8_t>meta_value

    cparams.prefilter = NULL
    cparams.preparams = NULL
    tuner = kwargs.get('tuner', blosc2.cparams_dflts['tuner'])
    cparams.tuner_id = tuner.value
    cparams.tuner_params = NULL
    cparams.instr_codec = False
    cparams.codec_params = NULL
    for i in range(len(filters)):
        cparams.filter_params[i] = NULL

    _check_cparams(cparams)


def compress2(src, **kwargs):
    cdef blosc2_cparams cparams
    create_cparams_from_kwargs(&cparams, kwargs)

    cdef blosc2_context *cctx
    cdef Py_buffer buf
    PyObject_GetBuffer(src, &buf, PyBUF_SIMPLE)
    cdef int size
    cdef int32_t len_dest = <int32_t> (buf.len + BLOSC2_MAX_OVERHEAD)
    dest = bytes(len_dest)
    _dest = <void*> <char *> dest
    cctx = blosc2_create_cctx(cparams)
    if cctx == NULL:
        raise RuntimeError("Could not create the compression context")
    if RELEASEGIL:
        with nogil:
            size = blosc2_compress_ctx(cctx, buf.buf, <int32_t> buf.len, _dest, len_dest)
    else:
        size = blosc2_compress_ctx(cctx, buf.buf, <int32_t> buf.len, _dest, len_dest)
    blosc2_free_ctx(cctx)
    PyBuffer_Release(&buf)
    if size < 0:
        raise RuntimeError("Could not compress the data")
    elif size == 0:
        del dest
        raise RuntimeError("The result could not fit ")
    return dest[:size]

cdef create_dparams_from_kwargs(blosc2_dparams *dparams, kwargs, blosc2_cparams* cparams=NULL):
    memcpy(dparams, &BLOSC2_DPARAMS_DEFAULTS, sizeof(BLOSC2_DPARAMS_DEFAULTS))
    dparams.nthreads = kwargs.get('nthreads', 1 if blosc2.IS_WASM else blosc2.nthreads)
    if blosc2.IS_WASM:
        dparams.nthreads = 1
    dparams.schunk = NULL
    dparams.postfilter = NULL
    dparams.postparams = NULL
    # TODO: support the next ones in the future
    #dparams.schunk = kwargs.get('schunk', blosc2.dparams_dflts['schunk'])
    #dparams.typesize = typesize = kwargs.get('typesize', blosc2.dparams_dflts['typesize'])
    _check_dparams(dparams, cparams)

def decompress2(src, dst=None, **kwargs):
    cdef blosc2_dparams dparams
    cdef char *dst_buf
    cdef void *view
    create_dparams_from_kwargs(&dparams, kwargs)

    cdef blosc2_context *dctx = blosc2_create_dctx(dparams)
    if dctx == NULL:
        raise RuntimeError("Could not create decompression context")
    cdef const uint8_t[:] typed_view_src
    mem_view_src = memoryview(src)
    typed_view_src = mem_view_src.cast('B')
    _check_comp_length('src', typed_view_src.nbytes)
    cdef int32_t nbytes
    cdef int32_t cbytes
    cdef int32_t blocksize
    cdef int32_t srcsize = <int32_t>typed_view_src.nbytes
    blosc2_cbuffer_sizes(<void*>&typed_view_src[0], &nbytes, &cbytes, &blocksize)
    cdef Py_buffer buf
    if dst is not None:
        PyObject_GetBuffer(dst, &buf, PyBUF_SIMPLE)
        if buf.len == 0:
            blosc2_free_ctx(dctx)
            raise ValueError("The dst length must be greater than 0")
        view = <void*>&typed_view_src[0]
        # For lazy chunks, blosc2_cbuffer_sizes() only reports the header cbytes.
        # The decode context needs the full source buffer length.
        if RELEASEGIL:
            with nogil:
                size = blosc2_decompress_ctx(dctx, view, srcsize, buf.buf, nbytes)
        else:
            size = blosc2_decompress_ctx(dctx, view, srcsize, buf.buf, nbytes)
        blosc2_free_ctx(dctx)
        PyBuffer_Release(&buf)
    else:
        dst = PyBytes_FromStringAndSize(NULL, nbytes)
        if dst is None:
            blosc2_free_ctx(dctx)
            raise RuntimeError("Could not get a bytes object")
        dst_buf = <char*>dst
        view = <void*>&typed_view_src[0]
        # For lazy chunks, blosc2_cbuffer_sizes() only reports the header cbytes.
        # The decode context needs the full source buffer length.
        if RELEASEGIL:
            with nogil:
                size = blosc2_decompress_ctx(dctx, view, srcsize, <void*>dst_buf, nbytes)
        else:
            size = blosc2_decompress_ctx(dctx, view, srcsize, <void*>dst_buf, nbytes)
        blosc2_free_ctx(dctx)
        if size >= 0:
            return dst
    if size < 0:
        raise ValueError("Error while decompressing, check the src data and/or the dparams")


def vlcompress(srcs, **kwargs):
    cdef blosc2_cparams cparams
    create_cparams_from_kwargs(&cparams, kwargs)

    cdef Py_ssize_t nblocks = len(srcs)
    if nblocks <= 0:
        raise ValueError("At least one block is required")

    cdef blosc2_context *cctx = NULL
    cdef Py_buffer *buffers = <Py_buffer*>calloc(nblocks, sizeof(Py_buffer))
    cdef const void **src_ptrs = <const void **>malloc(nblocks * sizeof(void *))
    cdef int32_t *srcsizes = <int32_t*>malloc(nblocks * sizeof(int32_t))
    cdef Py_ssize_t acquired = 0
    cdef Py_ssize_t i
    cdef int64_t total_nbytes = 0
    cdef int32_t len_dest
    cdef int size
    cdef Py_ssize_t release_i
    cdef void *_dest
    if buffers == NULL or src_ptrs == NULL or srcsizes == NULL:
        free(buffers)
        free(src_ptrs)
        free(srcsizes)
        raise MemoryError()

    try:
        for i in range(nblocks):
            PyObject_GetBuffer(srcs[i], &buffers[i], PyBUF_SIMPLE)
            acquired += 1
            if buffers[i].len <= 0:
                raise ValueError("Each VL block must have at least one byte")
            src_ptrs[i] = buffers[i].buf
            srcsizes[i] = <int32_t>buffers[i].len
            total_nbytes += buffers[i].len

        # VL blocks can carry enough per-block framing that the simple
        # total_nbytes + global_overhead estimate is too small for many tiny
        # buffers.  Budget one max-overhead chunk per block as a conservative
        # upper bound for the temporary destination.
        len_dest = <int32_t>(total_nbytes + BLOSC2_MAX_OVERHEAD * (nblocks + 1) + 64)
        dest = PyBytes_FromStringAndSize(NULL, len_dest)
        if dest is None:
            raise MemoryError()
        _dest = <void*><char *>dest
        cctx = blosc2_create_cctx(cparams)
        if cctx == NULL:
            raise RuntimeError("Could not create the compression context")
        if RELEASEGIL:
            with nogil:
                size = blosc2_vlcompress_ctx(cctx, src_ptrs, srcsizes, <int32_t>nblocks, _dest, len_dest)
        else:
            size = blosc2_vlcompress_ctx(cctx, src_ptrs, srcsizes, <int32_t>nblocks, _dest, len_dest)
    finally:
        if cctx != NULL:
            blosc2_free_ctx(cctx)
        for release_i in range(acquired):
            PyBuffer_Release(&buffers[release_i])
        free(buffers)
        free(src_ptrs)
        free(srcsizes)

    if size < 0:
        raise RuntimeError("Could not compress the data")
    elif size == 0:
        del dest
        raise RuntimeError("The result could not fit ")
    return dest[:size]


def vldecompress(src, **kwargs):
    cdef blosc2_dparams dparams
    create_dparams_from_kwargs(&dparams, kwargs)

    cdef blosc2_context *dctx = blosc2_create_dctx(dparams)
    if dctx == NULL:
        raise RuntimeError("Could not create decompression context")

    cdef const uint8_t[:] typed_view_src
    mem_view_src = memoryview(src)
    typed_view_src = mem_view_src.cast('B')
    _check_comp_length('src', typed_view_src.nbytes)
    cdef int32_t nbytes
    cdef int32_t cbytes
    cdef int32_t nblocks
    cdef int32_t srcsize = <int32_t>typed_view_src.nbytes
    blosc2_cbuffer_sizes(<void*>&typed_view_src[0], &nbytes, &cbytes, &nblocks)
    if nblocks <= 0:
        blosc2_free_ctx(dctx)
        raise ValueError("Chunk does not contain VL blocks")

    cdef void **dests = <void**>calloc(nblocks, sizeof(void *))
    cdef int32_t *destsizes = <int32_t*>malloc(nblocks * sizeof(int32_t))
    cdef int32_t rc
    cdef int32_t i
    cdef list out = []
    if dests == NULL or destsizes == NULL:
        blosc2_free_ctx(dctx)
        free(dests)
        free(destsizes)
        raise MemoryError()

    try:
        # For lazy chunks, blosc2_cbuffer_sizes() only reports the header cbytes.
        # The decode context needs the full source buffer length.
        rc = blosc2_vldecompress_ctx(dctx, <void*>&typed_view_src[0], srcsize, dests, destsizes, nblocks)
        if rc < 0:
            raise RuntimeError("Could not decompress the data")
        for i in range(rc):
            out.append(PyBytes_FromStringAndSize(<char*>dests[i], destsizes[i]))
            free(dests[i])
            dests[i] = NULL
        return out
    finally:
        for i in range(nblocks):
            if dests[i] != NULL:
                free(dests[i])
        free(dests)
        free(destsizes)
        blosc2_free_ctx(dctx)


def vldecompress_block(src, int32_t nblock, **kwargs):
    cdef blosc2_dparams dparams
    create_dparams_from_kwargs(&dparams, kwargs)

    cdef blosc2_context *dctx = blosc2_create_dctx(dparams)
    if dctx == NULL:
        raise RuntimeError("Could not create decompression context")

    cdef const uint8_t[:] typed_view_src
    mem_view_src = memoryview(src)
    typed_view_src = mem_view_src.cast('B')
    _check_comp_length('src', typed_view_src.nbytes)

    cdef uint8_t *dest = NULL
    cdef int32_t destsize = 0
    cdef int32_t rc
    try:
        rc = blosc2_vldecompress_block_ctx(
            dctx,
            <void*>&typed_view_src[0],
            <int32_t>typed_view_src.nbytes,
            nblock,
            &dest,
            &destsize,
        )
        if rc < 0:
            raise RuntimeError("Could not decompress the block")
        return PyBytes_FromStringAndSize(<char*>dest, destsize)
    finally:
        if dest != NULL:
            free(dest)
        blosc2_free_ctx(dctx)


# Immortal io object for the opt-in frame locking (default filesystem backend
# with blosc2_stdio_params.locking set).  The C side only reads the flag, so a
# single, module-lifetime instance can serve every locked schunk (and trivially
# satisfies the "params must outlive the schunk" contract).
cdef blosc2_stdio_params _locking_params
_locking_params.locking = True
cdef blosc2_io _locking_io
_locking_io.id = BLOSC2_IO_FILESYSTEM
_locking_io.name = "filesystem"
_locking_io.params = &_locking_params


cdef create_storage(blosc2_storage *storage, kwargs):
    contiguous = kwargs.get('contiguous', blosc2.storage_dflts['contiguous'])
    storage.contiguous = contiguous
    urlpath = kwargs.get('urlpath', blosc2.storage_dflts['urlpath'])
    if urlpath is None:
        storage.urlpath = NULL
    else:
        storage.urlpath = urlpath

    create_cparams_from_kwargs(storage.cparams, kwargs.get('cparams', {}))
    create_dparams_from_kwargs(storage.dparams, kwargs.get('dparams', {}), storage.cparams)

    cdef blosc2_io* io
    cdef blosc2_stdio_mmap* mmap_file
    mmap_mode = kwargs.get("mmap_mode")
    initial_mapping_size = kwargs.get("initial_mapping_size")
    locking = kwargs.get("locking")
    if mmap_mode is not None:
        if urlpath is None:
            raise ValueError("urlpath must be set when using mmap_mode")
        if not contiguous:
            raise ValueError("Only contiguous storage is supported for memory-mapped files")
        if locking:
            raise ValueError("locking is not supported together with mmap_mode")

        # sizeof(BLOSC2_STDIO_MMAP_DEFAULTS) yields the size of the full struct as defined in the C header
        mmap_file = <blosc2_stdio_mmap *>malloc(sizeof(BLOSC2_STDIO_MMAP_DEFAULTS))
        memcpy(mmap_file, &BLOSC2_STDIO_MMAP_DEFAULTS, sizeof(BLOSC2_STDIO_MMAP_DEFAULTS))

        # The storage for the bytes for the mmap_mode parameter need to be available even after this function
        kwargs["_mmap_mode_bytes"] = kwargs["mmap_mode"].encode("utf-8")
        mmap_file.mode = kwargs["_mmap_mode_bytes"]
        mmap_file.needs_free = True
        if initial_mapping_size is not None:
            mmap_file.initial_mapping_size = initial_mapping_size

        io = <blosc2_io *>malloc(sizeof(blosc2_io))
        io.id = BLOSC2_IO_FILESYSTEM_MMAP
        io.params = mmap_file
        storage.io = io
    elif locking:
        if urlpath is None:
            raise ValueError("urlpath must be set when using locking")
        storage.io = &_locking_io
    else:
        storage.io = NULL


cdef get_chunk_repeatval(blosc2_cparams cparams, const int32_t nbytes,
                        void *dest, int32_t destsize, Py_buffer *repeatval):
    if blosc2_chunk_repeatval(cparams, nbytes, dest, destsize, repeatval.buf) < 0:
        free(dest)
        PyBuffer_Release(repeatval)
        raise RuntimeError("Problems when creating the repeated values chunk")


cdef class SChunk:
    cdef blosc2_schunk *schunk
    cdef c_bool _is_view

    def __init__(self, _schunk=None, chunksize=2 ** 24, data=None, **kwargs):
        # hold on to a bytestring of urlpath for the lifetime of the instance
        # because its value is referenced via a C-pointer
        urlpath = kwargs.get("urlpath", None)
        if urlpath is not None:
            if isinstance(urlpath, pathlib.PurePath):
                urlpath = str(urlpath)
            self._urlpath = urlpath.encode() if isinstance(urlpath, str) else urlpath
            kwargs["urlpath"] = self._urlpath

        self.mode = blosc2.Storage().mode if kwargs.get("mode", None) is None else kwargs.get("mode")
        self.mmap_mode = kwargs.get("mmap_mode")
        self.initial_mapping_size = kwargs.get("initial_mapping_size")
        self.locking = bool(kwargs.get("locking", False))
        if self.locking and self.mmap_mode is not None:
            raise ValueError("locking is not supported together with mmap_mode")
        if self.mmap_mode is not None:
            self.mode = mode_from_mmap_mode(self.mmap_mode)
        if self.initial_mapping_size is not None:
            if self.mmap_mode is None:
                raise ValueError("initial_mapping_size can only be used with mmap_mode")

            if self.mmap_mode == "r":
                raise ValueError("initial_mapping_size can only be used with writing modes (r+, w+, c)")

        # `_is_view` indicates if a free should be done on this instance
        self._is_view = kwargs.get("_is_view", False)

        if _schunk is not None:
            self.schunk = <blosc2_schunk *> PyCapsule_GetPointer(_schunk, <char *> "blosc2_schunk*")
            if self.mode == "w" and urlpath is not None:
                blosc2.remove_urlpath(urlpath)
                self.schunk = blosc2_schunk_new(self.schunk.storage)
            return

        if kwargs is not None:
            if self.mode == "w":
                blosc2.remove_urlpath(urlpath)
            elif self.mode == "r":
                if urlpath is None:
                    raise ValueError("Cannot open the SChunk in reading mode (mode or mmap_mode is 'r') because you "
                                     "did not specify a urlpath pointing to an existing file on-disk")
                if not os.path.exists(urlpath):
                    raise ValueError("Cannot open the SChunk in reading mode (mode or mmap_mode is 'r') because the "
                                     f"file {urlpath} does not exist. Please use a writing mode if you want to create "
                                     "a new SChunk")

        cdef blosc2_storage storage
        # Create space for cparams and dparams in the stack
        cdef blosc2_cparams cparams
        cdef blosc2_dparams dparams
        storage.cparams = &cparams
        storage.dparams = &dparams
        if kwargs is None:
            storage = BLOSC2_STORAGE_DEFAULTS
        else:
            create_storage(&storage, kwargs)

        if self.mode == "r":
            offset = 0
            if storage.io != NULL:
                # mmap or locking: open through the user-defined io
                self.schunk = blosc2_schunk_open_offset_udio(storage.urlpath, offset, storage.io)
            else:
                self.schunk = blosc2_schunk_open_offset(storage.urlpath, offset)

            if kwargs is not None:
                check_schunk_params(self.schunk, kwargs)
            if schunk_is_ndarray(self.schunk):
                raise ValueError("Cannot open an NDArray as a SChunk. Please use blosc2.open instead")
        else:
            self.schunk = blosc2_schunk_new(&storage)

        if self.schunk == NULL:
            if self.mmap_mode is not None:
                free(storage.io)
            raise RuntimeError("Could not create the Schunk")

        # Add metalayers
        meta = kwargs.get("meta")
        if meta is not None:
            for (name, content) in meta.items():
                name = name.encode("utf-8") if isinstance(name, str) else name
                content = packb(content, default=encode_tuple, strict_types=True, use_bin_type=True)
                _check_rc(blosc2_meta_add(self.schunk, name, content, len(content)),
                          "Error while adding the metalayers")

        if chunksize > INT_MAX:
            raise ValueError("Maximum chunksize allowed is 2^31 - 1")
        self.schunk.chunksize = chunksize
        cdef const uint8_t[:] typed_view
        cdef int64_t index
        cdef Py_buffer buf
        cdef uint8_t *buf_ptr
        cdef int comp_size
        cdef int32_t csize
        cdef uint8_t* chunk
        cdef int32_t len_chunk
        if data is not None and len(data) > 0:
            PyObject_GetBuffer(data, &buf, PyBUF_SIMPLE)
            buf_ptr = <uint8_t *> buf.buf
            len_data = buf.len
            nchunks = len_data // chunksize + 1 if len_data % chunksize != 0 else len_data // chunksize
            len_chunk = chunksize
            for i in range(nchunks):
                if i == (nchunks - 1):
                    len_chunk = len_data - i * chunksize
                index = i * chunksize
                csize = <int32_t> (len_chunk + BLOSC2_MAX_OVERHEAD)
                chunk = <uint8_t*> malloc(csize)
                self.schunk.current_nchunk = i
                if RELEASEGIL:
                    with nogil:
                        comp_size = blosc2_compress_ctx(self.schunk.cctx, buf_ptr + index, len_chunk, chunk, csize)
                else:
                    comp_size = blosc2_compress_ctx(self.schunk.cctx, buf_ptr + index, len_chunk, chunk, csize)
                if comp_size < 0:
                    free(chunk)
                    PyBuffer_Release(&buf)
                    raise RuntimeError("Could not compress the data")
                elif comp_size == 0:
                    free(chunk)
                    PyBuffer_Release(&buf)
                    raise RuntimeError("The result could not fit")
                chunk = <uint8_t*> realloc(chunk, comp_size)
                _check_comp_length('chunk', comp_size)
                nchunks_ = blosc2_schunk_append_chunk(self.schunk, chunk, False)
                if nchunks_ != (i + 1):
                    free(chunk)
                    PyBuffer_Release(&buf)
                    raise RuntimeError("An error occurred while appending the chunks")
            PyBuffer_Release(&buf)

    @property
    def c_schunk(self):
        return <uintptr_t> self.schunk

    @property
    def chunksize(self):
        return self.schunk.chunksize

    @property
    def blocksize(self):
        return self.schunk.blocksize

    @property
    def nchunks(self):
        return self.schunk.nchunks

    @property
    def nbytes(self):
        return self.schunk.nbytes

    @property
    def cbytes(self):
        return self.schunk.cbytes

    @property
    def typesize(self):
        return self.schunk.typesize

    @property
    def urlpath(self):
        urlpath = self.schunk.storage.urlpath
        return urlpath.decode() if urlpath != NULL else None

    @property
    def contiguous(self):
        return self.schunk.storage.contiguous

    def get_cparams(self):
        if self.schunk.storage.cparams.compcode in blosc2.Codec._value2member_map_:
            codec = blosc2.Codec(self.schunk.storage.cparams.compcode)
        else:
            # User codec
            codec = self.schunk.storage.cparams.compcode

        filters = [0] * BLOSC2_MAX_FILTERS
        filters_meta = [0] * BLOSC2_MAX_FILTERS
        for i in range(BLOSC2_MAX_FILTERS):
            if self.schunk.filters[i] in blosc2.Filter._value2member_map_:
                filters[i] = blosc2.Filter(self.schunk.filters[i])
            else:
                # User filter
                filters[i] = self.schunk.filters[i]
            filters_meta[i] = self.schunk.filters_meta[i]

        cparams = blosc2.CParams(
            codec=codec,
            codec_meta=self.schunk.storage.cparams.compcode_meta,
            clevel=self.schunk.storage.cparams.clevel,
            use_dict=bool(self.schunk.storage.cparams.use_dict),
            typesize=self.schunk.storage.cparams.typesize,
            nthreads=self.schunk.storage.cparams.nthreads,
            blocksize=self.schunk.storage.cparams.blocksize,
            splitmode=blosc2.SplitMode(self.schunk.storage.cparams.splitmode),
            tuner=blosc2.Tuner(self.schunk.storage.cparams.tuner_id),
            filters=filters,
            filters_meta=filters_meta,
        )

        return cparams

    def update_cparams(self, new_cparams):
        cdef blosc2_cparams* cparams = self.schunk.storage.cparams
        codec = new_cparams.codec
        cparams.compcode = codec if not isinstance(codec, blosc2.Codec) else codec.value
        cparams.compcode_meta = new_cparams.codec_meta
        cparams.clevel = new_cparams.clevel
        cparams.use_dict = new_cparams.use_dict
        cparams.typesize = new_cparams.typesize
        cparams.nthreads = new_cparams.nthreads
        cparams.blocksize = new_cparams.blocksize
        cparams.splitmode = new_cparams.splitmode.value
        cparams.tuner_id = new_cparams.tuner.value

        filters = new_cparams.filters
        for i, filter in enumerate(filters):
            cparams.filters[i] = filter.value if isinstance(filter, Enum) else filter
        for i in range(len(filters), BLOSC2_MAX_FILTERS):
            cparams.filters[i] = 0

        filters_meta = new_cparams.filters_meta
        cdef int8_t meta_value
        for i, meta in enumerate(filters_meta):
            # We still may want to encode negative values
            meta_value = <int8_t> meta if meta < 0 else meta
            cparams.filters_meta[i] = <uint8_t> meta_value
        for i in range(len(filters_meta), BLOSC2_MAX_FILTERS):
            cparams.filters_meta[i] = 0

        _check_cparams(cparams)

        blosc2_free_ctx(self.schunk.cctx)
        self.schunk.cctx = blosc2_create_cctx(dereference(self.schunk.storage.cparams))
        if self.schunk.cctx == NULL:
            raise RuntimeError("Could not create compression context")
        self.schunk.compcode = self.schunk.storage.cparams.compcode
        self.schunk.compcode_meta = self.schunk.storage.cparams.compcode_meta
        self.schunk.clevel = self.schunk.storage.cparams.clevel
        self.schunk.splitmode = self.schunk.storage.cparams.splitmode
        self.schunk.typesize = self.schunk.storage.cparams.typesize
        self.schunk.blocksize = self.schunk.storage.cparams.blocksize
        self.schunk.filters = self.schunk.storage.cparams.filters
        self.schunk.filters_meta = self.schunk.storage.cparams.filters_meta

    def get_dparams(self):
        return blosc2.DParams(nthreads=self.schunk.storage.dparams.nthreads)

    def update_dparams(self, new_dparams):
        cdef blosc2_dparams* dparams = self.schunk.storage.dparams
        dparams.nthreads = new_dparams.nthreads

        _check_dparams(dparams, self.schunk.storage.cparams)

        blosc2_free_ctx(self.schunk.dctx)
        self.schunk.dctx = blosc2_create_dctx(dereference(self.schunk.storage.dparams))
        if self.schunk.dctx == NULL:
            raise RuntimeError("Could not create decompression context")

    def append_data(self, data):
        cdef Py_buffer buf
        PyObject_GetBuffer(data, &buf, PyBUF_SIMPLE)
        cdef int size
        cdef int32_t len_chunk = <int32_t> (buf.len + BLOSC2_MAX_OVERHEAD)
        cdef uint8_t* chunk = <uint8_t*> malloc(len_chunk)
        self.schunk.current_nchunk = self.schunk.nchunks
        if RELEASEGIL:
            with nogil:
                size = blosc2_compress_ctx(self.schunk.cctx, buf.buf, <int32_t> buf.len, chunk, len_chunk)
        else:
            size = blosc2_compress_ctx(self.schunk.cctx, buf.buf, <int32_t> buf.len, chunk, len_chunk)
        PyBuffer_Release(&buf)
        if size < 0:
            free(chunk)
            raise RuntimeError("Could not compress the data")
        elif size == 0:
            free(chunk)
            raise RuntimeError("The result could not fit")
        chunk = <uint8_t*> realloc(chunk, size)
        _check_comp_length('chunk', size)
        rc = blosc2_schunk_append_chunk(self.schunk, chunk, False)
        if rc < 0:
            free(chunk)
            raise RuntimeError("Could not append the chunk")
        return rc

    def fill_special(self, nitems, special_value, value):
        if value is None:
            return blosc2_schunk_fill_special(self.schunk, nitems, special_value, self.chunksize)

        if nitems == 0:
            return 0
        if nitems * self.typesize / self.chunksize > INT_MAX:
            raise RuntimeError("nitems is too large.  Try increasing the chunksize")
        if self.nbytes > 0 or self.cbytes > 0:
            raise RuntimeError("Filling with special values only works on empty SChunks")
        # Get a void pointer to the value
        array = np.array([value])
        if array.dtype.itemsize != self.typesize:
            if isinstance(value, int):
                dtype = np.dtype('i'+ str(self.typesize))
            elif isinstance(value, float):
                dtype = np.dtype('f' + str(self.typesize))
            else:
                raise ValueError("value size in bytes must match with typesize")
            array = np.array([value], dtype=dtype)
        cdef Py_buffer buf
        PyObject_GetBuffer(array, &buf, PyBUF_SIMPLE)
        # Create chunk with repeated values
        nchunks = nitems // self.chunkshape
        cdef blosc2_schunk *c_schunk = <blosc2_schunk *> self.c_schunk
        cdef blosc2_cparams *cparams = self.schunk.storage.cparams
        chunksize = BLOSC_EXTENDED_HEADER_LENGTH + self.typesize
        cdef void *chunk = malloc(chunksize)
        get_chunk_repeatval(dereference(cparams), self.chunksize, chunk, chunksize, &buf)

        for i in range(nchunks):
            if blosc2_schunk_append_chunk(self.schunk, <uint8_t *>chunk, True) < 0:
                free(chunk)
                PyBuffer_Release(&buf)
                raise RuntimeError("Error while appending the chunk")
        # Create and append last chunk if it is smaller than chunkshape
        remainder = nitems % self.chunkshape
        rc = 0
        if remainder != 0:
            get_chunk_repeatval(dereference(cparams), remainder * self.typesize, chunk, chunksize, &buf)
            rc = blosc2_schunk_append_chunk(self.schunk, <uint8_t *>chunk, True)
        free(chunk)
        PyBuffer_Release(&buf)
        if rc < 0:
            raise RuntimeError("Error while appending the chunk")

        return self.nchunks

    def decompress_chunk(self, nchunk, dst=None):
        cdef uint8_t *chunk
        cdef c_bool needs_free
        rc = blosc2_schunk_get_chunk(self.schunk, nchunk, &chunk, &needs_free)

        if rc < 0:
            raise RuntimeError("Error while getting the chunk")

        cdef int32_t nbytes
        cdef int32_t cbytes
        cdef int32_t blocksize
        blosc2_cbuffer_sizes(chunk, &nbytes, &cbytes, &blocksize)
        if needs_free:
            free(chunk)

        cdef Py_buffer buf
        if dst is not None:
            PyObject_GetBuffer(dst, &buf, PyBUF_SIMPLE)
            if buf.len == 0:
                raise ValueError("The dst length must be greater than 0")
            size = blosc2_schunk_decompress_chunk(self.schunk, nchunk, buf.buf, <int32_t>buf.len)
            PyBuffer_Release(&buf)
        else:
            dst = PyBytes_FromStringAndSize(NULL, nbytes)
            if dst is None:
                raise RuntimeError("Could not get a bytes object")
            size = blosc2_schunk_decompress_chunk(self.schunk, nchunk, <void*><char *>dst, nbytes)
            if size >= 0:
                return dst

        if size < 0:
            raise RuntimeError(f"Error while decompressing the specified chunk, error code: {size}")

    def get_chunk(self, nchunk):
        cdef uint8_t *chunk
        cdef c_bool needs_free
        cbytes = blosc2_schunk_get_chunk(self.schunk, nchunk, &chunk, &needs_free)
        if cbytes < 0:
           raise RuntimeError("Error while getting the chunk")
        ret_chunk = PyBytes_FromStringAndSize(<char*>chunk, cbytes)
        if needs_free:
            free(chunk)
        return ret_chunk

    def get_lazychunk(self, nchunk):
        cdef uint8_t *chunk
        cdef c_bool needs_free
        cbytes = blosc2_schunk_get_lazychunk(self.schunk, nchunk, &chunk, &needs_free)
        if cbytes < 0:
           raise RuntimeError("Error while getting the lazychunk")
        # The next does not always work (bug)
        # cdef uint8_t is_lazy = chunk[BLOSC2_MAX_OVERHEAD - 1] & 0x08
        # Workaround
        cdef uint8_t is_lazy = chunk[BLOSC2_MAX_OVERHEAD - 1] & 0x70
        if not is_lazy:
            # Put a cap on the buffer size for the non-lazy chunk
            cbytes = MAX_OVERHEAD
        ret_chunk = PyBytes_FromStringAndSize(<char*>chunk, cbytes)
        if needs_free:
            free(chunk)
        return ret_chunk

    def get_vlblock(self, nchunk, nblock):
        cdef uint8_t *block
        cdef int32_t destsize
        cbytes = blosc2_schunk_get_vlblock(self.schunk, nchunk, nblock, &block, &destsize)
        if cbytes < 0:
            raise RuntimeError("Error while getting the vlblock")
        ret_block = PyBytes_FromStringAndSize(<char*>block, destsize)
        free(block)
        return ret_block

    def delete_chunk(self, nchunk):
        rc = blosc2_schunk_delete_chunk(self.schunk, nchunk)
        if rc < 0:
            raise RuntimeError("Could not delete the desired chunk")
        return rc

    def append_chunk(self, chunk):
        cdef const uint8_t[:] typed_view_chunk
        mem_view_chunk = memoryview(chunk)
        typed_view_chunk = mem_view_chunk.cast('B')
        _check_comp_length('chunk', len(typed_view_chunk))
        rc = blosc2_schunk_append_chunk(self.schunk, &typed_view_chunk[0], True)
        if rc < 0:
            raise RuntimeError("Could not append the desired chunk")
        return rc

    def insert_chunk(self, nchunk, chunk):
        cdef const uint8_t[:] typed_view_chunk
        mem_view_chunk = memoryview(chunk)
        typed_view_chunk = mem_view_chunk.cast('B')
        _check_comp_length('chunk', len(typed_view_chunk))
        rc = blosc2_schunk_insert_chunk(self.schunk, nchunk, &typed_view_chunk[0], True)
        if rc < 0:
            raise RuntimeError("Could not insert the desired chunk")
        return rc

    def insert_data(self, nchunk, data, copy):
        cdef blosc2_context *cctx
        cdef Py_buffer buf
        PyObject_GetBuffer(data, &buf, PyBUF_SIMPLE)
        cdef int size
        cdef int32_t len_chunk = <int32_t> (buf.len + BLOSC2_MAX_OVERHEAD)
        cdef uint8_t* chunk = <uint8_t*> malloc(len_chunk)
        self.schunk.current_nchunk = nchunk  # prefilter needs this value to be set
        if RELEASEGIL:
            with nogil:
                # No need to create another cctx
                size = blosc2_compress_ctx(self.schunk.cctx, buf.buf, <int32_t> buf.len, chunk, len_chunk)
        else:
            size = blosc2_compress_ctx(self.schunk.cctx, buf.buf, <int32_t> buf.len, chunk, len_chunk)
        PyBuffer_Release(&buf)
        if size < 0:
            raise RuntimeError("Could not compress the data")
        elif size == 0:
            free(chunk)
            raise RuntimeError("The result could not fit ")

        chunk = <uint8_t*> realloc(chunk, size)
        _check_comp_length('chunk', size)
        rc = blosc2_schunk_insert_chunk(self.schunk, nchunk, chunk, copy)
        if copy:
            free(chunk)
        if rc < 0:
            raise RuntimeError("Could not insert the desired chunk")
        return rc

    def update_chunk(self, nchunk, chunk):
        cdef const uint8_t[:] typed_view_chunk
        mem_view_chunk = memoryview(chunk)
        typed_view_chunk = mem_view_chunk.cast('B')
        _check_comp_length('chunk', len(typed_view_chunk))
        rc = blosc2_schunk_update_chunk(self.schunk, nchunk, &typed_view_chunk[0], True)
        if rc < 0:
            raise RuntimeError("Could not update the desired chunk")
        return rc

    def reorder_offsets(self, order):
        cdef np.ndarray[np.int64_t, ndim=1] offsets_order = np.ascontiguousarray(order, dtype=np.int64)
        rc = blosc2_schunk_reorder_offsets(self.schunk, <int64_t*> offsets_order.data)
        if rc < 0:
            raise RuntimeError("Could not reorder the chunk offsets")
        return None

    def update_data(self, nchunk, data, copy):
        cdef Py_buffer buf
        PyObject_GetBuffer(data, &buf, PyBUF_SIMPLE)
        cdef int size
        cdef int32_t len_chunk = <int32_t> (buf.len + BLOSC2_MAX_OVERHEAD)
        cdef uint8_t* chunk = <uint8_t*> malloc(len_chunk)
        self.schunk.current_nchunk = nchunk  # prefilter needs this value to be set
        if RELEASEGIL:
            with nogil:
                size = blosc2_compress_ctx(self.schunk.cctx, buf.buf, <int32_t> buf.len, chunk, len_chunk)
        else:
            size = blosc2_compress_ctx(self.schunk.cctx, buf.buf, <int32_t> buf.len, chunk, len_chunk)

        PyBuffer_Release(&buf)
        if size < 0:
            raise RuntimeError("Could not compress the data")
        elif size == 0:
            free(chunk)
            raise RuntimeError("The result could not fit ")

        chunk = <uint8_t*> realloc(chunk, size)
        _check_comp_length('chunk', size)
        rc = blosc2_schunk_update_chunk(self.schunk, nchunk, chunk, copy)
        if copy:
            free(chunk)
        if rc < 0:
            raise RuntimeError("Could not update the desired chunk")
        return rc

    # This is used internally for prefiltering
    def _prefilter_data(self, nchunk, data, chunk_data):
        cdef Py_buffer buf
        PyObject_GetBuffer(data, &buf, PyBUF_SIMPLE)
        cdef Py_buffer chunk_buf
        PyObject_GetBuffer(chunk_data, &chunk_buf, PyBUF_SIMPLE)
        self.schunk.current_nchunk = nchunk  # prefilter needs this value to be set
        cdef int size = blosc2_compress_ctx(self.schunk.cctx, buf.buf, <int32_t> buf.len, chunk_buf.buf, chunk_buf.len)
        PyBuffer_Release(&buf)
        PyBuffer_Release(&chunk_buf)
        if size < 0:
            raise RuntimeError("Could not compress the data")
        elif size == 0:
            raise RuntimeError("The result could not fit ")
        return size

    def get_slice(self, start=0, stop=None, out=None):
        cdef int64_t nitems = self.schunk.nbytes // self.schunk.typesize
        start, stop, _ = slice(start, stop, 1).indices(nitems)
        if start >= stop:
            return b''

        cdef Py_ssize_t nbytes = (stop - start) * self.schunk.typesize
        cdef Py_buffer buf
        if out is not None:
            PyObject_GetBuffer(out, &buf, PyBUF_SIMPLE)
            if buf.len < nbytes:
                raise ValueError("Not enough space for writing the slice in out")
            rc = blosc2_schunk_get_slice_buffer(self.schunk, start, stop, buf.buf)
            PyBuffer_Release(&buf)
        else:
            out = PyBytes_FromStringAndSize(NULL, nbytes)
            if out is None:
                raise RuntimeError("Could not get a bytes object")
            rc = blosc2_schunk_get_slice_buffer(self.schunk, start, stop, <void*><char *> out)
            if rc >= 0:
                return out
        if rc < 0:
            raise RuntimeError("Error while getting the slice")

    def set_slice(self, value, start=0, stop=None):
        cdef int64_t nitems = self.schunk.nbytes // self.schunk.typesize
        start, stop = self._massage_key(start, stop, nitems)
        if start > nitems:
            raise ValueError("`start` cannot be greater than the SChunk nitems")

        cdef int64_t nbytes = (stop - start) * self.schunk.typesize

        cdef Py_buffer buf
        PyObject_GetBuffer(value, &buf, PyBUF_SIMPLE)
        cdef uint8_t *buf_ptr = <uint8_t *> buf.buf
        cdef int64_t buf_pos = 0
        cdef int64_t nbytes_copy = min(nbytes, buf.len - buf_pos)
        cdef int64_t data_start
        cdef uint8_t *data
        cdef uint8_t *chunk
        cdef int32_t alloc_len
        cdef int32_t chunk_nbytes
        cdef int32_t chunksize
        cdef int comp_rc
        if buf.len < nbytes:
            raise ValueError("Not enough data for writing the slice")

        if stop > nitems:
            # Increase SChunk's size
            if start < nitems:
                rc = blosc2_schunk_set_slice_buffer(self.schunk, start, nitems, buf.buf)
                buf_pos = (nitems - start) * self.schunk.typesize
            if self.schunk.nbytes % self.schunk.chunksize != 0:
                # Update last chunk before appending any other
                if stop * self.schunk.typesize >= self.schunk.chunksize * self.schunk.nchunks:
                    chunk_nbytes = self.schunk.chunksize
                    nbytes_copy = min(nbytes_copy, self.schunk.chunksize * self.schunk.nchunks - nitems * self.schunk.typesize)
                else:
                    chunk_nbytes = (stop * self.schunk.typesize) % self.schunk.chunksize
                data  = <uint8_t *> malloc(chunk_nbytes)
                rc = blosc2_schunk_decompress_chunk(self.schunk, self.schunk.nchunks - 1, data, chunk_nbytes)
                if rc < 0:
                    free(data)
                    PyBuffer_Release(&buf)
                    raise RuntimeError("Error while decompressing the chunk")
                data_start = self.schunk.nbytes - (self.schunk.nchunks - 1) * self.schunk.chunksize
                memcpy(data + data_start, buf_ptr + buf_pos, nbytes_copy)
                chunk = <uint8_t *> malloc(chunk_nbytes + BLOSC2_MAX_OVERHEAD)
                self.schunk.current_nchunk = self.schunk.nchunks - 1
                if RELEASEGIL:
                    with nogil:
                        comp_rc = blosc2_compress_ctx(self.schunk.cctx, data, chunk_nbytes, chunk, chunk_nbytes + BLOSC2_MAX_OVERHEAD)
                else:
                    comp_rc = blosc2_compress_ctx(self.schunk.cctx, data, chunk_nbytes, chunk, chunk_nbytes + BLOSC2_MAX_OVERHEAD)
                free(data)
                if comp_rc < 0:
                    free(chunk)
                    PyBuffer_Release(&buf)
                    raise RuntimeError("Error while compressing the data")
                elif comp_rc == 0:
                    free(chunk)
                    PyBuffer_Release(&buf)
                    raise RuntimeError("The result could not fit")
                rc = blosc2_schunk_update_chunk(self.schunk, self.schunk.nchunks - 1, chunk, True)
                free(chunk)
                if rc < 0:
                    PyBuffer_Release(&buf)
                    raise RuntimeError("Error while updating the chunk")
                buf_pos += nbytes_copy
            # Append data if needed
            if buf_pos < buf.len:
                nappends = int(stop * self.schunk.typesize / self.schunk.chunksize - self.schunk.nchunks)
                if (stop * self.schunk.typesize) % self.schunk.chunksize != 0:
                    nappends += 1
                for i in range(nappends):
                    if (self.schunk.nchunks + 1) * self.schunk.chunksize <= stop * self.schunk.typesize:
                        chunksize = self.schunk.chunksize
                    else:
                        chunksize = (stop * self.schunk.typesize) % self.schunk.chunksize
                    alloc_len = <int32_t> (chunksize + BLOSC2_MAX_OVERHEAD)
                    chunk = <uint8_t*> malloc(alloc_len)
                    self.schunk.current_nchunk = self.schunk.nchunks
                    if RELEASEGIL:
                        with nogil:
                            comp_rc = blosc2_compress_ctx(self.schunk.cctx, buf_ptr + buf_pos, chunksize, chunk, alloc_len)
                    else:
                        comp_rc = blosc2_compress_ctx(self.schunk.cctx, buf_ptr + buf_pos, chunksize, chunk, alloc_len)
                    if comp_rc < 0:
                        free(chunk)
                        PyBuffer_Release(&buf)
                        raise RuntimeError("Error while compressing the chunk")
                    elif comp_rc == 0:
                        free(chunk)
                        PyBuffer_Release(&buf)
                        raise RuntimeError("The result could not fit")
                    chunk = <uint8_t*> realloc(chunk, comp_rc)
                    _check_comp_length('chunk', comp_rc)
                    rc = blosc2_schunk_append_chunk(self.schunk, chunk, False)
                    if rc < 0:
                        free(chunk)
                        PyBuffer_Release(&buf)
                        raise RuntimeError("Error while appending the chunk")
                    buf_pos += chunksize
        else:
            rc = blosc2_schunk_set_slice_buffer(self.schunk, start, stop, buf.buf)
        PyBuffer_Release(&buf)
        if rc < 0:
            raise RuntimeError("Error while setting the slice")

    def to_cframe(self):
        cdef c_bool needs_free
        cdef uint8_t *cframe
        cdef int64_t cframe_len
        if RELEASEGIL:
            with nogil:
                cframe_len = blosc2_schunk_to_buffer(self.schunk, &cframe, &needs_free)
        else:
            cframe_len = blosc2_schunk_to_buffer(self.schunk, &cframe, &needs_free)
        if cframe_len < 0:
            raise RuntimeError("Error while getting the cframe")
        out = PyBytes_FromStringAndSize(<char*>cframe, cframe_len)
        if needs_free:
            free(cframe)

        return out

    def _avoid_cframe_free(self, avoid_cframe_free):
        blosc2_schunk_avoid_cframe_free(self.schunk, avoid_cframe_free)

    def lock(self):
        """Take the exclusive frame lock and hold it across several operations.

        A no-op when file locking is not enabled on this handle.  Use the
        `holding_lock()` context manager instead of calling this directly.
        """
        cdef int rc
        with nogil:  # may block waiting for other processes
            rc = blosc2_schunk_lock(self.schunk)
        if rc < 0:
            raise RuntimeError("Could not lock the super-chunk frame")

    def unlock(self):
        """Release the frame lock taken with `lock()`."""
        cdef int rc = blosc2_schunk_unlock(self.schunk)
        if rc < 0:
            raise RuntimeError("Could not unlock the super-chunk frame")

    def refresh(self):
        cdef int rc = blosc2_schunk_refresh(self.schunk)
        _check_rc(rc, "Error while refreshing the schunk")
        return rc == 1

    @property
    def change_tick(self):
        """Counter bumped whenever the (vl)metalayers change, locally or via a
        re-sync after a mutation made through another handle."""
        return self.schunk.change_tick

    def _massage_key(self, start, stop, nitems):
        if stop is None:
            stop = nitems
        elif stop < 0:
            stop += nitems
        if start is None:
            start = 0
        elif start < 0:
            start += nitems
        if stop - start <= 0:
            raise ValueError("`stop` mut be greater than `start`")

        return start, stop

    def _set_postfilter(self, func, dtype_input, dtype_output=None):
        # Get user data
        func_id = func.__name__
        blosc2.postfilter_funcs[func_id] = func
        func_id = func_id.encode("utf-8") if isinstance(func_id, str) else func_id

        dtype_output = dtype_input if dtype_output is None else dtype_output
        dtype_input = np.dtype(dtype_input)
        dtype_output = np.dtype(dtype_output)
        if dtype_output.itemsize != dtype_input.itemsize:
            del blosc2.postfilter_funcs[func_id]
            raise ValueError("`dtype_input` and `dtype_output` must have the same size")

        # Set postfilter
        cdef blosc2_dparams* dparams = self.schunk.storage.dparams
        dparams.postfilter = <blosc2_postfilter_fn> general_postfilter
        # Fill postparams
        cdef blosc2_postfilter_params* postparams = <blosc2_postfilter_params *> malloc(sizeof(blosc2_postfilter_params))
        cdef user_filters_udata* postf_udata = <user_filters_udata * > malloc(sizeof(user_filters_udata))
        postf_udata.py_func = <char * > malloc(strlen(func_id) + 1)
        strcpy(postf_udata.py_func, func_id)
        postf_udata.input_cdtype = dtype_input.num
        postf_udata.output_cdtype = dtype_output.num
        postf_udata.chunkshape = self.schunk.chunksize // self.schunk.typesize

        postparams.user_data = postf_udata
        dparams.postparams = postparams
        _check_dparams(dparams, self.schunk.storage.cparams)

        blosc2_free_ctx(self.schunk.dctx)
        self.schunk.dctx = blosc2_create_dctx(dereference(dparams))
        if self.schunk.dctx == NULL:
            raise RuntimeError("Could not create decompression context")

    cpdef remove_postfilter(self, func_name, _new_ctx=True):
        if func_name is not None:
            del blosc2.postfilter_funcs[func_name]

        cdef user_filters_udata* udata = <user_filters_udata * > self.schunk.storage.dparams.postparams.user_data
        free(udata.py_func)
        free(self.schunk.storage.dparams.postparams.user_data)
        free(self.schunk.storage.dparams.postparams)
        self.schunk.storage.dparams.postparams = NULL
        self.schunk.storage.dparams.postfilter = NULL

        blosc2_free_ctx(self.schunk.dctx)
        if _new_ctx:
            self.schunk.dctx = blosc2_create_dctx(dereference(self.schunk.storage.dparams))
            if self.schunk.dctx == NULL:
                raise RuntimeError("Could not create decompression context")
        else:
            # Avoid creating new dctx when calling this from the __dealloc__
            self.schunk.dctx = NULL

    def _set_filler(self, func, inputs_id, dtype_output):
        if self.schunk.storage.cparams.nthreads > 1:
            raise AttributeError("compress `nthreads` must be 1 when assigning a prefilter")

        func_id = func.__name__
        blosc2.prefilter_funcs[func_id] = func
        func_id = func_id.encode("utf-8") if isinstance(func_id, str) else func_id

        # Set prefilter
        cdef blosc2_cparams* cparams = self.schunk.storage.cparams
        cparams.prefilter = <blosc2_prefilter_fn> general_filler

        cdef blosc2_prefilter_params* preparams = <blosc2_prefilter_params *> calloc(1, sizeof(blosc2_prefilter_params))
        cdef filler_udata* fill_udata = <filler_udata *> malloc(sizeof(filler_udata))
        fill_udata.py_func = <char *> malloc(strlen(func_id) + 1)
        strcpy(fill_udata.py_func, func_id)
        fill_udata.inputs_id = inputs_id
        fill_udata.output_cdtype = np.dtype(dtype_output).num
        fill_udata.chunkshape = self.schunk.chunksize // self.schunk.typesize

        preparams.user_data = fill_udata
        cparams.preparams = preparams
        _check_cparams(cparams)

        blosc2_free_ctx(self.schunk.cctx)
        self.schunk.cctx = blosc2_create_cctx(dereference(cparams))
        if self.schunk.cctx == NULL:
            raise RuntimeError("Could not create compression context")

    def _set_prefilter(self, func, dtype_input, dtype_output=None):
        if self.schunk.storage.cparams.nthreads > 1:
            raise AttributeError("compress `nthreads` must be 1 when assigning a prefilter")
        func_id = func.__name__
        blosc2.prefilter_funcs[func_id] = func
        func_id = func_id.encode("utf-8") if isinstance(func_id, str) else func_id

        dtype_output = dtype_input if dtype_output is None else dtype_output
        dtype_input = np.dtype(dtype_input)
        dtype_output = np.dtype(dtype_output)
        if dtype_output.itemsize != dtype_input.itemsize:
            del blosc2.prefilter_funcs[func_id]
            raise ValueError("`dtype_input` and `dtype_output` must have the same size")

        cdef blosc2_cparams* cparams = self.schunk.storage.cparams
        cparams.prefilter = <blosc2_prefilter_fn> general_prefilter
        cdef blosc2_prefilter_params* preparams = <blosc2_prefilter_params *> calloc(1, sizeof(blosc2_prefilter_params))
        cdef user_filters_udata* pref_udata = <user_filters_udata*> malloc(sizeof(user_filters_udata))
        pref_udata.py_func = <char *> malloc(strlen(func_id) + 1)
        strcpy(pref_udata.py_func, func_id)
        pref_udata.input_cdtype = dtype_input.num
        pref_udata.output_cdtype = dtype_output.num
        pref_udata.chunkshape = self.schunk.chunksize // self.schunk.typesize

        preparams.user_data = pref_udata
        cparams.preparams = preparams
        _check_cparams(cparams)

        if self.schunk.cctx != NULL:
            # Freeing NULL context can lead to segmentation fault
            blosc2_free_ctx(self.schunk.cctx)
        self.schunk.cctx = blosc2_create_cctx(dereference(cparams))
        if self.schunk.cctx == NULL:
            raise RuntimeError("Could not create compression context")

    cpdef remove_prefilter(self, func_name, _new_ctx=True):
        cdef udf_udata* udf_data
        cdef user_filters_udata* udata
        cdef mm_udata* mm_data
        cdef me_udata* me_data
        cdef me_input_cache_s* input_cache
        cdef int i

        if func_name is not None and func_name in blosc2.prefilter_funcs:
            del blosc2.prefilter_funcs[func_name]

        # Clean up the miniexpr handle if this is a miniexpr_prefilter
        if self.schunk.storage.cparams.prefilter == <blosc2_prefilter_fn>miniexpr_prefilter:
            if self.schunk.storage.cparams.preparams != NULL:
                me_data = <me_udata*>self.schunk.storage.cparams.preparams.user_data
                if me_data != NULL:
                    if me_data.input_chunk_caches != NULL:
                        for i in range(me_data.ninputs):
                            input_cache = &me_data.input_chunk_caches[i]
                            if input_cache.data != NULL:
                                free(input_cache.data)
                                input_cache.data = NULL
                            input_cache.nchunk = -1
                            input_cache.state = ME_CACHE_EMPTY
                            if input_cache.state_lock != NULL:
                                PyThread_free_lock(input_cache.state_lock)
                                input_cache.state_lock = NULL
                            if input_cache.ready_lock != NULL:
                                PyThread_free_lock(input_cache.ready_lock)
                                input_cache.ready_lock = NULL
                        free(me_data.input_chunk_caches)
                    if me_data.inputs != NULL:
                        free(me_data.inputs)
                    if me_data.miniexpr_handle != NULL:  # XXX do we really need the conditional?
                        me_free(me_data.miniexpr_handle)
                    if me_data.eval_params != NULL:
                        free(me_data.eval_params)
                    free(me_data)
        elif self.schunk.storage.cparams.prefilter == <blosc2_prefilter_fn>matmul_prefilter:
            if self.schunk.storage.cparams.preparams != NULL:
                mm_data = <mm_udata*>self.schunk.storage.cparams.preparams.user_data
                if mm_data != NULL:
                    if mm_data.inputs != NULL:
                        free(mm_data.inputs)
                    free(mm_data)
        elif self.schunk.storage.cparams.prefilter != NULL:
            # From Python the preparams->udata with always have the field py_func
            if self.schunk.storage.cparams.preparams != NULL:
                udata = <user_filters_udata*>self.schunk.storage.cparams.preparams.user_data
                if udata != NULL:
                    if udata.py_func != NULL:
                        free(udata.py_func)
                    free(udata)

        if self.schunk.storage.cparams.preparams != NULL:
            free(self.schunk.storage.cparams.preparams)
        self.schunk.storage.cparams.preparams = NULL
        self.schunk.storage.cparams.prefilter = NULL

        if self.schunk.cctx != NULL:
            # Freeing NULL context can lead to segmentation fault
            blosc2_free_ctx(self.schunk.cctx)
        if _new_ctx:
            self.schunk.cctx = blosc2_create_cctx(dereference(self.schunk.storage.cparams))
            if self.schunk.cctx == NULL:
                raise RuntimeError("Could not create compression context")
        else:
            # Avoid creating new cctx when calling this from the __dealloc__
            self.schunk.cctx = NULL

    def __dealloc__(self):
        cdef blosc2_schunk *schunk_ptr
        if self.schunk != NULL and not self._is_view:
            # Free prefilters and postfilters params
            if self.schunk.storage.cparams.prefilter != NULL:
                self.remove_prefilter(func_name=None, _new_ctx=False)
            if self.schunk.storage.dparams.postfilter != NULL:
                self.remove_postfilter(func_name=None, _new_ctx=False)

            # Free the C-Blosc2 super-chunk with the GIL held so threadpool
            # teardown cannot race with active miniexpr workers.
            schunk_ptr = self.schunk
            self.schunk = NULL
            blosc2_schunk_free(schunk_ptr)


# postfilter
cdef int general_postfilter(blosc2_postfilter_params *params):
    cdef user_filters_udata *udata = <user_filters_udata *> params.user_data
    cdef int nd = 1
    cdef np.npy_intp dims = params.size // params.typesize
    input = np.PyArray_SimpleNewFromData(nd, &dims, udata.input_cdtype, <void*>params.input)
    output = np.PyArray_SimpleNewFromData(nd, &dims, udata.output_cdtype, <void*>params.output)
    offset = params.nchunk * udata.chunkshape + params.offset // params.typesize
    func_id = udata.py_func.decode("utf-8")
    blosc2.postfilter_funcs[func_id](input, output, offset)
    return 0


# filler
cdef int general_filler(blosc2_prefilter_params *params):
    cdef filler_udata *udata = <filler_udata *> params.user_data
    cdef int nd = 1
    cdef np.npy_intp dims = params.output_size // params.output_typesize

    inputs_tuple = _ctypes.PyObj_FromPtr(udata.inputs_id)

    output = np.PyArray_SimpleNewFromData(nd, &dims, udata.output_cdtype, <void*>params.output)
    offset = params.nchunk * udata.chunkshape + params.output_offset // params.output_typesize

    inputs = []
    for obj, dtype in inputs_tuple:
        if isinstance(obj, blosc2.SChunk):
            out = np.empty(dims, dtype=dtype)
            obj.get_slice(start=offset, stop=offset + dims, out=out)
            inputs.append(out)
        elif isinstance(obj, np.ndarray):
            inputs.append(obj[offset : offset + dims])
        elif isinstance(obj, (int, float, bool, complex)):
            inputs.append(np.full(dims, obj, dtype=dtype))
        else:
            raise ValueError("Unsupported operand")

    func_id = udata.py_func.decode("utf-8")
    blosc2.prefilter_funcs[func_id](tuple(inputs), output, offset)

    return 0


# Auxiliary function for miniexpr as a prefilter
# Only meant for (input and output) arrays that are blosc2.NDArray objects.
cdef int aux_miniexpr(me_udata *udata, int64_t nchunk, int32_t nblock,
                      c_bool is_postfilter, uint8_t *params_output, int32_t typesize) nogil:
    # Declare all C variables at the beginning
    cdef int64_t chunk_ndim[B2ND_MAX_DIM]
    cdef int64_t block_ndim[B2ND_MAX_DIM]
    cdef int64_t start_ndim[B2ND_MAX_DIM]
    cdef int64_t stop_ndim[B2ND_MAX_DIM]
    cdef int64_t buffershape[B2ND_MAX_DIM]

    cdef b2nd_array_t* ndarr
    cdef int rc
    cdef void** input_buffers = <void**> malloc(udata.ninputs * sizeof(uint8_t*))
    cdef uint8_t* src
    cdef uint8_t* chunk
    cdef c_bool needs_free
    cdef uint8_t* loaded_chunk
    cdef me_input_cache_s* input_cache
    cdef int32_t chunk_nbytes, chunk_cbytes, block_nbytes
    cdef int start, blocknitems, expected_blocknitems
    cdef int64_t valid_nitems
    cdef int64_t global_block
    cdef int32_t input_typesize
    cdef blosc2_context* dctx
    expected_blocknitems = -1
    valid_nitems = 0

    cdef me_expr* miniexpr_handle = udata.miniexpr_handle
    cdef void* aux_reduc_ptr

    if miniexpr_handle == NULL:
        raise ValueError("miniexpr: handle not assigned")
    if input_buffers == NULL:
        raise MemoryError("miniexpr: cannot allocate input buffer table")
    memset(input_buffers, 0, udata.ninputs * sizeof(uint8_t*))

    # Query valid (unpadded) items for this block
    rc = me_nd_valid_nitems(miniexpr_handle, nchunk, nblock, &valid_nitems)
    if rc != 0:
        raise RuntimeError(f"miniexpr: invalid block; error code: {rc}")
    if valid_nitems <= 0:
        # Nothing to compute for this block.
        # For reductions, keep aux_reduc neutral values untouched.
        if udata.aux_reduc_ptr == NULL:
            memset(params_output, 0, udata.array.blocknitems * typesize)
        free(input_buffers)
        return 0

    # Candidate-block pruning (1-D arrays only): when a bitmap is supplied and
    # this block is not a candidate, write a zero (false) result and return
    # without decompressing inputs or running miniexpr.  Disabled for reductions
    # (aux_reduc_ptr) to avoid perturbing accumulator state.
    if udata.candidate_blocks != NULL and udata.aux_reduc_ptr == NULL:
        global_block = nchunk * udata.blocks_in_chunk[0] + nblock
        if global_block < udata.candidate_blocks_len and udata.candidate_blocks[global_block] == 0:
            memset(params_output, 0, udata.array.blocknitems * typesize)
            free(input_buffers)
            return 0

    for i in range(udata.ninputs):
        ndarr = udata.inputs[i]
        if ndarr.sc.storage.urlpath == NULL:
            src = ndarr.sc.data[nchunk]
        else:
            input_cache = &udata.input_chunk_caches[i]
            if input_cache.state_lock == NULL or input_cache.ready_lock == NULL:
                raise MemoryError("miniexpr: cache locks not assigned")
            while True:
                PyThread_acquire_lock(input_cache.state_lock, 1)
                if input_cache.state == ME_CACHE_READY and input_cache.nchunk == nchunk and input_cache.data != NULL:
                    src = input_cache.data
                    PyThread_release_lock(input_cache.state_lock)
                    break
                if input_cache.state == ME_CACHE_ERROR and input_cache.nchunk == nchunk:
                    PyThread_release_lock(input_cache.state_lock)
                    raise ValueError("miniexpr: error getting chunk")
                if input_cache.state == ME_CACHE_LOADING:
                    PyThread_release_lock(input_cache.state_lock)
                    PyThread_acquire_lock(input_cache.ready_lock, 1)
                    PyThread_release_lock(input_cache.ready_lock)
                    continue
                PyThread_acquire_lock(input_cache.ready_lock, 1)
                input_cache.state = ME_CACHE_LOADING
                input_cache.nchunk = nchunk
                PyThread_release_lock(input_cache.state_lock)

                rc = blosc2_schunk_get_chunk(ndarr.sc, nchunk, &chunk, &needs_free)
                if rc < 0:
                    PyThread_acquire_lock(input_cache.state_lock, 1)
                    input_cache.state = ME_CACHE_ERROR
                    PyThread_release_lock(input_cache.state_lock)
                    PyThread_release_lock(input_cache.ready_lock)
                    raise ValueError("miniexpr: error getting chunk")

                if not needs_free:
                    loaded_chunk = <uint8_t*> malloc(rc)
                    if loaded_chunk == NULL:
                        PyThread_acquire_lock(input_cache.state_lock, 1)
                        input_cache.state = ME_CACHE_ERROR
                        PyThread_release_lock(input_cache.state_lock)
                        PyThread_release_lock(input_cache.ready_lock)
                        raise MemoryError("miniexpr: cannot allocate chunk copy")
                    memcpy(loaded_chunk, chunk, rc)
                else:
                    loaded_chunk = chunk

                PyThread_acquire_lock(input_cache.state_lock, 1)
                if input_cache.data != NULL:
                    free(input_cache.data)
                input_cache.data = loaded_chunk
                input_cache.nchunk = nchunk
                input_cache.state = ME_CACHE_READY
                src = input_cache.data
                PyThread_release_lock(input_cache.state_lock)
                PyThread_release_lock(input_cache.ready_lock)
                break
        rc = blosc2_cbuffer_sizes(src, &chunk_nbytes, &chunk_cbytes, &block_nbytes)
        if rc < 0:
            raise ValueError("miniexpr: error getting cbuffer sizes")
        if block_nbytes <= 0:
            raise ValueError("miniexpr: invalid block size")
        input_buffers[i] = malloc(block_nbytes)
        if input_buffers[i] == NULL:
            raise MemoryError("miniexpr: cannot allocate input block buffer")
        input_typesize = ndarr.sc.typesize
        blocknitems = block_nbytes // input_typesize
        if expected_blocknitems == -1:
            expected_blocknitems = blocknitems
        elif blocknitems != expected_blocknitems:
            raise ValueError("miniexpr: inconsistent block element counts across inputs")
        start = nblock * blocknitems
        # This is needed for thread safety, but adds a pretty low overhead (< 400ns on a modern CPU)
        # In the future, perhaps one can create a specific (serial) context just for
        # blosc2_getitem_ctx, but this is probably never going to be necessary.
        dctx = blosc2_create_dctx(BLOSC2_DPARAMS_DEFAULTS)
        # Unsafe, but it works for special arrays (e.g. blosc2.ones), and can be used for profiling
        # dctx = ndarr.sc.dctx
        if valid_nitems > blocknitems:
            raise ValueError("miniexpr: valid items exceed padded block size")
        rc = blosc2_getitem_ctx(dctx, src, chunk_cbytes, start, blocknitems,
                                input_buffers[i], block_nbytes)
        blosc2_free_ctx(dctx)
        if rc < 0:
            raise ValueError("miniexpr: error decompressing the chunk")
    # For reduction operations, we need to track which block we're processing
    # The linear_block_index should be based on the INPUT array structure, not the output array
    # Get the first input array's chunk and block structure
    cdef b2nd_array_t* first_input = udata.inputs[0]
    cdef int nblocks_per_chunk = 1
    for i in range(first_input.ndim):
        nblocks_per_chunk *= <int>udata.blocks_in_chunk[i]
    # Calculate the global linear block index: nchunk * blocks_per_chunk + nblock
    # This works because blocks never span chunks (chunks are padded to block boundaries)
    cdef int64_t linear_block_index = nchunk * nblocks_per_chunk + nblock
    cdef uintptr_t offset_bytes = typesize * linear_block_index

    # Call thread-safe miniexpr C API
    # NOTE: me_eval_nd expects the OUTPUT block size (in items), not the input block size.
    # For element-wise operations with same dtypes, they're equal, but for type-changing
    # operations (e.g., arccos(int32) -> float64), we must use the output's block item count.
    cdef int output_blocknitems = udata.array.blocknitems

    if udata.aux_reduc_ptr == NULL:
        aux_reduc_ptr = <void *> params_output
    else:
        # Reduction operation: evaluate only valid items into a single output element.
        # NOTE: miniexpr handles scalar outputs in me_eval_nd without touching tail bytes.
        aux_reduc_ptr = <void *> (<uintptr_t> udata.aux_reduc_ptr + offset_bytes)
    rc = me_eval_nd(miniexpr_handle, <const void**> input_buffers, udata.ninputs,
                    aux_reduc_ptr, output_blocknitems, nchunk, nblock, udata.eval_params)
    if rc != 0:
        raise RuntimeError(f"miniexpr: issues during evaluation; error code: {rc}")

    # Free resources
    for i in range(udata.ninputs):
        free(input_buffers[i])
    free(input_buffers)

    return 0

cdef int matmul_block_kernel(T* A, T* B, T* C, int M, int K, int N) nogil:
    cdef int r, c, k
    cdef T a
    cdef int rowA, rowC, rowB
    for r in range(M):
        rowA = r * K
        rowC = r * N
        for k in range(K):
            a = A[rowA + k]
            rowB = k * N
            for c in range(N):
                C[rowC + c] += <T>(a * B[rowB + c])
    return 0

cdef int aux_matmul(mm_udata *udata, int64_t nchunk, int32_t nblock, void *params_output, int32_t typesize, int typecode) nogil:
    # Declare all C variables at the beginning
    cdef b2nd_array_t* out_arr
    cdef b2nd_array_t* ndarr
    cdef c_bool first_run
    cdef int rc, p, q, r
    cdef void** input_buffers = <void**> malloc(2 * sizeof(uint8_t*))
    cdef uint8_t** src = <uint8_t**> malloc(2 * sizeof(uint8_t*))
    cdef int32_t chunk_nbytes[2]
    cdef int32_t chunk_cbytes[2]
    cdef int32_t block_nbytes[2]
    cdef int blocknitems[2]
    cdef int startA, startB, expected_blocknitems
    cdef blosc2_context* dctx
    cdef int i, j, block_i, block_j, chunk_i, chunk_j, ncols, block_ncols, Bblock_ncols, Bncols, Ablock_ncols, Ancols
    cdef int nchunkA = 0, nchunkB = 0, nblockA = 0, nblockB = 0, offsetA = 0, offsetB = 0, offset = 0
    out_arr = udata.array
    cdef int ndim = out_arr.ndim
    cdef int nchunk_ = nchunk
    cdef int coord, batch, batch_, batches = 1
    cdef int out_chunk_nrows, out_chunk_ncols, out_block_nrows, out_block_ncols
    cdef int selected_backend = b2_get_selected_matmul_backend()

    # batches = sum(strides[i]*elcoords[i])
    for i in range(ndim - 2):
        batches *= out_arr.blockshape[i]

    # nchunk = sum(strides[i]*chunkcoords[i])
    for i in range(ndim - 2):
        coord = nchunk_ // udata.chunks_strides[0][i]
        nchunk_ = nchunk_ % udata.chunks_strides[0][i]
        nchunkA += coord * udata.chunks_strides[1][i]
        nchunkB += coord * udata.chunks_strides[2][i]

    ncols = udata.chunks_strides[0][ndim - 2]
    Ancols = udata.chunks_strides[1][ndim - 2]
    Bncols = udata.chunks_strides[2][ndim - 2]
    out_chunk_nrows = out_arr.chunkshape[ndim - 2]
    out_chunk_ncols = out_arr.chunkshape[ndim - 1]

    # nblock = sum(strides[i]*blockcoords[i])
    cdef int nblock_ = nblock
    for i in range(ndim - 2):
        coord = nblock_ // udata.blocks_strides[0][i]
        nblock_ = nblock_ % udata.blocks_strides[0][i]
        nblockA += coord * udata.blocks_strides[1][i]
        nblockB += coord * udata.blocks_strides[2][i]

    block_ncols = udata.blocks_strides[0][ndim - 2]
    Ablock_ncols = udata.blocks_strides[1][ndim - 2]
    Bblock_ncols = udata.blocks_strides[2][ndim - 2]
    out_block_nrows = out_arr.blockshape[ndim - 2]
    out_block_ncols = out_arr.blockshape[ndim - 1]

    memset(params_output, 0, out_arr.blocknitems * typesize)

    dctx = blosc2_create_dctx(BLOSC2_DPARAMS_DEFAULTS)

    first_run = True
    while True: # chunk loop
        for i in range(2):
            chunk_idx = nchunkA if i == 0 else nchunkB
            ndarr = udata.inputs[i]
            ndim = ndarr.ndim
            src[i] = ndarr.sc.data[chunk_idx]
            rc = blosc2_cbuffer_sizes(src[i], &chunk_nbytes[i], &chunk_cbytes[i], &block_nbytes[i])
            if rc < 0:
                raise ValueError("miniexpr: error getting cbuffer sizes")
            if block_nbytes[i] <= 0:
                raise ValueError("miniexpr: invalid block size")
            if first_run:
                if i == 0:
                    q = ndarr.blockshape[ndim - 1]
                    p = ndarr.blockshape[ndim - 2]
                    # nchunk_ = chunks_in_row * chunk_row + chunk_col
                    # convert from chunk_idx to element idx chunk_i (row)
                    chunk_i = nchunk_ // ncols * out_chunk_nrows
                    chunk_startA = nchunkA + chunk_i // ndarr.chunkshape[ndim - 2] * Ancols
                    nchunkA = chunk_startA
                    # nblock_ = blocks_in_chunkrow * block_row + block_col
                    # convert from block_idx to element idx block_i (row)
                    block_i = nblock_ // block_ncols * out_block_nrows
                    block_startA = nblockA + block_i // p * Ablock_ncols
                else: # i = 1
                    r = ndarr.blockshape[ndim - 1]
                    # convert from chunk_idx to element idx chunk_j (col)
                    chunk_j = nchunk_ % ncols * out_chunk_ncols
                    chunk_startB = nchunkB + chunk_j // ndarr.chunkshape[ndim - 1]
                    nchunkB = chunk_startB
                    # convert from block_idx to element idx block_j (col)
                    block_j = nblock_ % block_ncols * out_block_ncols
                    block_startB = nblockB + block_j // r
                input_buffers[i] = malloc(block_nbytes[i])
            if input_buffers[i] == NULL:
                raise MemoryError("miniexpr: cannot allocate input block buffer")
            blocknitems[i] = block_nbytes[i] // <int> ndarr.sc.typesize

        first_run = False
        nblockA = block_startA
        nblockB = block_startB
        while True: # block loop
            startA = nblockA * blocknitems[0]
            startB = nblockB * blocknitems[1]
            rc = blosc2_getitem_ctx(dctx, src[0], chunk_cbytes[0], startA, blocknitems[0],
                                    input_buffers[0], block_nbytes[0])
            if rc < 0:
                raise ValueError("matmul: error decompressing the A chunk")
            rc = blosc2_getitem_ctx(dctx, src[1], chunk_cbytes[1], startB, blocknitems[1],
                                    input_buffers[1], block_nbytes[1])
            if rc < 0:
                raise ValueError("matmul: error decompressing the B chunk")
            batch = 0
            while batch < batches:
                batch_ = batch
                offsetA = 0
                offsetB = 0
                offset = 0
                for i in range(ndim - 2):
                    coord = batch_ // udata.el_strides[0][i]
                    batch_ = batch_ % udata.el_strides[0][i]
                    offsetA += coord * udata.el_strides[1][i]
                    offsetB += coord * udata.el_strides[2][i]
                    offset += coord * udata.el_strides[0][i]
                if typecode == 0:
                    if typesize == 4:
                        if selected_backend == B2_MATMUL_BACKEND_ACCELERATE:
                            rc = b2_gemm_accelerate_f32(
                                <float*>input_buffers[0] + offsetA,
                                <float*>input_buffers[1] + offsetB,
                                <float*>params_output + offset,
                                p,
                                q,
                                r,
                            )
                        elif selected_backend == B2_MATMUL_BACKEND_CBLAS:
                            rc = b2_gemm_cblas_f32(
                                <float*>input_buffers[0] + offsetA,
                                <float*>input_buffers[1] + offsetB,
                                <float*>params_output + offset,
                                p,
                                q,
                                r,
                            )
                        else:
                            rc = matmul_block_kernel[float](
                                <float*>input_buffers[0] + offsetA,
                                <float*>input_buffers[1] + offsetB,
                                <float*>params_output + offset,
                                p,
                                q,
                                r,
                            )
                    else:
                        if selected_backend == B2_MATMUL_BACKEND_ACCELERATE:
                            rc = b2_gemm_accelerate_f64(
                                <double*>input_buffers[0] + offsetA,
                                <double*>input_buffers[1] + offsetB,
                                <double*>params_output + offset,
                                p,
                                q,
                                r,
                            )
                        elif selected_backend == B2_MATMUL_BACKEND_CBLAS:
                            rc = b2_gemm_cblas_f64(
                                <double*>input_buffers[0] + offsetA,
                                <double*>input_buffers[1] + offsetB,
                                <double*>params_output + offset,
                                p,
                                q,
                                r,
                            )
                        else:
                            rc = matmul_block_kernel[double](
                                <double*>input_buffers[0] + offsetA,
                                <double*>input_buffers[1] + offsetB,
                                <double*>params_output + offset,
                                p,
                                q,
                                r,
                            )
                elif typecode == 1:
                    if typesize == 4:
                        rc = matmul_block_kernel[int32_t](<int32_t*>input_buffers[0] + offsetA, <int32_t*>input_buffers[1] + offsetB, <int32_t*>params_output + offset, p, q, r)
                    else:
                        rc = matmul_block_kernel[int64_t](<int64_t*>input_buffers[0] + offsetA, <int64_t*>input_buffers[1] + offsetB, <int64_t*>params_output + offset, p, q, r)
                else:
                    with gil:
                        raise ValueError("Unsupported dtype")
                batch += 1
            nblockA += 1
            nblockB += Bblock_ncols
            if (nblockA % Ablock_ncols == 0):
                break
        nchunkA += 1
        nchunkB += Bncols
        if (nchunkA % Ancols == 0):
            break


    blosc2_free_ctx(dctx)
    # Free resources
    for i in range(2):
        free(input_buffers[i])
    free(input_buffers)
    free(src)

    return 0

# Aux function for prefilter and postfilter udf
cdef int aux_udf(udf_udata *udata, int64_t nchunk, int32_t nblock,
                 c_bool is_postfilter, uint8_t *params_output, int32_t typesize):
    cdef int64_t chunk_ndim[B2ND_MAX_DIM]
    blosc2_unidim_to_multidim(udata.array.ndim, udata.chunks_in_array, nchunk, chunk_ndim)
    cdef int64_t block_ndim[B2ND_MAX_DIM]
    blosc2_unidim_to_multidim(udata.array.ndim, udata.blocks_in_chunk, nblock, block_ndim)
    cdef int64_t start_ndim[B2ND_MAX_DIM]
    for i in range(udata.array.ndim):
        start_ndim[i] = chunk_ndim[i] * udata.array.chunkshape[i] + block_ndim[i] * udata.array.blockshape[i]

    padding = False
    blockshape = []
    for i in range(udata.array.ndim):
        if start_ndim[i] + udata.array.blockshape[i] > udata.array.shape[i]:
            padding = True
            blockshape.append(udata.array.shape[i] - start_ndim[i])
            if blockshape[i] <= 0:
                # This block contains only padding, skip it
                return 0
        else:
            blockshape.append(udata.array.blockshape[i])
    cdef np.npy_intp dims[B2ND_MAX_DIM]
    for i in range(udata.array.ndim):
        dims[i] = blockshape[i]

    if padding:
        output = np.empty(blockshape, udata.array.dtype)
    else:
        output = np.PyArray_SimpleNewFromData(udata.array.ndim, dims, udata.output_cdtype, <void*>params_output)

    inputs_tuple = _ctypes.PyObj_FromPtr(udata.inputs_id)
    inputs_slice = []
    # Get slice of each operand
    l = []
    for i in range(udata.array.ndim):
        l.append(slice(start_ndim[i], start_ndim[i] + blockshape[i]))
    slices = tuple(l)
    for obj in inputs_tuple:
        if isinstance(obj, blosc2.NDArray | np.ndarray | blosc2.C2Array):
            inputs_slice.append(obj[slices])
        elif np.isscalar(obj):
            inputs_slice.append(obj)
        else:
            raise ValueError("Unsupported operand")

    # Call udf function
    func_id = udata.py_func.decode("utf-8")
    offset = tuple(start_ndim[i] for i in range(udata.array.ndim))
    if is_postfilter:
        blosc2.postfilter_funcs[func_id](tuple(inputs_slice), output, offset)
    else:
        blosc2.prefilter_funcs[func_id](tuple(inputs_slice), output, offset)

    cdef int64_t start[B2ND_MAX_DIM]
    cdef int64_t slice_shape[B2ND_MAX_DIM]
    cdef int64_t blockshape_int64[B2ND_MAX_DIM]
    cdef Py_buffer buf
    if padding:
        for i in range(udata.array.ndim):
            start[i] = 0
            slice_shape[i] = blockshape[i]
            blockshape_int64[i] = udata.array.blockshape[i]
        PyObject_GetBuffer(output, &buf, PyBUF_SIMPLE)
        rc = b2nd_copy_buffer2(udata.array.ndim, typesize,
                               buf.buf, slice_shape, start, slice_shape,
                               params_output, blockshape_int64, start)
        PyBuffer_Release(&buf)
        _check_rc(rc, "Could not copy the result into the buffer")

    return 0


cdef int miniexpr_prefilter(blosc2_prefilter_params *params):
    return aux_miniexpr(<me_udata *> params.user_data, params.nchunk, params.nblock, False,
                        params.output, params.output_typesize)

cdef int matmul_prefilter(blosc2_prefilter_params *params):
    cdef int typecode

    cdef mm_udata* udata = <mm_udata *> params.user_data
    cdef b2nd_array_t* out_arr = udata.array
    cdef char dtype_kind = out_arr.dtype[1]
    if dtype_kind == 'f':
        typecode = 0
    elif dtype_kind == 'i':
        typecode = 1
    else:
        raise ValueError("Unsupported dtype")
    return aux_matmul(udata, params.nchunk, params.nblock, params.output, params.output_typesize, typecode)

cdef int general_udf_prefilter(blosc2_prefilter_params *params):
    cdef udf_udata *udata = <udf_udata *> params.user_data
    return aux_udf(udata, params.nchunk, params.nblock, False, params.output, params.output_typesize)


cdef int general_udf_postfilter(blosc2_postfilter_params *params):
    cdef udf_udata *udata = <udf_udata *> params.user_data
    return aux_udf(udata, params.nchunk, params.nblock, True, params.output, params.typesize)


def nelem_from_inputs(inputs_tuple, nelem=None):
    for obj, dtype in inputs_tuple:
        if isinstance(obj, blosc2.SChunk):
            if nelem is not None and nelem != (obj.nbytes / obj.typesize):
                raise ValueError("operands must have same nelems")
            nelem = obj.nbytes / obj.typesize
        elif isinstance(obj, np.ndarray):
            if nelem is not None and nelem != obj.size:
                raise ValueError("operands must have same nelems")
            nelem = obj.size
    if nelem is None:
        raise ValueError("`nelem` must be set if none of the operands is a SChunk or a np.ndarray")
    return nelem

# prefilter
cdef int general_prefilter(blosc2_prefilter_params *params):
    cdef user_filters_udata *udata = <user_filters_udata *> params.user_data
    cdef int nd = 1
    cdef np.npy_intp dims = params.output_size // params.output_typesize


    input = np.PyArray_SimpleNewFromData(nd, &dims, udata.input_cdtype, <void*>params.input)
    output = np.PyArray_SimpleNewFromData(nd, &dims, udata.output_cdtype, <void*>params.output)
    offset = params.nchunk * udata.chunkshape + params.output_offset // params.output_typesize

    func_id = udata.py_func.decode("utf-8")
    blosc2.prefilter_funcs[func_id](input, output, offset)

    return 0


def remove_urlpath(path):
    blosc2_remove_urlpath(path)


# See https://github.com/dask/distributed/issues/3716#issuecomment-632913789
def encode_tuple(obj):
    if isinstance(obj, tuple):
        obj = ["__tuple__", *obj]
    return obj


def decode_tuple(obj):
    if obj[0] == "__tuple__":
        obj = tuple(obj[1:])
    return obj


cdef class vlmeta:
    cdef blosc2_schunk* schunk
    def __init__(self, schunk):
        self.schunk = <blosc2_schunk*> <uintptr_t>schunk

    def set_vlmeta(self, name, content, **cparams):
        cdef blosc2_cparams ccparams
        create_cparams_from_kwargs(&ccparams, cparams)
        name = name.encode("utf-8") if isinstance(name, str) else name
        content = content.encode("utf-8") if isinstance(content, str) else content
        cdef uint32_t len_content = <uint32_t> len(content)
        rc = blosc2_vlmeta_exists(self.schunk, name)
        if rc >= 0:
            rc = blosc2_vlmeta_update(self.schunk, name, <uint8_t*> content, len_content, &ccparams)
        else:
            rc = blosc2_vlmeta_add(self.schunk, name,  <uint8_t*> content, len_content, &ccparams)

        if rc < 0:
            raise RuntimeError

    def get_vlmeta(self, name):
        name = name.encode("utf-8") if isinstance(name, str) else name
        rc = blosc2_vlmeta_exists(self.schunk, name)
        cdef uint8_t* content
        cdef int32_t content_len
        if rc < 0:
            raise KeyError
        if rc >= 0:
            rc = blosc2_vlmeta_get(self.schunk, name, &content, &content_len)
        if rc < 0:
            raise RuntimeError
        return content[:content_len]

    def del_vlmeta(self, name):
        name = name.encode("utf-8") if isinstance(name, str) else name
        rc = blosc2_vlmeta_delete(self.schunk, name)
        if rc < 0:
            raise RuntimeError("Could not delete the vlmeta")

    def nvlmetalayers(self):
        return self.schunk.nvlmetalayers

    def get_names(self):
        cdef char** names = <char **> malloc(self.schunk.nvlmetalayers * sizeof (char *))
        rc = blosc2_vlmeta_get_names(self.schunk, names)
        if rc != self.schunk.nvlmetalayers:
            raise RuntimeError
        res = [names[i].decode("utf-8") for i in range(rc)]
        return res

    def to_dict(self):
        cdef char** names = <char **> malloc(self.schunk.nvlmetalayers * sizeof (char*))
        rc = blosc2_vlmeta_get_names(self.schunk, names)
        if rc != self.schunk.nvlmetalayers:
            raise RuntimeError
        res = {}
        for i in range(rc):
            res[names[i]] = unpackb(self.get_vlmeta(names[i]), list_hook=decode_tuple)
        return res


def meta__contains__(self, name):
    cdef blosc2_schunk *schunk = <blosc2_schunk *><uintptr_t> self.c_schunk
    name = name.encode("utf-8") if isinstance(name, str) else name
    n = blosc2_meta_exists(schunk, name)
    return False if n < 0 else True

def meta__getitem__(self, name):
    cdef blosc2_schunk *schunk = <blosc2_schunk *><uintptr_t> self.c_schunk
    name = name.encode("utf-8") if isinstance(name, str) else name
    cdef uint8_t *content
    cdef int32_t content_len
    n = blosc2_meta_get(schunk, name, &content, &content_len)
    res = PyBytes_FromStringAndSize(<char *> content, content_len)
    free(content)

    return res

def meta__setitem__(self, name, content):
    cdef blosc2_schunk *schunk = <blosc2_schunk *><uintptr_t> self.c_schunk
    name = name.encode("utf-8") if isinstance(name, str) else name
    old_content = meta__getitem__(self, name)
    if len(old_content) != len(content):
        raise ValueError("The length of the content in a metalayer cannot change.")
    blosc2_meta_update(schunk, name, content, len(content))

def meta__len__(self):
    cdef blosc2_schunk *schunk = <blosc2_schunk *><uintptr_t> self.c_schunk
    return schunk.nmetalayers

def meta_keys(self):
    cdef blosc2_schunk *schunk = <blosc2_schunk *><uintptr_t> self.c_schunk
    keys = []
    for i in range(meta__len__(self)):
        name = schunk.metalayers[i].name.decode("utf-8")
        keys.append(name)
    return keys


def open(urlpath, mode, offset, **kwargs):
    urlpath_ = urlpath.encode("utf-8") if isinstance(urlpath, str) else urlpath
    cdef blosc2_schunk* schunk
    cdef blosc2_stdio_mmap* mmap_file
    cdef blosc2_io* io

    mmap_mode = kwargs.get("mmap_mode")
    locking = kwargs.get("locking")
    if locking and mmap_mode is not None:
        raise ValueError("locking is not supported together with mmap_mode")
    if mmap_mode is not None:
        if mmap_mode == "w+":
            raise ValueError("w+ mmap_mode cannot be used to open an existing file")
        else:
            mode = mode_from_mmap_mode(mmap_mode)

    initial_mapping_size = kwargs.get("initial_mapping_size")
    if initial_mapping_size is not None:
        if mmap_mode is None:
            raise ValueError("initial_mapping_size can only be used with mmap_mode")

        if mmap_mode == "r":
            raise ValueError("initial_mapping_size can only be used with writing modes (r+, c)")

    if mmap_mode is None:
        if locking:
            schunk = blosc2_schunk_open_offset_udio(urlpath_, offset, &_locking_io)
        else:
            schunk = blosc2_schunk_open_offset(urlpath_, offset)
    else:
        mmap_file = <blosc2_stdio_mmap *>malloc(sizeof(BLOSC2_STDIO_MMAP_DEFAULTS))
        memcpy(mmap_file, &BLOSC2_STDIO_MMAP_DEFAULTS, sizeof(BLOSC2_STDIO_MMAP_DEFAULTS))

        mmap_mode_ = mmap_mode.encode("utf-8")
        mmap_file.mode = mmap_mode_
        mmap_file.needs_free = True
        if initial_mapping_size is not None:
            mmap_file.initial_mapping_size = initial_mapping_size

        io = <blosc2_io *>malloc(sizeof(blosc2_io))
        io.id = BLOSC2_IO_FILESYSTEM_MMAP
        io.params = mmap_file
        schunk = blosc2_schunk_open_offset_udio(urlpath_, offset, io)

    if schunk == NULL:
        if mmap_mode is not None:
            free(io)
        raise RuntimeError(f'blosc2_schunk_open_offset({urlpath!r}, {offset!r}) returned NULL')

    is_ndarray = schunk_is_ndarray(schunk)

    cdef b2nd_array_t *array
    if is_ndarray:
        _check_rc(b2nd_from_schunk(schunk, &array),
                  "Could not create array from schunk")

    kwargs["urlpath"] = urlpath
    kwargs["contiguous"] = schunk.storage.contiguous
    if mode != "w" and kwargs is not None:
        check_schunk_params(schunk, kwargs)
    cparams = kwargs.get("cparams")
    # nthreads is not stored in the frame; apply the live global when the caller
    # did not supply an explicit cparams — symmetric with the DParams default below.
    dparams = kwargs.get("dparams", blosc2.DParams())

    if is_ndarray:
        res = blosc2.NDArray(_schunk=PyCapsule_New(array.sc, <char *> "blosc2_schunk*", NULL),
                             _array=PyCapsule_New(array, <char *> "b2nd_array_t*", NULL), mode=mode)
        if cparams is not None:
            res._schunk.cparams = cparams if isinstance(cparams, blosc2.CParams) else blosc2.CParams(**cparams)
        else:
            res._schunk.cparams = dataclasses.replace(
                res._schunk.cparams, nthreads=(1 if blosc2.IS_WASM else blosc2.nthreads)
            )
        if dparams is not None:
            res._schunk.dparams = dparams if isinstance(dparams, blosc2.DParams) else blosc2.DParams(**dparams)
        res._schunk.mode = mode
    else:
        res = blosc2.SChunk(_schunk=PyCapsule_New(schunk, <char *> "blosc2_schunk*", NULL),
                            mode=mode, **kwargs)
        if cparams is not None:
            res.cparams = cparams if isinstance(cparams, blosc2.CParams) else blosc2.CParams(**cparams)
        else:
            res.cparams = dataclasses.replace(res.cparams, nthreads=(1 if blosc2.IS_WASM else blosc2.nthreads))
        if dparams is not None:
            res.dparams = dparams if isinstance(dparams, blosc2.DParams) else blosc2.DParams(**dparams)

    return res


def check_access_mode(urlpath, mode):
    if urlpath is not None and mode == "r":
        raise ValueError("Cannot do this action with reading mode")


def mode_from_mmap_mode(mmap_mode):
    # We ignore the user-supplied mode with mmap files and use a fixed mapping instead
    if mmap_mode == "r":
        mode = "r"
    elif mmap_mode == "r+":
        mode = "a"
    elif mmap_mode == "w+":
        mode = "w"
    elif mmap_mode == "c":
        # In terms of (internal) blosc, it is allowed to modify the file contents
        # The actual file is opened in read-only mode
        mode = "a"
    else:
        raise ValueError(f"Invalid mmap_mode: {mmap_mode}")

    return mode


cdef check_schunk_params(blosc2_schunk* schunk, kwargs):
    cparams = kwargs.get("cparams", None)
    if cparams is not None:
        blocksize = kwargs.get("blocksize", schunk.blocksize)
        if blocksize not in [0, schunk.blocksize]:
            raise ValueError("Cannot change blocksize with this mode")
        typesize = kwargs.get("typesize", schunk.typesize)
        if typesize != schunk.typesize:
            raise ValueError("Cannot change typesize with this mode")


cdef schunk_is_ndarray(blosc2_schunk* schunk):
    meta = "b2nd"
    meta = meta.encode("utf-8") if isinstance(meta, str) else meta
    return blosc2_meta_exists(schunk, meta) >= 0


def schunk_from_cframe(cframe, copy=False):
    cdef Py_buffer buf
    PyObject_GetBuffer(cframe, &buf, PyBUF_SIMPLE)
    cdef blosc2_schunk *schunk_ = blosc2_schunk_from_buffer(<uint8_t *>buf.buf, buf.len, copy)
    if schunk_ == NULL:
        raise RuntimeError("Could not get the schunk from the cframe")
    schunk = blosc2.SChunk(_schunk=PyCapsule_New(schunk_, <char *> "blosc2_schunk*", NULL))
    PyBuffer_Release(&buf)
    if not copy:
        schunk._avoid_cframe_free(True)
    return schunk


cdef int general_encoder(const uint8_t* input_buffer, int32_t input_len,
                        uint8_t* output_buffer, int32_t output_len,
                        uint8_t meta,
                        blosc2_cparams* cparams, const void* chunk):
    cdef int nd = 1
    cdef np.npy_intp input_dims = input_len
    cdef np.npy_intp output_dims = output_len
    input = np.PyArray_SimpleNewFromData(nd, &input_dims, np.NPY_UINT8, <void*>input_buffer)
    output = np.PyArray_SimpleNewFromData(nd, &output_dims, np.NPY_UINT8, <void*>output_buffer)

    cdef blosc2_schunk *sc = <blosc2_schunk *> cparams.schunk
    if sc != NULL:
        schunk = blosc2.SChunk(_schunk=PyCapsule_New(sc, <char *> "blosc2_schunk*", NULL), _is_view=True)
    else:
        raise RuntimeError("Cannot apply user codec without an SChunk")
    rc = blosc2.ucodecs_registry[cparams.compcode][1](input, output, meta, schunk)
    if rc is None:
        raise RuntimeError("encoder must return the number of compressed bytes")

    return rc


cdef int general_decoder(const uint8_t* input_buffer, int32_t input_len,
                        uint8_t* output_buffer, int32_t output_len,
                        uint8_t meta,
                        blosc2_dparams *dparams, const void* chunk):
    cdef int nd = 1
    cdef np.npy_intp input_dims = input_len
    cdef np.npy_intp output_dims = output_len
    input = np.PyArray_SimpleNewFromData(nd, &input_dims, np.NPY_UINT8, <void*>input_buffer)
    output = np.PyArray_SimpleNewFromData(nd, &output_dims, np.NPY_UINT8, <void*>output_buffer)

    cdef blosc2_schunk *sc = <blosc2_schunk *> dparams.schunk
    if sc != NULL:
        schunk = blosc2.SChunk(_schunk=PyCapsule_New(sc, <char *> "blosc2_schunk*", NULL), _is_view=True)
    else:
        raise RuntimeError("Cannot apply user codec without an SChunk")

    rc = blosc2.ucodecs_registry[sc.compcode][2](input, output, meta, schunk)
    if rc is None:
        raise RuntimeError("decoder must return the number of decompressed bytes")

    return rc


def register_codec(codec_name, id, encoder=None, decoder=None, version=1):
    if id < BLOSC2_USER_REGISTERED_CODECS_START or id > BLOSC2_USER_REGISTERED_CODECS_STOP:
        raise ValueError("`id` must be between ", BLOSC2_USER_REGISTERED_CODECS_START,
                         " and ", BLOSC2_USER_REGISTERED_CODECS_STOP)

    if (encoder is None and decoder is not None) or (encoder is not None and decoder is None):
        raise ValueError("both encoder and decoder must be given, or none")

    cdef blosc2_codec codec
    codec.compcode = id
    codec.version = version
    codec.complib = id
    codec_name_ = codec_name.encode() if isinstance(codec_name, str) else codec_name
    codec.compname = <char *> malloc(strlen(codec_name_) + 1)
    strcpy(codec.compname, codec_name_)
    if encoder is None:
        codec.encoder = NULL
    else:
        codec.encoder = <blosc2_codec_encoder_cb> general_encoder
    if decoder is None:
        codec.decoder = NULL
    else:
        codec.decoder = <blosc2_codec_decoder_cb> general_decoder

    rc = blosc2_register_codec(&codec)
    if rc < 0:
        raise RuntimeError("Error while registering codec")

    if encoder and decoder:
        blosc2.ucodecs_registry[id] = (codec_name, encoder, decoder)


cdef int general_forward(const uint8_t* input_buffer, uint8_t* output_buffer, int32_t size,
                        uint8_t meta, blosc2_cparams* cparams, uint8_t id):
    cdef int nd = 1
    cdef np.npy_intp dims = size
    input = np.PyArray_SimpleNewFromData(nd, &dims, np.NPY_UINT8, <void*>input_buffer)
    output = np.PyArray_SimpleNewFromData(nd, &dims, np.NPY_UINT8, <void*>output_buffer)

    cdef blosc2_schunk *sc = <blosc2_schunk *> cparams.schunk
    if sc != NULL:
        schunk = blosc2.SChunk(_schunk=PyCapsule_New(sc, <char *> "blosc2_schunk*", NULL), _is_view=True)
    else:
        raise RuntimeError("Cannot apply user codec without an SChunk")
    blosc2.ufilters_registry[id][0](input, output, meta, schunk)

    return BLOSC2_ERROR_SUCCESS


cdef int general_backward(const uint8_t* input_buffer, uint8_t* output_buffer, int32_t size,
                            uint8_t meta, blosc2_dparams* dparams, uint8_t id):
    cdef int nd = 1
    cdef np.npy_intp dims = size
    input = np.PyArray_SimpleNewFromData(nd, &dims, np.NPY_UINT8, <void*>input_buffer)
    output = np.PyArray_SimpleNewFromData(nd, &dims, np.NPY_UINT8, <void*>output_buffer)

    cdef blosc2_schunk *sc = <blosc2_schunk *> dparams.schunk
    if sc != NULL:
        schunk = blosc2.SChunk(_schunk=PyCapsule_New(sc, <char *> "blosc2_schunk*", NULL), _is_view=True)
    else:
        raise RuntimeError("Cannot apply user filter without an SChunk")

    blosc2.ufilters_registry[id][1](input, output, meta, schunk)

    return BLOSC2_ERROR_SUCCESS


def register_filter(id, forward, backward, filter_name):
    if id < BLOSC2_USER_REGISTERED_FILTERS_START or id > BLOSC2_USER_REGISTERED_FILTERS_STOP:
        raise ValueError("`id` must be between ", BLOSC2_USER_REGISTERED_FILTERS_START,
                         " and ", BLOSC2_USER_REGISTERED_FILTERS_STOP)
    if (forward is None and backward is not None) or (forward is not None and backward is None):
        raise ValueError("both encoder and decoder must be given, or none")

    cdef blosc2_filter filter
    filter.id = id
    if forward is None:
        filter.forward = NULL
    else:
        filter.forward = <blosc2_filter_forward_cb> general_forward
    if backward is None:
        filter.backward = NULL
    else:
        filter.backward = <blosc2_filter_backward_cb> general_backward
    if filter_name is None and not forward and not backward:
        raise ValueError("You need to pass the filter name or the forward and backward functions")
    if filter_name:
        filter_name_ = filter_name.encode() if isinstance(filter_name, str) else filter_name
        filter.name = <char *> malloc(strlen(filter_name_) + 1)
        strcpy(filter.name, filter_name_)

    rc = blosc2_register_filter(&filter)
    if rc < 0:
        raise RuntimeError("Error while registering filter")
    if forward and backward:
        blosc2.ufilters_registry[id] = (forward, backward)

cdef _check_rc(rc, message):
    if rc < 0:
        raise RuntimeError(message)


cdef class slice_flatter:
    cdef long long ndim
    cdef int done
    cdef long long[:] shape
    cdef long long[:] start
    cdef long long[:] stop
    cdef long long[:] strides
    cdef long long[:] indices
    cdef long long current_slice_start
    cdef long long current_slice_end
    cdef long long current_flat_idx  # Track the current flat index

    def __cinit__(self, long long[:] start not None, long long[:] stop not None, long long[:] strides not None):
        self.ndim = start.shape[0]
        self.done = 0
        self.start = start
        self.stop = stop
        self.strides = strides
        self.current_slice_start = -1
        self.current_slice_end = -1
        shape = tuple(stop[i] - start[i] for i in range(self.ndim))
        self.shape = np.array(shape, dtype=np.int64)
        self.indices = np.zeros(self.ndim, dtype=np.int64)
        # Initialize the flat index
        self.current_flat_idx = 0
        for j in range(self.ndim):
            self.current_flat_idx += self.start[j] * self.strides[j]

    def __iter__(self):
        return self

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __next__(self):
        cdef long long j, next_flat_idx
        cdef int extended_slice = 0

        # Check if we're done
        if self.done:
            if self.current_slice_start != -1:
                result = slice(self.current_slice_start, self.current_slice_end + 1)
                self.current_slice_start = -1
                return result
            raise StopIteration

        # Initialize first slice point if needed
        if self.current_slice_start == -1:
            next_flat_idx = 0
            for j in range(self.ndim):
                next_flat_idx += (self.start[j] + self.indices[j]) * self.strides[j]
            self.current_slice_start = next_flat_idx
            self.current_slice_end = next_flat_idx
            self.current_flat_idx = next_flat_idx
            self.incr_indices()

            # If we're done after the first element, return it
            if self.done:
                result = slice(self.current_slice_start, self.current_slice_end + 1)
                self.current_slice_start = -1
                return result

        # Extend slice as long as indices remain contiguous
        while not self.done:
            # Calculate next flat index
            next_flat_idx = 0
            for j in range(self.ndim):
                next_flat_idx += (self.start[j] + self.indices[j]) * self.strides[j]

            # If indices are contiguous, extend current slice
            if next_flat_idx == self.current_slice_end + 1:
                self.current_slice_end = next_flat_idx
                self.current_flat_idx = next_flat_idx
                self.incr_indices()
                extended_slice = 1
            else:
                # Non-contiguous index found, return current slice
                result = slice(self.current_slice_start, self.current_slice_end + 1)
                self.current_slice_start = next_flat_idx
                self.current_slice_end = next_flat_idx
                self.current_flat_idx = next_flat_idx
                self.incr_indices()
                return result

        # If we've reached the end after extending the slice
        if extended_slice:
            result = slice(self.current_slice_start, self.current_slice_end + 1)
            self.current_slice_start = -1
            return result

        # Should never reach here
        raise StopIteration

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void incr_indices(self) nogil:
        cdef long long i
        for i in range(self.ndim - 1, -1, -1):
            self.indices[i] += 1
            if self.indices[i] < self.shape[i]:
                break
            self.indices[i] = 0
            if i == 0:
                self.done = 1


cdef class NDArray:
    cdef b2nd_array_t* array

    def __init__(self, array, base=None):
        self._dtype = None
        self.array = <b2nd_array_t *> PyCapsule_GetPointer(array, <char *> "b2nd_array_t*")
        self.base = base # add reference to base if NDArray is a view

    @property
    def c_array(self):
        return <uintptr_t> self.array

    @property
    def shape(self) -> tuple[int]:
        return tuple([self.array.shape[i] for i in range(self.array.ndim)])

    @property
    def ext_shape(self):
        return tuple([self.array.extshape[i] for i in range(self.array.ndim)])

    @property
    def chunks(self):
        return tuple([self.array.chunkshape[i] for i in range(self.array.ndim)])

    @property
    def ext_chunks(self):
        return tuple([self.array.extchunkshape[i] for i in range(self.array.ndim)])

    @property
    def blocks(self):
        return tuple([self.array.blockshape[i] for i in range(self.array.ndim)])

    @property
    def ndim(self):
        return self.array.ndim

    @property
    def size(self):
        return self.array.nitems

    @property
    def chunksize(self):
        return self.array.chunknitems * self.array.sc.typesize

    @property
    def dtype(self):
        if self._dtype is not None:
            return self._dtype

        # Not in cache yet
        if self.array.dtype == NULL:
            return np.dtype(f"S{self.array.sc.typesize}")
        if self.array.dtype_format != B2ND_DEFAULT_DTYPE_FORMAT:
            raise ValueError("Only NumPy dtypes are supported")
        cdef char *bytes_dtype = self.array.dtype
        str_dtype = bytes_dtype.decode("utf-8")
        try:
            dtype = np.dtype(str_dtype)
        except (ValueError, TypeError):
            dtype = np.dtype(ast.literal_eval(str_dtype))
        self._dtype = dtype
        return dtype

    def get_slice_numpy(self, arr, key):
        start, stop = key

        cdef int64_t[B2ND_MAX_DIM] start_, stop_
        cdef int64_t[B2ND_MAX_DIM] buffershape_
        for i in range(self.ndim):
            start_[i] = start[i]
            stop_[i] = stop[i]
            buffershape_[i] = stop_[i] - start_[i]

        cdef Py_buffer view
        PyObject_GetBuffer(arr, &view, PyBUF_SIMPLE)
        _check_rc(b2nd_get_slice_cbuffer(self.array, start_, stop_,
                                         <void *> view.buf, buffershape_, view.len),
                  "Error while getting the buffer")
        PyBuffer_Release(&view)

        return arr

    def get_1d_span_numpy(self, arr, int64_t nchunk, int32_t start, int32_t nitems):
        if self.ndim != 1:
            raise ValueError("get_1d_span_numpy is only supported for 1-D arrays")
        if nchunk < 0 or nchunk >= self.array.sc.nchunks:
            raise IndexError("chunk index out of range")
        if start < 0 or nitems < 0:
            raise ValueError("start and nitems must be >= 0")
        if start + nitems > self.array.chunknitems:
            raise ValueError("requested span exceeds chunk size")

        cdef uint8_t *chunk = NULL
        cdef c_bool needs_free
        cdef int32_t chunk_nbytes
        cdef int32_t chunk_cbytes
        cdef int32_t block_nbytes
        cdef blosc2_context *dctx = self.array.sc.dctx
        cdef Py_buffer view
        cdef int rc
        cdef int32_t lazychunk_cbytes
        cdef c_bool owns_dctx = False

        lazychunk_cbytes = blosc2_schunk_get_lazychunk(self.array.sc, nchunk, &chunk, &needs_free)
        if lazychunk_cbytes < 0:
            raise RuntimeError("Error while getting the lazy chunk")

        rc = blosc2_cbuffer_sizes(chunk, &chunk_nbytes, &chunk_cbytes, &block_nbytes)
        if rc < 0:
            if needs_free:
                free(chunk)
            raise RuntimeError("Error while getting compressed buffer sizes")
        if start + nitems > chunk_nbytes // self.array.sc.typesize:
            if needs_free:
                free(chunk)
            raise ValueError("requested span exceeds decoded chunk size")

        PyObject_GetBuffer(arr, &view, PyBUF_SIMPLE)
        if view.len < nitems * self.array.sc.typesize:
            PyBuffer_Release(&view)
            if needs_free:
                free(chunk)
            raise ValueError("destination buffer is smaller than the requested decoded span")

        if dctx == NULL:
            dctx = blosc2_create_dctx(BLOSC2_DPARAMS_DEFAULTS)
            owns_dctx = True
        if dctx == NULL:
            PyBuffer_Release(&view)
            if needs_free:
                free(chunk)
            raise RuntimeError("Could not create decompression context")
        # For lazy chunks, blosc2_cbuffer_sizes() only reports the header cbytes.
        # blosc2_getitem_ctx() needs the full lazy chunk size returned by
        # blosc2_schunk_get_lazychunk().
        rc = blosc2_getitem_ctx(dctx, chunk, lazychunk_cbytes, start, nitems, view.buf, view.len)
        if owns_dctx:
            blosc2_free_ctx(dctx)
        PyBuffer_Release(&view)
        if needs_free:
            free(chunk)
        if rc < 0:
            raise RuntimeError("Error while decoding the requested span")

        return arr

    def get_sparse_numpy(self, arr, coords):
        cdef np.ndarray[np.int64_t, ndim=1, mode="c"] coords_ = np.ascontiguousarray(coords, dtype=np.int64)
        cdef Py_buffer view
        cdef int64_t ncoords = coords_.shape[0]
        cdef int rc

        PyObject_GetBuffer(arr, &view, PyBUF_SIMPLE)
        if view.len < ncoords * self.array.sc.typesize:
            PyBuffer_Release(&view)
            raise ValueError("destination buffer is smaller than the requested sparse selection")

        rc = b2nd_get_sparse_cbuffer(self.array, ncoords, <const int64_t *> coords_.data,
                                     <void *> view.buf, view.len)
        PyBuffer_Release(&view)
        _check_rc(rc, "Error while getting the sparse selection")
        return arr

    def get_oindex_numpy(self, arr, key):
        """
        Orthogonal indexing. Key is a tuple of lists of integer indices.
        """
        if len(key) != self.array.ndim:
            raise ValueError(f"Key must have {self.array.ndim} dimensions, got {len(key)}.")
        cdef int64_t[B2ND_MAX_DIM] buffershape_
        cdef int64_t** key_
        cdef int64_t buffersize_ = self.array.sc.typesize
        cdef int64_t[B2ND_MAX_DIM] sel_size

        key_ = <int64_t**> malloc(len(key) * sizeof(int64_t *))

        for i in range(self.array.ndim):
            buffershape_[i] = len(key[i])
            buffersize_ *= buffershape_[i]
            sel_size[i] = len(key[i])
            key_[i] = <int64_t *> malloc(sel_size[i] * sizeof(int64_t))
            for j in range(len(key[i])):
                key_[i][j] = key[i][j]

        cdef Py_buffer buf
        PyObject_GetBuffer(arr, &buf, PyBUF_SIMPLE)

        _check_rc(b2nd_get_orthogonal_selection(self.array, key_, sel_size, buf.buf,
                                                buffershape_, buffersize_), "Error while getting orthogonal selection")
        PyBuffer_Release(&buf)
        for i in range(len(key)):
            free(key_[i])  # Free the allocated memory for each key
        free(key_)
        return arr

    def set_oindex_numpy(self, key, arr):
        """
        Orthogonal indexing. Set elements of self with arr using key.
        """
        if len(key) != self.array.ndim:
            raise ValueError(f"Key must have {self.array.ndim} dimensions, got {len(key)}.")
        cdef int64_t[B2ND_MAX_DIM] buffershape_
        cdef int64_t** key_
        cdef int64_t buffersize_ = self.array.sc.typesize
        cdef int64_t[B2ND_MAX_DIM] sel_size

        key_ = <int64_t**> malloc(len(key) * sizeof(int64_t *))

        for i in range(self.array.ndim):
            buffershape_[i] = len(key[i])
            buffersize_ *= buffershape_[i]
            sel_size[i] = len(key[i])
            key_[i] = <int64_t *> malloc(sel_size[i] * sizeof(int64_t))
            for j in range(len(key[i])):
                key_[i][j] = key[i][j]

        cdef Py_buffer buf
        PyObject_GetBuffer(arr, &buf, PyBUF_SIMPLE)

        _check_rc(b2nd_set_orthogonal_selection(self.array, key_, sel_size, buf.buf,
                                                buffershape_, buffersize_), "Error while getting orthogonal selection")
        PyBuffer_Release(&buf)
        for i in range(len(key)):
            free(key_[i])  # Free the allocated memory for each key
        free(key_)
        return arr


    def get_slice(self, key, mask, **kwargs):
        start, stop = key
        shape = tuple(sp - st for sp, st in zip(stop, start))
        chunks = kwargs.pop("chunks", None)
        blocks = kwargs.pop("blocks", None)
        if blocks and len(shape) != len(blocks):
            for i in range(len(shape)):
                if shape[i] == 1:
                    blocks.insert(i, 1)
        if chunks and len(shape) != len(chunks):
            for i in range(len(shape)):
                if shape[i] == 1:
                    chunks.insert(i, 1)
        chunks, blocks = blosc2.compute_chunks_blocks(shape, chunks, blocks, self.dtype)

        # shape will be overwritten by get_slice
        cdef b2nd_context_t *ctx = create_b2nd_context(shape, chunks, blocks,
                                                       self.dtype, kwargs)
        if ctx == NULL:
            raise RuntimeError("Error while creating the context")
        ndim = self.ndim
        cdef int64_t[B2ND_MAX_DIM] start_, stop_
        for i in range(ndim):
            start_[i] = start[i]
            stop_[i] = stop[i]

        cdef b2nd_array_t *array
        _check_rc(b2nd_get_slice(ctx, &array, self.array, start_, stop_),
                  "Error while getting the slice")
        _check_rc(b2nd_free_ctx(ctx), "Error while freeing the context")

        cdef c_bool mask_[B2ND_MAX_DIM]
        for i in range(ndim):
            mask_[i] = mask[i]
        _check_rc(b2nd_squeeze_index(array, &array, mask_), "Error while squeezing sliced array")
        ndarray = blosc2.NDArray(_schunk=PyCapsule_New(array.sc, <char *> "blosc2_schunk*", NULL),
                                 _array=PyCapsule_New(array, <char *> "b2nd_array_t*", NULL))


        return ndarray

    def set_slice(self, key, ndarray):
        ndim = self.ndim
        start, stop = key
        cdef Py_buffer buf
        PyObject_GetBuffer(ndarray, &buf, PyBUF_SIMPLE)

        cdef int64_t[B2ND_MAX_DIM] buffershape_, start_, stop_
        for i in range(ndim):
            start_[i] = start[i]
            stop_[i] = stop[i]
            buffershape_[i] = stop[i] - start[i]

        _check_rc(b2nd_set_slice_cbuffer(buf.buf, buffershape_, buf.len, start_, stop_, self.array),
                  "Error while setting the slice")
        PyBuffer_Release(&buf)

        return self

    def tobytes(self):
        buffersize = self.size * self.array.sc.typesize
        buffer = bytes(buffersize)
        _check_rc(b2nd_to_cbuffer(self.array, <void *> <char *> buffer, buffersize),
                  "Error while filling the buffer")

        return buffer

    def to_cframe(self):
        cdef c_bool needs_free
        cdef uint8_t *cframe
        cdef int64_t cframe_len;
        cdef int rc;
        rc = b2nd_to_cframe(self.array, &cframe, &cframe_len, &needs_free)
        if rc < 0:
            raise RuntimeError("Error while getting the cframe")
        out = PyBytes_FromStringAndSize(<char*>cframe, cframe_len)
        if needs_free:
            free(cframe)

        return out

    def copy(self, dtype, **kwargs):
        chunks = kwargs.pop("chunks", self.chunks)
        blocks = kwargs.pop("blocks", self.blocks)
        kwargs["contiguous"] =  kwargs.get("contiguous", self.array.sc.storage.contiguous)

        chunks, blocks = blosc2.compute_chunks_blocks(self.shape, chunks, blocks, dtype, **kwargs)
        cdef b2nd_context_t *ctx = create_b2nd_context(self.shape, chunks, blocks, dtype, kwargs)
        if ctx == NULL:
            raise RuntimeError("Error while creating the context")

        cdef b2nd_array_t *array
        _check_rc(b2nd_copy(ctx, self.array, &array),
                  "Error while copying the array")

        ndarray = blosc2.NDArray(_schunk=PyCapsule_New(array.sc, <char *> "blosc2_schunk*", NULL),
                                 _array=PyCapsule_New(array, <char *> "b2nd_array_t*", NULL))
        _check_rc(b2nd_free_ctx(ctx), "Error while freeing the context")

        return ndarray

    def resize(self, new_shape):
        cdef int64_t new_shape_[B2ND_MAX_DIM]
        for i, s in enumerate(new_shape):
            new_shape_[i] = s
        _check_rc(b2nd_resize(self.array, new_shape_, NULL),
                  "Error while resizing the array")

    def refresh(self):
        cdef int rc = b2nd_refresh(self.array)
        _check_rc(rc, "Error while refreshing the array")
        return rc == 1

    def as_ffi_ptr(self):
        return PyCapsule_New(self.array, <char *> "b2nd_array_t*", NULL)

    cdef udf_udata *_fill_udf_udata(self, func_id, inputs):
        cdef udf_udata *udata = <udf_udata *> malloc(sizeof(udf_udata))
        udata.py_func = <char *> malloc(strlen(func_id) + 1)
        strcpy(udata.py_func, func_id)
        udata.inputs_id = id(inputs)
        udata.output_cdtype = np.dtype(self.dtype).num
        udata.array = self.array
        # Save these in udf_udata to avoid computing them for each block
        for i in range(self.array.ndim):
            udata.chunks_in_array[i] = udata.array.extshape[i] // udata.array.chunkshape[i]
            udata.blocks_in_chunk[i] = udata.array.extchunkshape[i] // udata.array.blockshape[i]

        return udata

    cdef me_udata *_fill_me_udata(self, inputs, fp_accuracy, aux_reduc, jit=None):
        cdef me_udata *udata = <me_udata *> calloc(1, sizeof(me_udata))
        cdef me_eval_params* eval_params
        cdef b2nd_array_t** inputs_
        cdef me_input_cache_s* input_chunk_caches
        cdef void* aux_reduc_ptr = NULL
        cdef int i
        if aux_reduc is not None:
            if not isinstance(aux_reduc, np.ndarray):
                raise TypeError("aux_reduc must be a NumPy array")
            aux_reduc_ptr = <void *> np.PyArray_DATA(<np.ndarray> aux_reduc)
        operands = list(inputs.values())
        ninputs = len(operands)
        if udata == NULL:
            raise MemoryError("Cannot allocate miniexpr user data")
        inputs_ = NULL
        if ninputs > 0:
            inputs_ = <b2nd_array_t**> malloc(ninputs * sizeof(b2nd_array_t*))
            if inputs_ == NULL:
                free(udata)
                raise MemoryError("Cannot allocate miniexpr input table")
        for i, operand in enumerate(operands):
            inputs_[i] = <b2nd_array_t*><uintptr_t>operand.c_array
        udata.inputs = inputs_
        udata.ninputs = ninputs
        input_chunk_caches = NULL
        if ninputs > 0:
            input_chunk_caches = <me_input_cache_s*> calloc(ninputs, sizeof(me_input_cache_s))
            if input_chunk_caches == NULL:
                free(inputs_)
                free(udata)
                raise MemoryError("Cannot allocate miniexpr chunk caches")
            for i in range(ninputs):
                input_chunk_caches[i].nchunk = -1
                input_chunk_caches[i].state = ME_CACHE_EMPTY
                input_chunk_caches[i].state_lock = PyThread_allocate_lock()
                if input_chunk_caches[i].state_lock == NULL:
                    while i > 0:
                        i -= 1
                        if input_chunk_caches[i].state_lock != NULL:
                            PyThread_free_lock(input_chunk_caches[i].state_lock)
                        if input_chunk_caches[i].ready_lock != NULL:
                            PyThread_free_lock(input_chunk_caches[i].ready_lock)
                    free(input_chunk_caches)
                    free(inputs_)
                    free(udata)
                    raise MemoryError("Cannot allocate miniexpr chunk cache state lock")
                input_chunk_caches[i].ready_lock = PyThread_allocate_lock()
                if input_chunk_caches[i].ready_lock == NULL:
                    PyThread_free_lock(input_chunk_caches[i].state_lock)
                    input_chunk_caches[i].state_lock = NULL
                    while i > 0:
                        i -= 1
                        if input_chunk_caches[i].state_lock != NULL:
                            PyThread_free_lock(input_chunk_caches[i].state_lock)
                        if input_chunk_caches[i].ready_lock != NULL:
                            PyThread_free_lock(input_chunk_caches[i].ready_lock)
                    free(input_chunk_caches)
                    free(inputs_)
                    free(udata)
                    raise MemoryError("Cannot allocate miniexpr chunk cache ready lock")
        udata.input_chunk_caches = input_chunk_caches
        eval_params = <me_eval_params*> malloc(sizeof(me_eval_params))
        if eval_params == NULL:
            for i in range(ninputs):
                if input_chunk_caches[i].state_lock != NULL:
                    PyThread_free_lock(input_chunk_caches[i].state_lock)
                if input_chunk_caches[i].ready_lock != NULL:
                    PyThread_free_lock(input_chunk_caches[i].ready_lock)
            free(input_chunk_caches)
            free(inputs_)
            free(udata)
            raise MemoryError("Cannot allocate miniexpr eval params")
        eval_params.disable_simd = False
        eval_params.simd_ulp_mode = ME_SIMD_ULP_3_5 if fp_accuracy == blosc2.FPAccuracy.MEDIUM else ME_SIMD_ULP_1
        if jit is None:
            eval_params.jit_mode = ME_JIT_DEFAULT
        elif jit:
            eval_params.jit_mode = ME_JIT_ON
        else:
            eval_params.jit_mode = ME_JIT_OFF
        udata.eval_params = eval_params
        udata.array = self.array
        udata.aux_reduc_ptr = aux_reduc_ptr
        # Save these in udf_udata to avoid computing them for each block
        for i in range(self.array.ndim):
            udata.chunks_in_array[i] = udata.array.extshape[i] // udata.array.chunkshape[i]
            udata.blocks_in_chunk[i] = udata.array.extchunkshape[i] // udata.array.blockshape[i]

        return udata

    cdef mm_udata *_fill_mm_udata(self, inputs):
        cdef mm_udata *udata = <mm_udata *> malloc(sizeof(mm_udata))
        cdef int cstrides, bstrides, estrides
        cdef b2nd_array_t* inp
        cdef b2nd_array_t** inputs_ = <b2nd_array_t**> malloc(2 * sizeof(b2nd_array_t*))
        for i in range(2):
            operand = inputs['x1'] if i == 0 else inputs['x2']
            inputs_[i] = <b2nd_array_t*><uintptr_t>operand.c_array
            inputs_[i].chunk_cache.nchunk = -1
            inputs_[i].chunk_cache.data = NULL
        udata.inputs = inputs_
        udata.array = self.array

        # Save these in udf_udata to avoid computing them for each block
        for i in range(3):
            udata.chunks_strides[i][self.array.ndim - 1] = 1
            udata.blocks_strides[i][self.array.ndim - 1] = 1
            udata.el_strides[i][self.array.ndim - 1] = 1
        for idx in range(2, self.array.ndim + 1):
            i = self.array.ndim - idx
            udata.chunks_strides[0][i] = udata.chunks_strides[0][i + 1] * udata.array.extshape[i + 1] // udata.array.chunkshape[i + 1]
            udata.blocks_strides[0][i] = udata.blocks_strides[0][i + 1] * udata.array.extchunkshape[i + 1] // udata.array.blockshape[i + 1]
            udata.el_strides[0][i] = udata.el_strides[0][i + 1] * udata.array.blockshape[i + 1]

        for j in range(1, 3):
            inp = inputs_[j - 1]
            cstrides = bstrides = estrides = 1
            for idx in range(2, self.array.ndim + 1):
                i = inp.ndim - idx
                if (inp.shape[i + 1] == 1 and i < inp.ndim - 3) or i < 0:
                    udata.chunks_strides[j][i] = 0
                    udata.blocks_strides[j][i] = 0
                    udata.el_strides[j][i] = 0
                else:
                    bstrides *= inp.extchunkshape[i + 1] // inp.blockshape[i + 1]
                    cstrides *= inp.extshape[i + 1] // inp.chunkshape[i + 1]
                    estrides *= inp.blockshape[i + 1]
                    udata.chunks_strides[j][i] = cstrides
                    udata.blocks_strides[j][i] = bstrides
                    udata.el_strides[j][i] = estrides

        return udata

    def _set_pref_expr(self, expression, inputs, fp_accuracy, aux_reduc=None, jit=None,
                       candidate_blocks=None):
        # Set prefilter for miniexpr
        cdef blosc2_cparams* cparams = self.array.sc.storage.cparams
        cparams.prefilter = <blosc2_prefilter_fn> miniexpr_prefilter

        cdef int jit_mode = ME_JIT_DEFAULT
        if jit is True:
            jit_mode = ME_JIT_ON
        elif jit is False:
            jit_mode = ME_JIT_OFF

        cdef me_udata* udata = self._fill_me_udata(inputs, fp_accuracy, aux_reduc, jit=jit)

        # Optional candidate-block bitmap (1-D only): a contiguous uint8 array,
        # one byte per global miniexpr block (0 = skip, non-zero = evaluate).
        # The contiguous buffer is returned so the caller can keep it alive until
        # the prefilter is removed (the prefilter runs synchronously during the
        # update_data loop); letting it be collected would dangle the pointer.
        cdef np.ndarray cb = None
        if candidate_blocks is not None and self.array.ndim == 1:
            cb = np.ascontiguousarray(candidate_blocks, dtype=np.uint8)
            udata.candidate_blocks = <const uint8_t*> np.PyArray_DATA(cb)
            udata.candidate_blocks_len = <int64_t> cb.shape[0]

        # Get the compiled expression handle for multi-threading
        cdef Py_ssize_t n = len(inputs)
        cdef me_variable* variables = <me_variable *> malloc(sizeof(me_variable) * n)
        if variables == NULL:
            raise MemoryError()
        cdef me_variable *var
        cdef np.dtype out_np_dtype = np.dtype(self.dtype)
        cdef me_dtype me_output_dtype = _me_dtype_from_numpy_dtype(out_np_dtype)
        if <int>me_output_dtype < 0:
            raise TypeError(f"miniexpr does not support output dtype: {out_np_dtype}")

        cdef np.dtype operand_dtype
        for i, (k, v) in enumerate(inputs.items()):
            var = &variables[i]
            var_name = k.encode("utf-8") if isinstance(k, str) else k
            var.name = <char *> malloc(strlen(var_name) + 1)
            strcpy(var.name, var_name)
            operand_dtype = np.dtype(v.dtype)
            var.dtype = _me_dtype_from_numpy_dtype(operand_dtype)
            if <int>var.dtype < 0:
                raise TypeError(f"miniexpr does not support operand dtype '{operand_dtype}' for input '{k}'")
            var.address = NULL  # chunked compile: addresses provided later
            var.type = 0  # auto-set to ME_VARIABLE inside compiler
            var.context = NULL
            var.itemsize = v.dtype.itemsize if v.dtype.num == 19 else 0 # only store item type if string

        cdef int error = 0
        cdef bytes expression_bytes
        cdef str expression_display
        if isinstance(expression, str):
            expression_display = expression
            expression_bytes = (<str>expression).encode("utf-8")
        elif isinstance(expression, bytes):
            expression_bytes = expression
            expression_display = (<bytes>expression).decode("utf-8", "replace")
        else:
            expression_display = str(expression)
            expression_bytes = expression_display.encode("utf-8")
        cdef me_dtype = me_output_dtype
        cdef me_expr *out_expr
        cdef int ndims = self.array.ndim
        cdef int64_t* shape = &self.array.shape[0]
        cdef int32_t* chunkshape = &self.array.chunkshape[0]
        cdef int32_t* blockshape = &self.array.blockshape[0]
        cdef int rc = me_compile_nd_jit(expression_bytes, variables, n, me_dtype, ndims,
                                        shape, chunkshape, blockshape, jit_mode,
                                        &error, &out_expr)
        cdef str me_error_msg = _me_compile_error_details(rc, error)
        if rc == ME_COMPILE_ERR_INVALID_ARG_TYPE:
            raise TypeError(f"miniexpr does not support operand or output dtype: {expression_display}; details: {me_error_msg}")
        if rc != ME_COMPILE_SUCCESS:
            raise NotImplementedError(f"Cannot compile expression: {expression_display}; details: {me_error_msg}")
        udata.miniexpr_handle = out_expr

        # Free resources
        for i in range(len(inputs)):
            free(variables[i].name)
        free(variables)

        cdef blosc2_prefilter_params* preparams = <blosc2_prefilter_params *> calloc(1, sizeof(blosc2_prefilter_params))
        preparams.user_data = udata
        preparams.output_is_disposable = False if aux_reduc is None else True
        cparams.preparams = preparams
        _check_cparams(cparams)

        if self.array.sc.cctx != NULL:
            # Freeing NULL context can lead to segmentation fault
            blosc2_free_ctx(self.array.sc.cctx)
        self.array.sc.cctx = blosc2_create_cctx(dereference(cparams))
        if self.array.sc.cctx == NULL:
            raise RuntimeError("Could not create compression context")

        # Return the anchored bitmap buffer (or None) so the caller keeps it
        # alive for the duration of the prefilter-driven update_data loop.
        return cb

    def _dsl_jit_status(self, expression, inputs):
        """Compile *expression* with JIT against this array's grid and report whether
        miniexpr produced a runtime JIT kernel (vs an interpreter fallback).

        Compiles and immediately frees the program; it does not set a prefilter or
        run the kernel.  Returns a dict with ``compiled`` (bool), ``jit`` (bool) and
        ``status`` (miniexpr compile-status name).  ``inputs`` is a name -> NDArray
        mapping; all operands must share this array's shape/chunks/blocks.
        """
        cdef Py_ssize_t n = len(inputs)
        cdef me_variable* variables = NULL
        if n > 0:
            variables = <me_variable *> malloc(sizeof(me_variable) * n)
            if variables == NULL:
                raise MemoryError()
        cdef me_variable *var
        cdef np.dtype out_np_dtype = np.dtype(self.dtype)
        cdef me_dtype me_output_dtype = _me_dtype_from_numpy_dtype(out_np_dtype)
        if <int>me_output_dtype < 0:
            free(variables)
            raise TypeError(f"miniexpr does not support output dtype: {out_np_dtype}")

        cdef np.dtype operand_dtype
        cdef bytes var_name
        for i, (k, v) in enumerate(inputs.items()):
            var = &variables[i]
            var_name = k.encode("utf-8") if isinstance(k, str) else k
            var.name = <char *> malloc(strlen(var_name) + 1)
            strcpy(var.name, var_name)
            operand_dtype = np.dtype(v.dtype)
            var.dtype = _me_dtype_from_numpy_dtype(operand_dtype)
            if <int>var.dtype < 0:
                for j in range(i + 1):
                    free(variables[j].name)
                free(variables)
                raise TypeError(f"miniexpr does not support operand dtype '{operand_dtype}' for input '{k}'")
            var.address = NULL
            var.type = 0
            var.context = NULL
            var.itemsize = v.dtype.itemsize if v.dtype.num == 19 else 0

        cdef bytes expression_bytes = (
            (<str>expression).encode("utf-8") if isinstance(expression, str) else expression
        )
        cdef int error = 0
        cdef me_expr *out_expr = NULL
        cdef int ndims = self.array.ndim
        cdef int64_t* shape = &self.array.shape[0]
        cdef int32_t* chunkshape = &self.array.chunkshape[0]
        cdef int32_t* blockshape = &self.array.blockshape[0]
        cdef int rc = me_compile_nd_jit(expression_bytes, variables, n, me_output_dtype, ndims,
                                        shape, chunkshape, blockshape, ME_JIT_ON,
                                        &error, &out_expr)
        for i in range(n):
            free(variables[i].name)
        free(variables)

        cdef bint jit = False
        if rc == ME_COMPILE_SUCCESS and out_expr != NULL:
            jit = me_expr_has_jit_kernel(out_expr)
        if out_expr != NULL:
            me_free(out_expr)
        return {
            "compiled": rc == ME_COMPILE_SUCCESS,
            "jit": bool(jit),
            "status": _me_compile_status_name(rc),
        }

    def _set_pref_matmul(self, inputs, fp_accuracy):
        # Set prefilter for miniexpr
        cdef blosc2_cparams* cparams = self.array.sc.storage.cparams
        cparams.prefilter = <blosc2_prefilter_fn> matmul_prefilter

        cdef mm_udata* udata = self._fill_mm_udata(inputs)
        cdef b2nd_array_t* out_arr = udata.array
        cdef blosc2_prefilter_params* preparams = <blosc2_prefilter_params *> calloc(1, sizeof(blosc2_prefilter_params))
        preparams.user_data = udata
        preparams.output_is_disposable = False
        cparams.preparams = preparams
        _check_cparams(cparams)

        if self.array.sc.cctx != NULL:
            # Freeing NULL context can lead to segmentation fault
            blosc2_free_ctx(self.array.sc.cctx)
        self.array.sc.cctx = blosc2_create_cctx(dereference(cparams))
        if self.array.sc.cctx == NULL:
            raise RuntimeError("Could not create compression context")

    def _set_pref_udf(self, func, inputs_id):
        if self.array.sc.storage.cparams.nthreads > 1:
            raise AttributeError("compress `nthreads` must be 1 when assigning a prefilter")

        func_id = func.__name__
        blosc2.prefilter_funcs[func_id] = func
        func_id = func_id.encode("utf-8") if isinstance(func_id, str) else func_id

        # Set prefilter
        cdef blosc2_cparams* cparams = self.array.sc.storage.cparams
        cparams.prefilter = <blosc2_prefilter_fn> general_udf_prefilter

        cdef blosc2_prefilter_params* preparams = <blosc2_prefilter_params *> calloc(1, sizeof(blosc2_prefilter_params))
        preparams.user_data = self._fill_udf_udata(func_id, inputs_id)
        cparams.preparams = preparams
        _check_cparams(cparams)

        blosc2_free_ctx(self.array.sc.cctx)
        self.array.sc.cctx = blosc2_create_cctx(dereference(cparams))
        if self.array.sc.cctx == NULL:
            raise RuntimeError("Could not create compression context")

    def _set_postf_udf(self, func, inputs_id):
        if self.array.sc.storage.dparams.nthreads > 1:
            raise AttributeError("decompress `nthreads` must be 1 when assigning a postfilter")

        func_id = func.__name__
        blosc2.postfilter_funcs[func_id] = func
        func_id = func_id.encode("utf-8") if isinstance(func_id, str) else func_id

        # Set postfilter
        cdef blosc2_dparams *dparams = self.array.sc.storage.dparams
        dparams.postfilter = <blosc2_postfilter_fn> general_udf_postfilter
        # Fill postparams
        cdef blosc2_postfilter_params *postparams = <blosc2_postfilter_params *> malloc(
            sizeof(blosc2_postfilter_params))
        postparams.user_data = self._fill_udf_udata(func_id,inputs_id)
        dparams.postparams = postparams
        _check_dparams(dparams, self.array.sc.storage.cparams)

        if self.array.sc.dctx != NULL:
            # Freeing NULL context can lead to segmentation fault
            blosc2_free_ctx(self.array.sc.dctx)
        self.array.sc.dctx = blosc2_create_dctx(dereference(dparams))
        if self.array.sc.dctx == NULL:
            raise RuntimeError("Could not create decompression context")

    def __dealloc__(self):
        if self.array != NULL:
            _check_rc(b2nd_free(self.array), "Error while freeing the array")


cdef b2nd_context_t* create_b2nd_context(shape, chunks, blocks, dtype, kwargs):
    if isinstance(dtype, list) and len(dtype) > 0 and isinstance(dtype[0], tuple):
        # Extract just the field names and basic dtype info
        fields = []
        for field in dtype:
            name = field[0]
            field_dtype = field[1]

            # Handle different field formats:
            # 1. ('name', ('|S10', {'h5py_encoding': 'ascii'})) - h5py style
            # 2. ('name', '<i4') - standard NumPy style
            # 3. ('name', '<i4', (shape,)) - NumPy with shape info

            if isinstance(field_dtype, tuple) and len(field_dtype) > 0:
                # h5py nested representation with metadata dict
                field_dtype = field_dtype[0]

            # Check if we have shape information as third element
            if len(field) > 2 and field[2] is not None:
                # Include the shape information
                fields.append((name, field_dtype, field[2]))
            else:
                fields.append((name, field_dtype))

        dtype = np.dtype(fields)
    else:
        dtype = np.dtype(dtype)

    typesize = dtype.itemsize
    if 'cparams' in kwargs:
        kwargs['cparams']['typesize'] = typesize
    else:
        kwargs['cparams'] = {'typesize': typesize} # last filter is shuffle
        if isinstance(dtype, np.dtypes.StrDType) or dtype == np.str_:
            kwargs['cparams']['filters'] = [blosc2.Filter.NOFILTER] * 5 + [blosc2.Filter.SHUFFLE]
            kwargs['cparams']['filters_meta'] = [0] * 5 + [4] # unicode char bytesize
    if dtype.kind == 'V':
        str_dtype = str(dtype)
    else:
        str_dtype = dtype.str
    str_dtype = str_dtype.encode("utf-8") if isinstance(str_dtype, str) else str_dtype

    urlpath = kwargs.get("urlpath")
    if 'contiguous' not in kwargs:
        # Make contiguous true for disk, else sparse (for in-memory performance)
        kwargs['contiguous'] = False if urlpath is None else True

    if urlpath is not None:
        if isinstance(urlpath, pathlib.PurePath):
            urlpath = str(urlpath)
        _urlpath = urlpath.encode() if isinstance(urlpath, str) else urlpath
        kwargs["urlpath"] = _urlpath

    if kwargs.get("mmap_mode") is not None:
        kwargs["mode"] = mode_from_mmap_mode(kwargs["mmap_mode"])

    mode = kwargs.get("mode", "a")
    if kwargs is not None:
        if mode == "w":
            blosc2.remove_urlpath(urlpath)
        elif mode == "r" and urlpath is not None:
            raise ValueError("NDArray must already exist")

    # Create storage
    cdef blosc2_storage storage
    cdef blosc2_cparams *cparams = <blosc2_cparams *>malloc(sizeof(blosc2_cparams))
    cdef blosc2_dparams *dparams = <blosc2_dparams *>malloc(sizeof(blosc2_dparams))
    storage.cparams = cparams
    storage.dparams = dparams
    create_storage(&storage, kwargs)

    # Shapes
    ndim = len(shape)
    cdef int64_t[B2ND_MAX_DIM] shape_
    cdef int32_t[B2ND_MAX_DIM] chunkshape
    cdef int32_t[B2ND_MAX_DIM] blockshape
    for i in range(ndim):
        chunkshape[i] = chunks[i]
        blockshape[i] = blocks[i]
        shape_[i] = shape[i]

    # Metalayers
    meta = kwargs.get('meta', None)
    cdef blosc2_metalayer[B2ND_MAX_METALAYERS] metalayers

    if meta is None:
        return b2nd_create_ctx(&storage, len(shape), shape_, chunkshape, blockshape, str_dtype,
                              B2ND_DEFAULT_DTYPE_FORMAT, NULL, 0)
    else:
        nmetalayers = len(meta)
        for i, (name, content) in enumerate(meta.items()):
            name2 = name.encode("utf-8") if isinstance(name, str) else name # do a copy
            metalayers[i].name = strdup(name2)
            content = packb(content, default=encode_tuple, strict_types=True, use_bin_type=True)
            metalayers[i].content = <uint8_t *> malloc(len(content))
            memcpy(metalayers[i].content, <uint8_t *> content, len(content))
            metalayers[i].content_len = len(content)

        return b2nd_create_ctx(&storage, len(shape), shape_, chunkshape, blockshape, str_dtype,
                              B2ND_DEFAULT_DTYPE_FORMAT, metalayers, nmetalayers)


def uninit(shape, chunks, blocks, dtype, **kwargs):
    cdef b2nd_context_t *ctx = create_b2nd_context(shape, chunks, blocks, dtype, kwargs)
    if ctx == NULL:
        raise RuntimeError("Error while creating the context")

    cdef b2nd_array_t *array
    _check_rc(b2nd_uninit(ctx, &array), "Could not build uninit array")
    _check_rc(b2nd_free_ctx(ctx), "Error while freeing the context")
    ndarray = blosc2.NDArray(_schunk=PyCapsule_New(array.sc, <char *> "blosc2_schunk*", NULL),
                             _array=PyCapsule_New(array, <char *> "b2nd_array_t*", NULL))
    ndarray._schunk.mode = kwargs.get("mode", "a")

    return ndarray


def nans(shape, chunks, blocks, dtype, **kwargs):
    cdef b2nd_context_t *ctx = create_b2nd_context(shape, chunks, blocks, dtype, kwargs)
    if ctx == NULL:
        raise RuntimeError("Error while creating the context")

    cdef b2nd_array_t *array
    _check_rc(b2nd_nans(ctx, &array), "Could not build nans array")
    _check_rc(b2nd_free_ctx(ctx), "Error while freeing the context")
    ndarray = blosc2.NDArray(_schunk=PyCapsule_New(array.sc, <char *> "blosc2_schunk*", NULL),
                             _array=PyCapsule_New(array, <char *> "b2nd_array_t*", NULL))
    ndarray._schunk.mode = kwargs.get("mode", "a")

    return ndarray


def empty(shape, chunks, blocks, dtype, **kwargs):
    cdef b2nd_context_t *ctx = create_b2nd_context(shape, chunks, blocks, dtype, kwargs)
    if ctx == NULL:
        raise RuntimeError("Error while creating the context")

    cdef b2nd_array_t *array
    _check_rc(b2nd_empty(ctx, &array), "Could not build empty array")
    _check_rc(b2nd_free_ctx(ctx), "Error while freeing the context")
    ndarray = blosc2.NDArray(_schunk=PyCapsule_New(array.sc, <char *> "blosc2_schunk*", NULL),
                             _array=PyCapsule_New(array, <char *> "b2nd_array_t*", NULL))
    ndarray._schunk.mode = kwargs.get("mode", "a")

    return ndarray


def zeros(shape, chunks, blocks, dtype, **kwargs):
    cdef b2nd_context_t *ctx = create_b2nd_context(shape, chunks, blocks, dtype, kwargs)
    if ctx == NULL:
        raise RuntimeError("Error while creating the context")

    cdef b2nd_array_t *array
    _check_rc(b2nd_zeros(ctx, &array), "Could not build zeros array")
    ndarray = blosc2.NDArray(_schunk=PyCapsule_New(array.sc, <char *> "blosc2_schunk*", NULL),
                             _array=PyCapsule_New(array, <char *> "b2nd_array_t*", NULL))
    _check_rc(b2nd_free_ctx(ctx), "Error while freeing the context")
    ndarray._schunk.mode = kwargs.get("mode", "a")

    return ndarray


def full(shape, chunks, blocks, fill_value, dtype, **kwargs):
    cdef b2nd_context_t *ctx = create_b2nd_context(shape, chunks, blocks, dtype, kwargs)
    if ctx == NULL:
        raise RuntimeError("Error while creating the context")

    dtype = np.dtype(dtype)
    nparr = np.array([fill_value], dtype=dtype)
    cdef Py_buffer val
    PyObject_GetBuffer(nparr, &val, PyBUF_SIMPLE)

    cdef b2nd_array_t *array
    _check_rc(b2nd_full(ctx, &array, val.buf), "Could not create full array")
    PyBuffer_Release(&val)

    ndarray = blosc2.NDArray(_schunk=PyCapsule_New(array.sc, <char *> "blosc2_schunk*", NULL),
                             _array=PyCapsule_New(array, <char *> "b2nd_array_t*", NULL))
    _check_rc(b2nd_free_ctx(ctx), "Error while freeing the context")
    ndarray._schunk.mode = kwargs.get("mode", "a")

    return ndarray


def from_buffer(buf, shape, chunks, blocks, dtype, **kwargs):
    cdef b2nd_context_t *ctx = create_b2nd_context(shape, chunks, blocks, dtype, kwargs)
    if ctx == NULL:
        raise RuntimeError("Error while creating the context")

    cdef b2nd_array_t *array
    _check_rc(b2nd_from_cbuffer(ctx, &array,  <void*> <char *> buf, len(buf)),
              "Error while creating the NDArray")
    ndarray = blosc2.NDArray(_schunk=PyCapsule_New(array.sc, <char *> "blosc2_schunk*", NULL),
                             _array=PyCapsule_New(array, <char *> "b2nd_array_t*", NULL))
    _check_rc(b2nd_free_ctx(ctx), "Error while freeing the context")
    ndarray._schunk.mode = kwargs.get("mode", "a")

    return ndarray


def asarray(ndarray, chunks, blocks, **kwargs):
    interface = ndarray.__array_interface__
    cdef Py_buffer buf
    PyObject_GetBuffer(ndarray, &buf, PyBUF_SIMPLE)

    shape = interface["shape"]
    dtype = interface["typestr"]
    if dtype.startswith("|V") and "descr" in interface:
        # Structured dtype
        dtype = interface["descr"]
    cdef b2nd_context_t *ctx = create_b2nd_context(shape, chunks, blocks, dtype, kwargs)
    if ctx == NULL:
        raise RuntimeError("Error while creating the context")

    cdef b2nd_array_t *array
    _check_rc(b2nd_from_cbuffer(ctx, &array, <void *> <char *> buf.buf, buf.len),
              "Error while creating the NDArray")
    PyBuffer_Release(&buf)
    ndarray = blosc2.NDArray(_schunk=PyCapsule_New(array.sc, <char *> "blosc2_schunk*", NULL),
                             _array=PyCapsule_New(array, <char *> "b2nd_array_t*", NULL))
    _check_rc(b2nd_free_ctx(ctx), "Error while freeing the context")
    ndarray._schunk.mode = kwargs.get("mode", "a")

    return ndarray

def array_from_ffi_ptr(array_ptr):
    array = <b2nd_array_t *> PyCapsule_GetPointer(array_ptr, <char *> "b2nd_array_t*")
    return blosc2.NDArray(_schunk=PyCapsule_New(array.sc, <char *> "blosc2_schunk*", NULL),
                             _array=array_ptr)

def ndarray_from_cframe(cframe, copy=False):
    cdef Py_buffer buf
    PyObject_GetBuffer(cframe, &buf, PyBUF_SIMPLE)
    cdef b2nd_array_t *array
    cdef int rc
    rc = b2nd_from_cframe(<uint8_t *>buf.buf, buf.len, copy, &array)
    if rc < 0:
        raise RuntimeError("Could not get the NDArray from the cframe")
    ndarray = blosc2.NDArray(_schunk=PyCapsule_New(array.sc, <char *> "blosc2_schunk*", NULL),
                             _array=PyCapsule_New(array, <char *> "b2nd_array_t*", NULL))

    PyBuffer_Release(&buf)
    if not copy:
        ndarray._schunk._avoid_cframe_free(True)
    return ndarray


def array_get_slice_nchunks(array: NDArray, key):
    start, stop = key
    cdef int64_t[B2ND_MAX_DIM] start_, stop_
    for i in range(array.ndim):
        start_[i] = start[i]
        stop_[i] = stop[i]
    cdef int64_t *chunks_idx
    rc = blosc2_get_slice_nchunks(array.array.sc, start_, stop_, &chunks_idx)
    _check_rc(rc, "Error while getting the chunk indexes")
    res = np.empty(rc, dtype=np.int64)
    for i in range(rc):
        res[i] = chunks_idx[i]
    free(chunks_idx)
    return res


def schunk_get_slice_nchunks(schunk: SChunk, key):
    start, stop = key
    nitems = schunk.nbytes // schunk.typesize
    start, stop, _ = slice(start, stop, 1).indices(nitems)

    cdef int64_t start_, stop_
    start_ = start
    stop_ = stop
    cdef int64_t *chunks_idx
    rc = blosc2_get_slice_nchunks(schunk.schunk, &start_, &stop_, &chunks_idx)
    _check_rc(rc, "Error while getting the chunk indexes")

    res = np.empty(rc, dtype=np.int64)
    for i in range(rc):
        res[i] = chunks_idx[i]
    free(chunks_idx)
    return res


def concat(arr1: NDArray, arr2: NDArray, axis: int, **kwargs):
    """
    Concatenate two NDArray objects along a specified axis.
    """
    cdef c_bool copy = kwargs.pop("copy", True)
    cdef b2nd_context_t *ctx = create_b2nd_context(arr1.shape, arr1.chunks, arr1.blocks, arr1.dtype, kwargs)
    if ctx == NULL:
        raise RuntimeError("Error while creating the context for concatenation")

    cdef b2nd_array_t *array
    _check_rc(b2nd_concatenate(ctx, arr1.array, arr2.array, axis, copy, &array),
              "Error while concatenating the arrays")
    _check_rc(b2nd_free_ctx(ctx), "Error while freeing the context")

    if copy:
        # We have copied the concatenated data into a new array
        return blosc2.NDArray(_schunk=PyCapsule_New(array.sc, <char *> "blosc2_schunk*", NULL),
                              _array=PyCapsule_New(array, <char *> "b2nd_array_t*", NULL))
    else:
        # Return the first array, which now contains the concatenated data
        return arr1

def expand_dims(arr1: NDArray, axis_mask: list[bool], final_dims: int) -> blosc2.NDArray:
    """
    Add new dummy axis to NDArray object at specified dimension.
    """
    cdef b2nd_array_t *view
    cdef c_bool mask_[B2ND_MAX_DIM]
    if final_dims > B2ND_MAX_DIM:
        raise ValueError(f"Cannot expand dimensions beyond {B2ND_MAX_DIM} dimensions")
    for i in range(final_dims):
        mask_[i] = axis_mask[i]
    _check_rc(b2nd_expand_dims(arr1.array, &view, mask_, final_dims),"Error while expanding the arrays")

    # create view with reference to arr1 to hold onto
    new_base = arr1 if arr1.base is None else arr1.base
    return blosc2.NDArray(_schunk=PyCapsule_New(view.sc, <char *> "blosc2_schunk*", NULL),
                          _array=PyCapsule_New(view, <char *> "b2nd_array_t*", NULL), _base=new_base)

def squeeze(arr1: NDArray, axis_mask: list[bool]) -> blosc2.NDArray:
    """
    Remove axis from NDArray object at specified dimensions.
    """
    cdef b2nd_array_t *view
    cdef c_bool mask_[B2ND_MAX_DIM]
    for i in range(arr1.ndim):
        mask_[i] = axis_mask[i]
    _check_rc(b2nd_squeeze_index(arr1.array, &view, mask_), "Error while squeezing array")

    # this squeezes even if not asked for by mask - may have to use in future though
    # if arr1.array.shape[0] == 1 and arr1.ndim == 1:
    #     arr1.array.ndim = 0

    # create view with reference to self to hold onto
    new_base = arr1 if arr1.base is None else arr1.base
    return blosc2.NDArray(_schunk=PyCapsule_New(view.sc, <char *> "blosc2_schunk*", NULL),
                        _array=PyCapsule_New(view, <char *> "b2nd_array_t*", NULL), _base=new_base)


def set_matmul_block_backend(mode):
    if mode == "auto":
        b2_set_matmul_backend(B2_MATMUL_BACKEND_AUTO)
    elif mode == "naive":
        b2_set_matmul_backend(B2_MATMUL_BACKEND_NAIVE)
    elif mode == "accelerate":
        b2_set_matmul_backend(B2_MATMUL_BACKEND_ACCELERATE)
    elif mode == "cblas":
        b2_set_matmul_backend(B2_MATMUL_BACKEND_CBLAS)
    else:
        raise ValueError("mode must be 'auto', 'naive', 'accelerate', or 'cblas'")


def get_matmul_block_backend():
    return b2_get_matmul_backend_name().decode("utf-8")


def get_selected_matmul_block_backend():
    return b2_get_selected_matmul_backend_name().decode("utf-8")


def get_loaded_matmul_cblas_library():
    cdef const char *path = b2_get_loaded_cblas_path()
    if path == NULL:
        return None
    return path.decode("utf-8")
