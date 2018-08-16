# -*- coding: utf-8 -*-
from cpython cimport PyObject
from numpy cimport int64_t, uint64_t, int32_t, uint32_t, float64_t

cdef extern from "khash_python.h":
    ctypedef uint32_t khint_t
    ctypedef khint_t khiter_t

    ctypedef struct kh_pymap_t:
        khint_t n_buckets, size, n_occupied, upper_bound
        uint32_t *flags
        PyObject **keys
        size_t *vals

    kh_pymap_t* kh_init_pymap()
    void kh_destroy_pymap(kh_pymap_t*)
    void kh_clear_pymap(kh_pymap_t*)
    khint_t kh_get_pymap(kh_pymap_t*, PyObject*)
    void kh_resize_pymap(kh_pymap_t*, khint_t)
    khint_t kh_put_pymap(kh_pymap_t*, PyObject*, int*)
    void kh_del_pymap(kh_pymap_t*, khint_t)

    bint kh_exist_pymap(kh_pymap_t*, khiter_t)

    ctypedef struct kh_pyset_t:
        khint_t n_buckets, size, n_occupied, upper_bound
        uint32_t *flags
        PyObject **keys
        size_t *vals

    kh_pyset_t* kh_init_pyset()
    void kh_destroy_pyset(kh_pyset_t*)
    void kh_clear_pyset(kh_pyset_t*)
    khint_t kh_get_pyset(kh_pyset_t*, PyObject*)
    void kh_resize_pyset(kh_pyset_t*, khint_t)
    khint_t kh_put_pyset(kh_pyset_t*, PyObject*, int*)
    void kh_del_pyset(kh_pyset_t*, khint_t)

    bint kh_exist_pyset(kh_pyset_t*, khiter_t)

    ctypedef char* kh_cstr_t

    ctypedef struct kh_str_t:
        khint_t n_buckets, size, n_occupied, upper_bound
        uint32_t *flags
        kh_cstr_t *keys
        size_t *vals

    kh_str_t* kh_init_str() nogil
    void kh_destroy_str(kh_str_t*) nogil
    void kh_clear_str(kh_str_t*) nogil
    khint_t kh_get_str(kh_str_t*, kh_cstr_t) nogil
    void kh_resize_str(kh_str_t*, khint_t) nogil
    khint_t kh_put_str(kh_str_t*, kh_cstr_t, int*) nogil
    void kh_del_str(kh_str_t*, khint_t) nogil

    bint kh_exist_str(kh_str_t*, khiter_t) nogil

    ctypedef struct kh_int64_t:
        khint_t n_buckets, size, n_occupied, upper_bound
        uint32_t *flags
        int64_t *keys
        size_t *vals

    kh_int64_t* kh_init_int64() nogil
    void kh_destroy_int64(kh_int64_t*) nogil
    void kh_clear_int64(kh_int64_t*) nogil
    khint_t kh_get_int64(kh_int64_t*, int64_t) nogil
    void kh_resize_int64(kh_int64_t*, khint_t) nogil
    khint_t kh_put_int64(kh_int64_t*, int64_t, int*) nogil
    void kh_del_int64(kh_int64_t*, khint_t) nogil

    bint kh_exist_int64(kh_int64_t*, khiter_t) nogil

    ctypedef uint64_t khuint64_t

    ctypedef struct kh_uint64_t:
        khint_t n_buckets, size, n_occupied, upper_bound
        uint32_t *flags
        khuint64_t *keys
        size_t *vals

    kh_uint64_t* kh_init_uint64() nogil
    void kh_destroy_uint64(kh_uint64_t*) nogil
    void kh_clear_uint64(kh_uint64_t*) nogil
    khint_t kh_get_uint64(kh_uint64_t*, uint64_t) nogil
    void kh_resize_uint64(kh_uint64_t*, khint_t) nogil
    khint_t kh_put_uint64(kh_uint64_t*, uint64_t, int*) nogil
    void kh_del_uint64(kh_uint64_t*, khint_t) nogil

    bint kh_exist_uint64(kh_uint64_t*, khiter_t) nogil

    ctypedef struct kh_float64_t:
        khint_t n_buckets, size, n_occupied, upper_bound
        uint32_t *flags
        float64_t *keys
        size_t *vals

    kh_float64_t* kh_init_float64() nogil
    void kh_destroy_float64(kh_float64_t*) nogil
    void kh_clear_float64(kh_float64_t*) nogil
    khint_t kh_get_float64(kh_float64_t*, float64_t) nogil
    void kh_resize_float64(kh_float64_t*, khint_t) nogil
    khint_t kh_put_float64(kh_float64_t*, float64_t, int*) nogil
    void kh_del_float64(kh_float64_t*, khint_t) nogil

    bint kh_exist_float64(kh_float64_t*, khiter_t) nogil

    ctypedef struct kh_int32_t:
        khint_t n_buckets, size, n_occupied, upper_bound
        uint32_t *flags
        int32_t *keys
        size_t *vals

    kh_int32_t* kh_init_int32() nogil
    void kh_destroy_int32(kh_int32_t*) nogil
    void kh_clear_int32(kh_int32_t*) nogil
    khint_t kh_get_int32(kh_int32_t*, int32_t) nogil
    void kh_resize_int32(kh_int32_t*, khint_t) nogil
    khint_t kh_put_int32(kh_int32_t*, int32_t, int*) nogil
    void kh_del_int32(kh_int32_t*, khint_t) nogil

    bint kh_exist_int32(kh_int32_t*, khiter_t) nogil

    # sweep factorize

    ctypedef struct kh_strbox_t:
        khint_t n_buckets, size, n_occupied, upper_bound
        uint32_t *flags
        kh_cstr_t *keys
        PyObject **vals

    kh_strbox_t* kh_init_strbox() nogil
    void kh_destroy_strbox(kh_strbox_t*) nogil
    void kh_clear_strbox(kh_strbox_t*) nogil
    khint_t kh_get_strbox(kh_strbox_t*, kh_cstr_t) nogil
    void kh_resize_strbox(kh_strbox_t*, khint_t) nogil
    khint_t kh_put_strbox(kh_strbox_t*, kh_cstr_t, int*) nogil
    void kh_del_strbox(kh_strbox_t*, khint_t) nogil

    bint kh_exist_strbox(kh_strbox_t*, khiter_t) nogil
