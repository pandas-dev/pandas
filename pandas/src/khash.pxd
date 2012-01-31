from cpython cimport PyObject
from numpy cimport int64_t, int32_t, uint32_t, float64_t

cdef extern from "khash.h":
    ctypedef uint32_t khint_t
    ctypedef khint_t khiter_t

    ctypedef struct kh_pymap_t:
        khint_t n_buckets, size, n_occupied, upper_bound
        uint32_t *flags
        PyObject **keys
        Py_ssize_t *vals

    inline kh_pymap_t* kh_init_pymap()
    inline void kh_destroy_pymap(kh_pymap_t*)
    inline void kh_clear_pymap(kh_pymap_t*)
    inline khint_t kh_get_pymap(kh_pymap_t*, PyObject*)
    inline void kh_resize_pymap(kh_pymap_t*, khint_t)
    inline khint_t kh_put_pymap(kh_pymap_t*, PyObject*, int*)
    inline void kh_del_pymap(kh_pymap_t*, khint_t)

    bint kh_exist_pymap(kh_pymap_t*, khiter_t)

    ctypedef struct kh_pyset_t:
        khint_t n_buckets, size, n_occupied, upper_bound
        uint32_t *flags
        PyObject **keys
        Py_ssize_t *vals

    inline kh_pyset_t* kh_init_pyset()
    inline void kh_destroy_pyset(kh_pyset_t*)
    inline void kh_clear_pyset(kh_pyset_t*)
    inline khint_t kh_get_pyset(kh_pyset_t*, PyObject*)
    inline void kh_resize_pyset(kh_pyset_t*, khint_t)
    inline khint_t kh_put_pyset(kh_pyset_t*, PyObject*, int*)
    inline void kh_del_pyset(kh_pyset_t*, khint_t)

    bint kh_exist_pyset(kh_pyset_t*, khiter_t)

    ctypedef char* kh_cstr_t

    ctypedef struct kh_str_t:
        khint_t n_buckets, size, n_occupied, upper_bound
        uint32_t *flags
        kh_cstr_t *keys
        Py_ssize_t *vals

    inline kh_str_t* kh_init_str()
    inline void kh_destroy_str(kh_str_t*)
    inline void kh_clear_str(kh_str_t*)
    inline khint_t kh_get_str(kh_str_t*, kh_cstr_t)
    inline void kh_resize_str(kh_str_t*, khint_t)
    inline khint_t kh_put_str(kh_str_t*, kh_cstr_t, int*)
    inline void kh_del_str(kh_str_t*, khint_t)

    bint kh_exist_str(kh_str_t*, khiter_t)

    ctypedef struct kh_int64_t:
        khint_t n_buckets, size, n_occupied, upper_bound
        uint32_t *flags
        int64_t *keys
        Py_ssize_t *vals

    inline kh_int64_t* kh_init_int64()
    inline void kh_destroy_int64(kh_int64_t*)
    inline void kh_clear_int64(kh_int64_t*)
    inline khint_t kh_get_int64(kh_int64_t*, int64_t)
    inline void kh_resize_int64(kh_int64_t*, khint_t)
    inline khint_t kh_put_int64(kh_int64_t*, int64_t, int*)
    inline void kh_del_int64(kh_int64_t*, khint_t)

    bint kh_exist_int64(kh_int64_t*, khiter_t)

    ctypedef struct kh_float64_t:
        khint_t n_buckets, size, n_occupied, upper_bound
        uint32_t *flags
        float64_t *keys
        Py_ssize_t *vals

    inline kh_float64_t* kh_init_float64()
    inline void kh_destroy_float64(kh_float64_t*)
    inline void kh_clear_float64(kh_float64_t*)
    inline khint_t kh_get_float64(kh_float64_t*, float64_t)
    inline void kh_resize_float64(kh_float64_t*, khint_t)
    inline khint_t kh_put_float64(kh_float64_t*, float64_t, int*)
    inline void kh_del_float64(kh_float64_t*, khint_t)

    bint kh_exist_float64(kh_float64_t*, khiter_t)

    ctypedef struct kh_int32_t:
        khint_t n_buckets, size, n_occupied, upper_bound
        uint32_t *flags
        int32_t *keys
        Py_ssize_t *vals

    inline kh_int32_t* kh_init_int32()
    inline void kh_destroy_int32(kh_int32_t*)
    inline void kh_clear_int32(kh_int32_t*)
    inline khint_t kh_get_int32(kh_int32_t*, int32_t)
    inline void kh_resize_int32(kh_int32_t*, khint_t)
    inline khint_t kh_put_int32(kh_int32_t*, int32_t, int*)
    inline void kh_del_int32(kh_int32_t*, khint_t)

    bint kh_exist_int32(kh_int32_t*, khiter_t)

