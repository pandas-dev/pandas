////////// MemviewSliceStruct.proto //////////
//@proto_block: utility_code_proto_before_types
//@substitute: naming
//@requires: Synchronization.c::Atomics

/* memoryview slice struct */
struct $memview_objstruct_cname;

typedef struct {
  struct $memview_objstruct_cname *memview;
  char *data;
  Py_ssize_t shape[{{max_dims}}];
  Py_ssize_t strides[{{max_dims}}];
  Py_ssize_t suboffsets[{{max_dims}}];
} $memviewslice_cname;

// used for "len(memviewslice)"
#define __Pyx_MemoryView_Len(m)  (m.shape[0])

// Metadata flags for slice creation and validation.
#define __Pyx_MEMVIEW_DIRECT   1
#define __Pyx_MEMVIEW_PTR      2
#define __Pyx_MEMVIEW_FULL     4
#define __Pyx_MEMVIEW_CONTIG   8
#define __Pyx_MEMVIEW_STRIDED  16
#define __Pyx_MEMVIEW_FOLLOW   32

#define __Pyx_IS_C_CONTIG 1
#define __Pyx_IS_F_CONTIG 2

#define __Pyx_MEMSLICE_INIT  {{memslice_init}}

#if CYTHON_ATOMICS
    #define __pyx_add_acquisition_count(memview) \
             __pyx_atomic_incr_relaxed(__pyx_get_slice_count_pointer(memview))
    #define __pyx_sub_acquisition_count(memview) \
            __pyx_atomic_decr_acq_rel(__pyx_get_slice_count_pointer(memview))
#else
    #define __pyx_add_acquisition_count(memview) \
            __pyx_add_acquisition_count_locked(__pyx_get_slice_count_pointer(memview), memview->lock)
    #define __pyx_sub_acquisition_count(memview) \
            __pyx_sub_acquisition_count_locked(__pyx_get_slice_count_pointer(memview), memview->lock)
#endif

/////////////// ObjectToMemviewSlice.proto ///////////////
//@substitute: naming

static CYTHON_INLINE $memviewslice_cname {{funcname}}(PyObject *, int writable_flag);


/////////////// ObjectToMemviewSlice ///////////////
//@requires: Buffer.c::BufferFormatStructs
//@requires: MemviewSliceValidateAndInit
//@substitute: naming

static CYTHON_INLINE $memviewslice_cname {{funcname}}(PyObject *obj, int writable_flag) {
    $memviewslice_cname result = __Pyx_MEMSLICE_INIT;
    __Pyx_BufFmt_StackElem stack[{{struct_nesting_depth}}];
    int axes_specs[] = { {{axes_specs}} };
    int retcode;

    if (obj == Py_None) {
        /* We don't bother to refcount None */
        result.memview = (struct $memview_objstruct_cname *) Py_None;
        return result;
    }

    retcode = __Pyx_ValidateAndInit_memviewslice(axes_specs, {{c_or_f_flag}},
                                                 {{buf_flag}} | writable_flag, {{ndim}},
                                                 &{{dtype_typeinfo}}, stack,
                                                 &result, obj);

    if (unlikely(retcode == -1))
        goto __pyx_fail;

    return result;
__pyx_fail:
    result.memview = NULL;
    result.data = NULL;
    return result;
}


/////////////// MemviewSliceValidateAndInit.export ///////////////
//@requires: Buffer.c::BufferFormatStructs
//@substitute: naming

static int __Pyx_ValidateAndInit_memviewslice(
                int *axes_specs,
                int c_or_f_flag,
                int buf_flags,
                int ndim,
                const __Pyx_TypeInfo *dtype,
                __Pyx_BufFmt_StackElem stack[],
                $memviewslice_cname *memviewslice,
                PyObject *original_obj);


/////////////// MemviewSliceValidateAndInit ///////////////
//@requires: Buffer.c::TypeInfoCompare
//@requires: Buffer.c::BufferFormatStructs
//@requires: Buffer.c::BufferFormatCheck
//@requires: MemviewSliceInit
//@substitute: naming

static int
__pyx_check_strides(Py_buffer *buf, int dim, int ndim, int spec)
{
    if (buf->shape[dim] <= 1)
        return 1;

    if (buf->strides) {
        if (spec & __Pyx_MEMVIEW_CONTIG) {
            if (spec & (__Pyx_MEMVIEW_PTR|__Pyx_MEMVIEW_FULL)) {
                if (unlikely(buf->strides[dim] != sizeof(void *))) {
                    PyErr_Format(PyExc_ValueError,
                                 "Buffer is not indirectly contiguous "
                                 "in dimension %d.", dim);
                    goto fail;
                }
            } else if (unlikely(buf->strides[dim] != buf->itemsize)) {
                PyErr_SetString(PyExc_ValueError,
                                "Buffer and memoryview are not contiguous "
                                "in the same dimension.");
                goto fail;
            }
        }

        if (spec & __Pyx_MEMVIEW_FOLLOW) {
            Py_ssize_t stride = buf->strides[dim];
            if (stride < 0)
                stride = -stride;
            if (unlikely(stride < buf->itemsize)) {
                PyErr_SetString(PyExc_ValueError,
                                "Buffer and memoryview are not contiguous "
                                "in the same dimension.");
                goto fail;
            }
        }
    } else {
        if (unlikely(spec & __Pyx_MEMVIEW_CONTIG && dim != ndim - 1)) {
            PyErr_Format(PyExc_ValueError,
                         "C-contiguous buffer is not contiguous in "
                         "dimension %d", dim);
            goto fail;
        } else if (unlikely(spec & (__Pyx_MEMVIEW_PTR))) {
            PyErr_Format(PyExc_ValueError,
                         "C-contiguous buffer is not indirect in "
                         "dimension %d", dim);
            goto fail;
        } else if (unlikely(buf->suboffsets)) {
            PyErr_SetString(PyExc_ValueError,
                            "Buffer exposes suboffsets but no strides");
            goto fail;
        }
    }

    return 1;
fail:
    return 0;
}

static int
__pyx_check_suboffsets(Py_buffer *buf, int dim, int ndim, int spec)
{
    CYTHON_UNUSED_VAR(ndim);
    // Todo: without PyBUF_INDIRECT we may not have suboffset information, i.e., the
    //       ptr may not be set to NULL but may be uninitialized?
    if (spec & __Pyx_MEMVIEW_DIRECT) {
        if (unlikely(buf->suboffsets && buf->suboffsets[dim] >= 0)) {
            PyErr_Format(PyExc_ValueError,
                         "Buffer not compatible with direct access "
                         "in dimension %d.", dim);
            goto fail;
        }
    }

    if (spec & __Pyx_MEMVIEW_PTR) {
        if (unlikely(!buf->suboffsets || (buf->suboffsets[dim] < 0))) {
            PyErr_Format(PyExc_ValueError,
                         "Buffer is not indirectly accessible "
                         "in dimension %d.", dim);
            goto fail;
        }
    }

    return 1;
fail:
    return 0;
}

static int
__pyx_verify_contig(Py_buffer *buf, int ndim, int c_or_f_flag)
{
    int i;

    if (c_or_f_flag & __Pyx_IS_F_CONTIG) {
        Py_ssize_t stride = 1;
        for (i = 0; i < ndim; i++) {
            if (unlikely(stride * buf->itemsize != buf->strides[i]  &&  buf->shape[i] > 1)) {
                PyErr_SetString(PyExc_ValueError,
                    "Buffer not fortran contiguous.");
                goto fail;
            }
            stride = stride * buf->shape[i];
        }
    } else if (c_or_f_flag & __Pyx_IS_C_CONTIG) {
        Py_ssize_t stride = 1;
        for (i = ndim - 1; i >- 1; i--) {
            if (unlikely(stride * buf->itemsize != buf->strides[i]  &&  buf->shape[i] > 1)) {
                PyErr_SetString(PyExc_ValueError,
                    "Buffer not C contiguous.");
                goto fail;
            }
            stride = stride * buf->shape[i];
        }
    }

    return 1;
fail:
    return 0;
}

static int __Pyx_ValidateAndInit_memviewslice(
                int *axes_specs,
                int c_or_f_flag,
                int buf_flags,
                int ndim,
                const __Pyx_TypeInfo *dtype,
                __Pyx_BufFmt_StackElem stack[],
                $memviewslice_cname *memviewslice,
                PyObject *original_obj)
{
    struct $memview_objstruct_cname *memview, *new_memview;
    __Pyx_RefNannyDeclarations
    Py_buffer *buf;
    int i, spec = 0, retval = -1;
    __Pyx_BufFmt_Context ctx;
    int from_memoryview = __pyx_memoryview_check(original_obj);

    __Pyx_RefNannySetupContext("ValidateAndInit_memviewslice", 0);

    if (from_memoryview && __pyx_typeinfo_cmp(dtype, ((struct $memview_objstruct_cname *)
                                                            original_obj)->typeinfo)) {
        /* We have a matching dtype, skip format parsing */
        memview = (struct $memview_objstruct_cname *) original_obj;
        new_memview = NULL;
    } else {
        memview = (struct $memview_objstruct_cname *) __pyx_memoryview_new(
                                            original_obj, buf_flags, 0, dtype);
        new_memview = memview;
        if (unlikely(!memview))
            goto fail;
    }

    buf = &memview->view;
    if (unlikely(buf->ndim != ndim)) {
        PyErr_Format(PyExc_ValueError,
                "Buffer has wrong number of dimensions (expected %d, got %d)",
                ndim, buf->ndim);
        goto fail;
    }

    if (new_memview) {
        __Pyx_BufFmt_Init(&ctx, stack, dtype);
        if (unlikely(!__Pyx_BufFmt_CheckString(&ctx, buf->format))) goto fail;
    }

    if (unlikely((unsigned) buf->itemsize != dtype->size)) {
        PyErr_Format(PyExc_ValueError,
                     "Item size of buffer (%" CYTHON_FORMAT_SSIZE_T "u byte%s) "
                     "does not match size of '%s' (%" CYTHON_FORMAT_SSIZE_T "u byte%s)",
                     buf->itemsize,
                     (buf->itemsize > 1) ? "s" : "",
                     dtype->name,
                     dtype->size,
                     (dtype->size > 1) ? "s" : "");
        goto fail;
    }

    /* Check axes */
    if (buf->len > 0) {
        // 0-sized arrays do not undergo these checks since their strides are
        // irrelevant and they are always both C- and F-contiguous.
        for (i = 0; i < ndim; i++) {
            spec = axes_specs[i];
            if (unlikely(!__pyx_check_strides(buf, i, ndim, spec)))
                goto fail;
            if (unlikely(!__pyx_check_suboffsets(buf, i, ndim, spec)))
                goto fail;
        }

        /* Check contiguity */
        if (unlikely(buf->strides && !__pyx_verify_contig(buf, ndim, c_or_f_flag)))
            goto fail;
    }

    /* Initialize */
    if (unlikely(__Pyx_init_memviewslice(memview, ndim, memviewslice,
                                         new_memview != NULL) == -1)) {
        goto fail;
    }

    retval = 0;
    goto no_fail;

fail:
    Py_XDECREF((PyObject*)new_memview);
    retval = -1;

no_fail:
    __Pyx_RefNannyFinishContext();
    return retval;
}


////////// MemviewSliceInit.proto //////////
//@substitute: naming

static int __Pyx_init_memviewslice(
                struct $memview_objstruct_cname *memview,
                int ndim,
                $memviewslice_cname *memviewslice,
                int memview_is_new_reference);


////////// MemviewSliceInit //////////
//@substitute: naming

static int
__Pyx_init_memviewslice(struct $memview_objstruct_cname *memview,
                        int ndim,
                        $memviewslice_cname *memviewslice,
                        int memview_is_new_reference)
{
    __Pyx_RefNannyDeclarations
    int i, retval=-1;
    Py_buffer *buf = &memview->view;
    __Pyx_RefNannySetupContext("init_memviewslice", 0);

    if (unlikely(memviewslice->memview || memviewslice->data)) {
        PyErr_SetString(PyExc_ValueError,
            "memviewslice is already initialized!");
        goto fail;
    }

    if (buf->strides) {
        for (i = 0; i < ndim; i++) {
            memviewslice->strides[i] = buf->strides[i];
        }
    } else {
        Py_ssize_t stride = buf->itemsize;
        for (i = ndim - 1; i >= 0; i--) {
            memviewslice->strides[i] = stride;
            stride *= buf->shape[i];
        }
    }

    for (i = 0; i < ndim; i++) {
        memviewslice->shape[i]   = buf->shape[i];
        if (buf->suboffsets) {
            memviewslice->suboffsets[i] = buf->suboffsets[i];
        } else {
            memviewslice->suboffsets[i] = -1;
        }
    }

    memviewslice->memview = memview;
    memviewslice->data = (char *)buf->buf;
    if (__pyx_add_acquisition_count(memview) == 0 && !memview_is_new_reference) {
        Py_INCREF((PyObject*)memview);
    }
    retval = 0;
    goto no_fail;

fail:
    /* Don't decref, the memoryview may be borrowed. Let the caller do the cleanup */
    /* __Pyx_XDECREF(memviewslice->memview); */
    memviewslice->memview = 0;
    memviewslice->data = 0;
    retval = -1;
no_fail:
    __Pyx_RefNannyFinishContext();
    return retval;
}


////////// MemviewRefcount.proto //////////
//@substitute: naming

static CYTHON_INLINE int __pyx_add_acquisition_count_locked(
    __pyx_atomic_int_type *acquisition_count, PyThread_type_lock lock);
static CYTHON_INLINE int __pyx_sub_acquisition_count_locked(
    __pyx_atomic_int_type *acquisition_count, PyThread_type_lock lock);

#define __pyx_get_slice_count_pointer(memview) (&memview->acquisition_count)
#define __PYX_INC_MEMVIEW(slice, have_gil) __Pyx_INC_MEMVIEW(slice, have_gil, __LINE__)
#define __PYX_XCLEAR_MEMVIEW(slice, have_gil) __Pyx_XCLEAR_MEMVIEW(slice, have_gil, __LINE__)
static CYTHON_INLINE void __Pyx_INC_MEMVIEW($memviewslice_cname *, int, int);
static CYTHON_INLINE void __Pyx_XCLEAR_MEMVIEW($memviewslice_cname *, int, int);


////////// MemviewRefcount //////////
//@substitute: naming

// vsnprintf
#include <stdio.h>

#ifndef _Py_NO_RETURN
// It's non-public and just a hint, so allow it to be missing.
#define _Py_NO_RETURN
#endif

_Py_NO_RETURN
static void __pyx_fatalerror(const char *fmt, ...) {
    va_list vargs;
    char msg[200];

#if PY_VERSION_HEX >= 0x030A0000 || defined(HAVE_STDARG_PROTOTYPES)
    va_start(vargs, fmt);
#else
    va_start(vargs);
#endif
    vsnprintf(msg, 200, fmt, vargs);
    va_end(vargs);

    Py_FatalError(msg);
}

static CYTHON_INLINE int
__pyx_add_acquisition_count_locked(__pyx_atomic_int_type *acquisition_count,
                                   PyThread_type_lock lock)
{
    int result;
    PyThread_acquire_lock(lock, 1);
    result = (*acquisition_count)++;
    PyThread_release_lock(lock);
    return result;
}

static CYTHON_INLINE int
__pyx_sub_acquisition_count_locked(__pyx_atomic_int_type *acquisition_count,
                                   PyThread_type_lock lock)
{
    int result;
    PyThread_acquire_lock(lock, 1);
    result = (*acquisition_count)--;
    PyThread_release_lock(lock);
    return result;
}


static CYTHON_INLINE void
__Pyx_INC_MEMVIEW($memviewslice_cname *memslice, int have_gil, int lineno)
{
    __pyx_nonatomic_int_type old_acquisition_count;
    struct $memview_objstruct_cname *memview = memslice->memview;
    if (unlikely(!memview || (PyObject *) memview == Py_None)) {
        // Allow uninitialized memoryview assignment and do not ref-count None.
        return;
    }

    old_acquisition_count = __pyx_add_acquisition_count(memview);
    if (unlikely(old_acquisition_count <= 0)) {
        if (likely(old_acquisition_count == 0)) {
            // First acquisition => keep the memoryview object alive.
            if (have_gil) {
                Py_INCREF((PyObject *) memview);
            } else {
                PyGILState_STATE _gilstate = PyGILState_Ensure();
                Py_INCREF((PyObject *) memview);
                PyGILState_Release(_gilstate);
            }
        } else {
            __pyx_fatalerror("Acquisition count is %d (line %d)",
                             old_acquisition_count+1, lineno);
        }
    }
}

static CYTHON_INLINE void __Pyx_XCLEAR_MEMVIEW($memviewslice_cname *memslice,
                                             int have_gil, int lineno) {
    __pyx_nonatomic_int_type old_acquisition_count;
    struct $memview_objstruct_cname *memview = memslice->memview;

    if (unlikely(!memview || (PyObject *) memview == Py_None)) {
        // Do not ref-count None.
        memslice->memview = NULL;
        return;
    }

    old_acquisition_count = __pyx_sub_acquisition_count(memview);
    memslice->data = NULL;
    if (likely(old_acquisition_count > 1)) {
        // Still other slices out there => we do not own the reference.
        memslice->memview = NULL;
    } else if (likely(old_acquisition_count == 1)) {
        // Last slice => discard owned Python reference to memoryview object.
        if (have_gil) {
            Py_CLEAR(memslice->memview);
        } else {
            PyGILState_STATE _gilstate = PyGILState_Ensure();
            Py_CLEAR(memslice->memview);
            PyGILState_Release(_gilstate);
        }
    } else {
        __pyx_fatalerror("Acquisition count is %d (line %d)",
                         old_acquisition_count-1, lineno);
    }
}


////////// MemviewSliceCopy.proto //////////
//@substitute: naming

static $memviewslice_cname
__pyx_memoryview_copy_new_contig(const $memviewslice_cname *from_mvs,
                                 const char *mode, int ndim,
                                 size_t sizeof_dtype, int contig_flag,
                                 int dtype_is_object);


////////// MemviewSliceCopy //////////
//@requires: MemviewSliceInit
//@substitute: naming

static $memviewslice_cname
__pyx_memoryview_copy_new_contig(const $memviewslice_cname *from_mvs,
                                 const char *mode, int ndim,
                                 size_t sizeof_dtype, int contig_flag,
                                 int dtype_is_object)
{
    __Pyx_RefNannyDeclarations
    int i;
    $memviewslice_cname new_mvs = __Pyx_MEMSLICE_INIT;
    struct $memview_objstruct_cname *from_memview = from_mvs->memview;
    Py_buffer *buf = &from_memview->view;
    PyObject *shape_tuple = NULL;
    struct __pyx_array_obj *array_obj = NULL;
    struct $memview_objstruct_cname *memview_obj = NULL;

    __Pyx_RefNannySetupContext("__pyx_memoryview_copy_new_contig", 0);

    for (i = 0; i < ndim; i++) {
        if (unlikely(from_mvs->suboffsets[i] >= 0)) {
            PyErr_Format(PyExc_ValueError, "Cannot copy memoryview slice with "
                                           "indirect dimensions (axis %d)", i);
            goto fail;
        }
    }

    shape_tuple = PyTuple_New(ndim);
    if (unlikely(!shape_tuple)) {
        goto fail;
    }
    __Pyx_GOTREF(shape_tuple);


    for(i = 0; i < ndim; i++) {
        PyObject *temp_int = PyLong_FromSsize_t(from_mvs->shape[i]);
        if(unlikely(!temp_int)) {
            goto fail;
        } else {
#if CYTHON_ASSUME_SAFE_MACROS
            PyTuple_SET_ITEM(shape_tuple, i, temp_int);
#else
            if (PyTuple_SetItem(shape_tuple, i, temp_int) < 0) {
                goto fail;
            }
#endif
        }
    }

    array_obj = __pyx_array_new(shape_tuple, sizeof_dtype, buf->format, mode, NULL);
    if (unlikely(!array_obj)) {
        goto fail;
    }
    __Pyx_GOTREF(array_obj);

    memview_obj = (struct $memview_objstruct_cname *) __pyx_memoryview_new(
                                    (PyObject *) array_obj, contig_flag,
                                    dtype_is_object,
                                    from_mvs->memview->typeinfo);
    if (unlikely(!memview_obj))
        goto fail;

    /* initialize new_mvs */
    if (unlikely(__Pyx_init_memviewslice(memview_obj, ndim, &new_mvs, 1) < 0))
        goto fail;

    if (unlikely(__pyx_memoryview_copy_contents(*from_mvs, new_mvs, ndim, ndim,
                                                dtype_is_object) < 0))
        goto fail;

    goto no_fail;

fail:
    __Pyx_XDECREF((PyObject *) new_mvs.memview);
    new_mvs.memview = NULL;
    new_mvs.data = NULL;
no_fail:
    __Pyx_XDECREF(shape_tuple);
    __Pyx_XDECREF((PyObject *) array_obj);
    __Pyx_RefNannyFinishContext();
    return new_mvs;
}


////////// CopyContentsUtility.proto /////////
//@requires: MemviewSliceCopy

#define {{func_cname}}(slice) \
        __pyx_memoryview_copy_new_contig(&slice, "{{mode}}", {{ndim}},            \
                                         sizeof({{dtype_decl}}), {{contig_flag}}, \
                                         {{dtype_is_object}})


////////// OverlappingSlices.proto //////////
//@substitute: naming

static int __pyx_slices_overlap($memviewslice_cname *slice1,
                                $memviewslice_cname *slice2,
                                int ndim, size_t itemsize);


////////// OverlappingSlices //////////
//@substitute: naming

/* Based on numpy's core/src/multiarray/array_assign.c */

/* Gets a half-open range [start, end) which contains the array data */
static void
__pyx_get_array_memory_extents($memviewslice_cname *slice,
                               void **out_start, void **out_end,
                               int ndim, size_t itemsize)
{
    char *start, *end;
    int i;

    start = end = slice->data;

    for (i = 0; i < ndim; i++) {
        Py_ssize_t stride = slice->strides[i];
        Py_ssize_t extent = slice->shape[i];

        if (extent == 0) {
            *out_start = *out_end = start;
            return;
        } else {
            if (stride > 0)
                end += stride * (extent - 1);
            else
                start += stride * (extent - 1);
        }
    }

    /* Return a half-open range */
    *out_start = start;
    *out_end = end + itemsize;
}

/* Returns 1 if the arrays have overlapping data, 0 otherwise */
static int
__pyx_slices_overlap($memviewslice_cname *slice1,
                     $memviewslice_cname *slice2,
                     int ndim, size_t itemsize)
{
    void *start1, *end1, *start2, *end2;

    __pyx_get_array_memory_extents(slice1, &start1, &end1, ndim, itemsize);
    __pyx_get_array_memory_extents(slice2, &start2, &end2, ndim, itemsize);

    return (start1 < end2) && (start2 < end1);
}


////////// MemviewSliceCheckContig.proto //////////

#define __pyx_memviewslice_is_contig_{{contig_type}}{{ndim}}(slice) \
    __pyx_memviewslice_is_contig(slice, '{{contig_type}}', {{ndim}})


////////// MemviewSliceIsContig.proto //////////
//@substitute: naming

static int __pyx_memviewslice_is_contig(const $memviewslice_cname mvs, char order, int ndim);/*proto*/


////////// MemviewSliceIsContig //////////
//@substitute: naming

static int
__pyx_memviewslice_is_contig(const $memviewslice_cname mvs, char order, int ndim)
{
    int i, index, step, start;
    Py_ssize_t itemsize = mvs.memview->view.itemsize;

    if (order == 'F') {
        step = 1;
        start = 0;
    } else {
        step = -1;
        start = ndim - 1;
    }

    for (i = 0; i < ndim; i++) {
        index = start + step * i;
        if (mvs.suboffsets[index] >= 0 || mvs.strides[index] != itemsize)
            return 0;

        itemsize *= mvs.shape[index];
    }

    return 1;
}


/////////////// MemviewSliceIndex.proto ///////////////

static CYTHON_INLINE char *
__pyx_memviewslice_index_full(const char *bufp, Py_ssize_t idx, Py_ssize_t stride, Py_ssize_t suboffset); /*proto*/


/////////////// MemviewSliceIndex ///////////////

static CYTHON_INLINE char *
__pyx_memviewslice_index_full(const char *bufp, Py_ssize_t idx,
                              Py_ssize_t stride, Py_ssize_t suboffset)
{
    bufp = bufp + idx * stride;
    if (suboffset >= 0) {
        bufp = *((char **) bufp) + suboffset;
    }
    return (char *) bufp;
}


/////////////// MemviewDtypeToObject.proto ///////////////

{{if to_py_function}}
static CYTHON_INLINE PyObject *{{get_function}}(const char *itemp); /* proto */
{{endif}}

{{if from_py_function}}
static CYTHON_INLINE int {{set_function}}(char *itemp, PyObject *obj); /* proto */
{{endif}}

/////////////// MemviewDtypeToObject ///////////////

{{#__pyx_memview_<dtype_name>_to_object}}

/* Convert a dtype to or from a Python object */

{{if to_py_function}}
static CYTHON_INLINE PyObject *{{get_function}}(const char *itemp) {
    return (PyObject *) {{to_py_function}}(*({{dtype}} {{'' if dtype_is_const else 'const'}} *) itemp);
}
{{endif}}

{{if from_py_function}}
static CYTHON_INLINE int {{set_function}}(char *itemp, PyObject *obj) {
    {{dtype}} value = {{from_py_function}}(obj);
    if (unlikely({{error_condition}}))
        return 0;
    *({{dtype}} *) itemp = value;
    return 1;
}
{{endif}}


/////////////// MemviewObjectToObject.proto ///////////////

/* Function callbacks (for memoryview object) for dtype object */
static PyObject *{{get_function}}(const char *itemp); /* proto */
static int {{set_function}}(const char *itemp, PyObject *obj); /* proto */


/////////////// MemviewObjectToObject ///////////////

static PyObject *{{get_function}}(const char *itemp) {
    PyObject *result = *(PyObject **) itemp;
    Py_INCREF(result);
    return result;
}

static int {{set_function}}(const char *itemp, PyObject *obj) {
    Py_INCREF(obj);
    Py_DECREF(*(PyObject **) itemp);
    *(PyObject **) itemp = obj;
    return 1;
}

/////////// ToughSlice //////////

/* Dimension is indexed with 'start:stop:step' */

if (unlikely(__pyx_memoryview_slice_memviewslice(
    &{{dst}},
    {{src}}.shape[{{dim}}], {{src}}.strides[{{dim}}], {{src}}.suboffsets[{{dim}}],
    {{dim}},
    {{new_ndim}},
    &{{get_suboffset_dim()}},
    {{start}},
    {{stop}},
    {{step}},
    {{int(have_start)}},
    {{int(have_stop)}},
    {{int(have_step)}},
    1) < 0))
{
    {{error_goto}}
}


////////// SimpleSlice //////////

/* Dimension is indexed with ':' only */

{{dst}}.shape[{{new_ndim}}] = {{src}}.shape[{{dim}}];
{{dst}}.strides[{{new_ndim}}] = {{src}}.strides[{{dim}}];

{{if access == 'direct'}}
    {{dst}}.suboffsets[{{new_ndim}}] = -1;
{{else}}
    {{dst}}.suboffsets[{{new_ndim}}] = {{src}}.suboffsets[{{dim}}];
    if ({{src}}.suboffsets[{{dim}}] >= 0)
        {{get_suboffset_dim()}} = {{new_ndim}};
{{endif}}


////////// SliceIndex //////////

// Dimension is indexed with an integer, we could use the ToughSlice
// approach, but this is faster

{
    Py_ssize_t __pyx_tmp_idx = {{idx}};

    {{if wraparound or boundscheck}}
        Py_ssize_t __pyx_tmp_shape = {{src}}.shape[{{dim}}];
    {{endif}}

    Py_ssize_t __pyx_tmp_stride = {{src}}.strides[{{dim}}];
    {{if wraparound}}
        if (__pyx_tmp_idx < 0)
            __pyx_tmp_idx += __pyx_tmp_shape;
    {{endif}}

    {{if boundscheck}}
        if (unlikely(!__Pyx_is_valid_index(__pyx_tmp_idx, __pyx_tmp_shape))) {
            {{if not have_gil}}
                PyGILState_STATE __pyx_gilstate_save = PyGILState_Ensure();
            {{endif}}

            PyErr_SetString(PyExc_IndexError,
                            "Index out of bounds (axis {{dim}})");

            {{if not have_gil}}
                PyGILState_Release(__pyx_gilstate_save);
            {{endif}}

            {{error_goto}}
        }
    {{endif}}

    {{if all_dimensions_direct}}
        {{dst}}.data += __pyx_tmp_idx * __pyx_tmp_stride;
    {{else}}
        if ({{get_suboffset_dim()}} < 0) {
            {{dst}}.data += __pyx_tmp_idx * __pyx_tmp_stride;

            /* This dimension is the first dimension, or is preceded by    */
            /* direct or indirect dimensions that are indexed away.        */
            /* Hence suboffset_dim must be less than zero, and we can have */
            /* our data pointer refer to another block by dereferencing.   */
            /*   slice.data -> B -> C     becomes     slice.data -> C      */

            {{if indirect}}
              {
                Py_ssize_t __pyx_tmp_suboffset = {{src}}.suboffsets[{{dim}}];

                {{if generic}}
                    if (__pyx_tmp_suboffset >= 0)
                {{endif}}

                    {{dst}}.data = *((char **) {{dst}}.data) + __pyx_tmp_suboffset;
              }
            {{endif}}

        } else {
            {{dst}}.suboffsets[{{get_suboffset_dim()}}] += __pyx_tmp_idx * __pyx_tmp_stride;

            /* Note: dimension can not be indirect, the compiler will have */
            /*       issued an error */
        }

    {{endif}}
}


////////// FillStrided1DScalar.proto //////////

static void
__pyx_fill_slice_{{dtype_name}}({{type_decl}} *p, Py_ssize_t extent, Py_ssize_t stride,
                                size_t itemsize, void *itemp);

////////// FillStrided1DScalar //////////

/* Fill a slice with a scalar value. The dimension is direct and strided or contiguous */
/* This can be used as a callback for the memoryview object to efficiently assign a scalar */
/* Currently unused */
static void
__pyx_fill_slice_{{dtype_name}}({{type_decl}} *p, Py_ssize_t extent, Py_ssize_t stride,
                                size_t itemsize, void *itemp)
{
    Py_ssize_t i;
    {{type_decl}} item = *(({{type_decl}} *) itemp);
    {{type_decl}} *endp;

    stride /= sizeof({{type_decl}});
    endp = p + stride * extent;

    while (p < endp) {
        *p = item;
        p += stride;
    }
}

/////////////////// SliceMemoryviewSlice.proto //////////////
//@substitute: naming

// Slicing in a single dimension of a memoryviewslice.

// Relies on being inlined to optimize out unused paths so
// not a good candidate for the shared library.
//
// Create a new slice dst given slice src.
//
//    dim             - the current src dimension (indexing will make dimensions
//                                                 disappear)
//    new_dim         - the new dst dimension
//    suboffset_dim   - pointer to a single int initialized to -1 to keep track of
//                      where slicing offsets should be added

static CYTHON_INLINE int __pyx_memoryview_slice_memviewslice(
        $memviewslice_cname *dst,
        Py_ssize_t shape, Py_ssize_t stride, Py_ssize_t suboffset,
        int dim, int new_ndim, int *suboffset_dim,
        Py_ssize_t start, Py_ssize_t stop, Py_ssize_t step,
        int have_start, int have_stop, int have_step,
        int is_slice); /* proto */

/////////////////// SliceMemoryviewSlice //////////////
//@substitute: naming

static void __pyx_memoryview_slice_memviewslice_err_dim(PyObject *error, const char* msg, int dim) {
    PyGILState_STATE gilstate = PyGILState_Ensure();
    PyErr_Format(error, msg, dim);
    PyGILState_Release(gilstate);
}

static CYTHON_INLINE int __pyx_memoryview_slice_memviewslice(
        $memviewslice_cname *dst,
        Py_ssize_t shape, Py_ssize_t stride, Py_ssize_t suboffset,
        int dim, int new_ndim, int *suboffset_dim,
        Py_ssize_t start, Py_ssize_t stop, Py_ssize_t step,
        int have_start, int have_stop, int have_step,
        int is_slice) {
    if (!is_slice) {
        // index is a normal integer-like index
        if (start < 0) {
            start += shape;
        }
        if (unlikely(!(0 <= start && start < shape))) {
            __pyx_memoryview_slice_memviewslice_err_dim(PyExc_IndexError, "Index out of bounds (axis %d)", dim);
            return -1;
        }
    } else {
        // index is a slice
        int negative_step;
        if (have_step) {
            negative_step = step < 0;
            if (unlikely(step == 0)) {
                __pyx_memoryview_slice_memviewslice_err_dim(PyExc_ValueError, "Step may not be zero (axis %d)", dim);
                return -1;
            }
        } else {
            negative_step = 0;
            step = 1;
        }

        // check our bounds and set defaults
        if (have_start) {
            if (start < 0) {
                start += shape;
                if (start < 0) {
                    start = 0;
                }
            } else if (start >= shape) {
                start = negative_step ? (shape - 1) : shape;
            }
        } else {
            start = negative_step ? (shape - 1) : 0;
        }
        if (have_stop) {
            if (stop < 0) {
                stop += shape;
                if (stop < 0) {
                    stop = 0;
                }
            } else if (stop > shape) {
                stop = shape;
            }
        } else {
            stop = negative_step ? -1 : shape;
        }

        // len = ceil( (stop - start) / step )
        Py_ssize_t new_shape = (stop - start) / step;
        if ((stop - start) - step * new_shape) {
            ++new_shape;
        }
        if (new_shape < 0) {
            new_shape = 0;
        }

        // shape/strides/suboffsets
        dst->strides[new_ndim] = stride * step;
        dst->shape[new_ndim] = new_shape;
        dst->suboffsets[new_ndim] = suboffset;
    }
    
    // Add the slicing or indexing offsets to the right suboffset or base data *
    if (suboffset_dim[0] < 0) {
        dst->data += start * stride;
    } else {
        dst->suboffsets[suboffset_dim[0]] += start * stride;
    }

    if (suboffset >= 0) {
        if (!is_slice) {
            if (likely(new_ndim == 0)) {
                dst->data = ((char **)(dst->data))[0] + suboffset;
            } else {
                __pyx_memoryview_slice_memviewslice_err_dim(
                    PyExc_IndexError,
                    "All dimensions preceding dimension %d must be indexed and not sliced",
                    dim);
                return -1;
            }
        } else {
            suboffset_dim[0] = new_ndim;
        }
    }
    return 0;
}
