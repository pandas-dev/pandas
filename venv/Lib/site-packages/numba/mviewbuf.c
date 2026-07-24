#include "_pymodule.h"

static int get_writable_buffer(PyObject* obj, Py_buffer *buf, int force)
{
    Py_buffer read_buf;
    int flags = PyBUF_ND|PyBUF_STRIDES|PyBUF_FORMAT;
    int ret;

    /* Attempt to get a writable buffer */
    if (!PyObject_GetBuffer(obj, buf, flags|PyBUF_WRITABLE))
        return 0;
    if (!force)
        return -1;

    /* Make a writable buffer from a read-only buffer */
    PyErr_Clear();
    if(-1 == PyObject_GetBuffer(obj, &read_buf, flags))
        return -1;
    ret = PyBuffer_FillInfo(buf, NULL, read_buf.buf, read_buf.len, 0,
                             flags|PyBUF_WRITABLE);
    PyBuffer_Release(&read_buf);
    return ret;
}

static int get_readonly_buffer(PyObject* obj, Py_buffer *buf)
{
    int flags = PyBUF_ND|PyBUF_STRIDES|PyBUF_FORMAT;

    return PyObject_GetBuffer(obj, buf, flags);
}


static void free_buffer(Py_buffer * buf)
{
    PyBuffer_Release(buf);
}

/**
 * Return a pointer to the data of a writable buffer from obj. If only a
 * read-only buffer is available and force is True, a read-write buffer based on
 * the read-only buffer is obtained. Note that this may have some surprising
 * effects on buffers which expect the data from their read-only buffer not to
 * be modified.
 */
static PyObject*
memoryview_get_buffer(PyObject *self, PyObject *args){
    PyObject *obj = NULL;
    int force = 0;
    int readonly = 0;
    PyObject *ret = NULL;
    Py_buffer buf;

    if (!PyArg_ParseTuple(args, "O|ii", &obj, &force, &readonly))
        return NULL;

    if (readonly) {
        if (get_readonly_buffer(obj, &buf))
            return NULL;
    } else {
        if (get_writable_buffer(obj, &buf, force))
            return NULL;
    }

    ret = PyLong_FromVoidPtr(buf.buf);
    free_buffer(&buf);
    return ret;
}

/**
 * Gets a half-open range [start, end) which contains the array data
 * Modified from numpy/core/src/multiarray/array_assign.c
 */
static PyObject*
get_extents(Py_ssize_t *shape, Py_ssize_t *strides, int ndim,
            Py_ssize_t itemsize, Py_ssize_t ptr)
{
    Py_ssize_t start, end;
    int idim;
    Py_ssize_t *dimensions = shape;
    PyObject *ret = NULL;

    if (ndim < 0 ){
        PyErr_SetString(PyExc_ValueError, "buffer ndim < 0");
        return NULL;
    }

    if (!dimensions) {
        if (ndim == 0) {
            start = end = ptr;
            end += itemsize;
            return Py_BuildValue("nn", start, end);
        }
        PyErr_SetString(PyExc_ValueError, "buffer shape is not defined");
        return NULL;
    }

    if (!strides) {
        PyErr_SetString(PyExc_ValueError, "buffer strides is not defined");
        return NULL;
    }

    /* Calculate with a closed range [start, end] */
    start = end = ptr;
    for (idim = 0; idim < ndim; ++idim) {
        Py_ssize_t stride = strides[idim], dim = dimensions[idim];
        /* If the array size is zero, return an empty range */
        if (dim == 0) {
            start = end = ptr;
            ret = Py_BuildValue("nn", start, end);
            break;
        }
        /* Expand either upwards or downwards depending on stride */
        else {
            if (stride > 0) {
                end += stride * (dim - 1);
            }
            else if (stride < 0) {
                start += stride * (dim - 1);
            }
        }
    }

    if (!ret) {
        /* Return a half-open range */
        Py_ssize_t out_start = start;
        Py_ssize_t out_end = end + itemsize;

        ret = Py_BuildValue("nn", out_start, out_end);
    }

    return ret;
}

static PyObject*
memoryview_get_extents(PyObject *self, PyObject *args)
{
    PyObject *obj = NULL;
    PyObject *ret = NULL;
    Py_buffer b;
    if (!PyArg_ParseTuple(args, "O", &obj))
        return NULL;

    if (get_readonly_buffer(obj, &b))
        return NULL;

    ret = get_extents(b.shape, b.strides, b.ndim, b.itemsize,
                      (Py_ssize_t)b.buf);
    free_buffer(&b);
    return ret;
}

static PyObject*
memoryview_get_extents_info(PyObject *self, PyObject *args)
{
    int i;
    Py_ssize_t *shape_ary = NULL;
    Py_ssize_t *strides_ary = NULL;
    PyObject *shape_tuple = NULL;
    PyObject *strides_tuple = NULL;
    PyObject *shape = NULL, *strides = NULL;
    Py_ssize_t itemsize = 0;
    int ndim = 0;
    PyObject* res = NULL;

    if (!PyArg_ParseTuple(args, "OOin", &shape, &strides, &ndim, &itemsize))
        goto cleanup;

    if (ndim < 0) {
        PyErr_SetString(PyExc_ValueError, "ndim is negative");
        goto cleanup;
    }

    if (itemsize <= 0) {
        PyErr_SetString(PyExc_ValueError, "ndim <= 0");
        goto cleanup;
    }

    shape_ary = malloc(sizeof(Py_ssize_t) * ndim + 1);
    strides_ary = malloc(sizeof(Py_ssize_t) * ndim + 1);

    shape_tuple = PySequence_Fast(shape, "shape is not a sequence");
    if (!shape_tuple) goto cleanup;

    for (i = 0; i < ndim; ++i) {
        shape_ary[i] = PyNumber_AsSsize_t(
                           PySequence_Fast_GET_ITEM(shape_tuple, i),
                           PyExc_OverflowError);
    }

    strides_tuple = PySequence_Fast(strides, "strides is not a sequence");
    if (!strides_tuple) goto cleanup;

    for (i = 0; i < ndim; ++i) {
        strides_ary[i] = PyNumber_AsSsize_t(
                           PySequence_Fast_GET_ITEM(strides_tuple, i),
                           PyExc_OverflowError);
    }

    res = get_extents(shape_ary, strides_ary, ndim, itemsize, 0);
cleanup:
    free(shape_ary);
    free(strides_ary);
    Py_XDECREF(shape_tuple);
    Py_XDECREF(strides_tuple);
    return res;
}


/* new type to expose buffer interface */
typedef struct {
    PyObject_HEAD
    /* Type-specific fields go here. */
} MemAllocObject;


static int
get_bufinfo(PyObject *self, Py_ssize_t *psize, void **pptr)
{
    PyObject *buflen = NULL;
    PyObject *bufptr = NULL;
    Py_ssize_t size = 0;
    void* ptr = NULL;
    int ret = -1;

    buflen = PyObject_GetAttrString(self, "_buflen_");
    if (!buflen) goto cleanup;

    bufptr = PyObject_GetAttrString(self, "_bufptr_");
    if (!bufptr) goto cleanup;

    size = PyNumber_AsSsize_t(buflen, PyExc_OverflowError);
    if (size == -1 && PyErr_Occurred()) goto cleanup;
    else if (size < 0) {
        PyErr_SetString(PyExc_ValueError, "negative buffer size");
        goto cleanup;
    }

    ptr = PyLong_AsVoidPtr(PyNumber_Long(bufptr));
    if (PyErr_Occurred())
        goto cleanup;
    else if (!ptr) {
        PyErr_SetString(PyExc_ValueError, "null buffer pointer");
        goto cleanup;
    }

    *psize = size;
    *pptr = ptr;
    ret = 0;
cleanup:
    Py_XDECREF(buflen);
    Py_XDECREF(bufptr);
    return ret;
}


static int
MemAllocObject_getbuffer(PyObject *self, Py_buffer *view, int flags)
{
    Py_ssize_t size = 0;
    void *ptr = 0;
    int readonly;

    if(-1 == get_bufinfo(self, &size, &ptr))
        return -1;

    readonly = (PyBUF_WRITABLE & flags) != PyBUF_WRITABLE;

    /* fill buffer */
    if (-1 == PyBuffer_FillInfo(view, self, (void*)ptr, size, readonly, flags))
        return -1;

    return 0;
}

static void
MemAllocObject_releasebuffer(PyObject *self, Py_buffer *view)
{
    /* Do nothing */
}

static PyBufferProcs MemAlloc_as_buffer = {
    MemAllocObject_getbuffer,
    MemAllocObject_releasebuffer,
};


static PyTypeObject MemAllocType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "mviewbuf.MemAlloc",                        /* tp_name */
    sizeof(MemAllocObject),                     /* tp_basicsize */
    0,                                          /* tp_itemsize */
    0,                                          /* tp_dealloc */
    0,                                          /* tp_vectorcall_offset */
    0,                                          /* tp_getattr */
    0,                                          /* tp_setattr */
    0,                                          /* tp_as_async */
    0,                                          /* tp_repr */
    0,                                          /* tp_as_number */
    0,                                          /* tp_as_sequence */
    0,                                          /* tp_as_mapping */
    0,                                          /* tp_hash */
    0,                                          /* tp_call */
    0,                                          /* tp_str */
    0,                                          /* tp_getattro */
    0,                                          /* tp_setattro */
    &MemAlloc_as_buffer,                        /* tp_as_buffer */
    (Py_TPFLAGS_DEFAULT| Py_TPFLAGS_BASETYPE),  /* tp_flags */
    0,                                          /* tp_doc */
    0,                                          /* tp_traverse */
    0,                                          /* tp_clear */
    0,                                          /* tp_richcompare */
    0,                                          /* tp_weaklistoffset */
    0,                                          /* tp_iter */
    0,                                          /* tp_iternext */
    0,                                          /* tp_methods */
    0,                                          /* tp_members */
    0,                                          /* tp_getset */
    0,                                          /* tp_base */
    0,                                          /* tp_dict */
    0,                                          /* tp_descr_get */
    0,                                          /* tp_descr_set */
    0,                                          /* tp_dictoffset */
    0,                                          /* tp_init */
    0,                                          /* tp_alloc */
    0,                                          /* tp_new */
    0,                                          /* tp_free */
    0,                                          /* tp_is_gc */
    0,                                          /* tp_bases */
    0,                                          /* tp_mro */
    0,                                          /* tp_cache */
    0,                                          /* tp_subclasses */
    0,                                          /* tp_weaklist */
    0,                                          /* tp_del */
    0,                                          /* tp_version_tag */
    0,                                          /* tp_finalize */
    0,                                          /* tp_vectorcall */
#if (PY_MAJOR_VERSION == 3) && (PY_MINOR_VERSION >= 12)
/* This was introduced first in 3.12
 * https://github.com/python/cpython/issues/91051
 */
    0,                                           /* tp_watched */
#endif
#if (PY_MAJOR_VERSION == 3) && (PY_MINOR_VERSION >= 13)
/* This was introduced in 3.13
 * https://github.com/python/cpython/pull/114900
 */
    0,                                           /* tp_versions_used */
#endif

/* WARNING: Do not remove this, only modify it! It is a version guard to
 * act as a reminder to update this struct on Python version update! */
#if (PY_MAJOR_VERSION == 3)
#if ! (NB_SUPPORTED_PYTHON_MINOR)
#error "Python minor version is not supported."
#endif
#else
#error "Python major version is not supported."
#endif
/* END WARNING*/
};


static PyMethodDef core_methods[] = {
#define declmethod(func) { #func , ( PyCFunction )func , METH_VARARGS , NULL }
    declmethod(memoryview_get_buffer),
    declmethod(memoryview_get_extents),
    declmethod(memoryview_get_extents_info),
    { NULL },
#undef declmethod
};


MOD_INIT(mviewbuf) {
    PyObject *module;
    MOD_DEF(module, "mviewbuf", "No docs", core_methods)
    if (module == NULL)
        return MOD_ERROR_VAL;

    MemAllocType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&MemAllocType) < 0){
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&MemAllocType);
    PyModule_AddObject(module, "MemAlloc", (PyObject*)&MemAllocType);

    return MOD_SUCCESS_VAL(module);
}

