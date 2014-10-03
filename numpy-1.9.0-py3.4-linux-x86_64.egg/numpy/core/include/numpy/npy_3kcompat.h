/*
 * This is a convenience header file providing compatibility utilities
 * for supporting Python 2 and Python 3 in the same code base.
 *
 * If you want to use this for your own projects, it's recommended to make a
 * copy of it. Although the stuff below is unlikely to change, we don't provide
 * strong backwards compatibility guarantees at the moment.
 */

#ifndef _NPY_3KCOMPAT_H_
#define _NPY_3KCOMPAT_H_

#include <Python.h>
#include <stdio.h>

#if PY_VERSION_HEX >= 0x03000000
#ifndef NPY_PY3K
#define NPY_PY3K 1
#endif
#endif

#include "numpy/npy_common.h"
#include "numpy/ndarrayobject.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * PyInt -> PyLong
 */

#if defined(NPY_PY3K)
/* Return True only if the long fits in a C long */
static NPY_INLINE int PyInt_Check(PyObject *op) {
    int overflow = 0;
    if (!PyLong_Check(op)) {
        return 0;
    }
    PyLong_AsLongAndOverflow(op, &overflow);
    return (overflow == 0);
}

#define PyInt_FromLong PyLong_FromLong
#define PyInt_AsLong PyLong_AsLong
#define PyInt_AS_LONG PyLong_AsLong
#define PyInt_AsSsize_t PyLong_AsSsize_t

/* NOTE:
 *
 * Since the PyLong type is very different from the fixed-range PyInt,
 * we don't define PyInt_Type -> PyLong_Type.
 */
#endif /* NPY_PY3K */

/*
 * PyString -> PyBytes
 */

#if defined(NPY_PY3K)

#define PyString_Type PyBytes_Type
#define PyString_Check PyBytes_Check
#define PyStringObject PyBytesObject
#define PyString_FromString PyBytes_FromString
#define PyString_FromStringAndSize PyBytes_FromStringAndSize
#define PyString_AS_STRING PyBytes_AS_STRING
#define PyString_AsStringAndSize PyBytes_AsStringAndSize
#define PyString_FromFormat PyBytes_FromFormat
#define PyString_Concat PyBytes_Concat
#define PyString_ConcatAndDel PyBytes_ConcatAndDel
#define PyString_AsString PyBytes_AsString
#define PyString_GET_SIZE PyBytes_GET_SIZE
#define PyString_Size PyBytes_Size

#define PyUString_Type PyUnicode_Type
#define PyUString_Check PyUnicode_Check
#define PyUStringObject PyUnicodeObject
#define PyUString_FromString PyUnicode_FromString
#define PyUString_FromStringAndSize PyUnicode_FromStringAndSize
#define PyUString_FromFormat PyUnicode_FromFormat
#define PyUString_Concat PyUnicode_Concat2
#define PyUString_ConcatAndDel PyUnicode_ConcatAndDel
#define PyUString_GET_SIZE PyUnicode_GET_SIZE
#define PyUString_Size PyUnicode_Size
#define PyUString_InternFromString PyUnicode_InternFromString
#define PyUString_Format PyUnicode_Format

#else

#define PyBytes_Type PyString_Type
#define PyBytes_Check PyString_Check
#define PyBytesObject PyStringObject
#define PyBytes_FromString PyString_FromString
#define PyBytes_FromStringAndSize PyString_FromStringAndSize
#define PyBytes_AS_STRING PyString_AS_STRING
#define PyBytes_AsStringAndSize PyString_AsStringAndSize
#define PyBytes_FromFormat PyString_FromFormat
#define PyBytes_Concat PyString_Concat
#define PyBytes_ConcatAndDel PyString_ConcatAndDel
#define PyBytes_AsString PyString_AsString
#define PyBytes_GET_SIZE PyString_GET_SIZE
#define PyBytes_Size PyString_Size

#define PyUString_Type PyString_Type
#define PyUString_Check PyString_Check
#define PyUStringObject PyStringObject
#define PyUString_FromString PyString_FromString
#define PyUString_FromStringAndSize PyString_FromStringAndSize
#define PyUString_FromFormat PyString_FromFormat
#define PyUString_Concat PyString_Concat
#define PyUString_ConcatAndDel PyString_ConcatAndDel
#define PyUString_GET_SIZE PyString_GET_SIZE
#define PyUString_Size PyString_Size
#define PyUString_InternFromString PyString_InternFromString
#define PyUString_Format PyString_Format

#endif /* NPY_PY3K */


static NPY_INLINE void
PyUnicode_ConcatAndDel(PyObject **left, PyObject *right)
{
    PyObject *newobj;
    newobj = PyUnicode_Concat(*left, right);
    Py_DECREF(*left);
    Py_DECREF(right);
    *left = newobj;
}

static NPY_INLINE void
PyUnicode_Concat2(PyObject **left, PyObject *right)
{
    PyObject *newobj;
    newobj = PyUnicode_Concat(*left, right);
    Py_DECREF(*left);
    *left = newobj;
}

/*
 * PyFile_* compatibility
 */
#if defined(NPY_PY3K)
/*
 * Get a FILE* handle to the file represented by the Python object
 */
static NPY_INLINE FILE*
npy_PyFile_Dup2(PyObject *file, char *mode, npy_off_t *orig_pos)
{
    int fd, fd2;
    PyObject *ret, *os;
    npy_off_t pos;
    FILE *handle;

    /* Flush first to ensure things end up in the file in the correct order */
    ret = PyObject_CallMethod(file, "flush", "");
    if (ret == NULL) {
        return NULL;
    }
    Py_DECREF(ret);
    fd = PyObject_AsFileDescriptor(file);
    if (fd == -1) {
        return NULL;
    }

    /*
     * The handle needs to be dup'd because we have to call fclose
     * at the end
     */
    os = PyImport_ImportModule("os");
    if (os == NULL) {
        return NULL;
    }
    ret = PyObject_CallMethod(os, "dup", "i", fd);
    Py_DECREF(os);
    if (ret == NULL) {
        return NULL;
    }
    fd2 = PyNumber_AsSsize_t(ret, NULL);
    Py_DECREF(ret);

    /* Convert to FILE* handle */
#ifdef _WIN32
    handle = _fdopen(fd2, mode);
#else
    handle = fdopen(fd2, mode);
#endif
    if (handle == NULL) {
        PyErr_SetString(PyExc_IOError,
                        "Getting a FILE* from a Python file object failed");
    }

    /* Record the original raw file handle position */
    *orig_pos = npy_ftell(handle);
    if (*orig_pos == -1) {
        PyErr_SetString(PyExc_IOError, "obtaining file position failed");
        fclose(handle);
        return NULL;
    }

    /* Seek raw handle to the Python-side position */
    ret = PyObject_CallMethod(file, "tell", "");
    if (ret == NULL) {
        fclose(handle);
        return NULL;
    }
    pos = PyLong_AsLongLong(ret);
    Py_DECREF(ret);
    if (PyErr_Occurred()) {
        fclose(handle);
        return NULL;
    }
    if (npy_fseek(handle, pos, SEEK_SET) == -1) {
        PyErr_SetString(PyExc_IOError, "seeking file failed");
        fclose(handle);
        return NULL;
    }
    return handle;
}

/*
 * Close the dup-ed file handle, and seek the Python one to the current position
 */
static NPY_INLINE int
npy_PyFile_DupClose2(PyObject *file, FILE* handle, npy_off_t orig_pos)
{
    int fd;
    PyObject *ret;
    npy_off_t position;

    position = npy_ftell(handle);

    /* Close the FILE* handle */
    fclose(handle);

    /*
     * Restore original file handle position, in order to not confuse
     * Python-side data structures
     */
    fd = PyObject_AsFileDescriptor(file);
    if (fd == -1) {
        return -1;
    }
    if (npy_lseek(fd, orig_pos, SEEK_SET) == -1) {
        PyErr_SetString(PyExc_IOError, "seeking file failed");
        return -1;
    }

    if (position == -1) {
        PyErr_SetString(PyExc_IOError, "obtaining file position failed");
        return -1;
    }

    /* Seek Python-side handle to the FILE* handle position */
    ret = PyObject_CallMethod(file, "seek", NPY_OFF_T_PYFMT "i", position, 0);
    if (ret == NULL) {
        return -1;
    }
    Py_DECREF(ret);
    return 0;
}

static NPY_INLINE int
npy_PyFile_Check(PyObject *file)
{
    int fd;
    fd = PyObject_AsFileDescriptor(file);
    if (fd == -1) {
        PyErr_Clear();
        return 0;
    }
    return 1;
}

/*
 * DEPRECATED DO NOT USE
 * use npy_PyFile_DupClose2 instead
 * this function will mess ups python3 internal file object buffering
 * Close the dup-ed file handle, and seek the Python one to the current position
 */
static NPY_INLINE int
npy_PyFile_DupClose(PyObject *file, FILE* handle)
{
    PyObject *ret;
    Py_ssize_t position;
    position = npy_ftell(handle);
    fclose(handle);

    ret = PyObject_CallMethod(file, "seek", NPY_SSIZE_T_PYFMT "i", position, 0);
    if (ret == NULL) {
        return -1;
    }
    Py_DECREF(ret);
    return 0;
}


#else

/* DEPRECATED, DO NOT USE */
#define npy_PyFile_DupClose(f, h, p) npy_PyFile_DupClose2((f), (h), (p))

/* use these */
static NPY_INLINE FILE *
npy_PyFile_Dup2(PyObject *file,
                const char *NPY_UNUSED(mode), npy_off_t *NPY_UNUSED(orig_pos))
{
    return PyFile_AsFile(file);
}

static NPY_INLINE int
npy_PyFile_DupClose2(PyObject *NPY_UNUSED(file), FILE* NPY_UNUSED(handle),
                     npy_off_t NPY_UNUSED(orig_pos))
{
    return 0;
}

#define npy_PyFile_Check PyFile_Check

#endif

/*
 * DEPRECATED, DO NOT USE
 * Use npy_PyFile_Dup2 instead.
 * This function will mess up python3 internal file object buffering.
 * Get a FILE* handle to the file represented by the Python object.
 */
static NPY_INLINE FILE*
npy_PyFile_Dup(PyObject *file, char *mode)
{
    npy_off_t orig;
    if (DEPRECATE("npy_PyFile_Dup is deprecated, use "
                  "npy_PyFile_Dup2") < 0) {
        return NULL;
    }

    return npy_PyFile_Dup2(file, mode, &orig);
}

static NPY_INLINE PyObject*
npy_PyFile_OpenFile(PyObject *filename, const char *mode)
{
    PyObject *open;
    open = PyDict_GetItemString(PyEval_GetBuiltins(), "open");
    if (open == NULL) {
        return NULL;
    }
    return PyObject_CallFunction(open, "Os", filename, mode);
}

static NPY_INLINE int
npy_PyFile_CloseFile(PyObject *file)
{
    PyObject *ret;

    ret = PyObject_CallMethod(file, "close", NULL);
    if (ret == NULL) {
        return -1;
    }
    Py_DECREF(ret);
    return 0;
}

/*
 * PyObject_Cmp
 */
#if defined(NPY_PY3K)
static NPY_INLINE int
PyObject_Cmp(PyObject *i1, PyObject *i2, int *cmp)
{
    int v;
    v = PyObject_RichCompareBool(i1, i2, Py_LT);
    if (v == 0) {
        *cmp = -1;
        return 1;
    }
    else if (v == -1) {
        return -1;
    }

    v = PyObject_RichCompareBool(i1, i2, Py_GT);
    if (v == 0) {
        *cmp = 1;
        return 1;
    }
    else if (v == -1) {
        return -1;
    }

    v = PyObject_RichCompareBool(i1, i2, Py_EQ);
    if (v == 0) {
        *cmp = 0;
        return 1;
    }
    else {
        *cmp = 0;
        return -1;
    }
}
#endif

/*
 * PyCObject functions adapted to PyCapsules.
 *
 * The main job here is to get rid of the improved error handling
 * of PyCapsules. It's a shame...
 */
#if PY_VERSION_HEX >= 0x03000000

static NPY_INLINE PyObject *
NpyCapsule_FromVoidPtr(void *ptr, void (*dtor)(PyObject *))
{
    PyObject *ret = PyCapsule_New(ptr, NULL, dtor);
    if (ret == NULL) {
        PyErr_Clear();
    }
    return ret;
}

static NPY_INLINE PyObject *
NpyCapsule_FromVoidPtrAndDesc(void *ptr, void* context, void (*dtor)(PyObject *))
{
    PyObject *ret = NpyCapsule_FromVoidPtr(ptr, dtor);
    if (ret != NULL && PyCapsule_SetContext(ret, context) != 0) {
        PyErr_Clear();
        Py_DECREF(ret);
        ret = NULL;
    }
    return ret;
}

static NPY_INLINE void *
NpyCapsule_AsVoidPtr(PyObject *obj)
{
    void *ret = PyCapsule_GetPointer(obj, NULL);
    if (ret == NULL) {
        PyErr_Clear();
    }
    return ret;
}

static NPY_INLINE void *
NpyCapsule_GetDesc(PyObject *obj)
{
    return PyCapsule_GetContext(obj);
}

static NPY_INLINE int
NpyCapsule_Check(PyObject *ptr)
{
    return PyCapsule_CheckExact(ptr);
}

#else

static NPY_INLINE PyObject *
NpyCapsule_FromVoidPtr(void *ptr, void (*dtor)(void *))
{
    return PyCObject_FromVoidPtr(ptr, dtor);
}

static NPY_INLINE PyObject *
NpyCapsule_FromVoidPtrAndDesc(void *ptr, void* context,
        void (*dtor)(void *, void *))
{
    return PyCObject_FromVoidPtrAndDesc(ptr, context, dtor);
}

static NPY_INLINE void *
NpyCapsule_AsVoidPtr(PyObject *ptr)
{
    return PyCObject_AsVoidPtr(ptr);
}

static NPY_INLINE void *
NpyCapsule_GetDesc(PyObject *obj)
{
    return PyCObject_GetDesc(obj);
}

static NPY_INLINE int
NpyCapsule_Check(PyObject *ptr)
{
    return PyCObject_Check(ptr);
}

#endif

/*
 * Hash value compatibility.
 * As of Python 3.2 hash values are of type Py_hash_t.
 * Previous versions use C long.
 */
#if PY_VERSION_HEX < 0x03020000
typedef long npy_hash_t;
#define NPY_SIZEOF_HASH_T NPY_SIZEOF_LONG
#else
typedef Py_hash_t npy_hash_t;
#define NPY_SIZEOF_HASH_T NPY_SIZEOF_INTP
#endif

#ifdef __cplusplus
}
#endif

#endif /* _NPY_3KCOMPAT_H_ */
