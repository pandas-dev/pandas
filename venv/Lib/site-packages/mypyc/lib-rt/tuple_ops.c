// Tuple primitive operations
//
// These are registered in mypyc.primitives.tuple_ops.

#include <Python.h>
#include "CPy.h"

PyObject *CPySequenceTuple_GetItem(PyObject *tuple, CPyTagged index) {
    if (CPyTagged_CheckShort(index)) {
        Py_ssize_t n = CPyTagged_ShortAsSsize_t(index);
        Py_ssize_t size = PyTuple_GET_SIZE(tuple);
        if (n >= 0) {
            if (n >= size) {
                PyErr_SetString(PyExc_IndexError, "tuple index out of range");
                return NULL;
            }
        } else {
            n += size;
            if (n < 0) {
                PyErr_SetString(PyExc_IndexError, "tuple index out of range");
                return NULL;
            }
        }
        PyObject *result = PyTuple_GET_ITEM(tuple, n);
        Py_INCREF(result);
        return result;
    } else {
        PyErr_SetString(PyExc_OverflowError, CPYTHON_LARGE_INT_ERRMSG);
        return NULL;
    }
}

PyObject *CPySequenceTuple_GetSlice(PyObject *obj, CPyTagged start, CPyTagged end) {
    if (likely(PyTuple_CheckExact(obj)
               && CPyTagged_CheckShort(start) && CPyTagged_CheckShort(end))) {
        Py_ssize_t startn = CPyTagged_ShortAsSsize_t(start);
        Py_ssize_t endn = CPyTagged_ShortAsSsize_t(end);
        if (startn < 0) {
            startn += PyTuple_GET_SIZE(obj);
        }
        if (endn < 0) {
            endn += PyTuple_GET_SIZE(obj);
        }
        return PyTuple_GetSlice(obj, startn, endn);
    }
    return CPyObject_GetSlice(obj, start, end);
}

// No error checking
PyObject *CPySequenceTuple_GetItemUnsafe(PyObject *tuple, Py_ssize_t index)
{
    PyObject *result = PyTuple_GET_ITEM(tuple, index);
    Py_INCREF(result);
    return result;
}

// PyTuple_SET_ITEM does no error checking,
// and should only be used to fill in brand new tuples.
void CPySequenceTuple_SetItemUnsafe(PyObject *tuple, Py_ssize_t index, PyObject *value)
{
    PyTuple_SET_ITEM(tuple, index, value);
}
