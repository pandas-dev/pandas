// Generic primitive operations
//
// These are registered in mypyc.primitives.generic_ops.

#include <Python.h>
#include "CPy.h"

CPyTagged CPyObject_Hash(PyObject *o) {
    Py_hash_t h = PyObject_Hash(o);
    if (h == -1) {
        return CPY_INT_TAG;
    } else {
        // This is tragically annoying. The range of hash values in
        // 64-bit python covers 64-bits, and our short integers only
        // cover 63. This means that half the time we are boxing the
        // result for basically no good reason. To add insult to
        // injury it is probably about to be immediately unboxed by a
        // tp_hash wrapper.
        return CPyTagged_FromSsize_t(h);
    }
}

PyObject *CPyObject_GetAttr3(PyObject *v, PyObject *name, PyObject *defl)
{
    PyObject *result = PyObject_GetAttr(v, name);
    if (!result && PyErr_ExceptionMatches(PyExc_AttributeError)) {
        PyErr_Clear();
        Py_INCREF(defl);
        result = defl;
    }
    return result;
}

PyObject *CPyIter_Next(PyObject *iter)
{
    return (*Py_TYPE(iter)->tp_iternext)(iter);
}

PyObject *CPyNumber_Power(PyObject *base, PyObject *index)
{
    return PyNumber_Power(base, index, Py_None);
}

PyObject *CPyNumber_InPlacePower(PyObject *base, PyObject *index)
{
    return PyNumber_InPlacePower(base, index, Py_None);
}

PyObject *CPyObject_GetSlice(PyObject *obj, CPyTagged start, CPyTagged end) {
    PyObject *start_obj = CPyTagged_AsObject(start);
    PyObject *end_obj = CPyTagged_AsObject(end);
    if (unlikely(start_obj == NULL || end_obj == NULL)) {
        return NULL;
    }
    PyObject *slice = PySlice_New(start_obj, end_obj, NULL);
    Py_DECREF(start_obj);
    Py_DECREF(end_obj);
    if (unlikely(slice == NULL)) {
        return NULL;
    }
    PyObject *result = PyObject_GetItem(obj, slice);
    Py_DECREF(slice);
    return result;
}

typedef PyObject *(*SetupFunction)(PyObject *);

PyObject *CPy_SetupObject(PyObject *type) {
    PyTypeObject *tp = (PyTypeObject *)type;
    PyMethodDef *def = NULL;
    for(; tp; tp = tp->tp_base) {
        def = tp->tp_methods;
        if (!def || !def->ml_name) {
            continue;
        }

        if (!strcmp(def->ml_name, "__internal_mypyc_setup")) {
            return ((SetupFunction)(void(*)(void))def->ml_meth)(type);
        }
    }

    PyErr_SetString(PyExc_RuntimeError, "Internal mypyc error: Unable to find object setup function");
    return NULL;
}
