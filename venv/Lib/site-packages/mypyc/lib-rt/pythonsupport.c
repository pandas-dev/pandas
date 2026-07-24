// Collects code that was copied in from cpython, for a couple of different reasons:
//  * We wanted to modify it to produce a more efficient version for our uses
//  * We needed to call it and it was static :(
//  * We wanted to call it and needed to backport it

#include "pythonsupport.h"

/////////////////////////////////////////
// Adapted from bltinmodule.c in Python 3.7.0
PyObject*
update_bases(PyObject *bases)
{
    Py_ssize_t i, j;
    PyObject *base, *meth, *new_base, *result, *new_bases = NULL;
    PyObject *stack[1] = {bases};
    assert(PyTuple_Check(bases));

    Py_ssize_t nargs = PyTuple_GET_SIZE(bases);
    for (i = 0; i < nargs; i++) {
        base = PyTuple_GET_ITEM(bases, i);
        if (PyType_Check(base)) {
            if (new_bases) {
                /* If we already have made a replacement, then we append every normal base,
                   otherwise just skip it. */
                if (PyList_Append(new_bases, base) < 0) {
                    goto error;
                }
            }
            continue;
        }
        if (PyObject_GetOptionalAttr(base, mypyc_interned_str.__mro_entries__, &meth) < 0) {
            goto error;
        }
        if (!meth) {
            if (new_bases) {
                if (PyList_Append(new_bases, base) < 0) {
                    goto error;
                }
            }
            continue;
        }
        new_base = PyObject_Vectorcall(meth, stack, 1, NULL);
        Py_DECREF(meth);
        if (!new_base) {
            goto error;
        }
        if (!PyTuple_Check(new_base)) {
            PyErr_SetString(PyExc_TypeError,
                            "__mro_entries__ must return a tuple");
            Py_DECREF(new_base);
            goto error;
        }
        if (!new_bases) {
            /* If this is a first successful replacement, create new_bases list and
               copy previously encountered bases. */
            if (!(new_bases = PyList_New(i))) {
                goto error;
            }
            for (j = 0; j < i; j++) {
                base = PyTuple_GET_ITEM(bases, j);
                PyList_SET_ITEM(new_bases, j, base);
                Py_INCREF(base);
            }
        }
        j = PyList_GET_SIZE(new_bases);
        if (PyList_SetSlice(new_bases, j, j, new_base) < 0) {
            goto error;
        }
        Py_DECREF(new_base);
    }
    if (!new_bases) {
        return bases;
    }
    result = PyList_AsTuple(new_bases);
    Py_DECREF(new_bases);
    return result;

error:
    Py_XDECREF(new_bases);
    return NULL;
}

// From Python 3.7's typeobject.c
int
init_subclass(PyTypeObject *type, PyObject *kwds)
{
    PyObject *super, *func, *result;
    PyObject *args[2] = {(PyObject *)type, (PyObject *)type};

    super = PyObject_Vectorcall((PyObject *)&PySuper_Type, args, 2, NULL);
    if (super == NULL) {
        return -1;
    }

    func = PyObject_GetAttr(super, mypyc_interned_str.__init_subclass__);
    Py_DECREF(super);
    if (func == NULL) {
        return -1;
    }

    result = _PyObject_FastCallDict(func, NULL, 0, kwds);
    Py_DECREF(func);
    if (result == NULL) {
        return -1;
    }

    Py_DECREF(result);
    return 0;
}

#if CPY_3_12_FEATURES

// Slow path of CPyLong_AsSsize_tAndOverflow (non-inlined)
Py_ssize_t
CPyLong_AsSsize_tAndOverflow_(PyObject *vv, int *overflow)
{
    PyLongObject *v = (PyLongObject *)vv;
    size_t x, prev;
    Py_ssize_t res;
    Py_ssize_t i;
    int sign;

    *overflow = 0;

    res = -1;
    i = CPY_LONG_TAG(v);

    sign = 1;
    x = 0;
    if (i & CPY_SIGN_NEGATIVE) {
        sign = -1;
    }
    i >>= CPY_NON_SIZE_BITS;
    while (--i >= 0) {
        prev = x;
        x = (x << PyLong_SHIFT) + CPY_LONG_DIGIT(v, i);
        if ((x >> PyLong_SHIFT) != prev) {
            *overflow = sign;
            goto exit;
        }
    }
    /* Haven't lost any bits, but casting to long requires extra
     * care.
     */
    if (x <= (size_t)CPY_TAGGED_MAX) {
        res = (Py_ssize_t)x * sign;
    }
    else if (sign < 0 && x == CPY_TAGGED_ABS_MIN) {
        res = CPY_TAGGED_MIN;
    }
    else {
        *overflow = sign;
        /* res is already set to -1 */
    }
  exit:
    return res;
}

#else

// Slow path of CPyLong_AsSsize_tAndOverflow (non-inlined, Python 3.11 and earlier)
Py_ssize_t
CPyLong_AsSsize_tAndOverflow_(PyObject *vv, int *overflow)
{
    /* This version by Tim Peters */
    PyLongObject *v = (PyLongObject *)vv;
    size_t x, prev;
    Py_ssize_t res;
    Py_ssize_t i;
    int sign;

    *overflow = 0;

    res = -1;
    i = Py_SIZE(v);

    sign = 1;
    x = 0;
    if (i < 0) {
        sign = -1;
        i = -(i);
    }
    while (--i >= 0) {
        prev = x;
        x = (x << PyLong_SHIFT) + CPY_LONG_DIGIT(v, i);
        if ((x >> PyLong_SHIFT) != prev) {
            *overflow = sign;
            goto exit;
        }
    }
    /* Haven't lost any bits, but casting to long requires extra
     * care.
     */
    if (x <= (size_t)CPY_TAGGED_MAX) {
        res = (Py_ssize_t)x * sign;
    }
    else if (sign < 0 && x == CPY_TAGGED_ABS_MIN) {
        res = CPY_TAGGED_MIN;
    }
    else {
        *overflow = sign;
        /* res is already set to -1 */
    }
  exit:
    return res;
}


#endif
