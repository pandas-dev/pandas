# -*- coding: utf-8 -*-
import operator

from cpython cimport (PyFloat_Check, PyBool_Check,
                      PyObject_RichCompareBool,
                      Py_EQ, Py_NE, Py_LT, Py_LE, Py_GT, Py_GE)

cimport cython
from cython cimport Py_ssize_t

import numpy as np
from numpy cimport ndarray, uint8_t, import_array
import_array()


from util cimport UINT8_MAX, is_nan

from missing cimport checknull


@cython.wraparound(False)
@cython.boundscheck(False)
def scalar_compare(object[:] values, object val, object op):
    """
    Compare each element of `values` array with the scalar `val`, with
    the comparison operation described by `op`.

    Parameters
    ----------
    values : ndarray[object]
    val : object
    op : {operator.eq, operator.ne,
          operator.le, operator.lt,
          operator.ge, operator.gt}

    Returns
    -------
    result : ndarray[bool]
    """
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[uint8_t, cast=True] result
        bint isnull_val
        int flag
        object x

    if op is operator.lt:
        flag = Py_LT
    elif op is operator.le:
        flag = Py_LE
    elif op is operator.gt:
        flag = Py_GT
    elif op is operator.ge:
        flag = Py_GE
    elif op is operator.eq:
        flag = Py_EQ
    elif op is operator.ne:
        flag = Py_NE
    else:
        raise ValueError('Unrecognized operator')

    result = np.empty(n, dtype=bool).view(np.uint8)
    isnull_val = checknull(val)

    if flag == Py_NE:
        for i in range(n):
            x = values[i]
            if checknull(x):
                result[i] = True
            elif isnull_val:
                result[i] = True
            else:
                try:
                    result[i] = PyObject_RichCompareBool(x, val, flag)
                except TypeError:
                    result[i] = True
    elif flag == Py_EQ:
        for i in range(n):
            x = values[i]
            if checknull(x):
                result[i] = False
            elif isnull_val:
                result[i] = False
            else:
                try:
                    result[i] = PyObject_RichCompareBool(x, val, flag)
                except TypeError:
                    result[i] = False

    else:
        for i in range(n):
            x = values[i]
            if checknull(x):
                result[i] = False
            elif isnull_val:
                result[i] = False
            else:
                result[i] = PyObject_RichCompareBool(x, val, flag)

    return result.view(bool)


@cython.wraparound(False)
@cython.boundscheck(False)
def vec_compare(object[:] left, object[:] right, object op):
    """
    Compare the elements of `left` with the elements of `right` pointwise,
    with the comparison operation described by `op`.

    Parameters
    ----------
    left : ndarray[object]
    right : ndarray[object]
    op : {operator.eq, operator.ne,
          operator.le, operator.lt,
          operator.ge, operator.gt}

    Returns
    -------
    result : ndarray[bool]
    """
    cdef:
        Py_ssize_t i, n = len(left)
        ndarray[uint8_t, cast=True] result
        int flag

    if n != len(right):
        raise ValueError('Arrays were different lengths: {n} vs {nright}'
                         .format(n=n, nright=len(right)))

    if op is operator.lt:
        flag = Py_LT
    elif op is operator.le:
        flag = Py_LE
    elif op is operator.gt:
        flag = Py_GT
    elif op is operator.ge:
        flag = Py_GE
    elif op is operator.eq:
        flag = Py_EQ
    elif op is operator.ne:
        flag = Py_NE
    else:
        raise ValueError('Unrecognized operator')

    result = np.empty(n, dtype=bool).view(np.uint8)

    if flag == Py_NE:
        for i in range(n):
            x = left[i]
            y = right[i]

            if checknull(x) or checknull(y):
                result[i] = True
            else:
                result[i] = PyObject_RichCompareBool(x, y, flag)
    else:
        for i in range(n):
            x = left[i]
            y = right[i]

            if checknull(x) or checknull(y):
                result[i] = False
            else:
                result[i] = PyObject_RichCompareBool(x, y, flag)

    return result.view(bool)


@cython.wraparound(False)
@cython.boundscheck(False)
def scalar_binop(object[:] values, object val, object op):
    """
    Apply the given binary operator `op` between each element of the array
    `values` and the scalar `val`.

    Parameters
    ----------
    values : ndarray[object]
    val : object
    op : binary operator

    Returns
    -------
    result : ndarray[object]
    """
    cdef:
        Py_ssize_t i, n = len(values)
        object[:] result
        object x

    result = np.empty(n, dtype=object)
    if val is None or is_nan(val):
        result[:] = val
        return result.base  # `.base` to access underlying np.ndarray

    for i in range(n):
        x = values[i]
        if x is None or is_nan(x):
            result[i] = x
        else:
            result[i] = op(x, val)

    return maybe_convert_bool(result.base)


@cython.wraparound(False)
@cython.boundscheck(False)
def vec_binop(object[:] left, object[:] right, object op):
    """
    Apply the given binary operator `op` pointwise to the elements of
    arrays `left` and `right`.

    Parameters
    ----------
    left : ndarray[object]
    right : ndarray[object]
    op : binary operator

    Returns
    -------
    result : ndarray[object]
    """
    cdef:
        Py_ssize_t i, n = len(left)
        object[:] result

    if n != len(right):
        raise ValueError('Arrays were different lengths: {n} vs {nright}'
                         .format(n=n, nright=len(right)))

    result = np.empty(n, dtype=object)

    for i in range(n):
        x = left[i]
        y = right[i]
        try:
            result[i] = op(x, y)
        except TypeError:
            if x is None or is_nan(x):
                result[i] = x
            elif y is None or is_nan(y):
                result[i] = y
            else:
                raise

    return maybe_convert_bool(result.base)  # `.base` to access np.ndarray


def maybe_convert_bool(ndarray[object] arr,
                       true_values=None, false_values=None):
    cdef:
        Py_ssize_t i, n
        ndarray[uint8_t] result
        object val
        set true_vals, false_vals
        int na_count = 0

    n = len(arr)
    result = np.empty(n, dtype=np.uint8)

    # the defaults
    true_vals = {'True', 'TRUE', 'true'}
    false_vals = {'False', 'FALSE', 'false'}

    if true_values is not None:
        true_vals = true_vals | set(true_values)

    if false_values is not None:
        false_vals = false_vals | set(false_values)

    for i in range(n):
        val = arr[i]

        if PyBool_Check(val):
            if val is True:
                result[i] = 1
            else:
                result[i] = 0
        elif val in true_vals:
            result[i] = 1
        elif val in false_vals:
            result[i] = 0
        elif PyFloat_Check(val):
            result[i] = UINT8_MAX
            na_count += 1
        else:
            return arr

    if na_count > 0:
        mask = result == UINT8_MAX
        arr = result.view(np.bool_).astype(object)
        np.putmask(arr, mask, np.nan)
        return arr
    else:
        return result.view(np.bool_)
