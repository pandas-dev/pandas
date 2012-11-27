### pytables extensions ###

from numpy cimport ndarray, int32_t, float64_t, int64_t
cimport numpy as np

cimport cython

import numpy as np
import operator
import sys

np.import_array()
np.import_ufunc()


from cpython cimport (PyDict_New, PyDict_GetItem, PyDict_SetItem,
                      PyDict_Contains, PyDict_Keys,
                      Py_INCREF, PyTuple_SET_ITEM,
                      PyTuple_SetItem,
                      PyTuple_New,
                      PyObject_SetAttrString)

@cython.boundscheck(False)
@cython.wraparound(False)
def create_hdf_rows_2d(ndarray index, ndarray[np.uint8_t, ndim=1] mask, list values):
    """ return a list of objects ready to be converted to rec-array format """

    cdef:
        unsigned int i, b, n_index, n_blocks, tup_size
        ndarray v
        list l
        object tup, val

    n_index   = index.shape[0]
    n_blocks  = len(values)
    tup_size  = n_blocks+1
    l = []
    for i from 0 <= i < n_index:
        
        if not mask[i]:

            tup = PyTuple_New(tup_size)
            val  = index[i]
            PyTuple_SET_ITEM(tup, 0, val)
            Py_INCREF(val)

            for b from 0 <= b < n_blocks:

                v   = values[b][:, i]
                PyTuple_SET_ITEM(tup, b+1, v)
                Py_INCREF(v)

            l.append(tup)

    return l

@cython.boundscheck(False)
@cython.wraparound(False)
def create_hdf_rows_3d(ndarray index, ndarray columns, ndarray[np.uint8_t, ndim=2] mask, list values):
    """ return a list of objects ready to be converted to rec-array format """

    cdef:
        unsigned int i, j, n_columns, n_index, n_blocks, tup_size
        ndarray v
        list l
        object tup, val

    n_index   = index.shape[0]
    n_columns = columns.shape[0]
    n_blocks  = len(values)
    tup_size  = n_blocks+2
    l = []
    for i from 0 <= i < n_index:

        for c from 0 <= c < n_columns:

            if not mask[i, c]:

                tup = PyTuple_New(tup_size)

                val  = columns[c]
                PyTuple_SET_ITEM(tup, 0, val)
                Py_INCREF(val)

                val  = index[i]
                PyTuple_SET_ITEM(tup, 1, val)
                Py_INCREF(val)

                for b from 0 <= b < n_blocks:

                    v   = values[b][:, i, c]
                    PyTuple_SET_ITEM(tup, b+2, v)
                    Py_INCREF(v)

                l.append(tup)

    return l
