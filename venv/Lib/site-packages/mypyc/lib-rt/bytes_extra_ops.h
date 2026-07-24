#ifndef MYPYC_BYTES_EXTRA_OPS_H
#define MYPYC_BYTES_EXTRA_OPS_H

#include <Python.h>
#include <stdint.h>
#include "CPy.h"

// Optimized bytes translate operation
PyObject *CPyBytes_Translate(PyObject *bytes, PyObject *table);

// Optimized bytes.__getitem__ operations

// If index is negative, convert to non-negative index (no range checking)
static inline int64_t CPyBytes_AdjustIndex(PyObject *obj, int64_t index) {
    if (index < 0) {
        return index + Py_SIZE(obj);
    }
    return index;
}

// Check if index is in valid range [0, len)
static inline bool CPyBytes_RangeCheck(PyObject *obj, int64_t index) {
    return index >= 0 && index < Py_SIZE(obj);
}

// Get byte at index (no bounds checking) - returns as CPyTagged
static inline CPyTagged CPyBytes_GetItemUnsafe(PyObject *obj, int64_t index) {
    return ((CPyTagged)(uint8_t)(PyBytes_AS_STRING(obj))[index]) << 1;
}

#endif
