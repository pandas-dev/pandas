#ifndef MYPYC_STR_EXTRA_OPS_H
#define MYPYC_STR_EXTRA_OPS_H

#include <Python.h>
#include <stdint.h>
#include "CPy.h"

// Optimized str indexing for ord(s[i])

// If index is negative, convert to non-negative index (no range checking)
static inline int64_t CPyStr_AdjustIndex(PyObject *obj, int64_t index) {
    if (index < 0) {
        return index + PyUnicode_GET_LENGTH(obj);
    }
    return index;
}

// Check if index is in valid range [0, len)
static inline bool CPyStr_RangeCheck(PyObject *obj, int64_t index) {
    return index >= 0 && index < PyUnicode_GET_LENGTH(obj);
}

// Get character at index as int (ord value) - no bounds checking, returns as CPyTagged
static inline CPyTagged CPyStr_GetItemUnsafeAsInt(PyObject *obj, int64_t index) {
    int kind = PyUnicode_KIND(obj);
    return PyUnicode_READ(kind, PyUnicode_DATA(obj), index) << 1;
}

#endif
