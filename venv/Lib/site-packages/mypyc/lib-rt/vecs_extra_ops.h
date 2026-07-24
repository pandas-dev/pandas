#ifndef VECS_EXTRA_OPS_H
#define VECS_EXTRA_OPS_H

#include "vecs/librt_vecs_api.h"

// Check if obj is an instance of vec (any vec type)
static inline int CPyVec_Check(PyObject *obj) {
    return PyObject_TypeCheck(obj, VecApi->get_vec_type());
}

static inline PyObject *CPyVecU8_ToBytes(VecU8 v) {
    if (v.len == 0) {
        return PyBytes_FromStringAndSize(NULL, 0);
    }
    return PyBytes_FromStringAndSize((const char *)v.items, v.len);
}

#endif
