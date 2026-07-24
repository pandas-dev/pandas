#include "bytearray_extra_ops.h"

PyObject *CPyByteArray_New(void) {
    return PyByteArray_FromStringAndSize(NULL, 0);
}
