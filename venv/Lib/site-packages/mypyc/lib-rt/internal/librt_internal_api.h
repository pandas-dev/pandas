#ifndef LIBRT_INTERNAL_API_H
#define LIBRT_INTERNAL_API_H

#include "librt_internal.h"

extern void *NativeInternal_API[LIBRT_INTERNAL_API_LEN];

#define ReadBuffer_internal (*(PyObject* (*)(PyObject *source)) NativeInternal_API[0])
#define WriteBuffer_internal (*(PyObject* (*)(void)) NativeInternal_API[1])
#define WriteBuffer_getvalue_internal (*(PyObject* (*)(PyObject *source)) NativeInternal_API[2])
#define write_bool_internal (*(char (*)(PyObject *source, char value)) NativeInternal_API[3])
#define read_bool_internal (*(char (*)(PyObject *source)) NativeInternal_API[4])
#define write_str_internal (*(char (*)(PyObject *source, PyObject *value)) NativeInternal_API[5])
#define read_str_internal (*(PyObject* (*)(PyObject *source)) NativeInternal_API[6])
#define write_float_internal (*(char (*)(PyObject *source, double value)) NativeInternal_API[7])
#define read_float_internal (*(double (*)(PyObject *source)) NativeInternal_API[8])
#define write_int_internal (*(char (*)(PyObject *source, CPyTagged value)) NativeInternal_API[9])
#define read_int_internal (*(CPyTagged (*)(PyObject *source)) NativeInternal_API[10])
#define write_tag_internal (*(char (*)(PyObject *source, uint8_t value)) NativeInternal_API[11])
#define read_tag_internal (*(uint8_t (*)(PyObject *source)) NativeInternal_API[12])
#define NativeInternal_ABI_Version (*(int (*)(void)) NativeInternal_API[13])
#define write_bytes_internal (*(char (*)(PyObject *source, PyObject *value)) NativeInternal_API[14])
#define read_bytes_internal (*(PyObject* (*)(PyObject *source)) NativeInternal_API[15])
#define cache_version_internal (*(uint8_t (*)(void)) NativeInternal_API[16])
#define ReadBuffer_type_internal (*(PyTypeObject* (*)(void)) NativeInternal_API[17])
#define WriteBuffer_type_internal (*(PyTypeObject* (*)(void)) NativeInternal_API[18])
#define NativeInternal_API_Version (*(int (*)(void)) NativeInternal_API[19])
#define extract_symbol_internal (*(PyObject* (*)(PyObject *source)) NativeInternal_API[20])

int
import_librt_internal(void);

static inline bool CPyReadBuffer_Check(PyObject *obj) {
    return Py_TYPE(obj) == ReadBuffer_type_internal();
}

static inline bool CPyWriteBuffer_Check(PyObject *obj) {
    return Py_TYPE(obj) == WriteBuffer_type_internal();
}

#endif  // LIBRT_INTERNAL_API_H
