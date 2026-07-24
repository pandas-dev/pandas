#ifndef LIBRT_STRINGS_API_H
#define LIBRT_STRINGS_API_H

int
import_librt_strings(void);

#include <stdbool.h>
#include <Python.h>
#include "librt_strings.h"

extern void *LibRTStrings_API[LIBRT_STRINGS_API_LEN];

#define LibRTStrings_ABIVersion (*(int (*)(void)) LibRTStrings_API[0])
#define LibRTStrings_APIVersion (*(int (*)(void)) LibRTStrings_API[1])
#define LibRTStrings_BytesWriter_internal (*(PyObject* (*)(void)) LibRTStrings_API[2])
#define LibRTStrings_BytesWriter_getvalue_internal (*(PyObject* (*)(PyObject *source)) LibRTStrings_API[3])
#define LibRTStrings_BytesWriter_append_internal (*(char (*)(PyObject *source, uint8_t value)) LibRTStrings_API[4])
#define LibRTStrings_ByteWriter_grow_buffer_internal (*(bool (*)(BytesWriterObject *obj, Py_ssize_t size)) LibRTStrings_API[5])
#define LibRTStrings_BytesWriter_type_internal (*(PyTypeObject* (*)(void)) LibRTStrings_API[6])
#define LibRTStrings_BytesWriter_truncate_internal (*(char (*)(PyObject *self, int64_t size)) LibRTStrings_API[7])
#define LibRTStrings_StringWriter_internal (*(PyObject* (*)(void)) LibRTStrings_API[8])
#define LibRTStrings_StringWriter_getvalue_internal (*(PyObject* (*)(PyObject *source)) LibRTStrings_API[9])
#define LibRTStrings_string_append_slow_path (*(char (*)(StringWriterObject *obj, int32_t value)) LibRTStrings_API[10])
#define LibRTStrings_StringWriter_type_internal (*(PyTypeObject* (*)(void)) LibRTStrings_API[11])
#define LibRTStrings_StringWriter_write_internal (*(char (*)(PyObject *source, PyObject *value)) LibRTStrings_API[12])
#define LibRTStrings_grow_string_buffer (*(bool (*)(StringWriterObject *obj, Py_ssize_t n)) LibRTStrings_API[13])


static inline bool CPyBytesWriter_Check(PyObject *obj) {
    return Py_TYPE(obj) == LibRTStrings_BytesWriter_type_internal();
}

static inline bool CPyStringWriter_Check(PyObject *obj) {
    return Py_TYPE(obj) == LibRTStrings_StringWriter_type_internal();
}

#endif  // LIBRT_STRINGS_API_H
