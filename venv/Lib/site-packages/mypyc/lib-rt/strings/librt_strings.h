#ifndef LIBRT_STRINGS_H
#define LIBRT_STRINGS_H

#include <stdbool.h>
#include <Python.h>
#include "librt_strings_common.h"

// ABI version -- only an exact match is compatible. This will only be changed in
// very exceptional cases (likely never) due to strict backward compatibility
// requirements.
#define LIBRT_STRINGS_ABI_VERSION 1

// API version -- more recent versions must maintain backward compatibility, i.e.
// we can add new features but not remove or change existing features (unless
// ABI version is changed, but see the comment above).
#define LIBRT_STRINGS_API_VERSION 4

// Number of functions in the capsule API. If you add a new function, also increase
// LIBRT_STRINGS_API_VERSION.
#define LIBRT_STRINGS_API_LEN 14

typedef struct {
    PyObject_HEAD
    char *buf;  // Beginning of the buffer
    char kind;  // Bytes per code point (1, 2 or 4)
    Py_ssize_t len;  // Current length (number of code points written)
    Py_ssize_t capacity;  // Total capacity of the buffer (number of code points)
    char data[WRITER_EMBEDDED_BUF_LEN];  // Default buffer
} StringWriterObject;

#endif  // LIBRT_STRINGS_H
