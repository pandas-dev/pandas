#ifndef LIBRT_INTERNAL_H
#define LIBRT_INTERNAL_H

#include <Python.h>
#include <stdbool.h>

// ABI version -- only an exact match is compatible. This will only be changed in
// very exceptional cases (likely never) due to strict backward compatibility
// requirements.
#define LIBRT_INTERNAL_ABI_VERSION 2

// API version -- more recent versions must maintain backward compatibility, i.e.
// we can add new features but not remove or change existing features (unless
// ABI version is changed, but see the comment above).
#define LIBRT_INTERNAL_API_VERSION 1

// Number of functions in the capsule API. If you add a new function, also increase
// LIBRT_INTERNAL_API_VERSION.
#define LIBRT_INTERNAL_API_LEN 21

#ifdef LIBRT_INTERNAL_MODULE

static PyObject *ReadBuffer_internal(PyObject *source);
static PyObject *WriteBuffer_internal(void);
static PyObject *WriteBuffer_getvalue_internal(PyObject *self);
static PyObject *ReadBuffer_internal(PyObject *source);
static PyObject *ReadBuffer_internal_empty(void);
static char write_bool_internal(PyObject *data, char value);
static char read_bool_internal(PyObject *data);
static char write_str_internal(PyObject *data, PyObject *value);
static PyObject *read_str_internal(PyObject *data);
static char write_float_internal(PyObject *data, double value);
static double read_float_internal(PyObject *data);
static char write_int_internal(PyObject *data, CPyTagged value);
static CPyTagged read_int_internal(PyObject *data);
static char write_tag_internal(PyObject *data, uint8_t value);
static uint8_t read_tag_internal(PyObject *data);
static int NativeInternal_ABI_Version(void);
static char write_bytes_internal(PyObject *data, PyObject *value);
static PyObject *read_bytes_internal(PyObject *data);
static uint8_t cache_version_internal(void);
static PyTypeObject *ReadBuffer_type_internal(void);
static PyTypeObject *WriteBuffer_type_internal(void);
static int NativeInternal_API_Version(void);
static PyObject *extract_symbol_internal(PyObject *data);

#endif

#endif  // LIBRT_INTERNAL_H
