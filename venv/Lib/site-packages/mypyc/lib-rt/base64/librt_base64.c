#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdbool.h>

#define BASE64_EXPORTS

#include "librt_base64.h"
#include "libbase64.h"
#include "pythoncapi_compat.h"

static PyObject *
b64decode_handle_invalid_input(
    PyObject *out_bytes, char *outbuf, size_t max_out, const char *src, size_t srclen, bool freesrc);

#define BASE64_MAXBIN ((PY_SSIZE_T_MAX - 3) / 2)

#define STACK_BUFFER_SIZE 1024

static void
convert_encoded_to_urlsafe(char *buf, size_t len) {
    // The loop is written to enable SIMD optimizations
    for (size_t i = 0; i < len; i++) {
        char ch = buf[i];
        if (ch == '+') {
            ch = '-';
        } else if (ch == '/') {
            ch = '_';
        }
        buf[i] = ch;
    }
}

static void
convert_urlsafe_to_encoded(const char *src, size_t len, char *buf) {
    // The loop is written to enable SIMD optimizations
    for (size_t i = 0; i < len; i++) {
        char ch = src[i];
        if (ch == '-') {
            ch = '+';
        } else if (ch == '_') {
            ch = '/';
        }
        buf[i] = ch;
    }
}

static PyObject *
b64encode_internal(PyObject *obj, bool urlsafe) {
    unsigned char *ascii_data;
    char *bin_data;
    int leftbits = 0;
    unsigned char this_ch;
    unsigned int leftchar = 0;
    Py_ssize_t bin_len, out_len;
    PyBytesWriter *writer;
    int newline = 0; // TODO

    if (!PyBytes_Check(obj)) {
        PyErr_SetString(PyExc_TypeError, "base64() expects a bytes object");
        return NULL;
    }

    bin_data = PyBytes_AS_STRING(obj);
    bin_len = PyBytes_GET_SIZE(obj);
    assert(bin_len >= 0);

    if (bin_len > BASE64_MAXBIN) {
        PyErr_SetString(PyExc_ValueError, "Too much data for base64 line");
        return NULL;
    }

    Py_ssize_t buflen = 4 * bin_len / 3 + 4;
    char *buf;
    char stack_buf[STACK_BUFFER_SIZE];
    if (buflen <= STACK_BUFFER_SIZE) {
        buf = stack_buf;
    } else {
        buf = PyMem_Malloc(buflen);
        if (buf == NULL) {
            return PyErr_NoMemory();
        }
    }
    size_t actual_len;
    base64_encode(bin_data, bin_len, buf, &actual_len, 0);

    if (urlsafe) {
        convert_encoded_to_urlsafe(buf, actual_len);
    }

    PyObject *res = PyBytes_FromStringAndSize(buf, actual_len);
    if (buflen > STACK_BUFFER_SIZE)
        PyMem_Free(buf);
    return res;
}

static PyObject*
b64encode(PyObject *self, PyObject *const *args, size_t nargs) {
    if (nargs != 1) {
        PyErr_SetString(PyExc_TypeError, "b64encode() takes exactly one argument");
        return 0;
    }
    return b64encode_internal(args[0], false);
}

static PyObject*
urlsafe_b64encode(PyObject *self, PyObject *const *args, size_t nargs) {
    if (nargs != 1) {
        PyErr_SetString(PyExc_TypeError, "urlsafe_b64encode() takes exactly one argument");
        return 0;
    }
    return b64encode_internal(args[0], true);
}

static inline int
is_valid_base64_char(char c, bool allow_padding) {
    return ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
            (c >= '0' && c <= '9') || (c == '+') || (c == '/') || (allow_padding && c == '='));
}

static PyObject *
b64decode_internal(PyObject *arg, bool urlsafe) {
    const char *src;
    Py_ssize_t srclen_ssz;

    // Get input pointer and length
    if (PyBytes_Check(arg)) {
        src = PyBytes_AS_STRING(arg);
        srclen_ssz = PyBytes_GET_SIZE(arg);
    } else if (PyUnicode_Check(arg)) {
        if (!PyUnicode_IS_ASCII(arg)) {
            PyErr_SetString(PyExc_ValueError,
                            "string argument should contain only ASCII characters");
            return NULL;
        }
        src = (const char *)PyUnicode_1BYTE_DATA(arg);
        srclen_ssz = PyUnicode_GET_LENGTH(arg);
    } else {
        PyErr_SetString(PyExc_TypeError,
                        "argument should be a bytes-like object or ASCII string");
        return NULL;
    }

    // Fast-path: empty input
    if (srclen_ssz == 0) {
        return PyBytes_FromStringAndSize(NULL, 0);
    }

    if (urlsafe) {
        char *new_src = PyMem_Malloc(srclen_ssz + 1);
        if (new_src == NULL) {
            return PyErr_NoMemory();
        }
        convert_urlsafe_to_encoded(src, srclen_ssz, new_src);
        src = new_src;
    }

    // Quickly ignore invalid characters at the end. Other invalid characters
    // are also accepted, but they need a slow path.
    while (srclen_ssz > 0 && !is_valid_base64_char(src[srclen_ssz - 1], true)) {
        srclen_ssz--;
    }

    // Compute an output capacity that's at least 3/4 of input, without overflow:
    // ceil(3/4 * N) == N - floor(N/4)
    size_t srclen = (size_t)srclen_ssz;
    size_t max_out = srclen - (srclen / 4);
    if (max_out == 0) {
        max_out = 1; // defensive (srclen > 0 implies >= 1 anyway)
    }
    if (max_out > (size_t)PY_SSIZE_T_MAX) {
        PyErr_SetString(PyExc_OverflowError, "input too large");
        return NULL;
    }

    // Allocate output bytes (uninitialized) of the max capacity
    PyObject *out_bytes = PyBytes_FromStringAndSize(NULL, (Py_ssize_t)max_out);
    if (out_bytes == NULL) {
        if (urlsafe) {
            PyMem_Free((void *)src);
        }
        return NULL; // Propagate memory error
    }

    char *outbuf = PyBytes_AS_STRING(out_bytes);
    size_t outlen = max_out;

    int ret = base64_decode(src, srclen, outbuf, &outlen, 0);

    if (ret != 1) {
        if (ret == 0) {
            // Slow path: handle non-base64 input
            return b64decode_handle_invalid_input(out_bytes, outbuf, max_out, src, srclen, urlsafe);
        }
        Py_DECREF(out_bytes);
        if (urlsafe) {
            PyMem_Free((void *)src);
        }
        if (ret == -1) {
            PyErr_SetString(PyExc_NotImplementedError, "base64 codec not available in this build");
        } else {
            PyErr_SetString(PyExc_RuntimeError, "base64_decode failed");
        }
        return NULL;
    }

    if (urlsafe) {
        PyMem_Free((void *)src);
    }

    // Sanity-check contract (decoder must not overflow our buffer)
    if (outlen > max_out) {
        Py_DECREF(out_bytes);
        PyErr_SetString(PyExc_RuntimeError, "decoder wrote past output buffer");
        return NULL;
    }

    // Shrink in place to the actual decoded length
    if (_PyBytes_Resize(&out_bytes, (Py_ssize_t)outlen) < 0) {
        // _PyBytes_Resize sets an exception and may free the old object
        return NULL;
    }
    return out_bytes;
}

// Process non-base64 input by ignoring non-base64 characters, for compatibility
// with stdlib b64decode.
static PyObject *
b64decode_handle_invalid_input(
    PyObject *out_bytes, char *outbuf, size_t max_out, const char *src, size_t srclen, bool freesrc)
{
    // Copy input to a temporary buffer, with non-base64 characters and extra suffix
    // characters removed
    size_t newbuf_len = 0;
    char *newbuf = PyMem_Malloc(srclen);
    if (newbuf == NULL) {
        Py_DECREF(out_bytes);
        if (freesrc) {
            PyMem_Free((void *)src);
        }
        return PyErr_NoMemory();
    }

    int pad_chars = 0;
    // Copy base64 characters to the new buffer. Ignore padding to conform to RFC 4648 section 3.3.
    for (size_t i = 0; i < srclen; i++) {
        char c = src[i];
        if (is_valid_base64_char(c, false)) {
            newbuf[newbuf_len++] = c;
            pad_chars = 0;
        } else if (c == '=') {
            pad_chars++;
        }
    }

    int quad_pos = newbuf_len % 4;
    // Stdlib always performs a non-strict padding check
    if (quad_pos != 0 && quad_pos + pad_chars < 4) {
        if (freesrc) {
            PyMem_Free((void *)src);
        }
        Py_DECREF(out_bytes);
        PyMem_Free(newbuf);
        PyErr_SetString(PyExc_ValueError, "Incorrect padding");
        return NULL;
    }

    if (quad_pos != 0) {
        // Add padding at the end to make the input length a multiple of 4. We know that this padding
        // is present in src because otherwise we would report the "Incorrect padding" error above.
        while (quad_pos < 4) {
            newbuf[newbuf_len++] = '=';
            quad_pos++;
        }
    }

    size_t outlen = max_out;
    int ret = base64_decode(newbuf, newbuf_len, outbuf, &outlen, 0);
    PyMem_Free(newbuf);
    if (freesrc) {
        PyMem_Free((void *)src);
    }

    if (ret != 1) {
        Py_DECREF(out_bytes);
        if (ret == 0) {
            PyErr_SetString(PyExc_ValueError, "Only base64 data is allowed");
        }
        if (ret == -1) {
            PyErr_SetString(PyExc_NotImplementedError, "base64 codec not available in this build");
        } else {
            PyErr_SetString(PyExc_RuntimeError, "base64_decode failed");
        }
        return NULL;
    }

    // Shrink in place to the actual decoded length
    if (_PyBytes_Resize(&out_bytes, (Py_ssize_t)outlen) < 0) {
        // _PyBytes_Resize sets an exception and may free the old object
        return NULL;
    }
    return out_bytes;
}

static PyObject*
b64decode(PyObject *self, PyObject *const *args, size_t nargs) {
    if (nargs != 1) {
        PyErr_SetString(PyExc_TypeError, "b64decode() takes exactly one argument");
        return 0;
    }
    return b64decode_internal(args[0], false);
}

static PyObject*
urlsafe_b64decode(PyObject *self, PyObject *const *args, size_t nargs) {
    if (nargs != 1) {
        PyErr_SetString(PyExc_TypeError, "urlsafe_b64decode() takes exactly one argument");
        return 0;
    }
    return b64decode_internal(args[0], true);
}

static PyMethodDef librt_base64_module_methods[] = {
    {"b64encode", (PyCFunction)b64encode, METH_FASTCALL, PyDoc_STR("Encode bytes object using Base64.")},
    {"b64decode", (PyCFunction)b64decode, METH_FASTCALL, PyDoc_STR("Decode a Base64 encoded bytes object or ASCII string.")},
    {"urlsafe_b64encode", (PyCFunction)urlsafe_b64encode, METH_FASTCALL, PyDoc_STR("Encode bytes object using URL and file system safe Base64 alphabet.")},
    {"urlsafe_b64decode", (PyCFunction)urlsafe_b64decode, METH_FASTCALL, PyDoc_STR("Decode bytes or ASCII string using URL and file system safe Base64 alphabet.")},
    {NULL, NULL, 0, NULL}
};

static int
base64_abi_version(void) {
    return LIBRT_BASE64_ABI_VERSION;
}

static int
base64_api_version(void) {
    return LIBRT_BASE64_API_VERSION;
}

static int
librt_base64_module_exec(PyObject *m)
{
    // Export mypy internal C API, be careful with the order!
    static void *base64_api[LIBRT_BASE64_API_LEN] = {
        (void *)base64_abi_version,
        (void *)base64_api_version,
        (void *)b64encode_internal,
        (void *)b64decode_internal,
    };
    PyObject *c_api_object = PyCapsule_New((void *)base64_api, "librt.base64._C_API", NULL);
    if (PyModule_Add(m, "_C_API", c_api_object) < 0) {
        return -1;
    }
    return 0;
}

static PyModuleDef_Slot librt_base64_module_slots[] = {
    {Py_mod_exec, librt_base64_module_exec},
#ifdef Py_MOD_GIL_NOT_USED
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL}
};

static PyModuleDef librt_base64_module = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "base64",
    .m_doc = "Fast base64 encoding and decoding optimized for mypyc",
    .m_size = 0,
    .m_methods = librt_base64_module_methods,
    .m_slots = librt_base64_module_slots,
};

PyMODINIT_FUNC
PyInit_base64(void)
{
    return PyModuleDef_Init(&librt_base64_module);
}
