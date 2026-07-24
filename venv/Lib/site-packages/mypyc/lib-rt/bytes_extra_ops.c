#include "bytes_extra_ops.h"

PyObject *CPyBytes_Translate(PyObject *bytes, PyObject *table) {
    // Fast path: exact bytes object with exact bytes table
    if (PyBytes_CheckExact(bytes) && PyBytes_CheckExact(table)) {
        Py_ssize_t table_len = PyBytes_GET_SIZE(table);
        if (table_len != 256) {
            PyErr_SetString(PyExc_ValueError,
                           "translation table must be 256 characters long");
            return NULL;
        }

        Py_ssize_t len = PyBytes_GET_SIZE(bytes);
        const char *input = PyBytes_AS_STRING(bytes);
        const char *trans_table = PyBytes_AS_STRING(table);

        PyObject *result = PyBytes_FromStringAndSize(NULL, len);
        if (result == NULL) {
            return NULL;
        }

        char *output = PyBytes_AS_STRING(result);
        bool changed = false;

        // Without a loop unrolling hint performance can be worse than CPython
        CPY_UNROLL_LOOP(4)
        for (Py_ssize_t i = len; --i >= 0;) {
            char c = *input++;
            if ((*output++ = trans_table[(unsigned char)c]) != c)
                changed = true;
        }

        // If nothing changed, discard result and return the original object
        if (!changed) {
            Py_DECREF(result);
            Py_INCREF(bytes);
            return bytes;
        }

        return result;
    }

    // Fallback to Python method call for non-exact types or non-standard tables
    return PyObject_CallMethodOneArg(bytes, mypyc_interned_str.translate, table);
}
