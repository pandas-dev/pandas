/*
Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#ifndef PANDAS__LIBS_SRC_PARSE_HELPER_H_
#define PANDAS__LIBS_SRC_PARSE_HELPER_H_

#include <float.h>
#include "inline_helper.h"
#include "headers/portable.h"
#include "parser/tokenizer.h"

int to_double(char *item, double *p_value, char sci, char decimal,
              int *maybe_int) {
    char *p_end = NULL;
    int error = 0;

    *p_value = xstrtod(item, &p_end, decimal, sci, '\0', 1, &error, maybe_int);

    return (error == 0) && (!*p_end);
}

#if PY_VERSION_HEX < 0x02060000
#define PyBytes_Check PyString_Check
#define PyBytes_AS_STRING PyString_AS_STRING
#endif  // PY_VERSION_HEX

int floatify(PyObject *str, double *result, int *maybe_int) {
    int status;
    char *data;
    PyObject *tmp = NULL;
    const char sci = 'E';
    const char dec = '.';

    if (PyBytes_Check(str)) {
        data = PyBytes_AS_STRING(str);
    } else if (PyUnicode_Check(str)) {
        tmp = PyUnicode_AsUTF8String(str);
        data = PyBytes_AS_STRING(tmp);
    } else {
        PyErr_SetString(PyExc_TypeError, "Invalid object type");
        return -1;
    }

    status = to_double(data, result, sci, dec, maybe_int);

    if (!status) {
        /* handle inf/-inf */
        if (strlen(data) == 3) {
            if (0 == strcasecmp(data, "inf")) {
                *result = HUGE_VAL;
                *maybe_int = 0;
            } else {
                goto parsingerror;
            }
        } else if (strlen(data) == 4) {
            if (0 == strcasecmp(data, "-inf")) {
                *result = -HUGE_VAL;
                *maybe_int = 0;
            } else if (0 == strcasecmp(data, "+inf")) {
                *result = HUGE_VAL;
                *maybe_int = 0;
            } else {
                goto parsingerror;
            }
        } else {
            goto parsingerror;
        }
    }

    Py_XDECREF(tmp);
    return 0;

parsingerror:
    PyErr_Format(PyExc_ValueError, "Unable to parse string \"%s\"", data);
    Py_XDECREF(tmp);
    return -1;
}

PANDAS_INLINE void lowercase(char *p) {
    for (; *p; ++p) *p = tolower_ascii(*p);
}

PANDAS_INLINE void uppercase(char *p) {
    for (; *p; ++p) *p = toupper_ascii(*p);
}

#endif  // PANDAS__LIBS_SRC_PARSE_HELPER_H_
