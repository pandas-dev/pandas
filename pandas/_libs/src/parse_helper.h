/*
Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#ifndef PANDAS__LIBS_SRC_PARSE_HELPER_H_
#define PANDAS__LIBS_SRC_PARSE_HELPER_H_

#include <errno.h>
#include <float.h>
#include "inline_helper.h"
#include "headers/portable.h"

static double xstrtod(const char *p, char **q, char decimal, char sci,
                      int skip_trailing, int *maybe_int);

int to_double(char *item, double *p_value, char sci, char decimal,
              int *maybe_int) {
    char *p_end = NULL;

    *p_value = xstrtod(item, &p_end, decimal, sci, 1, maybe_int);

    return (errno == 0) && (!*p_end);
}

#if PY_VERSION_HEX < 0x02060000
#define PyBytes_Check PyString_Check
#define PyBytes_AS_STRING PyString_AS_STRING
#endif

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

    /*
    #if PY_VERSION_HEX >= 0x03000000
      return PyFloat_FromString(str);
    #else
      return PyFloat_FromString(str, NULL);
    #endif
    */
}

// ---------------------------------------------------------------------------
// Implementation of xstrtod

//
// strtod.c
//
// Convert string to double
//
// Copyright (C) 2002 Michael Ringgaard. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. Neither the name of the project nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.
//
// -----------------------------------------------------------------------
// Modifications by Warren Weckesser, March 2011:
// * Rename strtod() to xstrtod().
// * Added decimal and sci arguments.
// * Skip trailing spaces.
// * Commented out the other functions.
//

PANDAS_INLINE void lowercase(char *p) {
    for (; *p; ++p) *p = tolower_ascii(*p);
}

PANDAS_INLINE void uppercase(char *p) {
    for (; *p; ++p) *p = toupper_ascii(*p);
}

static double xstrtod(const char *str, char **endptr, char decimal, char sci,
                      int skip_trailing, int *maybe_int) {
    double number;
    int exponent;
    int negative;
    char *p = (char *)str;
    double p10;
    int n;
    int num_digits;
    int num_decimals;

    errno = 0;
    *maybe_int = 1;

    // Skip leading whitespace
    while (isspace(*p)) p++;

    // Handle optional sign
    negative = 0;
    switch (*p) {
        case '-':
            negative = 1;  // Fall through to increment position
        case '+':
            p++;
    }

    number = 0.;
    exponent = 0;
    num_digits = 0;
    num_decimals = 0;

    // Process string of digits
    while (isdigit_ascii(*p)) {
        number = number * 10. + (*p - '0');
        p++;
        num_digits++;
    }

    // Process decimal part
    if (*p == decimal) {
        *maybe_int = 0;
        p++;

        while (isdigit_ascii(*p)) {
            number = number * 10. + (*p - '0');
            p++;
            num_digits++;
            num_decimals++;
        }

        exponent -= num_decimals;
    }

    if (num_digits == 0) {
        errno = ERANGE;
        return 0.0;
    }

    // Correct for sign
    if (negative) number = -number;

    // Process an exponent string
    if (toupper_ascii(*p) == toupper_ascii(sci)) {
        *maybe_int = 0;

        // Handle optional sign
        negative = 0;
        switch (*++p) {
            case '-':
                negative = 1;  // Fall through to increment pos
            case '+':
                p++;
        }

        // Process string of digits
        num_digits = 0;
        n = 0;
        while (isdigit_ascii(*p)) {
            n = n * 10 + (*p - '0');
            num_digits++;
            p++;
        }

        if (negative)
            exponent -= n;
        else
            exponent += n;

        // If no digits, after the 'e'/'E', un-consume it
        if (num_digits == 0) p--;
    }

    if (exponent < DBL_MIN_EXP || exponent > DBL_MAX_EXP) {
        errno = ERANGE;
        return HUGE_VAL;
    }

    // Scale the result
    p10 = 10.;
    n = exponent;
    if (n < 0) n = -n;
    while (n) {
        if (n & 1) {
            if (exponent < 0)
                number /= p10;
            else
                number *= p10;
        }
        n >>= 1;
        p10 *= p10;
    }

    if (number == HUGE_VAL) {
        errno = ERANGE;
    }

    if (skip_trailing) {
        // Skip trailing whitespace
        while (isspace_ascii(*p)) p++;
    }

    if (endptr) *endptr = p;

    return number;
}

#endif  // PANDAS__LIBS_SRC_PARSE_HELPER_H_
