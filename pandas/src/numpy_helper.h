#include "Python.h"
#include "numpy/arrayobject.h"
#include "numpy/arrayscalars.h"

#ifndef PANDAS_INLINE
  #if defined(__GNUC__)
    #define PANDAS_INLINE __inline__
  #elif defined(_MSC_VER)
    #define PANDAS_INLINE __inline
  #elif defined (__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
    #define PANDAS_INLINE inline
  #else
    #define PANDAS_INLINE
  #endif
#endif

#define PANDAS_FLOAT 0
#define PANDAS_INT 1
#define PANDAS_BOOL 2
#define PANDAS_STRING 3
#define PANDAS_OBJECT 4
#define PANDAS_DATETIME 5

PANDAS_INLINE int
infer_type(PyObject* obj) {
  if (PyBool_Check(obj)) {
    return PANDAS_BOOL;
  }
  else if (PyArray_IsIntegerScalar(obj)) {
    return PANDAS_INT;
  }
  else if (PyArray_IsScalar(obj, Datetime)) {
    return PANDAS_DATETIME;
  }
  else if (PyFloat_Check(obj) || PyArray_IsScalar(obj, Floating)) {
    return PANDAS_FLOAT;
  }
  else if (PyString_Check(obj) || PyUnicode_Check(obj)) {
    return PANDAS_STRING;
  }
  else {
    return PANDAS_OBJECT;
  }
}

PANDAS_INLINE npy_int64
get_nat(void) {
  return NPY_MIN_INT64;
}

PANDAS_INLINE npy_datetime
get_datetime64_value(PyObject* obj) {
  return ((PyDatetimeScalarObject*) obj)->obval;

}


PANDAS_INLINE int
is_integer_object(PyObject* obj) {
  return (!PyBool_Check(obj)) && PyArray_IsIntegerScalar(obj);
//  return PyArray_IsIntegerScalar(obj);
}

PANDAS_INLINE int
is_float_object(PyObject* obj) {
  return (PyFloat_Check(obj) || PyArray_IsScalar(obj, Floating));
}
PANDAS_INLINE int
is_complex_object(PyObject* obj) {
  return (PyComplex_Check(obj) || PyArray_IsScalar(obj, ComplexFloating));
}

PANDAS_INLINE int
is_bool_object(PyObject* obj) {
  return (PyBool_Check(obj) || PyArray_IsScalar(obj, Bool));
}

PANDAS_INLINE int
is_string_object(PyObject* obj) {
  return (PyString_Check(obj) || PyUnicode_Check(obj));
}

PANDAS_INLINE int
is_datetime64_object(PyObject *obj) {
  return PyArray_IsScalar(obj, Datetime);
}

PANDAS_INLINE int
assign_value_1d(PyArrayObject* ap, Py_ssize_t _i, PyObject* v) {
  npy_intp i = (npy_intp) _i;
  char *item = (char *) PyArray_DATA(ap) + i * PyArray_STRIDE(ap, 0);
  return PyArray_DESCR(ap)->f->setitem(v, item, ap);
}

PANDAS_INLINE PyObject*
get_value_1d(PyArrayObject* ap, Py_ssize_t i) {
  char *item = (char *) PyArray_DATA(ap) + i * PyArray_STRIDE(ap, 0);
  return PyArray_Scalar(item, PyArray_DESCR(ap), (PyObject*) ap);
}


PANDAS_INLINE char*
get_c_string(PyObject* obj) {
#if PY_VERSION_HEX >= 0x03000000
  PyObject* enc_str = PyUnicode_AsEncodedString(obj, "utf-8", "error");

  char *ret;
  ret = PyBytes_AS_STRING(enc_str);

  // TODO: memory leak here

  // Py_XDECREF(enc_str);
  return ret;
#else
  return PyString_AsString(obj);
#endif
}

PANDAS_INLINE PyObject*
char_to_string(char* data) {
#if PY_VERSION_HEX >= 0x03000000
    return PyUnicode_FromString(data);
#else
    return PyString_FromString(data);
#endif
}

// PANDAS_INLINE int
// is_string(PyObject* obj) {
// #if PY_VERSION_HEX >= 0x03000000
//   return PyUnicode_Check(obj);
// #else
//   return PyString_Check(obj);
// #endif

#include <errno.h>
#include <float.h>

PANDAS_INLINE double
xstrtod(const char *p, char **q, char decimal, char sci, int skip_trailing);

int to_double(char *item, double *p_value, char sci, char decimal)
{
    char *p_end;

    *p_value = xstrtod(item, &p_end, decimal, sci, 1);

    return (errno == 0) && (!*p_end);
}

#if PY_VERSION_HEX < 0x02060000
  #define PyBytes_Check                PyString_Check
  #define PyBytes_AS_STRING            PyString_AS_STRING
#endif

PANDAS_INLINE int floatify(PyObject* str, double *result) {
    int status;
    char *data;
    PyObject* tmp = NULL;
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

    status = to_double(data, result, sci, dec);

    if (!status) {
        /* handle inf/-inf */
        if (0 == strcmp(data, "-inf")) {
            *result = -HUGE_VAL;
        } else if (0 == strcmp(data, "inf")) {
            *result = HUGE_VAL;
        } else {
            PyErr_SetString(PyExc_ValueError, "Unable to parse string");
            Py_XDECREF(tmp);
            return -1;
        }
    }

    Py_XDECREF(tmp);
    return 0;

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
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
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
    for ( ; *p; ++p) *p = tolower(*p);
}

PANDAS_INLINE void uppercase(char *p) {
    for ( ; *p; ++p) *p = toupper(*p);
}


PANDAS_INLINE double xstrtod(const char *str, char **endptr, char decimal,
                             char sci, int skip_trailing)
{
  double number;
  int exponent;
  int negative;
  char *p = (char *) str;
  double p10;
  int n;
  int num_digits;
  int num_decimals;

  errno = 0;

  // Skip leading whitespace
  while (isspace(*p)) p++;

  // Handle optional sign
  negative = 0;
  switch (*p)
  {
    case '-': negative = 1; // Fall through to increment position
    case '+': p++;
  }

  number = 0.;
  exponent = 0;
  num_digits = 0;
  num_decimals = 0;

  // Process string of digits
  while (isdigit(*p))
  {
    number = number * 10. + (*p - '0');
    p++;
    num_digits++;
  }

  // Process decimal part
  if (*p == decimal)
  {
    p++;

    while (isdigit(*p))
    {
      number = number * 10. + (*p - '0');
      p++;
      num_digits++;
      num_decimals++;
    }

    exponent -= num_decimals;
  }

  if (num_digits == 0)
  {
    errno = ERANGE;
    return 0.0;
  }

  // Correct for sign
  if (negative) number = -number;

  // Process an exponent string
  if (toupper(*p) == toupper(sci))
  {
    // Handle optional sign
    negative = 0;
    switch (*++p)
    {
      case '-': negative = 1;   // Fall through to increment pos
      case '+': p++;
    }

    // Process string of digits
    n = 0;
    while (isdigit(*p))
    {
      n = n * 10 + (*p - '0');
      p++;
    }

    if (negative)
      exponent -= n;
    else
      exponent += n;
  }


  if (exponent < DBL_MIN_EXP  || exponent > DBL_MAX_EXP)
  {

    errno = ERANGE;
    return HUGE_VAL;
  }

  // Scale the result
  p10 = 10.;
  n = exponent;
  if (n < 0) n = -n;
  while (n)
  {
    if (n & 1)
    {
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
      while (isspace(*p)) p++;
  }

  if (endptr) *endptr = p;


  return number;
}

void set_array_owndata(PyArrayObject *ao) {
    ao->flags |= NPY_OWNDATA;
}


// PANDAS_INLINE PyObject*
// get_base_ndarray(PyObject* ap) {
//   // if (!ap || (NULL == ap)) {
//   //   Py_RETURN_NONE;
//   // }

//   while (!PyArray_CheckExact(ap)) {
//     ap = PyArray_BASE((PyArrayObject*) ap);
//     if (ap == Py_None) Py_RETURN_NONE;
//   }
//   // PyArray_BASE is a borrowed reference
//   if(ap) {
//     Py_INCREF(ap);
//   }
//   return ap;
// }

