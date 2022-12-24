/*
Copyright (c) 2022-, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#ifndef PANDAS__LIBS_PANDAS_TYPE_H_
#define PANDAS__LIBS_PANDAS_TYPE_H_
#include "pyerrors.h"
#ifdef __cplusplus
extern "C" {
#endif

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/ndarrayobject.h>

/*
Cython equivalent of `isinstance(val, np.timedelta64)`

Parameters
----------
val : object

Returns
-------
is_timedelta64 : bool
*/
int is_timedelta64_object(PyObject *obj) {
  return PyObject_TypeCheck(obj, &PyTimedeltaArrType_Type);
}

/*
Cython equivalent of

`isinstance(val, (int, long, np.integer)) and not isinstance(val, bool)`

Parameters
----------
val : object

Returns
-------
is_integer : bool

Notes
-----
This counts np.timedelta64 objects as integers.
*/
int is_integer_object(PyObject *obj) {
  return !PyBool_Check(obj) && PyArray_IsIntegerScalar(obj)
    && !is_timedelta64_object(obj);
}

/*
Cython equivalent of `isinstance(val, (float, np.complex_))`

Parameters
----------
val : object

Returns
-------
is_float : bool
*/
int is_float_object(PyObject *obj) {
  return PyFloat_Check(obj) || PyObject_TypeCheck(obj, &PyFloatingArrType_Type);
}

/*
Cython equivalent of `isinstance(val, (complex, np.complex_))`

Parameters
----------
val : object

Returns
-------
is_complex : bool
*/
int is_complex_object(PyObject *obj) {
  return PyComplex_Check(obj) ||
    PyObject_TypeCheck(obj, &PyComplexFloatingArrType_Type);
}

/*
Cython equivalent of `isinstance(val, (bool, np.bool_))`

Parameters
----------
val : object

Returns
-------
is_bool : bool
*/
int is_bool_object(PyObject *obj) {
  return PyBool_Check(obj) || PyObject_TypeCheck(obj, &PyBoolArrType_Type);
}

int is_real_number_object(PyObject *obj) {
  return is_bool_object(obj) || is_integer_object(obj) || is_float_object(obj);
}

/*
Cython equivalent of `isinstance(val, np.datetime64)`

Parameters
----------
val : object

Returns
-------
is_datetime64 : bool
*/
int is_datetime64_object(PyObject *obj) {
  return PyObject_TypeCheck(obj, &PyDatetimeArrType_Type);
}

/*
Cython equivalent of `isinstance(val, np.ndarray)`

Parameters
----------
val : object

Returns
-------
is_ndarray : bool
*/
int is_array(PyObject *obj) { return PyArray_Check(obj); }


/*
Check if val is a Not-A-Number float or complex, including
float('NaN') and np.nan.

Parameters
----------
val : object

Returns
-------
is_nan : bool
*/
int is_nan(PyObject *obj) {
  if (is_float_object(obj)) {
    double fobj = PyFloat_AsDouble(obj);
    if (fobj == -1.0 && PyErr_Occurred()) {
      // TODO(wayd): handle this error!
    }

    return fobj != fobj;
  }

  if (is_complex_object(obj)) {
    int ret = PyObject_RichCompareBool(obj, obj, Py_NE) == 1;
    if (ret == -1) {
      // TODO(wayd): handle this error!
    }
  }

  return 0;
}

#ifdef __cplusplus
}
#endif
#endif  // PANDAS__LIBS_PANDAS_TYPE_H_
