/*
Copyright (c) 2011-2013, ESN Social Software AB and Jonas Tarnstrom
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the ESN Social Software AB nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ESN SOCIAL SOFTWARE AB OR JONAS TARNSTROM BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Portions of code from MODP_ASCII - Ascii transformations (upper/lower, etc)
https://github.com/client9/stringencoders
Copyright (c) 2007  Nick Galbreath -- nickg [at] modp [dot] com. All rights
reserved.

Numeric decoder derived from TCL library
https://www.opensource.apple.com/source/tcl/tcl-14/tcl/license.terms
 * Copyright (c) 1988-1993 The Regents of the University of California.
 * Copyright (c) 1994 Sun Microsystems, Inc.
*/

// Licence at LICENSES/ULTRAJSON_LICENSE

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "pandas/vendored/ujson/lib/ultrajson.h"

static int Object_objectAddKey(void *Py_UNUSED(prv), JSOBJ obj, JSOBJ name,
                               JSOBJ value) {
  int ret = PyDict_SetItem(obj, name, value);
  Py_DECREF((PyObject *)name);
  Py_DECREF((PyObject *)value);
  return ret == 0 ? 1 : 0;
}

static int Object_arrayAddItem(void *Py_UNUSED(prv), JSOBJ obj, JSOBJ value) {
  int ret = PyList_Append(obj, value);
  Py_DECREF((PyObject *)value);
  return ret == 0 ? 1 : 0;
}

static JSOBJ Object_newString(void *Py_UNUSED(prv), wchar_t *start,
                              wchar_t *end) {
  return PyUnicode_FromWideChar(start, (end - start));
}

static JSOBJ Object_newTrue(void *Py_UNUSED(prv)) { Py_RETURN_TRUE; }

static JSOBJ Object_newFalse(void *Py_UNUSED(prv)) { Py_RETURN_FALSE; }

static JSOBJ Object_newNull(void *Py_UNUSED(prv)) { Py_RETURN_NONE; }

static JSOBJ Object_newPosInf(void *Py_UNUSED(prv)) {
  return PyFloat_FromDouble(Py_HUGE_VAL);
}

static JSOBJ Object_newNegInf(void *Py_UNUSED(prv)) {
  return PyFloat_FromDouble(-Py_HUGE_VAL);
}

static JSOBJ Object_newObject(void *Py_UNUSED(prv), void *Py_UNUSED(decoder)) {
  return PyDict_New();
}

static JSOBJ Object_endObject(void *Py_UNUSED(prv), JSOBJ obj) { return obj; }

static JSOBJ Object_newArray(void *Py_UNUSED(prv), void *Py_UNUSED(decoder)) {
  return PyList_New(0);
}

static JSOBJ Object_endArray(void *Py_UNUSED(prv), JSOBJ obj) { return obj; }

static JSOBJ Object_newInteger(void *Py_UNUSED(prv), JSINT32 value) {
  return PyLong_FromLong(value);
}

static JSOBJ Object_newLong(void *Py_UNUSED(prv), JSINT64 value) {
  return PyLong_FromLongLong(value);
}

static JSOBJ Object_newUnsignedLong(void *Py_UNUSED(prv), JSUINT64 value) {
  return PyLong_FromUnsignedLongLong(value);
}

static JSOBJ Object_newDouble(void *Py_UNUSED(prv), double value) {
  return PyFloat_FromDouble(value);
}

static void Object_releaseObject(void *Py_UNUSED(prv), JSOBJ obj,
                                 void *Py_UNUSED(decoder)) {
  Py_XDECREF(((PyObject *)obj));
}

PyObject *JSONToObj(PyObject *Py_UNUSED(self), PyObject *args,
                    PyObject *kwargs) {
  JSONObjectDecoder dec = {.newString = Object_newString,
                           .objectAddKey = Object_objectAddKey,
                           .arrayAddItem = Object_arrayAddItem,
                           .newTrue = Object_newTrue,
                           .newFalse = Object_newFalse,
                           .newNull = Object_newNull,
                           .newPosInf = Object_newPosInf,
                           .newNegInf = Object_newNegInf,
                           .newObject = Object_newObject,
                           .endObject = Object_endObject,
                           .newArray = Object_newArray,
                           .endArray = Object_endArray,
                           .newInt = Object_newInteger,
                           .newLong = Object_newLong,
                           .newUnsignedLong = Object_newUnsignedLong,
                           .newDouble = Object_newDouble,
                           .releaseObject = Object_releaseObject,
                           .malloc = PyObject_Malloc,
                           .free = PyObject_Free,
                           .realloc = PyObject_Realloc,
                           .errorStr = NULL,
                           .errorOffset = NULL,
                           .preciseFloat = 0,
                           .prv = NULL};

  char *kwlist[] = {"obj", "precise_float", NULL};
  char *buf;
  Py_ssize_t len;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s#|b", kwlist, &buf, &len,
                                   &dec.preciseFloat)) {
    return NULL;
  }

  PyObject *ret = JSON_DecodeObject(&dec, buf, len);

  if (PyErr_Occurred()) {
    if (ret) {
      Py_DECREF((PyObject *)ret);
    }
    return NULL;
  }

  if (dec.errorStr) {
    /*
    FIXME: It's possible to give a much nicer error message here with actual
    failing element in input etc*/

    PyErr_Format(PyExc_ValueError, "%s", dec.errorStr);

    if (ret) {
      Py_DECREF((PyObject *)ret);
    }

    return NULL;
  }

  return ret;
}
