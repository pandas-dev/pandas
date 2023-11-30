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

#define PY_ARRAY_UNIQUE_SYMBOL UJSON_NUMPY
#define NO_IMPORT_ARRAY
#define PY_SSIZE_T_CLEAN
#include "pandas/vendored/ujson/lib/ultrajson.h"
#include <Python.h>
#include <numpy/arrayobject.h>

typedef struct __PyObjectDecoder {
  JSONObjectDecoder dec;

  void *npyarr;      // Numpy context buffer
  void *npyarr_addr; // Ref to npyarr ptr to track DECREF calls
} PyObjectDecoder;

typedef struct __NpyArrContext {
  PyObject *ret;
  PyObject *labels[2];
  PyArray_Dims shape;

  PyObjectDecoder *dec;
} NpyArrContext;

// free the numpy context buffer
void Npy_releaseContext(NpyArrContext *npyarr) {
  if (npyarr) {
    if (npyarr->shape.ptr) {
      PyObject_Free(npyarr->shape.ptr);
    }
    if (npyarr->dec) {
      npyarr->dec->npyarr = NULL;
    }
    Py_XDECREF(npyarr->labels[0]);
    Py_XDECREF(npyarr->labels[1]);
    Py_XDECREF(npyarr->ret);
    PyObject_Free(npyarr);
  }
}

static int Object_objectAddKey(void *prv, JSOBJ obj, JSOBJ name, JSOBJ value) {
  int ret = PyDict_SetItem(obj, name, value);
  Py_DECREF((PyObject *)name);
  Py_DECREF((PyObject *)value);
  return ret == 0 ? 1 : 0;
}

static int Object_arrayAddItem(void *prv, JSOBJ obj, JSOBJ value) {
  int ret = PyList_Append(obj, value);
  Py_DECREF((PyObject *)value);
  return ret == 0 ? 1 : 0;
}

static JSOBJ Object_newString(void *prv, wchar_t *start, wchar_t *end) {
  return PyUnicode_FromWideChar(start, (end - start));
}

static JSOBJ Object_newTrue(void *prv) { Py_RETURN_TRUE; }

static JSOBJ Object_newFalse(void *prv) { Py_RETURN_FALSE; }

static JSOBJ Object_newNull(void *prv) { Py_RETURN_NONE; }

static JSOBJ Object_newPosInf(void *prv) {
  return PyFloat_FromDouble(Py_HUGE_VAL);
}

static JSOBJ Object_newNegInf(void *prv) {
  return PyFloat_FromDouble(-Py_HUGE_VAL);
}

static JSOBJ Object_newObject(void *prv, void *decoder) { return PyDict_New(); }

static JSOBJ Object_endObject(void *prv, JSOBJ obj) { return obj; }

static JSOBJ Object_newArray(void *prv, void *decoder) { return PyList_New(0); }

static JSOBJ Object_endArray(void *prv, JSOBJ obj) { return obj; }

static JSOBJ Object_newInteger(void *prv, JSINT32 value) {
  return PyLong_FromLong((long)value);
}

static JSOBJ Object_newLong(void *prv, JSINT64 value) {
  return PyLong_FromLongLong(value);
}

static JSOBJ Object_newUnsignedLong(void *prv, JSUINT64 value) {
  return PyLong_FromUnsignedLongLong(value);
}

static JSOBJ Object_newDouble(void *prv, double value) {
  return PyFloat_FromDouble(value);
}

static void Object_releaseObject(void *prv, JSOBJ obj, void *_decoder) {
  PyObjectDecoder *decoder = (PyObjectDecoder *)_decoder;
  if (obj != decoder->npyarr_addr) {
    Py_XDECREF(((PyObject *)obj));
  }
}

static char *g_kwlist[] = {"obj", "precise_float", "labelled", "dtype", NULL};

PyObject *JSONToObj(PyObject *self, PyObject *args, PyObject *kwargs) {
  PyObject *ret;
  PyObject *sarg;
  PyObject *arg;
  PyObject *opreciseFloat = NULL;
  JSONObjectDecoder *decoder;
  PyObjectDecoder pyDecoder;
  PyArray_Descr *dtype = NULL;
  int labelled = 0;

  JSONObjectDecoder dec = {
      Object_newString,  Object_objectAddKey,  Object_arrayAddItem,
      Object_newTrue,    Object_newFalse,      Object_newNull,
      Object_newPosInf,  Object_newNegInf,     Object_newObject,
      Object_endObject,  Object_newArray,      Object_endArray,
      Object_newInteger, Object_newLong,       Object_newUnsignedLong,
      Object_newDouble,  Object_releaseObject, PyObject_Malloc,
      PyObject_Free,     PyObject_Realloc};

  dec.preciseFloat = 0;
  dec.prv = NULL;

  pyDecoder.dec = dec;
  pyDecoder.npyarr = NULL;
  pyDecoder.npyarr_addr = NULL;

  decoder = (JSONObjectDecoder *)&pyDecoder;

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|OiiO&", g_kwlist, &arg,
                                   &opreciseFloat, &labelled,
                                   PyArray_DescrConverter2, &dtype)) {
    Npy_releaseContext(pyDecoder.npyarr);
    return NULL;
  }

  if (opreciseFloat && PyObject_IsTrue(opreciseFloat)) {
    decoder->preciseFloat = 1;
  }

  if (PyBytes_Check(arg)) {
    sarg = arg;
  } else if (PyUnicode_Check(arg)) {
    sarg = PyUnicode_AsUTF8String(arg);
    if (sarg == NULL) {
      // Exception raised above us by codec according to docs
      return NULL;
    }
  } else {
    PyErr_Format(PyExc_TypeError, "Expected 'str' or 'bytes'");
    return NULL;
  }

  decoder->errorStr = NULL;
  decoder->errorOffset = NULL;

  ret = JSON_DecodeObject(decoder, PyBytes_AS_STRING(sarg),
                          PyBytes_GET_SIZE(sarg));

  if (sarg != arg) {
    Py_DECREF(sarg);
  }

  if (PyErr_Occurred()) {
    if (ret) {
      Py_DECREF((PyObject *)ret);
    }
    Npy_releaseContext(pyDecoder.npyarr);
    return NULL;
  }

  if (decoder->errorStr) {
    /*
    FIXME: It's possible to give a much nicer error message here with actual
    failing element in input etc*/

    PyErr_Format(PyExc_ValueError, "%s", decoder->errorStr);

    if (ret) {
      Py_DECREF((PyObject *)ret);
    }
    Npy_releaseContext(pyDecoder.npyarr);

    return NULL;
  }

  return ret;
}
