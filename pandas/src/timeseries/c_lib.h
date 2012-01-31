#ifndef C_LIB_H
#define C_LIB_H

#include <Python.h>
#include <structmember.h>
#include "numpy/arrayobject.h"

/* c_lib defines generic functions that aren't inherently time series/date
specific but are needed in various parts of the module. */

#define INT_ERR_CODE -999

#define MEM_CHECK(item) if (item == NULL) { return PyErr_NoMemory(); }
#define ERR_CHECK(item) if (item == NULL) { return NULL; }

char *str_uppercase(char *);
char *str_replace(const char*, const char*, const char*);

PyObject *np_add(PyObject*, PyObject*);
PyObject *np_multiply(PyObject*, PyObject*);
PyObject *np_subtract(PyObject*, PyObject*);
PyObject *np_sqrt(PyObject*);
int np_greater(PyObject*, PyObject*);
int np_greater_equal(PyObject*, PyObject*);

PyObject *set_callback(PyObject*, PyObject**);

void import_c_lib(PyObject*);

#endif
