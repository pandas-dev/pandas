#ifndef NUMBA_PY_MODULE_H_
#define NUMBA_PY_MODULE_H_

#define PY_SSIZE_T_CLEAN

#include "Python.h"
#include "structmember.h"
#include "frameobject.h"

#define MOD_ERROR_VAL NULL
#define MOD_SUCCESS_VAL(val) val
#define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
#define MOD_DEF(ob, name, doc, methods) { \
        static struct PyModuleDef moduledef = { \
          PyModuleDef_HEAD_INIT, name, doc, -1, methods, NULL, NULL, NULL, NULL }; \
        ob = PyModule_Create(&moduledef); }
#define MOD_INIT_EXEC(name) PyInit_##name();

#define PyString_AsString PyUnicode_AsUTF8
#define PyString_Check PyUnicode_Check
#define PyString_FromFormat PyUnicode_FromFormat
#define PyString_FromString PyUnicode_FromString
#define PyString_InternFromString PyUnicode_InternFromString
#define PyInt_Type PyLong_Type
#define PyInt_Check PyLong_Check
#define PyInt_CheckExact PyLong_CheckExact
#define SetAttrStringFromVoidPointer(m, name) do { \
        PyObject *tmp = PyLong_FromVoidPtr((void *) &name); \
        PyObject_SetAttrString(m, #name, tmp); \
        Py_DECREF(tmp); } while (0)


#define NB_SUPPORTED_PYTHON_MINOR ((PY_MINOR_VERSION == 10) || (PY_MINOR_VERSION == 11) || (PY_MINOR_VERSION == 12) || (PY_MINOR_VERSION == 13))

#endif /* NUMBA_PY_MODULE_H_ */
