/*

Copyright (c) 2023, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

*/
#define _PANDAS_PARSER_IMPL

#include "pandas/parser/pd_parser.h"
#include "pandas/parser/io.h"

static int to_double(const char *item, const char *end, double *p_value,
                     char sci, char decimal, int *maybe_int) {
  char *p_end = NULL;
  int error = 0;

  /* Switch to precise xstrtod GH 31364 */
  *p_value = precise_xstrtod_with_end(item, &p_end, decimal, sci, '\0', 1,
                                      &error, maybe_int, end);

  // The value may contain an embedded NUL, so it is fully consumed only when
  // the parse reached `end`, not when it reached a NUL (GH#66524).
  return (error == 0) && (p_end == end);
}

static int floatify(PyObject *str, double *result, int *maybe_int) {
  const char *data;
  Py_ssize_t length;
  const char sci = 'E';
  const char dec = '.';

  if (PyBytes_Check(str)) {
    data = PyBytes_AS_STRING(str);
    length = PyBytes_GET_SIZE(str);
  } else if (PyUnicode_Check(str)) {
    // PyUnicode_AsUTF8AndSize returns a pointer to the string's internal UTF-8
    // buffer (caching it on the object), avoiding a heap allocation.
    data = PyUnicode_AsUTF8AndSize(str, &length);
    if (data == NULL) {
      return -1;
    }
  } else {
    PyErr_SetString(PyExc_TypeError, "Invalid object type");
    return -1;
  }

  const int status =
      to_double(data, data + length, result, sci, dec, maybe_int);

  if (!status) {
    /* handle inf/-inf infinity/-infinity */
    const int sign = infinity_sign(data, length);
    if (sign == 0) {
      goto parsingerror;
    }
    *result = sign > 0 ? HUGE_VAL : -HUGE_VAL;
    *maybe_int = 0;
  }

  return 0;

parsingerror:;
  // Report the value through its repr rather than with "%s", which stops at an
  // embedded NUL and so would name a truncated value that parses fine
  // (GH#66524). Rebuild an exact str/bytes rather than repr-ing `str` itself,
  // whose subclasses carry their type into the repr (numpy scalars render as
  // `np.str_('a')`); this stays lossless for bytes, which no decode would be.
  PyObject *display = PyBytes_Check(str)
                          ? PyBytes_FromStringAndSize(data, length)
                          : PyUnicode_FromStringAndSize(data, length);
  if (display == NULL) {
    return -1;
  }
  PyErr_Format(PyExc_ValueError, "Unable to parse string %R", display);
  Py_DECREF(display);
  return -1;
}

static void pandas_parser_destructor(PyObject *op) {
  void *ptr = PyCapsule_GetPointer(op, PandasParser_CAPSULE_NAME);
  PyMem_Free(ptr);
}

static int pandas_parser_exec(PyObject *Py_UNUSED(module)) {
  PandasParser_CAPI *capi = PyMem_Malloc(sizeof(PandasParser_CAPI));
  if (capi == NULL) {
    PyErr_NoMemory();
    return -1;
  }

  capi->to_double = to_double;
  capi->floatify = floatify;
  capi->new_rd_source = new_rd_source;
  capi->del_rd_source = del_rd_source;
  capi->buffer_rd_bytes = buffer_rd_bytes;
  capi->uint_state_init = uint_state_init;
  capi->uint64_conflict = uint64_conflict;
  capi->coliter_setup = coliter_setup;
  capi->parser_new = parser_new;
  capi->parser_init = parser_init;
  capi->parser_free = parser_free;
  capi->parser_del = parser_del;
  capi->parser_add_skiprow = parser_add_skiprow;
  capi->parser_set_skipfirstnrows = parser_set_skipfirstnrows;
  capi->parser_set_default_options = parser_set_default_options;
  capi->parser_consume_rows = parser_consume_rows;
  capi->parser_trim_buffers = parser_trim_buffers;
  capi->tokenize_all_rows = tokenize_all_rows;
  capi->tokenize_nrows = tokenize_nrows;
  capi->str_to_int64 = str_to_int64;
  capi->str_to_uint64 = str_to_uint64;
  capi->precise_xstrtod = precise_xstrtod;
  capi->to_boolean = to_boolean;

  PyObject *capsule =
      PyCapsule_New(capi, PandasParser_CAPSULE_NAME, pandas_parser_destructor);
  if (capsule == NULL) {
    PyMem_Free(capi);
    return -1;
  }

  // Monkeypatch the top level pandas module to have an attribute for the
  // C-API. This is required because Python capsules do not support setting
  // this attribute on anything but the top level package. Ideally not
  // done when cpython gh-6898 gets implemented
  PyObject *pandas = PyImport_ImportModule("pandas");
  if (!pandas) {
    PyErr_SetString(PyExc_ImportError,
                    "pd_parser.c could not import module pandas");
    Py_DECREF(capsule);
    return -1;
  }

  if (PyModule_AddObject(pandas, "_pandas_parser_CAPI", capsule) < 0) {
    Py_DECREF(capsule);
    return -1;
  }

  return 0;
}

static PyModuleDef_Slot pandas_parser_slots[] = {
    {Py_mod_exec, pandas_parser_exec},
#if PY_VERSION_HEX >= 0x030D0000
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL},
};

static struct PyModuleDef pandas_parsermodule = {
    PyModuleDef_HEAD_INIT,
    .m_name = "pandas._libs.pandas_parser",

    .m_doc = "Internal module with parser support for other extensions",
    .m_size = 0,
    .m_methods = NULL,
    .m_slots = pandas_parser_slots};

PyMODINIT_FUNC PyInit_pandas_parser(void) {
  return PyModuleDef_Init(&pandas_parsermodule);
}
