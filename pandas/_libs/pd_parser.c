/*

Copyright (c) 2023, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

*/
#define _PANDAS_PARSER_IMPL

#include "pd_parser.h"
#include "src/parser/io.h"

static int pandas_to_double(char *item, double *p_value, char sci, char decimal,
                            int *maybe_int) {
  char *p_end = NULL;
  int error = 0;

  /* Switch to precise xstrtod GH 31364 */
  *p_value =
      precise_xstrtod(item, &p_end, decimal, sci, '\0', 1, &error, maybe_int);

  return (error == 0) && (!*p_end);
}

static int pandas_floatify(PyObject *str, double *result, int *maybe_int) {
  int status;
  char *data;
  PyObject *tmp = NULL;
  const char sci = 'E';
  const char dec = '.';

  if (PyBytes_Check(str)) {
    data = PyBytes_AS_STRING(str);
  } else if (PyUnicode_Check(str)) {
    tmp = PyUnicode_AsUTF8String(str);
    if (tmp == NULL) {
      return -1;
    }
    data = PyBytes_AS_STRING(tmp);
  } else {
    PyErr_SetString(PyExc_TypeError, "Invalid object type");
    return -1;
  }

  status = pandas_to_double(data, result, sci, dec, maybe_int);

  if (!status) {
    /* handle inf/-inf infinity/-infinity */
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
    } else if (strlen(data) == 8) {
      if (0 == strcasecmp(data, "infinity")) {
        *result = HUGE_VAL;
        *maybe_int = 0;
      } else {
        goto parsingerror;
      }
    } else if (strlen(data) == 9) {
      if (0 == strcasecmp(data, "-infinity")) {
        *result = -HUGE_VAL;
        *maybe_int = 0;
      } else if (0 == strcasecmp(data, "+infinity")) {
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

static void *pandas_new_rd_source(PyObject *obj) { return new_rd_source(obj); }

static int pandas_del_rd_source(void *src) { return del_rd_source(src); }

static void *pandas_buffer_rd_bytes(void *source, size_t nbytes,
                                    size_t *bytes_read, int *status,
                                    const char *encoding_errors) {
  return buffer_rd_bytes(source, nbytes, bytes_read, status, encoding_errors);
}

static void pandas_uint_state_init(uint_state *self) { uint_state_init(self); }

static int pandas_uint64_conflict(uint_state *self) {
  return uint64_conflict(self);
}

static void pandas_coliter_setup(coliter_t *self, parser_t *parser, int64_t i,
                                 int64_t start) {
  coliter_setup(self, parser, i, start);
}

static parser_t *pandas_parser_new() { return parser_new(); }

static int pandas_parser_init(parser_t *self) { return parser_init(self); }

static void pandas_parser_free(parser_t *self) { parser_free(self); }

static void pandas_parser_del(parser_t *self) { parser_del(self); }

static int pandas_parser_add_skiprow(parser_t *self, int64_t row) {
  return parser_add_skiprow(self, row);
}

static int pandas_parser_set_skipfirstnrows(parser_t *self, int64_t nrows) {
  return parser_set_skipfirstnrows(self, nrows);
}

static void pandas_parser_set_default_options(parser_t *self) {
  parser_set_default_options(self);
}

static int pandas_parser_consume_rows(parser_t *self, size_t nrows) {
  return parser_consume_rows(self, nrows);
}

static int pandas_parser_trim_buffers(parser_t *self) {
  return parser_trim_buffers(self);
}

static int pandas_tokenize_all_rows(parser_t *self,
                                    const char *encoding_errors) {
  return tokenize_all_rows(self, encoding_errors);
}

static int pandas_tokenize_nrows(parser_t *self, size_t nrows,
                                 const char *encoding_errors) {
  return tokenize_nrows(self, nrows, encoding_errors);
}

static int64_t pandas_str_to_int64(const char *p_item, int64_t int_min,
                                   int64_t int_max, int *error, char tsep) {
  return str_to_int64(p_item, int_min, int_max, error, tsep);
}

static uint64_t pandas_str_to_uint64(uint_state *state, const char *p_item,
                                     int64_t int_max, uint64_t uint_max,
                                     int *error, char tsep) {
  return str_to_uint64(state, p_item, int_max, uint_max, error, tsep);
}

static double pandas_xstrtod(const char *p, char **q, char decimal, char sci,
                             char tsep, int skip_trailing, int *error,
                             int *maybe_int) {
  return xstrtod(p, q, decimal, sci, tsep, skip_trailing, error, maybe_int);
}

static double pandas_precise_xstrtod(const char *p, char **q, char decimal,
                                     char sci, char tsep, int skip_trailing,
                                     int *error, int *maybe_int) {
  return precise_xstrtod(p, q, decimal, sci, tsep, skip_trailing, error,
                         maybe_int);
}

static double pandas_round_trip(const char *p, char **q, char decimal, char sci,
                                char tsep, int skip_trailing, int *error,
                                int *maybe_int) {
  return round_trip(p, q, decimal, sci, tsep, skip_trailing, error, maybe_int);
}

static int pandas_to_boolean(const char *item, uint8_t *val) {
  return to_boolean(item, val);
}

static void pandas_parser_destructor(PyObject *op) {
  void *ptr = PyCapsule_GetPointer(op, PandasParser_CAPSULE_NAME);
  PyMem_Free(ptr);
}

static int pandas_parser_exec(PyObject *module) {
  PandasParser_CAPI *capi = PyMem_Malloc(sizeof(PandasParser_CAPI));
  if (capi == NULL) {
    PyErr_NoMemory();
    return -1;
  }

  capi->to_double = pandas_to_double;
  capi->floatify = pandas_floatify;
  capi->new_rd_source = pandas_new_rd_source;
  capi->del_rd_source = pandas_del_rd_source;
  capi->buffer_rd_bytes = pandas_buffer_rd_bytes;
  capi->uint_state_init = pandas_uint_state_init;
  capi->uint64_conflict = pandas_uint64_conflict;
  capi->coliter_setup = pandas_coliter_setup;
  capi->parser_new = pandas_parser_new;
  capi->parser_init = pandas_parser_init;
  capi->parser_free = pandas_parser_free;
  capi->parser_del = pandas_parser_del;
  capi->parser_add_skiprow = pandas_parser_add_skiprow;
  capi->parser_set_skipfirstnrows = pandas_parser_set_skipfirstnrows;
  capi->parser_set_default_options = pandas_parser_set_default_options;
  capi->parser_consume_rows = pandas_parser_consume_rows;
  capi->parser_trim_buffers = pandas_parser_trim_buffers;
  capi->tokenize_all_rows = pandas_tokenize_all_rows;
  capi->tokenize_nrows = pandas_tokenize_nrows;
  capi->str_to_int64 = pandas_str_to_int64;
  capi->str_to_uint64 = pandas_str_to_uint64;
  capi->xstrtod = pandas_xstrtod;
  capi->precise_xstrtod = pandas_precise_xstrtod;
  capi->round_trip = pandas_round_trip;
  capi->to_boolean = pandas_to_boolean;

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
    {Py_mod_exec, pandas_parser_exec}, {0, NULL}};

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
