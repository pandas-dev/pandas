/*
Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#pragma once

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "tokenizer.h"

#define FS(source) ((file_source *)source)

typedef struct _rd_source {
  PyObject *obj;
  PyObject *buffer;
  size_t position;
} rd_source;

#define RDS(source) ((rd_source *)source)

void *new_rd_source(PyObject *obj);

void del_rd_source(void *src);

char *buffer_rd_bytes(void *source, size_t nbytes, size_t *bytes_read,
                      int *status, const char *encoding_errors);
