/*
Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#ifndef PANDAS__LIBS_SRC_PARSER_IO_H_
#define PANDAS__LIBS_SRC_PARSER_IO_H_

#include "Python.h"
#include "tokenizer.h"

typedef struct _file_source {
    /* The file being read. */
    int fd;

    char *buffer;
    size_t size;
} file_source;

#define FS(source) ((file_source *)source)

#if !defined(_WIN32) && !defined(HAVE_MMAP)
#define HAVE_MMAP
#endif

typedef struct _memory_map {
    int fd;

    /* Size of the file, in bytes. */
    char *memmap;
    size_t size;

    size_t position;
} memory_map;

#define MM(src) ((memory_map *)src)

void *new_mmap(char *fname);

int del_mmap(void *src);

void *buffer_mmap_bytes(void *source, size_t nbytes, size_t *bytes_read,
                        int *status);

typedef struct _rd_source {
    PyObject *obj;
    PyObject *buffer;
    size_t position;
} rd_source;

#define RDS(source) ((rd_source *)source)

void *new_file_source(char *fname, size_t buffer_size);

void *new_rd_source(PyObject *obj);

int del_file_source(void *src);
int del_rd_source(void *src);

void *buffer_file_bytes(void *source, size_t nbytes, size_t *bytes_read,
                        int *status);

void *buffer_rd_bytes(void *source, size_t nbytes, size_t *bytes_read,
                      int *status);

#endif  // PANDAS__LIBS_SRC_PARSER_IO_H_
