/*
Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#include "io.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#ifndef O_BINARY
#define O_BINARY 0
#endif  // O_BINARY

#if PY_VERSION_HEX >= 0x03060000 && defined(_WIN32)
#define USE_WIN_UTF16
#include <Windows.h>
#endif

/*
  On-disk FILE, uncompressed
*/

void *new_file_source(char *fname, size_t buffer_size) {
    file_source *fs = (file_source *)malloc(sizeof(file_source));
    if (fs == NULL) {
        return NULL;
    }

#ifdef USE_WIN_UTF16
    // Fix gh-15086 properly - convert UTF8 to UTF16 that Windows widechar API
    // accepts. This is needed because UTF8 might _not_ be convertible to MBCS
    // for some conditions, as MBCS is locale-dependent, and not all unicode
    // symbols can be expressed in it.
    {
        wchar_t* wname = NULL;
        int required = MultiByteToWideChar(CP_UTF8, 0, fname, -1, NULL, 0);
        if (required == 0) {
            free(fs);
            return NULL;
        }
        wname = (wchar_t*)malloc(required * sizeof(wchar_t));
        if (wname == NULL) {
            free(fs);
            return NULL;
        }
        if (MultiByteToWideChar(CP_UTF8, 0, fname, -1, wname, required) <
                                                                required) {
            free(wname);
            free(fs);
            return NULL;
        }
        fs->fd = _wopen(wname, O_RDONLY | O_BINARY);
        free(wname);
    }
#else
    fs->fd = open(fname, O_RDONLY | O_BINARY);
#endif
    if (fs->fd == -1) {
        free(fs);
        return NULL;
    }

    // Only allocate this heap memory if we are not memory-mapping the file
    fs->buffer = (char *)malloc((buffer_size + 1) * sizeof(char));

    if (fs->buffer == NULL) {
        close(fs->fd);
        free(fs);
        return NULL;
    }

    memset(fs->buffer, '\0', buffer_size + 1);
    fs->size = buffer_size;

    return (void *)fs;
}

void *new_rd_source(PyObject *obj) {
    rd_source *rds = (rd_source *)malloc(sizeof(rd_source));

    /* hold on to this object */
    Py_INCREF(obj);
    rds->obj = obj;
    rds->buffer = NULL;
    rds->position = 0;

    return (void *)rds;
}

/*

  Cleanup callbacks

 */

int del_file_source(void *ptr) {
    file_source *fs = ptr;
    if (fs == NULL) return 0;

    free(fs->buffer);
    close(fs->fd);
    free(fs);

    return 0;
}

int del_rd_source(void *rds) {
    Py_XDECREF(RDS(rds)->obj);
    Py_XDECREF(RDS(rds)->buffer);
    free(rds);

    return 0;
}

/*

  IO callbacks

 */

void *buffer_file_bytes(void *source, size_t nbytes, size_t *bytes_read,
                        int *status) {
    file_source *fs = FS(source);
    ssize_t rv;

    if (nbytes > fs->size) {
        nbytes = fs->size;
    }

    rv = read(fs->fd, fs->buffer, nbytes);
    switch (rv) {
    case -1:
        *status = CALLING_READ_FAILED;
        *bytes_read = 0;
        return NULL;
    case 0:
        *status = REACHED_EOF;
        *bytes_read = 0;
        return NULL;
    default:
        *status = 0;
        *bytes_read = rv;
        fs->buffer[rv] = '\0';
        break;
    }

    return (void *)fs->buffer;
}

void *buffer_rd_bytes(void *source, size_t nbytes, size_t *bytes_read,
                      int *status) {
    PyGILState_STATE state;
    PyObject *result, *func, *args, *tmp;

    void *retval;

    size_t length;
    rd_source *src = RDS(source);
    state = PyGILState_Ensure();

    /* delete old object */
    Py_XDECREF(src->buffer);
    src->buffer = NULL;
    args = Py_BuildValue("(i)", nbytes);

    func = PyObject_GetAttrString(src->obj, "read");

    /* TODO: does this release the GIL? */
    result = PyObject_CallObject(func, args);
    Py_XDECREF(args);
    Py_XDECREF(func);

    if (result == NULL) {
        PyGILState_Release(state);
        *bytes_read = 0;
        *status = CALLING_READ_FAILED;
        return NULL;
    } else if (!PyBytes_Check(result)) {
        tmp = PyUnicode_AsUTF8String(result);
        Py_DECREF(result);
        if (tmp == NULL) {
            PyGILState_Release(state);
            return NULL;
        }
        result = tmp;
    }

    length = PySequence_Length(result);

    if (length == 0)
        *status = REACHED_EOF;
    else
        *status = 0;

    /* hang on to the Python object */
    src->buffer = result;
    retval = (void *)PyBytes_AsString(result);

    PyGILState_Release(state);

    /* TODO: more error handling */
    *bytes_read = length;

    return retval;
}

#ifdef HAVE_MMAP

#include <sys/mman.h>

void *new_mmap(char *fname) {
    memory_map *mm;
    struct stat stat;
    size_t filesize;

    mm = (memory_map *)malloc(sizeof(memory_map));
    if (mm == NULL) {
        fprintf(stderr, "new_file_buffer: malloc() failed.\n");
        return (NULL);
    }
    mm->fd = open(fname, O_RDONLY | O_BINARY);
    if (mm->fd == -1) {
        fprintf(stderr, "new_file_buffer: open(%s) failed. errno =%d\n",
          fname, errno);
        free(mm);
        return NULL;
    }

    if (fstat(mm->fd, &stat) == -1) {
        fprintf(stderr, "new_file_buffer: fstat() failed. errno =%d\n",
          errno);
        close(mm->fd);
        free(mm);
        return NULL;
    }
    filesize = stat.st_size; /* XXX This might be 32 bits. */

    mm->memmap = mmap(NULL, filesize, PROT_READ, MAP_SHARED, mm->fd, 0);
    if (mm->memmap == MAP_FAILED) {
        /* XXX Eventually remove this print statement. */
        fprintf(stderr, "new_file_buffer: mmap() failed.\n");
        close(mm->fd);
        free(mm);
        return NULL;
    }

    mm->size = (off_t)filesize;
    mm->position = 0;

    return mm;
}

int del_mmap(void *ptr) {
    memory_map *mm = ptr;

    if (mm == NULL) return 0;

    munmap(mm->memmap, mm->size);
    close(mm->fd);
    free(mm);

    return 0;
}

void *buffer_mmap_bytes(void *source, size_t nbytes, size_t *bytes_read,
                        int *status) {
    void *retval;
    memory_map *src = source;
    size_t remaining = src->size - src->position;

    if (remaining == 0) {
        *bytes_read = 0;
        *status = REACHED_EOF;
        return NULL;
    }

    if (nbytes > remaining) {
        nbytes = remaining;
    }

    retval = src->memmap + src->position;

    /* advance position in mmap data structure */
    src->position += nbytes;

    *bytes_read = nbytes;
    *status = 0;

    return retval;
}

#else

/* kludgy */

void *new_mmap(char *fname) { return NULL; }

int del_mmap(void *src) { return 0; }

/* don't use this! */

void *buffer_mmap_bytes(void *source, size_t nbytes, size_t *bytes_read,
                        int *status) {
    return NULL;
}

#endif  // HAVE_MMAP
