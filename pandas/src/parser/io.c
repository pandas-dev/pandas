#include "io.h"

 /*
   On-disk FILE, uncompressed
  */


void *new_file_source(char *fname, size_t buffer_size) {
    file_source *fs = (file_source *) malloc(sizeof(file_source));
    fs->fp = fopen(fname, "rb");

    if (fs->fp == NULL) {
        free(fs);
        return NULL;
    }
    setbuf(fs->fp, NULL);

    fs->initial_file_pos = ftell(fs->fp);

    // Only allocate this heap memory if we are not memory-mapping the file
    fs->buffer = (char*) malloc((buffer_size + 1) * sizeof(char));

    if (fs->buffer == NULL) {
        return NULL;
    }

    memset(fs->buffer, 0, buffer_size + 1);
    fs->buffer[buffer_size] = '\0';

    return (void *) fs;
}


// XXX handle on systems without the capability


/*
 *  void *new_file_buffer(FILE *f, int buffer_size)
 *
 *  Allocate a new file_buffer.
 *  Returns NULL if the memory allocation fails or if the call to mmap fails.
 *
 *  buffer_size is ignored.
 */


void* new_rd_source(PyObject *obj) {
    rd_source *rds = (rd_source *) malloc(sizeof(rd_source));

    /* hold on to this object */
    Py_INCREF(obj);
    rds->obj = obj;
    rds->buffer = NULL;
    rds->position = 0;

    return (void*) rds;
}

/*

  Cleanup callbacks

 */

int del_file_source(void *fs) {
    // fseek(FS(fs)->fp, FS(fs)->initial_file_pos, SEEK_SET);
    if (fs == NULL)
        return 0;

    /* allocated on the heap */
    free(FS(fs)->buffer);
    fclose(FS(fs)->fp);
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


void* buffer_file_bytes(void *source, size_t nbytes,
                        size_t *bytes_read, int *status) {
    file_source *src = FS(source);

    *bytes_read = fread((void*) src->buffer, sizeof(char), nbytes,
                        src->fp);

    if (*bytes_read == 0) {
        *status = REACHED_EOF;
    } else {
        *status = 0;
    }

    return (void*) src->buffer;

}


void* buffer_rd_bytes(void *source, size_t nbytes,
                      size_t *bytes_read, int *status) {
    PyGILState_STATE state;
    PyObject *result, *func, *args, *tmp;

    void *retval;

    size_t length;
    rd_source *src = RDS(source);

    /* delete old object */
    Py_XDECREF(src->buffer);
    args = Py_BuildValue("(i)", nbytes);

    state = PyGILState_Ensure();
    func = PyObject_GetAttrString(src->obj, "read");
    /* printf("%s\n", PyBytes_AsString(PyObject_Repr(func))); */

    /* TODO: does this release the GIL? */
    result = PyObject_CallObject(func, args);
    Py_XDECREF(args);
    Py_XDECREF(func);

    /* PyObject_Print(PyObject_Type(result), stdout, 0); */
    if (result == NULL) {
        PyGILState_Release(state);
        *bytes_read = 0;
        *status = CALLING_READ_FAILED;
        return NULL;
    }
    else if (!PyBytes_Check(result)) {
        tmp = PyUnicode_AsUTF8String(result);
        Py_XDECREF(result);
        result = tmp;
    }

    length = PySequence_Length(result);

    if (length == 0)
        *status = REACHED_EOF;
    else
        *status = 0;

    /* hang on to the Python object */
    src->buffer = result;
    retval = (void*) PyBytes_AsString(result);


    PyGILState_Release(state);

    /* TODO: more error handling */
    *bytes_read = length;

    return retval;
}


#ifdef HAVE_MMAP

#include <sys/stat.h>
#include <sys/mman.h>

void *new_mmap(char *fname)
{
    struct stat buf;
    int fd;
    memory_map *mm;
    /* off_t position; */
    off_t filesize;

    mm = (memory_map *) malloc(sizeof(memory_map));
    mm->fp = fopen(fname, "rb");

    fd = fileno(mm->fp);
    if (fstat(fd, &buf) == -1) {
        fprintf(stderr, "new_file_buffer: fstat() failed. errno =%d\n", errno);
        return NULL;
    }
    filesize = buf.st_size;  /* XXX This might be 32 bits. */


    if (mm == NULL) {
        /* XXX Eventually remove this print statement. */
        fprintf(stderr, "new_file_buffer: malloc() failed.\n");
        return NULL;
    }
    mm->size = (off_t) filesize;
    mm->line_number = 0;

    mm->fileno = fd;
    mm->position = ftell(mm->fp);
    mm->last_pos = (off_t) filesize;

    mm->memmap = mmap(NULL, filesize, PROT_READ, MAP_SHARED, fd, 0);
    if (mm->memmap == NULL) {
        /* XXX Eventually remove this print statement. */
        fprintf(stderr, "new_file_buffer: mmap() failed.\n");
        free(mm);
        mm = NULL;
    }

    return (void*) mm;
}


int del_mmap(void *src)
{
    munmap(MM(src)->memmap, MM(src)->size);

    fclose(MM(src)->fp);

    /*
     *  With a memory mapped file, there is no need to do
     *  anything if restore == RESTORE_INITIAL.
     */
    /* if (restore == RESTORE_FINAL) { */
    /*     fseek(FB(fb)->file, FB(fb)->current_pos, SEEK_SET); */
    /* } */
    free(src);

    return 0;
}

void* buffer_mmap_bytes(void *source, size_t nbytes,
                        size_t *bytes_read, int *status) {
    void *retval;
    memory_map *src = MM(source);

    if (src->position == src->last_pos) {
        *bytes_read = 0;
        *status = REACHED_EOF;
        return NULL;
    }

    retval = src->memmap + src->position;

    if (src->position + nbytes > src->last_pos) {
        // fewer than nbytes remaining
        *bytes_read = src->last_pos - src->position;
    } else {
        *bytes_read = nbytes;
    }

    *status = 0;

    /* advance position in mmap data structure */
    src->position += *bytes_read;

    return retval;
}

#else

/* kludgy */

void *new_mmap(char *fname) {
  return NULL;
}

int del_mmap(void *src) {
  return 0;
}

/* don't use this! */

void* buffer_mmap_bytes(void *source, size_t nbytes,
                        size_t *bytes_read, int *status) {
  return NULL;
}

#endif
