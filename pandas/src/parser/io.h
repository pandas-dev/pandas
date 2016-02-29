#include "Python.h"
#include "tokenizer.h"


typedef struct _file_source {
    /* The file being read. */
    FILE *fp;

    char *buffer;
    /* Size of the file, in bytes. */
    /* off_t size; */

    /* file position when the file_buffer was created. */
    off_t initial_file_pos;

    /* Offset in the file of the data currently in the buffer. */
    off_t buffer_file_pos;

    /* Actual number of bytes in the current buffer. (Can be less than buffer_size.) */
    off_t last_pos;

    /* Size (in bytes) of the buffer. */
    // off_t buffer_size;

    /* Pointer to the buffer. */
    // char *buffer;

} file_source;

#define FS(source) ((file_source *)source)

#if !defined(_WIN32) && !defined(HAVE_MMAP)
#define HAVE_MMAP
#endif

typedef struct _memory_map {

    FILE *fp;

    /* Size of the file, in bytes. */
    off_t size;

    /* file position when the file_buffer was created. */
    off_t initial_file_pos;

    int line_number;

    int fileno;
    off_t position;
    off_t last_pos;
    char *memmap;

} memory_map;

#define MM(src) ((memory_map*) src)

void *new_mmap(char *fname);

int del_mmap(void *src);

void* buffer_mmap_bytes(void *source, size_t nbytes,
                        size_t *bytes_read, int *status);


typedef struct _rd_source {
    PyObject* obj;
    PyObject* buffer;
    size_t position;
} rd_source;

#define RDS(source) ((rd_source *)source)

void *new_file_source(char *fname, size_t buffer_size);

void *new_rd_source(PyObject *obj);

int del_file_source(void *src);
int del_rd_source(void *src);

void* buffer_file_bytes(void *source, size_t nbytes,
                        size_t *bytes_read, int *status);

void* buffer_rd_bytes(void *source, size_t nbytes,
                      size_t *bytes_read, int *status);

