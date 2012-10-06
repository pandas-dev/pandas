/*

Copyright (c) 2012, Lambda Foundry, Inc.
All rights reserved.

Released subject to the terms of the GNU Lesser Public License

*/

 /*
   Low-level ascii-file processing for pandas. Combines some elements from
   Python's built-in csv module and Warren Weckesser's textreader project on
   GitHub. See Python Software Foundation License and BSD licenses for these.

  */


 #include "Python.h"
 #include "structmember.h"

 #include "parser.h"
 #include "conversions.h"

 #define READ_ERROR_OUT_OF_MEMORY   1

 #define REACHED_EOF 1

 #define HAVE_MEMMAP
 #define HAVE_GZIP

 #define FB_EOF   -1
 #define FB_ERROR -2

 /*
  * restore:
  *  RESTORE_NOT     (0):
  *      Free memory, but leave the file position wherever it
  *      happend to be.
  *  RESTORE_INITIAL (1):
  *      Restore the file position to the location at which
  *      the file_buffer was created.
  *  RESTORE_FINAL   (2):
  *      Put the file position at the next byte after the
  *      data read from the file_buffer.
  */
 #define RESTORE_NOT     0
 #define RESTORE_INITIAL 1
 #define RESTORE_FINAL   2




void *safe_realloc(void *buffer, size_t size) {
    void *result;
    // OS X is weird
    // http://stackoverflow.com/questions/9560609/
    // different-realloc-behaviour-in-linux-and-osx

    result = realloc(buffer, size);

    if (result != NULL) {
        // errno gets set to 12 on my OS Xmachine in some cases even when the
        // realloc succeeds. annoying
        errno = 0;
    } else {
        return buffer;
    }

    return result;
}


void coliter_setup(coliter_t *self, parser_t *parser, int i, int start) {
    // column i, starting at 0
    self->words = parser->words;
    self->col = i;
    self->line_start = parser->line_start + start;
}

coliter_t *coliter_new(parser_t *self, int i) {
    // column i, starting at 0
    coliter_t *iter = (coliter_t*) malloc(sizeof(coliter_t));

    if (NULL == iter) {
        return NULL;
    }

    coliter_setup(iter, self, i, 0);
    return iter;
}



 /* int64_t str_to_int64(const char *p_item, int64_t int_min, int64_t int_max, int *error); */
 /* uint64_t str_to_uint64(const char *p_item, uint64_t uint_max, int *error); */


 void free_if_not_null(void *ptr) {
     if (ptr != NULL) free(ptr);
 }


 /*

   On-disk FILE, uncompressed

  */

 typedef struct _file_source {
     /* The file being read. */
     FILE *fp;

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


 void *new_file_source(FILE *fp) {
     file_source *fs = (file_source *) malloc(sizeof(file_source));
     fs->fp = fp;

     fs->initial_file_pos = ftell(fp);

     return (void *) fs;
 }

 void del_file_source(void *fs) {
     // TODO: error codes?
     // fclose(FS(fs)->fp);

     // fseek(FS(fs)->fp, FS(fs)->initial_file_pos, SEEK_SET);
     // allocated on the heap
     free(fs);
 }


int parser_file_source_init(parser_t *self, FILE* fp) {
    self->sourcetype = 'F';
    self->source = new_file_source(fp);

    // Only allocate this heap memory if we are not memory-mapping the file
    self->data = (char*) malloc((self->chunksize + 1) * sizeof(char));
    memset(self->data, 0, self->chunksize + 1);
    self->data[self->chunksize] = '\0';

    if (self->data == NULL) {
        return PARSER_OUT_OF_MEMORY;
    }

    return 0;
}

 /*

   In-memory bytes

  */


 typedef struct _array_source {
     char *data;
     size_t position, length;
 } array_source;

 #define ARS(source) ((array_source *)source)


 void *new_array_source(char *data, size_t length) {
     array_source *ars = (array_source *) malloc(sizeof(array_source));

     // to be safe, copy the data from the Python string
     ars->data = malloc(sizeof(char) * (length + 1));
     strcpy(ars->data, data);

     ars->position = 0;
     ars->length = length;

     return (void *) ars;
 }

 void del_array_source(void *ars) {
     // I made a copy
     free(ARS(ars)->data);

     free(ars);
 }

 /*

   Parser / tokenizer

 */


 void *grow_buffer(void *buffer, int length, int *capacity,
                   int space, int elsize, int *error) {
     int cap = *capacity;

     // Can we fit potentially nbytes tokens (+ null terminators) in the stream?
     while (length + space > cap) {
         cap = cap? cap << 1 : 2;

         buffer = safe_realloc(buffer, elsize * cap);

         if (buffer == NULL) {
             // TODO: error codes
             *error = -1;
         }
     }

     // sigh, multiple return values
     *capacity = cap;
     *error = 0;
     return buffer;
 }

 typedef struct _typed_array {
     char type_code;
     int elsize;
     size_t length;
     void *data;
 } typed_array;

 #define array_t typed_array



 /*
  *  XXX Handle errors in any of the functions called by read_rows().
  *
  *  XXX Currently *nrows must be at least 1.
  */

 /* void *read_rows(FILE *f, */
 /*              int *nrows, */
 /*              char *fmt, */
 /*                 char delimiter, */
 /*              char quote, */
 /*              char comment, */
 /*                 char sci, */
 /*              char decimal, */
 /*                 int allow_embedded_newline, */
 /*                 char *datetime_fmt, */
 /*                 int tz_offset, */
 /*                 int32_t *usecols, */
 /*              int num_usecols, */
 /*                 int skiprows, */
 /*                 void *data_array, */
 /*                 int *p_error_type, */
 /*              int *p_error_lineno) */
 /* { */
 /*     void *fb; */
 /*     char *data_ptr; */
 /*     int num_fields, current_num_fields; */
 /*     char **result; */
 /*     int fmt_nfields; */
 /*     field_type *ftypes; */
 /*     int size; */
 /*     int row_count; */
 /*     int j; */
 /*     int *valid_usecols; */
 /*     char word_buffer[WORD_BUFFER_SIZE]; */
 /*     int tok_error_type; */

 /*     *p_error_type = 0; */
 /*     *p_error_lineno = 0; */

 /*     if (datetime_fmt == NULL || strlen(datetime_fmt) == 0) { */
 /*         datetime_fmt = "%Y-%m-%d %H:%M:%S"; */
 /*     } */

 /*     size = (*nrows) * calc_size(fmt, &fmt_nfields); */

 /*     ftypes = enumerate_fields(fmt);  /\* Must free this when finished. *\/ */
 /*     if (ftypes == NULL) { */
 /*         /\* Out of memory. *\/ */
 /*         *p_error_type = READ_ERROR_OUT_OF_MEMORY; */
 /*         return NULL; */
 /*     } */

 /*     /\* */
 /*     for (k = 0; k < fmt_nfields; ++k) { */
 /*         printf("k = %d  typechar = '%c'  size = %d\n", k, ftypes[k].typechar, ftypes[k].size); */
 /*     } */
 /*     printf("size = %d\n", size); */
 /*     printf("-----\n"); */
 /*     *\/ */

 /*     if (data_array == NULL) { */
 /*         /\* XXX The case where data_ptr is allocated here is untested. *\/ */
 /*         data_ptr = malloc(size); */
 /*     } */
 /*     else { */
 /*         data_ptr = data_array; */
 /*     } */

 /*     fb = new_file_buffer(f, -1); */
 /*     if (fb == NULL) { */
 /*         free(ftypes); */
 /*         *p_error_type = ERROR_OUT_OF_MEMORY; */
 /*         return NULL; */
 /*     } */

 /*     /\* XXX Check interaction of skiprows with comments. *\/ */
 /*     while ((skiprows > 0) && ((result = tokenize(fb, word_buffer, WORD_BUFFER_SIZE, */
 /*                               delimiter, quote, comment, &num_fields, TRUE, &tok_error_type)) != NULL)) { */
 /*         if (result == NULL) { */
 /*             break; */
 /*         } */
 /*         free(result); */
 /*         --skiprows; */
 /*     } */

 /*     if (skiprows > 0) { */
 /*         /\* There were fewer rows in the file than skiprows. *\/ */
 /*         /\* This is not treated as an error. The result should be an empty array. *\/ */
 /*         *nrows = 0; */
 /*         free(ftypes); */
 /*         del_file_buffer(fb, RESTORE_FINAL); */
 /*         return data_ptr; */
 /*     } */

 /*     /\* XXX Assume *nrows > 0! *\/ */
 /*     /\* */
 /*      *  Read the first row to get the number of fields in the file. */
 /*      *  We'll then use this to pre-validate the values in usecols. */
 /*      *  (It might be easier to do this in the Python wrapper, but that */
 /*      *  would require refactoring the C interface a bit to expose more */
 /*      *  to Python.) */
 /*      *\/ */
 /*     row_count = 0; */
 /*     result = tokenize(fb, word_buffer, WORD_BUFFER_SIZE, */
 /*                               delimiter, quote, comment, &num_fields, TRUE, &tok_error_type); */
 /*     if (result == NULL) { */
 /*         *p_error_type = tok_error_type; */
 /*         *p_error_lineno = 1; */
 /*         free(ftypes); */
 /*         del_file_buffer(fb, RESTORE_FINAL); */
 /*         return NULL; */
 /*     } */

 /*     valid_usecols = (int *) malloc(num_usecols * sizeof(int)); */
 /*     if (valid_usecols == NULL) { */
 /*         /\* Out of memory. *\/ */
 /*         *p_error_type = ERROR_OUT_OF_MEMORY; */
 /*         free(result); */
 /*         free(ftypes); */
 /*         del_file_buffer(fb, RESTORE_FINAL); */
 /*         return NULL; */
 /*     } */

 /*     /\* */
 /*      *  Validate the column indices in usecols, and put the validated */
 /*      *  column indices in valid_usecols. */
 /*      *\/ */
 /*     for (j = 0; j < num_usecols; ++j) { */

 /*         int32_t k; */
 /*         k = usecols[j]; */
 /*         if (k < -num_fields || k >= num_fields) { */
 /*             /\* Invalid column index. *\/ */
 /*             *p_error_type = ERROR_INVALID_COLUMN_INDEX; */
 /*             *p_error_lineno = j;  /\* Abuse 'lineno' and put the bad column index there. *\/ */
 /*             free(valid_usecols); */
 /*             free(result); */
 /*             free(ftypes); */
 /*             del_file_buffer(fb, RESTORE_FINAL); */
 /*             return NULL; */
 /*         } */
 /*         if (k < 0) { */
 /*             k += num_fields; */
 /*         } */
 /*         valid_usecols[j] = k; */
 /*     } */

 /*     current_num_fields = num_fields; */
 /*     row_count = 0; */
 /*     do { */
 /*         int j, k; */

 /*         if (current_num_fields != num_fields) { */
 /*             *p_error_type = ERROR_CHANGED_NUMBER_OF_FIELDS; */
 /*             *p_error_lineno = line_number(fb); */
 /*             break; */
 /*         } */

 /*         for (j = 0; j < num_usecols; ++j) { */

 /*             int error; */
 /*             char typ = ftypes[j].typechar; */
 /*             /\* k is the column index of the field in the file. *\/ */
 /*             k = valid_usecols[j]; */

 /*             /\* XXX Handle error != 0 in the following cases. *\/ */
 /*             if (typ == 'b') { */
 /*                 int8_t x = (int8_t) str_to_int64(result[k], INT8_MIN, INT8_MAX, &error); */
 /*                 *(int8_t *) data_ptr = x; */
 /*                 data_ptr += ftypes[j].size; */
 /*             } */
 /*             else if (typ == 'B') { */
 /*                 uint8_t x = (uint8_t) str_to_uint64(result[k], UINT8_MAX, &error); */
 /*                 *(uint8_t *) data_ptr = x; */
 /*                 data_ptr += ftypes[j].size; */
 /*             } */
 /*             else if (typ == 'h') { */
 /*                 int16_t x = (int16_t) str_to_int64(result[k], INT16_MIN, INT16_MAX, &error); */
 /*                 *(int16_t *) data_ptr = x; */
 /*                 data_ptr += ftypes[j].size; */
 /*             } */
 /*             else if (typ == 'H') { */
 /*                 uint16_t x = (uint16_t) str_to_uint64(result[k], UINT16_MAX, &error); */
 /*                 *(uint16_t *) data_ptr = x; */
 /*                 data_ptr += ftypes[j].size; */
 /*             } */
 /*             else if (typ == 'i') { */
 /*                 int32_t x = (int32_t) str_to_int64(result[k], INT32_MIN, INT32_MAX, &error); */
 /*                 *(int32_t *) data_ptr = x; */
 /*                 data_ptr += ftypes[j].size; */
 /*             } */
 /*             else if (typ == 'I') { */
 /*                 uint32_t x = (uint32_t) str_to_uint64(result[k], UINT32_MAX, &error); */
 /*                 *(uint32_t *) data_ptr = x; */
 /*                 data_ptr += ftypes[j].size; */
 /*             } */
 /*             else if (typ == 'q') { */
 /*                 int64_t x = (int64_t) str_to_int64(result[k], INT64_MIN, INT64_MAX, &error); */
 /*                 *(int64_t *) data_ptr = x; */
 /*                 data_ptr += ftypes[j].size; */
 /*             } */
 /*             else if (typ == 'Q') { */
 /*                 uint64_t x = (uint64_t) str_to_uint64(result[k], UINT64_MAX, &error); */
 /*                 *(uint64_t *) data_ptr = x; */
 /*                 data_ptr += ftypes[j].size; */
 /*             } */
 /*             else if (typ == 'f' || typ == 'd') { */
 /*                 // Convert to float. */
 /*                 double x; */
 /*                 if ((strlen(result[k]) == 0) || !to_double(result[k], &x, sci, decimal)) { */
 /*                     // XXX  Find the canonical platform-independent method to assign nan. */
 /*                     x = 0.0 / 0.0; */
 /*                 } */
 /*                 if (typ == 'f') */
 /*                     *(float *) data_ptr = (float) x; */
 /*                 else */
 /*                     *(double *) data_ptr = x; */
 /*                 data_ptr += ftypes[j].size; */
 /*             } */
 /*             else if (typ == 'c' || typ == 'z') { */
 /*                 // Convert to complex. */
 /*                 double x, y; */
 /*                 if ((strlen(result[k]) == 0) || !to_complex(result[k], &x, &y, sci, decimal)) { */
 /*                     // XXX  Find the canonical platform-independent method to assign nan. */
 /*                     x = 0.0 / 0.0; */
 /*                     y = x; */
 /*                 } */
 /*                 if (typ == 'c') { */
 /*                     *(float *) data_ptr = (float) x; */
 /*                     data_ptr += ftypes[j].size / 2; */
 /*                     *(float *) data_ptr = (float) y; */
 /*                 } */
 /*                 else { */
 /*                     *(double *) data_ptr = x; */
 /*                     data_ptr += ftypes[j].size / 2; */
 /*                     *(double *) data_ptr = y; */
 /*                 } */
 /*                 data_ptr += ftypes[j].size / 2; */
 /*             } */
 /*             else if (typ == 'U') { */
 /*                 // Datetime64, microseconds. */
 /*                 struct tm tm = {0,0,0,0,0,0,0,0,0}; */
 /*                 time_t t; */

 /*                 if (strptime(result[k], datetime_fmt, &tm) == NULL) { */
 /*                     memset(data_ptr, 0, 8); */
 /*                 } */
 /*                 else { */
 /*                     tm.tm_isdst = -1; */
 /*                     t = mktime(&tm); */
 /*                     if (t == -1) { */
 /*                         memset(data_ptr, 0, 8); */
 /*                     } */
 /*                     else { */
 /*                         *(uint64_t *) data_ptr = (long long) (t - tz_offset) * 1000000L; */
 /*                     } */
 /*                 } */
 /*                 data_ptr += 8; */
 /*             } */
 /*             else { */
 /*                 // String */
 /*                 strncpy(data_ptr, result[k], ftypes[j].size); */
 /*                 data_ptr += ftypes[j].size; */
 /*             } */
 /*         } */
 /*         free(result); */
 /*         ++row_count; */
 /*     } while ((row_count < *nrows) && (result = tokenize(fb, word_buffer, WORD_BUFFER_SIZE, */
 /*                               delimiter, quote, comment, &current_num_fields, TRUE, &tok_error_type)) != NULL); */

 /*     del_file_buffer(fb, RESTORE_FINAL); */

 /*     *nrows = row_count; */

 /*     free(valid_usecols); */

 /*     return (void *) data_ptr; */
 /* } */


int merge_chunks(parser_t *parser) {
    int i, j, ncols;

    // Get a consensus on number of columns and check types
    for (i = 0; i < parser->nchunks; ++i)
    {
        if (i == 0) {
            ncols = parser->chunks[i].ncols;
        } else {
            // XXX this should not happen
            if (ncols != parser->chunks[i].ncols) {
                return -1;
            }
        }
    }

    for (i = 0; i < parser->nchunks; ++i)
    {
        for (j = 0; j < ncols; ++j)
        {
            return -1;
        }
    }

    return 0;
}

int _try_int64(parser_t *parser, array_t *arr, char** strings, size_t length) {
    int i, error;
    int64_t *data;

    arr->data = malloc(length * sizeof(int64_t));
    data = (int64_t*) arr->data;

    for (i = 0; i < length; ++i)
    {
        *data++ = (int64_t) str_to_int64(strings[i], INT64_MIN,
                                         INT64_MAX, &error,
                                         parser->thousands);

        if (error != 0) {
            return -1;
        }
    }

    return 0;
}

int _try_float(parser_t *parser, array_t* arr, char** strings, size_t length) {
    int i, error;
    double *data;

    arr->data = malloc(length * sizeof(double));
    data = (double*) arr->data;

    for (i = 0; i < length; ++i)
    {
        error = to_double(strings[i], data, parser->sci, parser->decimal);

        if (error != 1) {
            return -1;
        }
    }

    return 0;
}


int _try_boolean(parser_t *parser, array_t* arr, char** strings, size_t length) {
    int i, error;
    uint8_t *data;

    arr->data = malloc(length * sizeof(uint8_t));
    data = (uint8_t*) arr->data;

    for (i = 0; i < length; ++i)
    {
        error = to_boolean(strings[i], data);

        if (error != 1) {
            return -1;
        }
    }

    return 0;
}

typedef int (*cast_func)(parser_t *parser, array_t* arr,
                         char** strings, size_t length);

static cast_func _inference_order[3] = {_try_int64, _try_float, _try_boolean};

int convert_infer(parser_t *parser, array_t* result,
                  char** strings, size_t length) {

    int i, status;
    /* array_t* result = (array_t*) malloc(sizeof(array_t*)); */

    for (i = 0; i < sizeof(_inference_order); ++i)
    {
        status = _inference_order[i](parser, result, strings, length);

        if (status == 0) {
            // success
            return 0;
        }
    }

    free(result);

    return 0;
}


void parser_set_default_options(parser_t *self) {
    // parsing, type inference
    self->infer_types = 1;
    self->decimal = '.';
    self->sci = 'E';

    // For tokenization
    self->state = START_RECORD;

    self->delimiter = ','; // XXX
    self->delim_whitespace = 0;

    self->doublequote = 0;
    self->quotechar = '"';
    self->escapechar = 0;
    self->skipinitialspace = 0;
    self->quoting = QUOTE_MINIMAL;
    self->allow_embedded_newline = 1;
    self->strict = 0;

    self->error_bad_lines = 0;
    self->warn_bad_lines = 0;

    self->commentchar = '#';
    self->thousands = '\0';

    self->skipset = NULL;
    self->skip_footer = 0;
}

int get_parser_memory_footprint(parser_t *self) {
    return 0;
}

parser_t* parser_new() {
    return (parser_t*) calloc(1, sizeof(parser_t));
}

// XXX handle on systems without the capability

#include <sys/stat.h>
#include <sys/mman.h>

typedef struct _memory_map {

    FILE *file;

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


/*
 *  void *new_file_buffer(FILE *f, int buffer_size)
 *
 *  Allocate a new file_buffer.
 *  Returns NULL if the memory allocation fails or if the call to mmap fails.
 *
 *  buffer_size is ignored.
 */

void *new_mmap(FILE *f)
{
    struct stat buf;
    int fd;
    memory_map *mm;
    off_t position;
    off_t filesize;

    fd = fileno(f);
    if (fstat(fd, &buf) == -1) {
        fprintf(stderr, "new_file_buffer: fstat() failed. errno =%d\n", errno);
        return NULL;
    }
    filesize = buf.st_size;  /* XXX This might be 32 bits. */

    mm = (memory_map *) malloc(sizeof(memory_map));
    if (mm == NULL) {
        /* XXX Eventually remove this print statement. */
        fprintf(stderr, "new_file_buffer: malloc() failed.\n");
        return NULL;
    }
    mm->file = f;
    mm->size = (off_t) filesize;
    mm->line_number = 0;

    mm->fileno = fd;
    mm->position = ftell(f);
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


void del_mmap(void *src)
{
    munmap(MM(src)->memmap, MM(src)->size);

    /*
     *  With a memory mapped file, there is no need to do
     *  anything if restore == RESTORE_INITIAL.
     */
    /* if (restore == RESTORE_FINAL) { */
    /*     fseek(FB(fb)->file, FB(fb)->current_pos, SEEK_SET); */
    /* } */
    free(src);
}

int _buffer_mmap_bytes(parser_t *self, size_t nbytes) {
    memory_map *src = MM(self->source);

    if (src->position == src->last_pos) {
        self->datalen = 0;
        return REACHED_EOF;
    }

    self->data = src->memmap + src->position;

    if (src->position + nbytes > src->last_pos) {
        // fewer than nbytes remaining
        self->datalen = src->last_pos - src->position;
    } else {
        self->datalen = nbytes;
    }

    src->position += self->datalen;
    return 0;
}

int parser_mmap_init(parser_t *self, FILE* fp) {
    self->sourcetype = 'M';
    self->source = new_mmap(fp);

    // TODO: better error message
    if (NULL == self->source)
        return -1;

    return 0;
}

int parser_gzip_source_init(parser_t *self, FILE* fp) {
    return 0;
}

int parser_array_source_init(parser_t *self, char *bytes, size_t length) {
    self->sourcetype = 'A';
    self->source = new_array_source(bytes, length);
    return 0;
}

int parser_init(parser_t *self) {
    int sz;

    /*
      Initialize data buffers
    */

    self->stream = NULL;
    self->words = NULL;
    self->word_starts = NULL;
    self->line_start = NULL;
    self->line_fields = NULL;

    // token stream
    self->stream = (char*) malloc(STREAM_INIT_SIZE * sizeof(char));
    if (self->stream == NULL) {
        return PARSER_OUT_OF_MEMORY;
    }
    self->stream_cap = STREAM_INIT_SIZE;
    self->stream_len = 0;

    // word pointers and metadata
    sz = STREAM_INIT_SIZE / 10;
    sz = sz? sz : 1;
    self->words = (char**) malloc(sz * sizeof(char*));
    self->word_starts = (int*) malloc(sz * sizeof(int));
    self->words_cap = sz;
    self->words_len = 0;

    // line pointers and metadata
    self->line_start = (int*) malloc(sz * sizeof(int));

    self->line_fields = (int*) malloc(sz * sizeof(int));

    self->lines_cap = sz;
    self->lines = 0;
    self->file_lines = 0;

    if (self->stream == NULL || self->words == NULL ||
        self->word_starts == NULL || self->line_start == NULL ||
        self->line_fields == NULL) {

        parser_cleanup(self);

        return PARSER_OUT_OF_MEMORY;
    }

    /* amount of bytes buffered */
    self->datalen = 0;
    self->datapos = 0;

    self->line_start[0] = 0;
    self->line_fields[0] = 0;

    self->pword_start = self->stream;
    self->word_start = 0;

    self->state = START_RECORD;

    self->error_msg = NULL;

    return 0;
}


int make_stream_space(parser_t *self, size_t nbytes) {
    int i, status, cap;
    void *orig_ptr;

    // Can we fit potentially nbytes tokens (+ null terminators) in the stream?

    /* TRACE(("maybe growing buffers\n")); */

    /*
      TOKEN STREAM
    */

    orig_ptr = (void *) self->stream;
    self->stream = (char*) grow_buffer((void *) self->stream,
                                       self->stream_len,
                                       &self->stream_cap, nbytes * 2,
                                       sizeof(char), &status);

    if (status != 0) {
        return PARSER_OUT_OF_MEMORY;
    }

    // realloc sets errno when moving buffer?
    if (self->stream != orig_ptr) {
        // uff
        /* TRACE(("Moving word pointers\n")) */

        self->pword_start = self->stream + self->word_start;

        for (i = 0; i < self->words_len; ++i)
        {
            self->words[i] = self->stream + self->word_starts[i];
        }
    }


    /*
      WORD VECTORS
    */

    cap = self->words_cap;
    self->words = (char**) grow_buffer((void *) self->words,
                                       self->words_len,
                                       &self->words_cap, nbytes,
                                       sizeof(char*), &status);
    if (status != 0) {
        return PARSER_OUT_OF_MEMORY;
    }


    // realloc took place
    if (cap != self->words_cap) {
        self->word_starts = (int*) safe_realloc((void *) self->word_starts,
                                                sizeof(int) * self->words_cap);
        if (self->word_starts == NULL) {
            return PARSER_OUT_OF_MEMORY;
        }
    }


    /*
      LINE VECTORS
    */

    cap = self->lines_cap;
    self->line_start = (int*) grow_buffer((void *) self->line_start,
                                          self->lines,
                                          &self->lines_cap, nbytes,
                                          sizeof(int), &status);
    if (status != 0) {
        return PARSER_OUT_OF_MEMORY;
    }

    // realloc took place
    if (cap != self->lines_cap) {
        self->line_fields = (int*) safe_realloc((void *) self->line_fields,
                                                sizeof(int) * self->lines_cap);

        if (self->line_fields == NULL) {
            return PARSER_OUT_OF_MEMORY;
        }
    }


    /* TRACE(("finished growing buffers\n")); */

    return 0;
}


int inline push_char(parser_t *self, char c) {
    /* TRACE(("pushing %c \n", c)) */
    self->stream[self->stream_len++] = c;
    return 0;
}

int inline end_field(parser_t *self) {
    int pos;

    // XXX cruft
    self->numeric_field = 0;

    // null terminate token
    push_char(self, '\0');

    // set pointer and metadata
    self->words[self->words_len] = self->pword_start;

    TRACE(("Saw word %s at: %d\n", self->pword_start, self->word_start))

    self->word_starts[self->words_len] = self->word_start;
    self->words_len++;

    // increment line field count
    self->line_fields[self->lines]++;

    // New field begin in stream
    self->pword_start = self->stream + self->stream_len;
    self->word_start = self->stream_len;

    return 0;
}

int inline end_line(parser_t *self) {
    int fields;
    khiter_t k;  /* for hash set detection */
    int ex_fields = -1;

    fields = self->line_fields[self->lines];

    if (self->lines > 0) {
        ex_fields = self->line_fields[self->lines - 1];
    }

    if (self->skipset != NULL) {
        k = kh_get_int64((kh_int64_t*) self->skipset, self->file_lines);

        if (k != ((kh_int64_t*)self->skipset)->n_buckets) {
            TRACE(("Skipping row %d\n", self->file_lines));
            // increment file line count
            self->file_lines++;

            // skip the tokens from this bad line
            self->line_start[self->lines] += fields;

            // reset field count
            self->line_fields[self->lines] = 0;
            return 0;
        }
    }

    if (!(self->lines <= self->header + 1) && fields != ex_fields) {
        // increment file line count
        self->file_lines++;

        // skip the tokens from this bad line
        self->line_start[self->lines] += fields;

        // reset field count
        self->line_fields[self->lines] = 0;

        // file_lines is now the _actual_ file line number (starting at 1)

        if (self->error_bad_lines) {
            self->error_msg = (char*) malloc(100);
            sprintf(self->error_msg, "Expected %d fields in line %d, saw %d\n",
                    ex_fields, self->file_lines, fields);

            TRACE(("Error at line %d, %d fields\n", self->file_lines, fields));

            return -1;
        } else {
            // simply skip bad lines
            if (self->warn_bad_lines) {
                // print error message
                printf("Skipping line %d: expected %d fields, saw %d\n",
                       self->file_lines, ex_fields, fields);
            }
        }
    } else {
        // increment both line counts
        self->file_lines++;

        self->lines++;

        // good line, set new start point
        self->line_start[self->lines] = (self->line_start[self->lines - 1] +
                                         fields);

        TRACE(("new line start: %d\n", self->line_start[self->lines]));

        // new line start with 0 fields
        self->line_fields[self->lines] = 0;
    }

    TRACE(("Finished line, at %d\n", self->lines));

    return 0;
}

int parser_clear_data_buffers(parser_t *self) {
    free_if_not_null(self->stream);
    free_if_not_null(self->words);
    free_if_not_null(self->word_starts);
    free_if_not_null(self->line_start);
    free_if_not_null(self->line_fields);

    return 0;
}

void parser_free(parser_t *self) {
    // opposite of parser_init
    parser_cleanup(self);
    free(self);
}



int parser_cleanup(parser_t *self) {
    if (parser_cleanup_filebuffers(self) < 0) {
        return -1;
    }

    if (parser_clear_data_buffers(self) < 0) {
        return -1;
    }

    // XXX where to put this
    free_if_not_null(self->error_msg);

    if (self->skipset != NULL)
        kh_destroy_int64((kh_int64_t*) self->skipset);

    return 0;
}

int parser_add_skiprow(parser_t *self, int64_t row) {
    khiter_t k;
    kh_int64_t *set;
    int ret = 0;

    if (self->skipset == NULL) {
        self->skipset = (void*) kh_init_int64();
    }

    set = (kh_int64_t*) self->skipset;

    k = kh_put_int64(set, row, &ret);
    set->keys[k] = row;

    return 0;
}

int parser_buffer_bytes(parser_t *self, size_t nbytes) {
    int status;
    size_t bytes;
    void *src = self->source;

    // This should probably end up as a method table

    status = 0;

    self->datapos = 0;

    switch(self->sourcetype) {
        case 'F': // basic FILE*

            bytes = fread((void *) self->data, sizeof(char), nbytes,
                          FS(src)->fp);
            self->datalen = bytes;

            TRACE(("Read %d bytes\n", (int) bytes));

            // printf("%s\n", self->data);

            if (bytes == 0) {
                status = REACHED_EOF;
            }
            break;

        case 'A': // in-memory bytes (e.g. from StringIO)
            // ew, side effects
            status = _buffer_array_bytes(self, nbytes);
            break;

#ifdef HAVE_MEMMAP
        case 'M': // memory map
            status = _buffer_mmap_bytes(self, nbytes);

            break;
#endif

#ifdef HAVE_GZIP
        case 'G': // gzip'd file

            break;
#endif

    }

    return status;
}

int _buffer_array_bytes(parser_t *self, size_t nbytes) {
    array_source *src = ARS(self->source);

    if (src->position == src->length) {
        self->datalen = 0;
        return REACHED_EOF;
    }

    self->data = src->data + src->position;

    if (src->position + nbytes > src->length) {
        // fewer than nbytes remaining
        self->datalen = src->length - src->position;
    } else {
        self->datalen = nbytes;
    }

    src->position += self->datalen;

    TRACE(("datalen: %d\n", self->datalen));

    TRACE(("pos: %d, length: %d", (int) src->position, (int) src->length));
    return 0;
}


int parser_cleanup_filebuffers(parser_t *self) {
    switch(self->sourcetype) {

        case 'F':
            free(self->data);
            del_file_source(self->source);
            break;

        case 'A': // in-memory bytes (e.g. from StringIO)
            del_array_source(self->source);
            break;

#ifdef HAVE_MEMMAP
        case 'M': // memory map
            del_mmap(self->source);
            break;
#endif


#ifdef HAVE_GZIP
        case 'G': // gzip'd file

            break;
#endif

    }

    return 0;
}

/*

  Tokenization macros and state machine code

*/

//    printf("pushing %c\n", c);                \

#define PUSH_CHAR(c)                           \
    *stream++ = c;                             \
    slen++;

// This is a little bit of a hack but works for now

#define END_FIELD()                                \
    self->stream_len = slen;                   \
    if (end_field(self) < 0) {                 \
        goto parsingerror;                     \
    }                                          \
    stream = self->stream + self->stream_len;  \
    slen = self->stream_len;

#define END_LINE()                                                      \
    self->stream_len = slen;                                            \
    if (end_line(self) < 0) {                                           \
        goto parsingerror;                                              \
    }                                                                   \
    self->state = START_RECORD;                                         \
    if (line_limit > 0 && self->lines == start_lines + line_limit) {    \
        goto linelimit;                                                 \
                                                                        \
    }                                                                   \
    stream = self->stream + self->stream_len;                           \
    slen = self->stream_len;

#define IS_WHITESPACE(c) ((c == ' ' || c == '\t'))

typedef int (*parser_op)(parser_t *self, size_t line_limit);

#define _TOKEN_CLEANUP()                                                \
    self->stream_len = slen;                                            \
    self->datapos = i;                                                  \
    TRACE(("datapos: %d, datalen: %d\n", self->datapos, self->datalen));



int tokenize_delimited(parser_t *self, size_t line_limit)
{
    int i, slen, start_lines;
    char c;
    char *stream;
    char *buf = self->data + self->datapos;

    start_lines = self->lines;

    if (make_stream_space(self, self->datalen - self->datapos) < 0) {
        self->error_msg = "out of memory";
        return -1;
    }

    stream = self->stream + self->stream_len;
    slen = self->stream_len;

    TRACE(("%s\n", buf));

    for (i = self->datapos; i < self->datalen; ++i)
    {
        // Next character in file
        c = *buf++;

        TRACE(("Iter: %d Char: %c Line %d field_count %d, state %d\n",
               i, c, self->file_lines + 1, self->line_fields[self->lines],
               self->state));

        switch(self->state) {
        case START_RECORD:
            // start of record

            if (c == '\n') {
                // \n\r possible?
                END_LINE();
                break;
            } else if (c == '\r') {
                self->state = EAT_CRNL;
                break;
            }

            /* normal character - handle as START_FIELD */
            self->state = START_FIELD;
            /* fallthru */
        case START_FIELD:
            /* expecting field */
            if (c == '\n') {
                END_FIELD();
                END_LINE();
                /* self->state = START_RECORD; */
            } else if (c == '\r') {
                END_FIELD();
                self->state = EAT_CRNL;
            }
            else if (c == self->quotechar &&
                     self->quoting != QUOTE_NONE) {
                /* start quoted field */
                self->state = IN_QUOTED_FIELD;
            }
            else if (c == self->escapechar) {
                /* possible escaped character */
                self->state = ESCAPED_CHAR;
            }
            else if (c == ' ' && self->skipinitialspace)
                /* ignore space at start of field */
                ;
            else if (c == self->delimiter) {
                /* save empty field */
                END_FIELD();
            }
            else {
                /* begin new unquoted field */
                if (self->quoting == QUOTE_NONNUMERIC)
                    self->numeric_field = 1;

                // TRACE(("pushing %c", c));
                PUSH_CHAR(c);
                self->state = IN_FIELD;
            }
            break;

        case ESCAPED_CHAR:
            /* if (c == '\0') */
            /*  c = '\n'; */

            PUSH_CHAR(c);
            self->state = IN_FIELD;
            break;

        case IN_FIELD:
            /* in unquoted field */
            if (c == '\n') {
                END_FIELD();
                END_LINE();
                /* self->state = START_RECORD; */
            } else if (c == '\r') {
                END_FIELD();
                self->state = EAT_CRNL;
            }
            else if (c == self->escapechar) {
                /* possible escaped character */
                self->state = ESCAPED_CHAR;
            }
            else if (c == self->delimiter) {
                // End of field. End of line not reached yet
                END_FIELD();
                self->state = START_FIELD;
            }
            else {
                /* normal character - save in field */
                PUSH_CHAR(c);
            }
            break;

        case IN_QUOTED_FIELD:
            /* in quoted field */
            if (c == self->escapechar) {
                /* Possible escape character */
                self->state = ESCAPE_IN_QUOTED_FIELD;
            }
            else if (c == self->quotechar &&
                     self->quoting != QUOTE_NONE) {
                if (self->doublequote) {
                    /* doublequote; " represented by "" */
                    self->state = QUOTE_IN_QUOTED_FIELD;
                }
                else {
                    /* end of quote part of field */
                    self->state = IN_FIELD;
                }
            }
            else {
                /* normal character - save in field */
                PUSH_CHAR(c);
            }
            break;

        case ESCAPE_IN_QUOTED_FIELD:
            /* if (c == '\0') */
            /*  c = '\n'; */

            PUSH_CHAR(c);
            self->state = IN_QUOTED_FIELD;
            break;

        case QUOTE_IN_QUOTED_FIELD:
            /* doublequote - seen a quote in an quoted field */
            if (self->quoting != QUOTE_NONE && c == self->quotechar) {
                /* save "" as " */

                PUSH_CHAR(c);
                self->state = IN_QUOTED_FIELD;
            }
            else if (c == self->delimiter) {
                // End of field. End of line not reached yet

                END_FIELD();
                self->state = START_FIELD;
            }
            else if (c == '\n') {
                END_FIELD();
                END_LINE();
                /* self->state = START_RECORD; */
            }
            else if (c == '\r') {
                END_FIELD();
                self->state = EAT_CRNL;
            }
            else if (!self->strict) {
                PUSH_CHAR(c);
                self->state = IN_FIELD;
            }
            else {
                self->error_msg = (char*) malloc(50);
                sprintf(self->error_msg, "'%c' expected after '%c'",
                        self->delimiter, self->quotechar);
                goto parsingerror;
            }
            break;

        case EAT_CRNL:
            if (c == '\n') {
                END_LINE();
                /* self->state = START_RECORD; */
            } else {
                /* self->error_msg = ("new-line character seen in" */
                /*                 " unquoted field - do you need" */
                /*                 " to open the file in " */
                /*                 "universal-newline mode?"); */
                goto parsingerror;
            }
            break;

        }

    }

    _TOKEN_CLEANUP();

    TRACE(("Finished tokenizing input\n"))

    return 0;

parsingerror:
    i++;
    _TOKEN_CLEANUP();

    return -1;

linelimit:
    i++;
    _TOKEN_CLEANUP();

    return 0;
}

int tokenize_whitespace(parser_t *self, size_t line_limit)
{
    int i, slen, start_lines;
    char c;
    char *stream;
    char *buf = self->data + self->datapos;

    start_lines = self->lines;

    if (make_stream_space(self, self->datalen - self->datapos) < 0) {
        self->error_msg = "out of memory";
        return -1;
    }

    stream = self->stream + self->stream_len;
    slen = self->stream_len;

    TRACE(("%s\n", buf));

    for (i = self->datapos; i < self->datalen; ++i)
    {
        // Next character in file
        c = *buf++;

        TRACE(("Iter: %d Char: %c Line %d field_count %d\n",
               i, c, self->file_lines + 1, self->line_fields[self->lines]));

        switch(self->state) {

        case EAT_WHITESPACE:
            if (!IS_WHITESPACE(c)) {
                // END_FIELD();
                self->state = START_FIELD;
                // Fall through to subsequent state
            } else {
                // if whitespace char, keep slurping
                break;
            }

        case START_RECORD:
            // start of record
            if (c == '\n') {
                // \n\r possible?
                END_LINE();
                break;
            } else if (c == '\r') {
                self->state = EAT_CRNL;
                break;
            }

            /* normal character - handle as START_FIELD */
            self->state = START_FIELD;
            /* fallthru */
        case START_FIELD:
            /* expecting field */
            if (c == '\n') {
                END_FIELD();
                END_LINE();
                /* self->state = START_RECORD; */
            } else if (c == '\r') {
                END_FIELD();
                self->state = EAT_CRNL;
            }
            else if (c == self->quotechar &&
                     self->quoting != QUOTE_NONE) {
                /* start quoted field */
                self->state = IN_QUOTED_FIELD;
            }
            else if (c == self->escapechar) {
                /* possible escaped character */
                self->state = ESCAPED_CHAR;
            }
            /* else if (c == ' ' && self->skipinitialspace) */
            /*     /\* ignore space at start of field *\/ */
            /*     ; */
            else if (IS_WHITESPACE(c)) {
                self->state = EAT_WHITESPACE;
            }
            else {
                /* begin new unquoted field */
                if (self->quoting == QUOTE_NONNUMERIC)
                    self->numeric_field = 1;

                // TRACE(("pushing %c", c));
                PUSH_CHAR(c);
                self->state = IN_FIELD;
            }
            break;

        case ESCAPED_CHAR:
            /* if (c == '\0') */
            /*  c = '\n'; */

            PUSH_CHAR(c);
            self->state = IN_FIELD;
            break;

        case IN_FIELD:
            /* in unquoted field */
            if (c == '\n') {
                END_FIELD();
                END_LINE();
                /* self->state = START_RECORD; */
            } else if (c == '\r') {
                END_FIELD();
                self->state = EAT_CRNL;
            }
            else if (c == self->escapechar) {
                /* possible escaped character */
                self->state = ESCAPED_CHAR;
            }
            else if (IS_WHITESPACE(c)) {
                // End of field. End of line not reached yet
                END_FIELD();
                self->state = EAT_WHITESPACE;
            }
            else {
                /* normal character - save in field */
                PUSH_CHAR(c);
            }
            break;

        case IN_QUOTED_FIELD:
            /* in quoted field */
            if (c == self->escapechar) {
                /* Possible escape character */
                self->state = ESCAPE_IN_QUOTED_FIELD;
            }
            else if (c == self->quotechar &&
                     self->quoting != QUOTE_NONE) {
                if (self->doublequote) {
                    /* doublequote; " represented by "" */
                    self->state = QUOTE_IN_QUOTED_FIELD;
                }
                else {
                    /* end of quote part of field */
                    self->state = IN_FIELD;
                }
            }
            else {
                /* normal character - save in field */
                PUSH_CHAR(c);
            }
            break;

        case ESCAPE_IN_QUOTED_FIELD:
            /* if (c == '\0') */
            /*  c = '\n'; */

            PUSH_CHAR(c);
            self->state = IN_QUOTED_FIELD;
            break;

        case QUOTE_IN_QUOTED_FIELD:
            /* doublequote - seen a quote in an quoted field */
            if (self->quoting != QUOTE_NONE && c == self->quotechar) {
                /* save "" as " */

                PUSH_CHAR(c);
                self->state = IN_QUOTED_FIELD;
            }
            else if (IS_WHITESPACE(c)) {
                // End of field. End of line not reached yet

                END_FIELD();
                self->state = EAT_WHITESPACE;
            }
            else if (c == '\n') {
                END_FIELD();
                END_LINE();
                /* self->state = START_RECORD; */
            }
            else if (c == '\r') {
                END_FIELD();
                self->state = EAT_CRNL;
            }
            else if (!self->strict) {
                PUSH_CHAR(c);
                self->state = IN_FIELD;
            }
            else {
                self->error_msg = (char*) malloc(50);
                sprintf(self->error_msg, "'%c' expected after '%c'",
                        self->delimiter, self->quotechar);
                goto parsingerror;
            }
            break;

        case EAT_CRNL:
            if (c == '\n') {
                END_LINE();
                /* self->state = START_RECORD; */
            } else {
                /* self->error_msg = ("new-line character seen in" */
                /*                 " unquoted field - do you need" */
                /*                 " to open the file in " */
                /*                 "universal-newline mode?"); */
                goto parsingerror;
            }
            break;

        }

    }

    _TOKEN_CLEANUP();

    TRACE(("Finished tokenizing input\n"))

    return 0;

parsingerror:
    i++;
    _TOKEN_CLEANUP();

    return -1;

linelimit:
    i++;
    _TOKEN_CLEANUP();

    return 0;
}


int parser_handle_eof(parser_t *self) {
    TRACE(("handling eof, datalen: %d\n", self->datalen))
    if (self->datalen == 0 && (self->state != START_RECORD)) {
        // test cases needed here
        // TODO: empty field at end of line
        TRACE(("handling eof\n"));

        if (self->state == IN_FIELD) {
            if (end_field(self) < 0)
                return -1;
        } else if (self->state == QUOTE_IN_QUOTED_FIELD) {
            if (end_field(self) < 0)
                return -1;
        } else if (self->state == IN_QUOTED_FIELD) {
            self->error_msg = (char*) malloc(100);
            sprintf(self->error_msg, "EOF inside string starting at line %d",
                    self->file_lines);
            return -1;
        }

        if (end_line(self) < 0)
            return -1;

        return 0;
    }
}

int parser_consume_rows(parser_t *self, size_t nrows) {
    int i, offset, word_deletions, char_count;

    if (nrows > self->lines) {
        nrows = self->lines;
    }

    /* do nothing */
    if (nrows == 0)
        return 0;

    /* cannot guarantee that nrows + 1 has been observed */
    word_deletions = self->line_start[nrows - 1] + self->line_fields[nrows - 1];
    char_count = (self->word_starts[word_deletions - 1] +
                  strlen(self->words[word_deletions - 1]) + 1);

    TRACE(("Deleting %d words, %d chars\n", word_deletions, char_count));

    /* move stream, only if something to move */
    if (char_count < self->stream_len) {
        memmove((void*) self->stream, (void*) (self->stream + char_count),
                self->stream_len - char_count);
    }
    /* buffer counts */
    self->stream_len -= char_count;

    /* move token metadata */
    for (i = 0; i < self->words_len - word_deletions; ++i) {
        offset = i + word_deletions;

        self->words[i] = self->words[offset] - char_count;
        self->word_starts[i] = self->word_starts[offset] - char_count;
    }
    self->words_len -= word_deletions;

    /* move current word pointer to stream */
    self->pword_start -= char_count;
    self->word_start -= char_count;

    /* printf("Line_start: "); */
    /* for (i = 0; i < self->lines; ++i) { */
    /*     printf("%d ", self->line_start[i]); */
    /* } */
    /* printf("\n"); */

    /* move line metadata */
    for (i = 0; i < self->lines - nrows + 1; ++i)
    {
        offset = i + nrows;
        self->line_start[i] = self->line_start[offset] - word_deletions;

        /* TRACE(("First word in line %d is now %s\n", i, */
        /*        self->words[self->line_start[i]])); */

        self->line_fields[i] = self->line_fields[offset];
    }
    self->lines -= nrows;
    /* self->line_fields[self->lines] = 0; */

    /* printf("Line_start: "); */
    /* for (i = 0; i < self->lines; ++i) { */
    /*     printf("%d ", self->line_start[i]); */
    /* } */
    /* printf("\n"); */

    return 0;
}

static size_t _next_pow2(size_t sz) {
    size_t result = 1;
    while (result < sz) result *= 2;
    return result;
}

int parser_trim_buffers(parser_t *self) {
    /*
      Free memory
     */
    size_t new_cap;

    /* trim stream */
    new_cap = _next_pow2(self->stream_len);
    if (new_cap < self->stream_cap) {
        self->stream = safe_realloc((void*) self->stream, new_cap);
        self->stream_cap = new_cap;
    }

    /* trim words, word_starts */
    new_cap = _next_pow2(self->words_len);
    if (new_cap < self->words_cap) {
        self->words = (char**) safe_realloc((void*) self->words,
                                            new_cap * sizeof(char*));
        self->word_starts = (int*) safe_realloc((void*) self->word_starts,
                                                new_cap * sizeof(int));
        self->words_cap = new_cap;
    }

    /* trim line_start, line_fields */
    new_cap = _next_pow2(self->lines);
    if (new_cap < self->lines_cap) {
        self->line_start = (int*) safe_realloc((void*) self->line_start,
                                               new_cap * sizeof(int));
        self->line_fields = (int*) safe_realloc((void*) self->line_fields,
                                                new_cap * sizeof(int));
        self->lines_cap = new_cap;
    }

    return 0;
}

void debug_print_parser(parser_t *self) {
    int i, j, line;
    char *token;

    for (line = 0; line < self->lines; ++line)
    {
        printf("(Parsed) Line %d: ", line);

        for (j = 0; j < self->line_fields[j]; ++j)
        {
            token = self->words[j + self->line_start[line]];
            printf("%s ", token);
        }
        printf("\n");
    }
}

int clear_parsed_lines(parser_t *self, size_t nlines) {
    // TODO. move data up in stream, shift relevant word pointers

    return 0;
}


/*
  nrows : number of rows to tokenize (or until reach EOF)
  all : tokenize all the data vs. certain number of rows
 */
int _tokenize_helper(parser_t *self, size_t nrows, int all) {
    parser_op tokenize_bytes;

    int status = 0;
    int start_lines = self->lines;

    if (self->delim_whitespace) {
        tokenize_bytes = tokenize_whitespace;
    } else {
        tokenize_bytes = tokenize_delimited;
    }

    if (self->state == FINISHED) {
        return 0;
    }

    TRACE(("Asked to tokenize %d rows\n", (int) nrows));

    while (1) {
        if (!all && self->lines - start_lines >= nrows)
            break;

        if (self->datapos == self->datalen) {
            status = parser_buffer_bytes(self, self->chunksize);

            if (status == REACHED_EOF) {
                // close out last line
                status = parser_handle_eof(self);
                self->state = FINISHED;
                break;
            }
        }

        TRACE(("Trying to process %d bytes\n", self->datalen - self->datapos));
        TRACE(("sourcetype: %c, status: %d\n", self->sourcetype, status));

        status = tokenize_bytes(self, nrows);

        if (status < 0) {
            // XXX
            TRACE(("Status %d returned from tokenize_bytes, breaking\n",
                   status));
            status = -1;
            break;
        }
    }
    TRACE(("leaving tokenize_helper\n"));
    return status;
}

int tokenize_nrows(parser_t *self, size_t nrows) {
    int status = _tokenize_helper(self, nrows, 0);
    return status;
}

int tokenize_all_rows(parser_t *self) {
    int status = _tokenize_helper(self, -1, 1);
    return status;
}

/*
  Iteration through ragged matrix structure
*/

int test_tokenize(char *fname) {
    parser_t *self;
    coliter_t citer;
    char *error_msg;
    int status = 0;
    int nbytes = CHUNKSIZE;

    clock_t start = clock();

    self = parser_new();
    self->chunksize = nbytes;

    // self->source = malloc(sizeof(file_source));

    FILE* fp = fopen(fname, "rb");
    parser_file_source_init(self, fp);

    parser_set_default_options(self);

    if (parser_init(self) < 0) {
        return -1;
    }

    self->header = 0;

    status = tokenize_all_rows(self);

    if (status != 0) {
        if (self->error_msg == NULL) {
            printf("PARSE_ERROR: no message\n");
        }
        else {
            printf("PARSE_ERROR: %s", self->error_msg);
        }
    }


    // debug_print_parser(parser);

    // return 0;
    /* if (status < 0) { */
    /*  return status; */
    /* } */

    /* int i, words = 0; */
    /* for (i = 0; i < parser.stream_len; ++i) */
    /* { */
    /*     if (parser.stream[i] == '\0') words++; */
    /* } */

    printf("Time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
    /* return 0; */

    int i, j, error, columns;
    char *word;
    double data;

    columns = self->line_fields[0];

    printf("columns: %d\n", columns);
    printf("errno is %d\n", errno);

    for (j = 0; j < columns; ++j)
    {
        coliter_setup(&citer, self, j, 0);

        for (i = 0; i < self->lines; ++i)
        {
            if (j >= self->line_fields[i]) continue;

            word = COLITER_NEXT(citer);
            error = to_double(word, &data, self->sci, self->decimal);
            if (error != 1) {
                printf("error at %d, errno: %d\n", i, errno);
                printf("failed: %s %d\n", word, (int) strlen(word));
                break;
            }
        }
    }

    /* for (j = 0; j < columns; ++j) */
    /* { */
    /*  // citer = coliter_new(&parser, j); */
    /*  for (i = 0; i < parser->lines; ++i) */
    /*  { */
    /*      if (j >= parser->line_fields[i]) continue; */
    /*      word = parser->words[parser->line_start[i] + j]; */
    /*      error = to_double(word, &data, parser->sci, parser->decimal); */
    /*      if (error != 1) { */
    /*          printf("error at %d, errno: %d\n", i, errno); */
    /*          printf("failed: %s %d\n", word, (int) strlen(word)); */
    /*          break; */
    /*      } */
    /*  } */
    /*  // free(citer); */
    /* } */

    printf("Time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);

    /* for (i = 0; i < parser.words_len; ++i) */
    /* { */
    /*     error = to_double(parser.words[i], &data, parser.sci, parser.decimal); */
    /*  if (error != 1) { */
    /*      ; */
    /*      printf("failed: %s %d\n", parser.words[i], */
    /*             (int) strlen(parser.words[i])); */
    /*  } else { */
    /*      ; */
    /*      /\* printf("success: %.4f\n", data); *\/ */
    /*  } */
    /* } */


    /* for (i = 0; i < parser.words_len; ++i) */
    /* { */
    /*     error = to_double(parser.words[i], &data, parser.sci, parser.decimal); */
    /*  if (error != 1) { */
    /*      ; */
            /* printf("failed: %s %d\n", parser.words[i], */
            /*     (int) strlen(parser.words[i])); */
    /*  } else { */
    /*      ; */
    /*      /\* printf("success: %.4f\n", data); *\/ */
    /*  } */
    /* } */


    /* printf("saw %d words\n", words); */

    // debug_print_parser(&parser);

    parser_free(self);

    fclose(fp);

    /* if (parser_cleanup(self) < 0) { */
    /*     return -1; */
    /* } */

    /* free(parser); */

    return status;
}


void test_count_lines(char *fname) {
    clock_t start = clock();

    char *buffer, *tmp;
    size_t bytes, lines = 0;
    int i;
    FILE *fp = fopen(fname, "rb");

    buffer = (char*) malloc(CHUNKSIZE * sizeof(char));

    while(1) {
        tmp = buffer;
        bytes = fread((void *) buffer, sizeof(char), CHUNKSIZE, fp);
        // printf("Read %d bytes\n", bytes);

        if (bytes == 0) {
            break;
        }

        for (i = 0; i < bytes; ++i)
        {
            if (*tmp++ == '\n') {
                lines++;
            }
        }
    }


    printf("Saw %d lines\n", (int) lines);

    free(buffer);
    fclose(fp);

    printf("Time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
}


int main(int argc, char *argv[])
{
    // import_array();

    int i;
    TRACE(("hello: %s\n", "Wes"));

    test_tokenize("/Users/wesm/code/pandas/pandas/io/tests/test1.csv");

    // char *msg = (char*) malloc(50);
    // sprintf(msg, "Hello: %s\n", "wes");
    // printf("%s", msg);

    /* for (i = 0; i < 10; ++i) */
    /* { */
    /*  test_count_lines("../foo.csv"); */
    /* } */

    return 0;
}
