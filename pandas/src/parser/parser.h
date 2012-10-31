/*

Copyright (c) 2012, Lambda Foundry, Inc.

See LICENSE for the license

*/

#ifndef _PARSER_COMMON_H_
#define _PARSER_COMMON_H_

#include "Python.h"

/* #include "structmember.h" */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <errno.h>

#if defined(_MSC_VER)
#include "ms_stdint.h"
#else
#include <stdint.h>
#endif

#include "khash.h"

#define CHUNKSIZE 1024*256
#define KB 1024
#define MB 1024 * KB
#define STREAM_INIT_SIZE 32


/*

  C flat file parsing low level code for pandas / NumPy

 */

#define FALSE 0
#define TRUE  1

/* Maximum number of columns in a file. */
#define MAX_NUM_COLUMNS    2000

/* Maximum number of characters in single field. */

#define FIELD_BUFFER_SIZE  2000


/*
 *  Common set of error types for the read_rows() and tokenize()
 *  functions.
 */

#define ERROR_OUT_OF_MEMORY             1
#define ERROR_INVALID_COLUMN_INDEX     10
#define ERROR_CHANGED_NUMBER_OF_FIELDS 12
#define ERROR_TOO_MANY_CHARS           21
#define ERROR_TOO_MANY_FIELDS          22
#define ERROR_NO_DATA                  23


// #define VERBOSE

#if defined(VERBOSE)
#define TRACE(X) printf X;
#else
#define TRACE(X)
#endif


#define PARSER_OUT_OF_MEMORY -1


/*
 *  XXX Might want to couple count_rows() with read_rows() to avoid duplication
 *      of some file I/O.
 */

/*
 *  WORD_BUFFER_SIZE determines the maximum amount of non-delimiter
 *  text in a row.
 */
#define WORD_BUFFER_SIZE 4000


typedef enum {
    START_RECORD,
    START_FIELD,
    ESCAPED_CHAR,
    IN_FIELD,
    IN_QUOTED_FIELD,
    ESCAPE_IN_QUOTED_FIELD,
    QUOTE_IN_QUOTED_FIELD,
    EAT_CRNL,
    EAT_WHITESPACE,
    FINISHED
} ParserState;

typedef enum {
    QUOTE_MINIMAL, QUOTE_ALL, QUOTE_NONNUMERIC, QUOTE_NONE
} QuoteStyle;


typedef struct parser_t {
    void *source;
    char sourcetype;   // 'M' for mmap, 'F' for FILE, 'A' for array

    int chunksize;  // Number of bytes to prepare for each chunk
    char *data;     // pointer to data to be processed
    int datalen;    // amount of data available
    int datapos;

    // where to write out tokenized data
    char *stream;
    int stream_len;
    int stream_cap;

    // Store words in (potentially ragged) matrix for now, hmm
    char **words;
    int *word_starts; // where we are in the stream
    int words_len;
    int words_cap;

    char *pword_start;    // pointer to stream start of current field
    int word_start;       // position start of current field

    int *line_start;      // position in words for start of line
    int *line_fields;     // Number of fields in each line
    int lines;            // Number of (good) lines observed
    int file_lines;       // Number of file lines observed (including bad or skipped)
    int lines_cap;        // Vector capacity

    // Tokenizing stuff
    ParserState state;
    int doublequote;            /* is " represented by ""? */
    char delimiter;             /* field separator */
    int delim_whitespace;       /* delimit by consuming space/tabs instead */
    char quotechar;             /* quote character */
    char escapechar;            /* escape character */
    int skipinitialspace;       /* ignore spaces following delimiter? */
    int quoting;                /* style of quoting to write */

    // krufty, hmm =/
    int numeric_field;

    char commentchar;
    int allow_embedded_newline;
    int strict;                 /* raise exception on bad CSV */

    int error_bad_lines;
    int warn_bad_lines;

    int infer_types;

    // floating point options
    char decimal;
    char sci;

    // thousands separator (comma, period)
    char thousands;

    int header; // Boolean: 1: has header, 0: no header

    void *skipset;
    int skip_footer;

    // error handling
    char *error_msg;
} parser_t;


void *safe_realloc(void *buffer, size_t size);


typedef struct coliter_t {
    char **words;
    int *line_start;
    int col;
} coliter_t;

void coliter_setup(coliter_t *self, parser_t *parser, int i, int start);
coliter_t *coliter_new(parser_t *self, int i);

/* #define COLITER_NEXT(iter) iter->words[iter->line_start[iter->line++] + iter->col] */
// #define COLITER_NEXT(iter) iter.words[iter.line_start[iter.line++] + iter.col]

#define COLITER_NEXT(iter) iter.words[*iter.line_start++ + iter.col]

parser_t* parser_new();

int parser_init(parser_t *self);

int parser_file_source_init(parser_t *self, FILE *fp);

int parser_mmap_init(parser_t *self, FILE *fp);

int parser_rd_source_init(parser_t *self, PyObject *source);

int parser_gzip_source_init(parser_t *self, FILE *fp);

int parser_consume_rows(parser_t *self, size_t nrows);

int parser_trim_buffers(parser_t *self);

int parser_add_skiprow(parser_t *self, int64_t row);

void parser_free(parser_t *self);

void parser_set_default_options(parser_t *self);

void debug_print_parser(parser_t *self);

int tokenize_nrows(parser_t *self, size_t nrows);

int tokenize_all_rows(parser_t *self);

/*

  Have parsed / type-converted a chunk of data and want to free memory from the
  token stream

 */
int clear_parsed_lines(parser_t *self, size_t nlines);

int64_t str_to_int64(const char *p_item, int64_t int_min,
                     int64_t int_max, int *error, char tsep);
uint64_t str_to_uint64(const char *p_item, uint64_t uint_max, int *error);

int inline to_double(char *item, double *p_value, char sci, char decimal);
int inline to_complex(char *item, double *p_real, double *p_imag, char sci, char decimal);
int inline to_longlong(char *item, long long *p_value);
int inline to_longlong_thousands(char *item, long long *p_value, char tsep);
int inline to_boolean(char *item, uint8_t *val);

#endif // _PARSER_COMMON_H_
