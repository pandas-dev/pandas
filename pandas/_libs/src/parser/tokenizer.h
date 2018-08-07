/*

Copyright (c) 2012, Lambda Foundry, Inc., except where noted

Incorporates components of WarrenWeckesser/textreader, licensed under 3-clause
BSD

See LICENSE for the license

*/

#ifndef PANDAS__LIBS_SRC_PARSER_TOKENIZER_H_
#define PANDAS__LIBS_SRC_PARSER_TOKENIZER_H_

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "Python.h"

#include <ctype.h>

#define ERROR_OK 0
#define ERROR_NO_DIGITS 1
#define ERROR_OVERFLOW 2
#define ERROR_INVALID_CHARS 3

#include "../headers/stdint.h"
#include "../helper.h"

#include "khash.h"

#define CHUNKSIZE 1024 * 256
#define KB 1024
#define MB 1024 * KB
#define STREAM_INIT_SIZE 32

#define REACHED_EOF 1
#define CALLING_READ_FAILED 2


#if defined(_MSC_VER)
#define strtoll _strtoi64
#endif

/*

  C flat file parsing low level code for pandas / NumPy

 */

#define FALSE 0
#define TRUE 1

// Maximum number of columns in a file.
#define MAX_NUM_COLUMNS 2000

// Maximum number of characters in single field.
#define FIELD_BUFFER_SIZE 2000

/*
 *  Common set of error types for the read_rows() and tokenize()
 *  functions.
 */
#define ERROR_OUT_OF_MEMORY 1
#define ERROR_INVALID_COLUMN_INDEX 10
#define ERROR_CHANGED_NUMBER_OF_FIELDS 12
#define ERROR_TOO_MANY_CHARS 21
#define ERROR_TOO_MANY_FIELDS 22
#define ERROR_NO_DATA 23

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
    EAT_CRNL_NOP,
    EAT_WHITESPACE,
    EAT_COMMENT,
    EAT_LINE_COMMENT,
    WHITESPACE_LINE,
    START_FIELD_IN_SKIP_LINE,
    IN_FIELD_IN_SKIP_LINE,
    IN_QUOTED_FIELD_IN_SKIP_LINE,
    QUOTE_IN_QUOTED_FIELD_IN_SKIP_LINE,
    FINISHED
} ParserState;

typedef enum {
    QUOTE_MINIMAL,
    QUOTE_ALL,
    QUOTE_NONNUMERIC,
    QUOTE_NONE
} QuoteStyle;

typedef void *(*io_callback)(void *src, size_t nbytes, size_t *bytes_read,
                             int *status);
typedef int (*io_cleanup)(void *src);

typedef struct parser_t {
    void *source;
    io_callback cb_io;
    io_cleanup cb_cleanup;

    int64_t chunksize;      // Number of bytes to prepare for each chunk
    char *data;             // pointer to data to be processed
    int64_t datalen;        // amount of data available
    int64_t datapos;

    // where to write out tokenized data
    char *stream;
    int64_t stream_len;
    int64_t stream_cap;

    // Store words in (potentially ragged) matrix for now, hmm
    char **words;
    int64_t *word_starts;   // where we are in the stream
    int64_t words_len;
    int64_t words_cap;

    char *pword_start;      // pointer to stream start of current field
    int64_t word_start;     // position start of current field

    int64_t *line_start;    // position in words for start of line
    int64_t *line_fields;   // Number of fields in each line
    int64_t lines;          // Number of (good) lines observed
    int64_t file_lines;     // Number of lines (including bad or skipped)
    int64_t lines_cap;      // Vector capacity

    // Tokenizing stuff
    ParserState state;
    int doublequote;      /* is " represented by ""? */
    char delimiter;       /* field separator */
    int delim_whitespace; /* delimit by consuming space/tabs instead */
    char quotechar;       /* quote character */
    char escapechar;      /* escape character */
    char lineterminator;
    int skipinitialspace; /* ignore spaces following delimiter? */
    int quoting;          /* style of quoting to write */

    // krufty, hmm =/
    int numeric_field;

    char commentchar;
    int allow_embedded_newline;
    int strict; /* raise exception on bad CSV */

    int usecols;  // Boolean: 1: usecols provided, 0: none provided

    int expected_fields;
    int error_bad_lines;
    int warn_bad_lines;

    // floating point options
    char decimal;
    char sci;

    // thousands separator (comma, period)
    char thousands;

    int header;            // Boolean: 1: has header, 0: no header
    int64_t header_start;  // header row start
    int64_t header_end;    // header row end

    void *skipset;
    PyObject *skipfunc;
    int64_t skip_first_N_rows;
    int skip_footer;
    // pick one, depending on whether the converter requires GIL
    double (*double_converter_nogil)(const char *, char **,
                                     char, char, char, int);
    double (*double_converter_withgil)(const char *, char **,
                                       char, char, char, int);

    // error handling
    char *warn_msg;
    char *error_msg;

    int skip_empty_lines;
} parser_t;

typedef struct coliter_t {
    char **words;
    int64_t *line_start;
    int col;
} coliter_t;

void coliter_setup(coliter_t *self, parser_t *parser, int i, int start);
coliter_t *coliter_new(parser_t *self, int i);

#define COLITER_NEXT(iter, word)                          \
    do {                                                  \
        const int64_t i = *iter.line_start++ + iter.col;      \
        word = i < *iter.line_start ? iter.words[i] : ""; \
    } while (0)

parser_t *parser_new(void);

int parser_init(parser_t *self);

int parser_consume_rows(parser_t *self, size_t nrows);

int parser_trim_buffers(parser_t *self);

int parser_add_skiprow(parser_t *self, int64_t row);

int parser_set_skipfirstnrows(parser_t *self, int64_t nrows);

void parser_free(parser_t *self);

void parser_del(parser_t *self);

void parser_set_default_options(parser_t *self);

int tokenize_nrows(parser_t *self, size_t nrows);

int tokenize_all_rows(parser_t *self);

// Have parsed / type-converted a chunk of data
// and want to free memory from the token stream

typedef struct uint_state {
    int seen_sint;
    int seen_uint;
    int seen_null;
} uint_state;

void uint_state_init(uint_state *self);

int uint64_conflict(uint_state *self);

uint64_t str_to_uint64(uint_state *state, const char *p_item, int64_t int_max,
                       uint64_t uint_max, int *error, char tsep);
int64_t str_to_int64(const char *p_item, int64_t int_min, int64_t int_max,
                     int *error, char tsep);
double xstrtod(const char *p, char **q, char decimal, char sci, char tsep,
               int skip_trailing);
double precise_xstrtod(const char *p, char **q, char decimal, char sci,
                       char tsep, int skip_trailing);
double round_trip(const char *p, char **q, char decimal, char sci, char tsep,
                  int skip_trailing);
int to_boolean(const char *item, uint8_t *val);

#endif  // PANDAS__LIBS_SRC_PARSER_TOKENIZER_H_
