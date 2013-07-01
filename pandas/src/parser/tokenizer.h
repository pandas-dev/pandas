/*

Copyright (c) 2012, Lambda Foundry, Inc., except where noted

Incorporates components of WarrenWeckesser/textreader, licensed under 3-clause
BSD

See LICENSE for the license

*/

#ifndef _PARSER_COMMON_H_
#define _PARSER_COMMON_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>

#include <ctype.h>

#define ERROR_OK             0
#define ERROR_NO_DIGITS      1
#define ERROR_OVERFLOW       2
#define ERROR_INVALID_CHARS  3
#define ERROR_MINUS_SIGN     4

#if defined(_MSC_VER)
#include "../headers/ms_stdint.h"
#else
#include <stdint.h>
#endif

#include "khash.h"

#define CHUNKSIZE 1024*256
#define KB 1024
#define MB 1024 * KB
#define STREAM_INIT_SIZE 32

#define REACHED_EOF 1
#define CALLING_READ_FAILED 2

#ifndef P_INLINE
  #if defined(__GNUC__)
    #define P_INLINE __inline__
  #elif defined(_MSC_VER)
    #define P_INLINE
  #elif defined (__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
    #define P_INLINE inline
  #else
    #define P_INLINE
  #endif
#endif

#if defined(_MSC_VER)
#define strtoll _strtoi64
#endif

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


/* #define VERBOSE */

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
    EAT_COMMENT,
    FINISHED
} ParserState;

typedef enum {
    QUOTE_MINIMAL, QUOTE_ALL, QUOTE_NONNUMERIC, QUOTE_NONE
} QuoteStyle;


typedef void* (*io_callback)(void *src, size_t nbytes, size_t *bytes_read,
                            int *status);
typedef int (*io_cleanup)(void *src);

typedef struct parser_t {
    void *source;
    io_callback cb_io;
    io_cleanup cb_cleanup;

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
    int lines;            // Number of (good) lines observedb
    int file_lines;       // Number of file lines observed (including bad or skipped)
    int lines_cap;        // Vector capacity

    // Tokenizing stuff
    ParserState state;
    int doublequote;            /* is " represented by ""? */
    char delimiter;             /* field separator */
    int delim_whitespace;       /* delimit by consuming space/tabs instead */
    char quotechar;             /* quote character */
    char escapechar;            /* escape character */
    char lineterminator;
    int skipinitialspace;       /* ignore spaces following delimiter? */
    int quoting;                /* style of quoting to write */

    // krufty, hmm =/
    int numeric_field;

    char commentchar;
    int allow_embedded_newline;
    int strict;                 /* raise exception on bad CSV */

    int expected_fields;
    int error_bad_lines;
    int warn_bad_lines;

    // floating point options
    char decimal;
    char sci;

    // thousands separator (comma, period)
    char thousands;

    int header; // Boolean: 1: has header, 0: no header
    int header_start; // header row start
    int header_end;   // header row end

    void *skipset;
    int skip_footer;

    // error handling
    char *warn_msg;
    char *error_msg;
} parser_t;




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

int P_INLINE to_double(char *item, double *p_value, char sci, char decimal);
int P_INLINE to_complex(char *item, double *p_real, double *p_imag, char sci, char decimal);
int P_INLINE to_longlong(char *item, long long *p_value);
int P_INLINE to_longlong_thousands(char *item, long long *p_value, char tsep);
int P_INLINE to_boolean(char *item, uint8_t *val);

#endif // _PARSER_COMMON_H_
