/*

Copyright (c) 2012, Lambda Foundry, Inc., except where noted

Incorporates components of WarrenWeckesser/textreader, licensed under 3-clause
BSD

See LICENSE for the license

*/

#pragma once

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#define ERROR_NO_DIGITS 1
#define ERROR_OVERFLOW 2
#define ERROR_INVALID_CHARS 3

#include <stdint.h>

#include "pandas/parser/pd_strtoi.h"

#ifdef __cplusplus
extern "C" {
#endif

#define STREAM_INIT_SIZE 32

#define REACHED_EOF 1
#define CALLING_READ_FAILED 2

/*

  C flat file parsing low level code for pandas / NumPy

 */

/*
 *  Common set of error types for the read_rows() and tokenize()
 *  functions.
 */

#define PARSER_OUT_OF_MEMORY -1

/*
 *  TODO: Might want to couple count_rows() with read_rows() to avoid
 *        duplication of some file I/O.
 */

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

typedef enum { BLHM_ERROR, BLHM_WARN, BLHM_SKIP } BadLineHandleMethod;

typedef char *(*io_callback)(void *src, size_t nbytes, size_t *bytes_read,
                             int *status, const char *encoding_errors);
typedef void (*io_cleanup)(void *src);

typedef struct parser_t {
  void *source;
  io_callback cb_io;
  io_cleanup cb_cleanup;

  int64_t chunksize; // Number of bytes to prepare for each chunk
  char *data;        // pointer to data to be processed
  int64_t datalen;   // amount of data available
  int64_t datapos;

  // where to write out tokenized data
  char *stream;
  uint64_t stream_len;
  uint64_t stream_cap;

  // Store words in (potentially ragged) matrix for now, hmm
  char **words;
  int64_t *word_starts; // where we are in the stream
  uint64_t words_len;
  uint64_t words_cap;
  uint64_t max_words_cap; // maximum word cap encountered

  char *pword_start;  // pointer to stream start of current field
  int64_t word_start; // position start of current field

  int64_t *line_start;  // position in words for start of line
  int64_t *line_fields; // Number of fields in each line
  uint64_t lines;       // Number of (good) lines observed
  uint64_t file_lines;  // Number of lines (including bad or skipped)
  uint64_t lines_cap;   // Vector capacity

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

  char commentchar;
  int allow_embedded_newline;

  int usecols; // Boolean: 1: usecols provided, 0: none provided

  Py_ssize_t expected_fields;
  BadLineHandleMethod on_bad_lines;

  // floating point options
  char decimal;
  char sci;

  // thousands separator (comma, period)
  char thousands;

  int header;           // Boolean: 1: has header, 0: no header
  int64_t header_start; // header row start
  uint64_t header_end;  // header row end

  void *skipset;
  PyObject *skipfunc;
  int64_t skip_first_N_rows;
  int64_t skip_footer;
  double (*double_converter)(const char *, char **, char, char, char, int,
                             int *, int *, const char *);

  // error handling
  char *warn_msg;
  char *error_msg;

  int skip_empty_lines;

  // Boolean: 1 when the input was pre-loaded via TextReader.load_buffer.
  // The buffer then never starts with a header row, so the first line gets
  // no special treatment: no BOM strip, no exemption from field-count checks.
  int preloaded;
} parser_t;

typedef struct coliter_t {
  char **words;
  int64_t *line_start;
  int64_t col;
} coliter_t;

void coliter_setup(coliter_t *self, parser_t *parser, int64_t i, int64_t start);

// Advance the column iterator and return the next field's token, emitting its
// resolved index via idx_out so callers needing the token length can compute
// it as word_starts[idx+1] - word_starts[idx] - 1 (using parser->stream_len
// for the last token where idx+1 == words_len). A missing field yields "" and
// idx_out = -1, which callers must treat as length 0 rather than indexing into
// word_starts.
static inline const char *coliter_next_with_idx(coliter_t *self,
                                                int64_t *idx_out) {
  const int64_t idx = *self->line_start++ + self->col;
  if (idx >= *self->line_start) {
    *idx_out = -1;
    return "";
  }
  *idx_out = idx;
  return self->words[idx];
}

static inline const char *coliter_next(coliter_t *self) {
  int64_t idx;
  return coliter_next_with_idx(self, &idx);
}

parser_t *parser_new(void);

int parser_init(parser_t *self);

int parser_consume_rows(parser_t *self, uint64_t nrows);

int parser_trim_buffers(parser_t *self);

int parser_add_skiprow(parser_t *self, int64_t row);

void parser_set_skipfirstnrows(parser_t *self, int64_t nrows);

void parser_free(parser_t *self);

void parser_del(parser_t *self);

void parser_set_default_options(parser_t *self);

int tokenize_nrows(parser_t *self, uint64_t nrows, const char *encoding_errors);

int tokenize_all_rows(parser_t *self, const char *encoding_errors);

// Have parsed / type-converted a chunk of data
// and want to free memory from the token stream

// uint_state itself is defined in pd_strtoi.h (included above) alongside
// str_to_uint64; these helpers remain implemented in tokenizer.c.
void uint_state_init(uint_state *self);

int uint64_conflict(uint_state *self);

double precise_xstrtod(const char *p, char **q, char decimal, char sci,
                       char tsep, int skip_trailing, int *error,
                       int *maybe_int);
// As precise_xstrtod, but takes the known end of the token (one past its
// last byte) to skip the end-of-token scan; pass NULL to locate it as usual.
double precise_xstrtod_with_end(const char *p, char **q, char decimal, char sci,
                                char tsep, int skip_trailing, int *error,
                                int *maybe_int, const char *end);
// Hot-path double parse (fast_float) for read_csv with default settings (no
// thousands separator, 'e'/'E' exponent). Returns 0 and writes *out only when
// a plain numeric token parses cleanly and consumes exactly [start, end);
// anything else — leading/trailing spaces, inf/nan spellings, junk — returns
// nonzero so the caller retries through the full converter with its legacy
// semantics. Overflow parses to +/-inf like precise_xstrtod. On nonzero
// return *out may still have been clobbered.
int try_parse_plain_double(const char *start, const char *end, char decimal,
                           double *out);
int to_boolean(const char *item, uint8_t *val);

#ifdef __cplusplus
} // extern "C"
#endif
