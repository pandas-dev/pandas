/*

Copyright (c) 2012, Lambda Foundry, Inc., except where noted

Incorporates components of WarrenWeckesser/textreader, licensed under 3-clause
BSD

See LICENSE for the license

*/

/*

Low-level ascii-file processing for pandas. Combines some elements from
Python's built-in csv module and Warren Weckesser's textreader project on
GitHub. See Python Software Foundation License and BSD licenses for these.

*/
#include "pandas/parser/tokenizer.h"
#include "pandas/parser/pd_strtoi.h"
#include "pandas/portable.h"

#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

// Arch selection comes from the build-time SIMD config (pandas/_libs/simd),
// which sets exactly one of PANDAS_HAVE_NEON / PANDAS_HAVE_SSE2 /
// PANDAS_HAVE_SCALAR from the host CPU family. Pull in the matching intrinsics
// header for whichever applies; the scalar case defines neither and falls back
// to the byte-at-a-time scan.
#include "pandas_simd_config.h"

#if defined(PANDAS_HAVE_NEON)
#  include <arm_neon.h>
#  ifdef _MSC_VER
#    include <intrin.h>
#  endif
#elif defined(PANDAS_HAVE_SSE2)
#  include <emmintrin.h>
#  ifdef _MSC_VER
#    include <intrin.h>
#  endif
#endif

// Portable count-trailing-zeros, defined only for the architecture that uses
// it: pandas_ctzll for NEON (64-bit lane), pandas_ctz for SSE2 (32-bit
// movemask). The __builtin_ctz* form is preferred where available (including
// clang-cl and mingw on Windows); MSVC falls back to _BitScanForward*.
#ifndef __has_builtin
#  define __has_builtin(x) 0
#endif

#if defined(PANDAS_HAVE_NEON)
static inline int pandas_ctzll(unsigned long long mask) {
#  if __has_builtin(__builtin_ctzll)
  return __builtin_ctzll(mask);
#  elif defined(_MSC_VER)
  unsigned long index;
  _BitScanForward64(&index, mask);
  return (int)index;
#  else
#    error "no count-trailing-zeros builtin available"
#  endif
}
#elif defined(PANDAS_HAVE_SSE2)
static inline int pandas_ctz(unsigned int mask) {
#  if __has_builtin(__builtin_ctz)
  return __builtin_ctz(mask);
#  elif defined(_MSC_VER)
  unsigned long index;
  _BitScanForward(&index, mask);
  return (int)index;
#  else
#    error "no count-trailing-zeros builtin available"
#  endif
}
#endif

#include "pandas/portable.h"
#include "pandas/vendored/klib/khash.h" // for kh_int64_t, kh_destroy_int64

// Arrow256 allows up to 76 decimal digits.
// We rounded up to the next power of 2.
#define PROCESSED_WORD_CAPACITY 128

void coliter_setup(coliter_t *self, parser_t *parser, int64_t i,
                   int64_t start) {
  // column i, starting at 0
  self->words = parser->words;
  self->col = i;
  self->line_start = parser->line_start + start;
}

static void free_if_not_null(void **ptr) {
  if (*ptr != NULL) {
    free(*ptr);
    *ptr = NULL;
  }
}

/*

  Parser / tokenizer

*/

static void *grow_buffer(void *buffer, uint64_t length, uint64_t *capacity,
                         int64_t space, int64_t elsize, int *error) {
  uint64_t cap = *capacity;
  void *newbuffer = buffer;

  // Can we fit potentially nbytes tokens (+ null terminators) in the stream?
  while ((length + space >= cap) && (newbuffer != NULL)) {
    cap = cap ? cap << 1 : 2;
    buffer = newbuffer;
    newbuffer = realloc(newbuffer, elsize * cap);
  }

  if (newbuffer == NULL) {
    // realloc failed so don't change *capacity, set *error to errno
    // and return the last good realloc'd buffer so it can be freed
    *error = errno;
    newbuffer = buffer;
  } else {
    // realloc worked, update *capacity and set *error to 0
    // sigh, multiple return values
    *capacity = cap;
    *error = 0;
  }
  return newbuffer;
}

void parser_set_default_options(parser_t *self) {
  self->decimal = '.';
  self->sci = 'E';

  // For tokenization
  self->state = START_RECORD;

  self->delimiter = ','; // XXX
  self->delim_whitespace = 0;

  self->doublequote = 0;
  self->quotechar = '"';
  self->escapechar = 0;

  self->lineterminator = '\0'; /* NUL->standard logic */

  self->skipinitialspace = 0;
  self->quoting = QUOTE_MINIMAL;
  self->allow_embedded_newline = 1;

  self->expected_fields = -1;
  self->on_bad_lines = BLHM_ERROR;
  self->preloaded = 0;

  self->commentchar = '#';
  self->thousands = '\0';

  self->skipset = NULL;
  self->skipfunc = NULL;
  self->skip_first_N_rows = -1;
  self->skip_footer = 0;
}

parser_t *parser_new(void) { return (parser_t *)calloc(1, sizeof(parser_t)); }

static void parser_clear_data_buffers(parser_t *self) {
  free_if_not_null((void *)&self->stream);
  free_if_not_null((void *)&self->words);
  free_if_not_null((void *)&self->word_starts);
  free_if_not_null((void *)&self->line_start);
  free_if_not_null((void *)&self->line_fields);
}

static void parser_cleanup(parser_t *self) {
  // self can be NULL when cleanup runs on a TextReader whose __cinit__
  // raised before parser_new() was called (GH#53131).
  if (self == NULL) {
    return;
  }

  // XXX where to put this
  free_if_not_null((void *)&self->error_msg);
  free_if_not_null((void *)&self->warn_msg);

  if (self->skipset != NULL) {
    kh_destroy_int64((kh_int64_t *)self->skipset);
    self->skipset = NULL;
  }

  parser_clear_data_buffers(self);
  if (self->cb_cleanup != NULL) {
    self->cb_cleanup(self->source);
    self->cb_cleanup = NULL;
  }
}

int parser_init(parser_t *self) {
  /*
    Initialize data buffers
  */

  self->stream = NULL;
  self->words = NULL;
  self->word_starts = NULL;
  self->line_start = NULL;
  self->line_fields = NULL;
  self->error_msg = NULL;
  self->warn_msg = NULL;

  // token stream
  self->stream = malloc(STREAM_INIT_SIZE);
  if (self->stream == NULL) {
    parser_cleanup(self);
    return PARSER_OUT_OF_MEMORY;
  }
  self->stream_cap = STREAM_INIT_SIZE;
  self->stream_len = 0;

  // word pointers and metadata
  _Static_assert(STREAM_INIT_SIZE / 10 > 0,
                 "STREAM_INIT_SIZE must be defined and >= 10");
  const int64_t sz = STREAM_INIT_SIZE / 10;
  self->words = malloc(sz * sizeof(char *));
  self->word_starts = malloc(sz * sizeof(int64_t));
  self->max_words_cap = sz;
  self->words_cap = sz;
  self->words_len = 0;

  // line pointers and metadata
  self->line_start = malloc(sz * sizeof(int64_t));

  self->line_fields = malloc(sz * sizeof(int64_t));

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
  self->warn_msg = NULL;

  self->commentchar = '\0';

  return 0;
}

void parser_free(parser_t *self) {
  // opposite of parser_init
  parser_cleanup(self);
}

void parser_del(parser_t *self) { free(self); }

static int make_stream_space(parser_t *self, size_t nbytes) {
  // Can we fit potentially nbytes tokens (+ null terminators) in the stream?

  /*
    TOKEN STREAM
  */

  int status;
  char *orig_ptr = (void *)self->stream;
  self->stream = (char *)grow_buffer((void *)self->stream, self->stream_len,
                                     &self->stream_cap, nbytes * 2, 1, &status);

  if (status != 0) {
    return PARSER_OUT_OF_MEMORY;
  }

  // realloc sets errno when moving buffer?
  if (self->stream != orig_ptr) {
    self->pword_start = self->stream + self->word_start;

    for (uint64_t i = 0; i < self->words_len; ++i) {
      self->words[i] = self->stream + self->word_starts[i];
    }
  }

  /*
    WORD VECTORS
  */

  const uint64_t words_cap = self->words_cap;

  /**
   * If we are reading in chunks, we need to be aware of the maximum number
   * of words we have seen in previous chunks (self->max_words_cap), so
   * that way, we can properly allocate when reading subsequent ones.
   *
   * Otherwise, we risk a buffer overflow if we mistakenly under-allocate
   * just because a recent chunk did not have as many words.
   */
  const uint64_t length = self->words_len + nbytes < self->max_words_cap
                              ? self->max_words_cap - nbytes - 1
                              : self->words_len;

  self->words =
      (char **)grow_buffer((void *)self->words, length, &self->words_cap,
                           nbytes, sizeof(char *), &status);

  if (status != 0) {
    return PARSER_OUT_OF_MEMORY;
  }

  // realloc took place
  if (words_cap != self->words_cap) {
    int64_t *newptr = (int64_t *)realloc(self->word_starts,
                                         sizeof(int64_t) * self->words_cap);
    if (newptr == NULL) {
      return PARSER_OUT_OF_MEMORY;
    } else {
      self->word_starts = newptr;
    }
  }

  /*
    LINE VECTORS
  */
  const uint64_t lines_cap = self->lines_cap;
  self->line_start = (int64_t *)grow_buffer((void *)self->line_start,
                                            self->lines + 1, &self->lines_cap,
                                            nbytes, sizeof(int64_t), &status);
  if (status != 0) {
    return PARSER_OUT_OF_MEMORY;
  }

  // realloc took place
  if (lines_cap != self->lines_cap) {
    int64_t *newptr = (int64_t *)realloc(self->line_fields,
                                         sizeof(int64_t) * self->lines_cap);
    if (newptr == NULL) {
      return PARSER_OUT_OF_MEMORY;
    } else {
      self->line_fields = newptr;
    }
  }

  return 0;
}

static int push_char(parser_t *self, char c) {
  if (self->stream_len >= self->stream_cap) {
    const size_t bufsize = 100;
    self->error_msg = malloc(bufsize);
    snprintf(self->error_msg, bufsize,
             "Buffer overflow caught - possible malformed input file.\n");
    return PARSER_OUT_OF_MEMORY;
  }
  self->stream[self->stream_len++] = c;
  return 0;
}

static inline int end_field(parser_t *self) {
  // XXX cruft
  if (self->words_len >= self->words_cap) {
    const size_t bufsize = 100;
    self->error_msg = malloc(bufsize);
    snprintf(self->error_msg, bufsize,
             "Buffer overflow caught - possible malformed input file.\n");
    return PARSER_OUT_OF_MEMORY;
  }

  // null terminate token
  push_char(self, '\0');

  // set pointer and metadata
  self->words[self->words_len] = self->pword_start;
  self->word_starts[self->words_len] = self->word_start;
  self->words_len++;

  // increment line field count
  self->line_fields[self->lines]++;

  // New field begin in stream
  self->pword_start = self->stream + self->stream_len;
  self->word_start = self->stream_len;

  return 0;
}

static void append_warning(parser_t *self, const char *msg) {
  const int64_t length = strlen(msg);

  if (self->warn_msg == NULL) {
    self->warn_msg = malloc(length + 1);
    snprintf(self->warn_msg, length + 1, "%s", msg);
  } else {
    const int64_t ex_length = strlen(self->warn_msg);
    char *newptr = (char *)realloc(self->warn_msg, ex_length + length + 1);
    if (newptr != NULL) {
      self->warn_msg = newptr;
      snprintf(self->warn_msg + ex_length, length + 1, "%s", msg);
    }
  }
}

static int end_line(parser_t *self) {
  int64_t ex_fields = self->expected_fields;
  int64_t fields = self->line_fields[self->lines];

  if (self->lines > 0) {
    if (self->expected_fields >= 0) {
      ex_fields = self->expected_fields;
    } else {
      ex_fields = self->line_fields[self->lines - 1];
    }
  }

  if (self->state == START_FIELD_IN_SKIP_LINE ||
      self->state == IN_FIELD_IN_SKIP_LINE ||
      self->state == IN_QUOTED_FIELD_IN_SKIP_LINE ||
      self->state == QUOTE_IN_QUOTED_FIELD_IN_SKIP_LINE) {
    // increment file line count
    self->file_lines++;

    // skip the tokens from this bad line
    self->line_start[self->lines] += fields;

    // reset field count
    self->line_fields[self->lines] = 0;
    return 0;
  }

  /* The first line is normally exempt from the field-count check: it is the
   * header, or it defines the table width (implicit-index inference).  In
   * preloaded mode (TextReader.load_buffer) the width was fixed up front and
   * the first line is plain data, so the exemption must not apply. */
  if ((self->preloaded || !(self->lines <= self->header_end + 1)) &&
      (fields > ex_fields) && !(self->usecols)) {
    // increment file line count
    self->file_lines++;

    // skip the tokens from this bad line
    self->line_start[self->lines] += fields;

    // reset field count
    self->line_fields[self->lines] = 0;

    // file_lines is now the actual file line number (starting at 1)
    if (self->on_bad_lines == BLHM_ERROR) {
      const size_t bufsize = 100;
      self->error_msg = malloc(bufsize);
      snprintf(self->error_msg, bufsize,
               "Expected %" PRId64 " fields in line %" PRIu64 ", saw %" PRId64
               "\n",
               ex_fields, self->file_lines, fields);
      return -1;
    } else {
      // simply skip bad lines
      if (self->on_bad_lines == BLHM_WARN) {
        // pass up error message
        const size_t bufsize = 100;
        char *msg = (char *)malloc(bufsize);
        snprintf(msg, bufsize,
                 "Skipping line %" PRIu64 ": expected %" PRId64
                 " fields, saw %" PRId64 "\n",
                 self->file_lines, ex_fields, fields);
        append_warning(self, msg);
        free(msg);
      }
    }
  } else {
    // missing trailing delimiters
    if ((self->lines >= self->header_end + 1) && fields < ex_fields) {
      // might overrun the buffer when closing fields
      if (make_stream_space(self, ex_fields - fields) < 0) {
        const size_t bufsize = 100;
        self->error_msg = malloc(bufsize);
        snprintf(self->error_msg, bufsize, "out of memory");
        return -1;
      }

      while (fields < ex_fields) {
        end_field(self);
        fields++;
      }
    }

    // increment both line counts
    self->file_lines++;
    self->lines++;

    // good line, set new start point
    if (self->lines >= self->lines_cap) {
      const size_t bufsize = 100;
      self->error_msg = malloc(bufsize);
      snprintf(self->error_msg, bufsize,
               "Buffer overflow caught - "
               "possible malformed input file.\n");
      return PARSER_OUT_OF_MEMORY;
    }
    self->line_start[self->lines] =
        (self->line_start[self->lines - 1] + fields);

    // new line start with 0 fields
    self->line_fields[self->lines] = 0;
  }
  return 0;
}

int parser_add_skiprow(parser_t *self, int64_t row) {
  khiter_t k;
  kh_int64_t *set;
  int ret = 0;

  if (self->skipset == NULL) {
    self->skipset = (void *)kh_init_int64();
  }

  set = (kh_int64_t *)self->skipset;

  k = kh_put_int64(set, row, &ret);
  set->keys[k] = row;

  return 0;
}

void parser_set_skipfirstnrows(parser_t *self, int64_t nrows) {
  // self->file_lines is zero based so subtract 1 from nrows
  if (nrows > 0) {
    self->skip_first_N_rows = nrows - 1;
  }
}

static int parser_buffer_bytes(parser_t *self, size_t nbytes,
                               const char *encoding_errors) {
  int status;
  size_t bytes_read;

  status = 0;
  self->datapos = 0;
  self->data =
      self->cb_io(self->source, nbytes, &bytes_read, &status, encoding_errors);
  self->datalen = bytes_read;

  if (status != REACHED_EOF && self->data == NULL) {
    const size_t bufsize = 200;
    self->error_msg = malloc(bufsize);

    if (status == CALLING_READ_FAILED) {
      snprintf(self->error_msg, bufsize,
               "Calling read(nbytes) on source failed. "
               "Try engine='python'.");
    } else {
      snprintf(self->error_msg, bufsize, "Unknown error in IO callback");
    }
    return -1;
  }
  return status;
}

/*

  Tokenization macros and state machine code

*/

#define PUSH_CHAR(c)                                                           \
  if (slen >= self->stream_cap) {                                              \
    const size_t bufsize = 100;                                                \
    self->error_msg = malloc(bufsize);                                         \
    snprintf(self->error_msg, bufsize,                                         \
             "Buffer overflow caught - possible malformed input file.\n");     \
    return PARSER_OUT_OF_MEMORY;                                               \
  }                                                                            \
  *stream++ = c;                                                               \
  slen++;

// This is a little bit of a hack but works for now

#define END_FIELD()                                                            \
  self->stream_len = slen;                                                     \
  if (end_field(self) < 0) {                                                   \
    goto parsingerror;                                                         \
  }                                                                            \
  stream = self->stream + self->stream_len;                                    \
  slen = self->stream_len;

#define END_LINE_STATE(STATE)                                                  \
  self->stream_len = slen;                                                     \
  if (end_line(self) < 0) {                                                    \
    goto parsingerror;                                                         \
  }                                                                            \
  stream = self->stream + self->stream_len;                                    \
  slen = self->stream_len;                                                     \
  self->state = STATE;                                                         \
  if (line_limit > 0 && self->lines == start_lines + line_limit) {             \
    goto linelimit;                                                            \
  }

#define END_LINE_AND_FIELD_STATE(STATE)                                        \
  self->stream_len = slen;                                                     \
  if (end_line(self) < 0) {                                                    \
    goto parsingerror;                                                         \
  }                                                                            \
  if (end_field(self) < 0) {                                                   \
    goto parsingerror;                                                         \
  }                                                                            \
  stream = self->stream + self->stream_len;                                    \
  slen = self->stream_len;                                                     \
  self->state = STATE;                                                         \
  if (line_limit > 0 && self->lines == start_lines + line_limit) {             \
    goto linelimit;                                                            \
  }

#define END_LINE() END_LINE_STATE(START_RECORD)

#define IS_TERMINATOR(c) (c == lineterminator)

#define IS_QUOTE(c) ((c == self->quotechar && self->quoting != QUOTE_NONE))

// don't parse '\r' with a custom line terminator
#define IS_CARRIAGE(c) (has_carriage && c == carriage_symbol)

#define IS_COMMENT_CHAR(c) (has_comment && c == comment_symbol)

#define IS_ESCAPE_CHAR(c) (has_escape && c == escape_symbol)

#define IS_SKIPPABLE_SPACE(c)                                                  \
  ((!self->delim_whitespace && c == ' ' && self->skipinitialspace))

// applied when in a field
#define IS_DELIMITER(c)                                                        \
  ((!delim_whitespace && c == delimiter) || (delim_whitespace && isblank(c)))

#define _TOKEN_CLEANUP()                                                       \
  self->stream_len = slen;                                                     \
  self->datapos = i;

#define CHECK_FOR_BOM()                                                        \
  if (self->datalen - self->datapos >= 3 && *buf == '\xef' &&                  \
      *(buf + 1) == '\xbb' && *(buf + 2) == '\xbf') {                          \
    buf += 3;                                                                  \
    self->datapos += 3;                                                        \
  }

static int skip_this_line(parser_t *self, int64_t rownum) {
  if (self->skipfunc != NULL) {
    PyGILState_STATE state = PyGILState_Ensure();
    PyObject *result = PyObject_CallFunction(self->skipfunc, "i", rownum);

    // Error occurred. It will be processed
    // and caught at the Cython level.
    const int should_skip = result == NULL ? -1 : PyObject_IsTrue(result);

    Py_XDECREF(result);
    PyGILState_Release(state);

    return should_skip;
  } else if (self->skipset != NULL) {
    return (kh_get_int64((kh_int64_t *)self->skipset, self->file_lines) !=
            ((kh_int64_t *)self->skipset)->n_buckets);
  } else {
    return (rownum <= self->skip_first_N_rows);
  }
}

#ifdef PANDAS_HAVE_NEON

// Scan data for any special character using NEON.
// Returns the byte offset of the first special character,
// or the number of bytes scanned (all clean) if none found.
static inline size_t fast_scan_simd(const char *data, size_t len,
                                    uint8x16_t vdelim, uint8x16_t vterm,
                                    uint8x16_t vcr, uint8x16_t vquote,
                                    uint8x16_t vescape, uint8x16_t vcomment) {
  size_t i = 0;
  for (; i + 15 < len; i += 16) {
    uint8x16_t chunk = vld1q_u8((const uint8_t *)&data[i]);
    uint8x16_t m = vceqq_u8(chunk, vdelim);
    m = vorrq_u8(m, vceqq_u8(chunk, vterm));
    m = vorrq_u8(m, vceqq_u8(chunk, vcr));
    m = vorrq_u8(m, vceqq_u8(chunk, vquote));
    m = vorrq_u8(m, vceqq_u8(chunk, vescape));
    m = vorrq_u8(m, vceqq_u8(chunk, vcomment));

    // Each matching byte is 0xFF, non-matching is 0x00. NEON has no movemask,
    // so narrow each byte to a nibble (shift-right-narrow by 4): byte k becomes
    // nibble k of a single u64, 0xF if matched else 0x0. ctzll/4 is the offset.
    const uint8x8_t m8 = vshrn_n_u16(vreinterpretq_u16_u8(m), 4);
    const uint64_t matches = vget_lane_u64(vreinterpret_u64_u8(m8), 0);
    if (matches)
      return i + pandas_ctzll(matches) / 4;
  }
  return i;
}

// Lighter scan for quoted fields: only check quote and escape chars.
static inline size_t fast_scan_quoted_simd(const char *data, size_t len,
                                           uint8x16_t vquote,
                                           uint8x16_t vescape) {
  size_t i = 0;
  for (; i + 15 < len; i += 16) {
    uint8x16_t chunk = vld1q_u8((const uint8_t *)&data[i]);
    uint8x16_t m = vceqq_u8(chunk, vquote);
    m = vorrq_u8(m, vceqq_u8(chunk, vescape));

    const uint8x8_t m8 = vshrn_n_u16(vreinterpretq_u16_u8(m), 4);
    const uint64_t matches = vget_lane_u64(vreinterpret_u64_u8(m8), 0);
    if (matches)
      return i + pandas_ctzll(matches) / 4;
  }
  return i;
}

#elif defined(PANDAS_HAVE_SSE2)

static inline size_t fast_scan_simd(const char *data, size_t len,
                                    __m128i vdelim, __m128i vterm, __m128i vcr,
                                    __m128i vquote, __m128i vescape,
                                    __m128i vcomment) {
  size_t i = 0;
  for (; i + 15 < len; i += 16) {
    __m128i chunk = _mm_loadu_si128((const __m128i *)&data[i]);
    __m128i m = _mm_cmpeq_epi8(chunk, vdelim);
    m = _mm_or_si128(m, _mm_cmpeq_epi8(chunk, vterm));
    m = _mm_or_si128(m, _mm_cmpeq_epi8(chunk, vcr));
    m = _mm_or_si128(m, _mm_cmpeq_epi8(chunk, vquote));
    m = _mm_or_si128(m, _mm_cmpeq_epi8(chunk, vescape));
    m = _mm_or_si128(m, _mm_cmpeq_epi8(chunk, vcomment));

    int mask = _mm_movemask_epi8(m);
    if (mask)
      return i + pandas_ctz(mask);
  }
  return i;
}

static inline size_t fast_scan_quoted_simd(const char *data, size_t len,
                                           __m128i vquote, __m128i vescape) {
  size_t i = 0;
  for (; i + 15 < len; i += 16) {
    __m128i chunk = _mm_loadu_si128((const __m128i *)&data[i]);
    __m128i m = _mm_cmpeq_epi8(chunk, vquote);
    m = _mm_or_si128(m, _mm_cmpeq_epi8(chunk, vescape));

    int mask = _mm_movemask_epi8(m);
    if (mask)
      return i + pandas_ctz(mask);
  }
  return i;
}

#endif

// Per-byte match masks for the block fast lane. On NEON the mask carries
// 4 bits (one nibble) per input byte; on SSE2 one bit per byte. The
// BLOCK_MASK_* macros hide that difference.
#ifdef PANDAS_HAVE_NEON
typedef uint64_t block_mask_t;
static inline block_mask_t block_movemask(uint8x16_t m) {
  return vget_lane_u64(
      vreinterpret_u64_u8(vshrn_n_u16(vreinterpretq_u16_u8(m), 4)), 0);
}
#  define BLOCK_MASK_POS(mask) (pandas_ctzll(mask) >> 2)
#  define BLOCK_MASK_CLEAR(mask, p) ((mask) &= ~(0xFULL << ((p) * 4)))
#  define BLOCK_MASK_TEST(mask, p) (((mask) >> ((p) * 4)) & 0xFULL)
#elif defined(PANDAS_HAVE_SSE2)
typedef uint32_t block_mask_t;
#  define BLOCK_MASK_POS(mask) (pandas_ctz(mask))
#  define BLOCK_MASK_CLEAR(mask, p) ((mask) &= ~(1u << (p)))
#  define BLOCK_MASK_TEST(mask, p) (((mask) >> (p)) & 1u)
#endif

static int tokenize_bytes(parser_t *self, uint64_t line_limit,
                          uint64_t start_lines) {
  char *buf = self->data + self->datapos;

  const char lineterminator =
      (self->lineterminator == '\0') ? '\n' : self->lineterminator;

  const int delim_whitespace = self->delim_whitespace;
  const char delimiter = self->delimiter;

  const char carriage_symbol = '\r';
  const bool has_carriage = (self->lineterminator == '\0');
  const char comment_symbol = self->commentchar;
  const bool has_comment = (self->commentchar != '\0');
  const char escape_symbol = self->escapechar;
  const bool has_escape = (self->escapechar != '\0');
  const bool has_skip = (self->skipfunc != NULL || self->skipset != NULL ||
                         self->skip_first_N_rows >= 0);

#if defined(PANDAS_HAVE_NEON) || defined(PANDAS_HAVE_SSE2)
  // The block fast lane models plain fields separated by single-char
  // delimiters/terminators; configurations that transform field content at
  // tokenize time fall back to the state machine entirely. Quotes, escapes,
  // carriage returns and comment chars are handled per block: any block
  // containing one is left to the state machine.
  const bool lane_ok =
      !delim_whitespace && !has_skip && !self->skipinitialspace;
#endif

#ifdef PANDAS_HAVE_NEON
  const uint8x16_t vdelim = vdupq_n_u8((uint8_t)delimiter);
  const uint8x16_t vterm = vdupq_n_u8((uint8_t)lineterminator);
  const uint8x16_t vcr =
      has_carriage ? vdupq_n_u8((uint8_t)carriage_symbol) : vterm;
  const uint8x16_t vquote = (self->quoting != QUOTE_NONE)
                                ? vdupq_n_u8((uint8_t)self->quotechar)
                                : vterm;
  // Fall back to vquote (not vterm) so the quoted-field scan is not
  // broken by line terminators inside quoted fields.
  const uint8x16_t vescape = (self->escapechar != '\0')
                                 ? vdupq_n_u8((uint8_t)self->escapechar)
                                 : vquote;
  const uint8x16_t vcomment = (self->commentchar != '\0')
                                  ? vdupq_n_u8((uint8_t)self->commentchar)
                                  : vterm;
#elif defined(PANDAS_HAVE_SSE2)
  const __m128i vdelim = _mm_set1_epi8(delimiter);
  const __m128i vterm = _mm_set1_epi8(lineterminator);
  const __m128i vcr = has_carriage ? _mm_set1_epi8(carriage_symbol) : vterm;
  const __m128i vquote =
      (self->quoting != QUOTE_NONE) ? _mm_set1_epi8(self->quotechar) : vterm;
  // Fall back to vquote (not vterm) so the quoted-field scan is not
  // broken by line terminators inside quoted fields.
  const __m128i vescape =
      (self->escapechar != '\0') ? _mm_set1_epi8(self->escapechar) : vquote;
  const __m128i vcomment =
      (self->commentchar != '\0') ? _mm_set1_epi8(self->commentchar) : vterm;
#endif

  // Reserve worst-case stream space up front: the bulk-scan copies below
  // (scalar and SIMD) write input data 1:1 with no per-copy capacity check, so
  // the stream must already hold all remaining input. Field null-terminators
  // can exceed 1:1 but go through PUSH_CHAR/END_FIELD, which re-check capacity.
  if (make_stream_space(self, self->datalen - self->datapos) < 0) {
    const size_t bufsize = 100;
    self->error_msg = malloc(bufsize);
    snprintf(self->error_msg, bufsize, "out of memory");
    return -1;
  }

  char *stream = self->stream + self->stream_len;
  uint64_t slen = self->stream_len;

  // Lookup table marking characters that force a state-machine transition
  // during bulk scanning in IN_FIELD and IN_QUOTED_FIELD.
  // Bit 0 (0x1): breaks scan in an unquoted field.
  // Bit 1 (0x2): breaks scan in a quoted field.
  uint8_t breaks_field_scan[256] = {0};
  uint8_t index;

  memcpy(&index, &lineterminator, sizeof(lineterminator));
  breaks_field_scan[index] |= 0x1;
  if (has_carriage) {
    memcpy(&index, &carriage_symbol, sizeof(carriage_symbol));
    breaks_field_scan[index] |= 0x1;
  }
  if (has_escape) {
    memcpy(&index, &escape_symbol, sizeof(escape_symbol));
    breaks_field_scan[index] |= 0x1 | 0x2;
  }
  if (!delim_whitespace) {
    memcpy(&index, &delimiter, sizeof(delimiter));
    breaks_field_scan[index] |= 0x1;
  } else {
    // Mirrors IS_DELIMITER's use of isblank(), which matches ' ' and '\t'.
    char space = ' ', tab = '\t';
    memcpy(&index, &space, sizeof(space));
    breaks_field_scan[index] |= 0x1;
    memcpy(&index, &tab, sizeof(tab));
    breaks_field_scan[index] |= 0x1;
  }
  if (has_comment) {
    memcpy(&index, &comment_symbol, sizeof(comment_symbol));
    breaks_field_scan[index] |= 0x1;
  }
  if (self->quoting != QUOTE_NONE) {
    memcpy(&index, &self->quotechar, sizeof(self->quotechar));
    breaks_field_scan[index] |= 0x2;
  }

  if (self->file_lines == 0 && !self->preloaded) {
    /* In preloaded mode the buffer usually starts mid-file, where a BOM
     * byte sequence is real data; load_buffer strips a leading BOM itself
     * when the buffer really is the start of the file. */
    CHECK_FOR_BOM();
  }

  char c;
  int64_t i;
  for (i = self->datapos; i < self->datalen; ++i) {
#if defined(PANDAS_HAVE_NEON) || defined(PANDAS_HAVE_SSE2)
    // Block fast lane: while the next 16 bytes contain only ordinary
    // characters, delimiters and line terminators, copy them to the stream
    // in one shot and emit every field/line boundary from the SIMD match
    // mask, bypassing the per-character state machine. The block is copied
    // verbatim; end_field()'s NUL then lands exactly on the copied
    // delimiter/terminator byte, so the stream layout is identical to the
    // state machine's.
    if (lane_ok && (self->state == START_FIELD || self->state == IN_FIELD ||
                    self->state == START_RECORD)) {
      const int64_t lane_start_i = i;
      self->stream_len = slen;
      // True when the next byte to scan is the first byte of a line (and
      // nothing of it is buffered yet). Tracked across blocks so the
      // whitespace-line probe below runs once per line, not once per block;
      // line ends that fall inside a block are probed in the specials loop.
      bool at_line_edge = (self->line_fields[self->lines] == 0 &&
                           self->stream_len == (uint64_t)self->word_start);
      // The 16-byte block copy and the <=15-byte tail re-copy below are the
      // lane's only unchecked stream writes; requiring 32 bytes of slack per
      // block keeps them within stream_cap no matter how end_line moves
      // stream_len (checked writes inside end_field/end_line error out on
      // their own).
      while (i + 16 <= self->datalen &&
             self->stream_len + 32 <= self->stream_cap) {
        if (at_line_edge && (*buf == ' ' || *buf == '\t') &&
            *buf != delimiter) {
          // A line starting with a blank may be a whitespace-only line
          // (WHITESPACE_LINE handling); defer to the state machine.
          break;
        }
        at_line_edge = false;
#  ifdef PANDAS_HAVE_NEON
        const uint8x16_t chunk = vld1q_u8((const uint8_t *)buf);
        const uint8x16_t eq_delim = vceqq_u8(chunk, vdelim);
        const uint8x16_t eq_term = vceqq_u8(chunk, vterm);
        uint8x16_t dirty = vdupq_n_u8(0);
        if (self->quoting != QUOTE_NONE)
          dirty = vorrq_u8(dirty, vceqq_u8(chunk, vquote));
        if (has_escape)
          dirty = vorrq_u8(dirty, vceqq_u8(chunk, vescape));
        if (has_comment)
          dirty = vorrq_u8(dirty, vceqq_u8(chunk, vcomment));
        if (block_movemask(dirty))
          break;
        const block_mask_t crs =
            has_carriage ? block_movemask(vceqq_u8(chunk, vcr)) : 0;
        block_mask_t specials =
            block_movemask(vorrq_u8(eq_delim, eq_term)) | crs;
        const block_mask_t terms = block_movemask(eq_term);
#  else // PANDAS_HAVE_SSE2
        const __m128i chunk = _mm_loadu_si128((const __m128i *)buf);
        const __m128i eq_delim = _mm_cmpeq_epi8(chunk, vdelim);
        const __m128i eq_term = _mm_cmpeq_epi8(chunk, vterm);
        __m128i dirty = _mm_setzero_si128();
        if (self->quoting != QUOTE_NONE)
          dirty = _mm_or_si128(dirty, _mm_cmpeq_epi8(chunk, vquote));
        if (has_escape)
          dirty = _mm_or_si128(dirty, _mm_cmpeq_epi8(chunk, vescape));
        if (has_comment)
          dirty = _mm_or_si128(dirty, _mm_cmpeq_epi8(chunk, vcomment));
        if (_mm_movemask_epi8(dirty))
          break;
        const block_mask_t crs =
            has_carriage
                ? (block_mask_t)_mm_movemask_epi8(_mm_cmpeq_epi8(chunk, vcr))
                : 0;
        block_mask_t specials =
            (block_mask_t)_mm_movemask_epi8(_mm_or_si128(eq_delim, eq_term)) |
            crs;
        const block_mask_t terms = (block_mask_t)_mm_movemask_epi8(eq_term);
#  endif
        memcpy(self->stream + self->stream_len, buf, 16);
        if (!specials) {
          self->stream_len += 16;
          buf += 16;
          i += 16;
          continue;
        }
        uint64_t block_base = self->stream_len;
        while (specials) {
          const int p = (int)BLOCK_MASK_POS(specials);
          int is_term = BLOCK_MASK_TEST(terms, p) != 0;
          // Extra input bytes the terminator consumes (1 for the \n of a
          // \r\n pair; the field's NUL replaces the \r).
          int term_extra = 0;
          if (has_carriage && BLOCK_MASK_TEST(crs, p)) {
            if (p >= 15 || buf[p + 1] != lineterminator) {
              // \r\n split across blocks, or a lone \r: defer to the
              // state machine's EAT_CRNL handling.
              self->stream_len = block_base + (uint64_t)p;
              buf += p;
              i += p;
              goto block_lane_exit;
            }
            is_term = 1;
            term_extra = 1;
          }
          if (is_term && self->line_fields[self->lines] == 0 &&
              block_base + (uint64_t)p == (uint64_t)self->word_start) {
            // Line terminator with no line content: defer to the state
            // machine, which knows about skip_empty_lines.
            self->stream_len = block_base + (uint64_t)p;
            buf += p;
            i += p;
            goto block_lane_exit;
          }
          self->stream_len = block_base + (uint64_t)p;
          if (end_field(self) < 0) {
            slen = self->stream_len;
            i += p;
            goto parsingerror;
          }
          if (is_term) {
            if (end_line(self) < 0) {
              slen = self->stream_len;
              i += p;
              goto parsingerror;
            }
            if (line_limit > 0 && self->lines == start_lines + line_limit) {
              // Must run before the tail re-copy bail below: end_line has
              // already counted this line, and the state machine's ==-based
              // limit checks would never match again.
              slen = self->stream_len;
              // index of the terminator's last byte; the label increments
              // one past it
              i += p + term_extra;
              self->state = START_RECORD;
              goto linelimit;
            }
            if (term_extra ||
                self->stream_len != block_base + (uint64_t)(p + 1)) {
              // The write position no longer tracks the input position:
              // the terminator consumed an extra input byte (\r\n), or
              // end_line moved stream_len (synthetic trailing fields for a
              // short row, or a skipped bad row). Re-copy the unconsumed
              // tail of the block after the new position and rebase so
              // block_base + q stays correct for later specials.
              if (self->stream_len + 16 > self->stream_cap) {
                // Not enough slack for the unchecked tail re-copy; hand
                // the rest of the input back to the state machine, whose
                // writes are capacity-checked.
                buf += p + 1 + term_extra;
                i += p + 1 + term_extra;
                goto block_lane_exit;
              }
              memcpy(self->stream + self->stream_len, buf + p + 1 + term_extra,
                     (size_t)(15 - p - term_extra));
              block_base = self->stream_len - (uint64_t)(p + 1 + term_extra);
            }
            if (line_limit > 0 && self->lines == start_lines + line_limit) {
              slen = self->stream_len;
              // index of the terminator's last byte; the label increments
              // one past it
              i += p + term_extra;
              self->state = START_RECORD;
              goto linelimit;
            }
            if (p + 1 + term_extra < 16) {
              const char nxt = buf[p + 1 + term_extra];
              if ((nxt == ' ' || nxt == '\t') && nxt != delimiter) {
                // Next line starts with a blank: defer to the state
                // machine (see the whitespace-line check at block top).
                buf += p + 1 + term_extra;
                i += p + 1 + term_extra;
                goto block_lane_exit;
              }
            } else {
              // Line ends exactly at the block edge: probe the next
              // line's first byte at the next block top.
              at_line_edge = true;
            }
            if (term_extra) {
              BLOCK_MASK_CLEAR(specials, p + 1); // the \n of the pair
            }
          }
          BLOCK_MASK_CLEAR(specials, p);
        }
        self->stream_len = block_base + 16;
        buf += 16;
        i += 16;
      }
    block_lane_exit:
      slen = self->stream_len;
      stream = self->stream + slen;
      if (i != lane_start_i) {
        // Only reset the state when the lane consumed input; a no-progress
        // exit must preserve states like the post-WHITESPACE_LINE
        // START_FIELD, or the machine and the lane ping-pong forever.
        if (slen != (uint64_t)self->word_start) {
          self->state = IN_FIELD;
        } else if (self->line_fields[self->lines] > 0) {
          self->state = START_FIELD;
        } else {
          self->state = START_RECORD;
        }
      }
      if (i >= self->datalen) {
        break;
      }
    }
#endif
    // next character in file
    c = *buf++;

    switch (self->state) {
    case START_FIELD_IN_SKIP_LINE:
      if (IS_TERMINATOR(c)) {
        END_LINE();
      } else if (IS_CARRIAGE(c)) {
        self->file_lines++;
        self->state = EAT_CRNL_NOP;
      } else if (IS_QUOTE(c)) {
        self->state = IN_QUOTED_FIELD_IN_SKIP_LINE;
      } else if (IS_DELIMITER(c)) {
        // Do nothing, we're starting a new field again.
      } else {
        self->state = IN_FIELD_IN_SKIP_LINE;
      }
      break;

    case IN_FIELD_IN_SKIP_LINE:
      if (IS_TERMINATOR(c)) {
        END_LINE();
      } else if (IS_CARRIAGE(c)) {
        self->file_lines++;
        self->state = EAT_CRNL_NOP;
      } else if (IS_DELIMITER(c)) {
        self->state = START_FIELD_IN_SKIP_LINE;
      }
      break;

    case IN_QUOTED_FIELD_IN_SKIP_LINE:
      if (IS_QUOTE(c)) {
        if (self->doublequote) {
          self->state = QUOTE_IN_QUOTED_FIELD_IN_SKIP_LINE;
        } else {
          self->state = IN_FIELD_IN_SKIP_LINE;
        }
      }
      break;

    case QUOTE_IN_QUOTED_FIELD_IN_SKIP_LINE:
      if (IS_QUOTE(c)) {
        self->state = IN_QUOTED_FIELD_IN_SKIP_LINE;
      } else if (IS_TERMINATOR(c)) {
        END_LINE();
      } else if (IS_CARRIAGE(c)) {
        self->file_lines++;
        self->state = EAT_CRNL_NOP;
      } else if (IS_DELIMITER(c)) {
        self->state = START_FIELD_IN_SKIP_LINE;
      } else {
        self->state = IN_FIELD_IN_SKIP_LINE;
      }
      break;

    case WHITESPACE_LINE:
      if (IS_TERMINATOR(c)) {
        self->file_lines++;
        self->state = START_RECORD;
        break;
      } else if (IS_CARRIAGE(c)) {
        self->file_lines++;
        self->state = EAT_CRNL_NOP;
        break;
      } else if (!self->delim_whitespace) {
        if (isblank(c) && c != self->delimiter) {
        } else { // backtrack
          // use i + 1 because buf has been incremented but not i
          do {
            --buf;
            --i;
          } while (i + 1 > self->datapos && !IS_TERMINATOR(*buf) &&
                   !IS_CARRIAGE(*buf));

          // reached a newline/carriage return rather than the beginning
          if (IS_TERMINATOR(*buf) || IS_CARRIAGE(*buf)) {
            ++buf; // move pointer to first char after newline
            ++i;
          }
          self->state = START_FIELD;
        }
        break;
      }
      // fall through

    case EAT_WHITESPACE:
      if (IS_TERMINATOR(c)) {
        END_LINE();
        self->state = START_RECORD;
        break;
      } else if (IS_CARRIAGE(c)) {
        self->state = EAT_CRNL;
        break;
      } else if (IS_COMMENT_CHAR(c)) {
        self->state = EAT_COMMENT;
        break;
      } else if (!isblank(c)) {
        self->state = START_FIELD;
        PD_FALLTHROUGH; // fall through to subsequent state
      } else {
        // if whitespace char, keep slurping
        break;
      }

    case START_RECORD: {
      // start of record
      const int should_skip =
          has_skip ? skip_this_line(self, self->file_lines) : 0;

      if (should_skip == -1) {
        goto parsingerror;
      } else if (should_skip) {
        if (IS_QUOTE(c)) {
          self->state = IN_QUOTED_FIELD_IN_SKIP_LINE;
        } else {
          self->state = IN_FIELD_IN_SKIP_LINE;

          if (IS_TERMINATOR(c)) {
            END_LINE();
          }
        }
        break;
      } else if (IS_TERMINATOR(c)) {
        // \n\r possible?
        if (self->skip_empty_lines) {
          self->file_lines++;
        } else {
          END_LINE();
        }
        break;
      } else if (IS_CARRIAGE(c)) {
        if (self->skip_empty_lines) {
          self->file_lines++;
          self->state = EAT_CRNL_NOP;
        } else {
          self->state = EAT_CRNL;
        }
        break;
      } else if (IS_COMMENT_CHAR(c)) {
        self->state = EAT_LINE_COMMENT;
        break;
      } else if (isblank(c)) {
        if (self->delim_whitespace) {
          if (self->skip_empty_lines) {
            self->state = WHITESPACE_LINE;
          } else {
            self->state = EAT_WHITESPACE;
          }
          break;
        } else if (c != self->delimiter && self->skip_empty_lines) {
          self->state = WHITESPACE_LINE;
          break;
        }
      }

      // normal character - fall through
      // to handle as START_FIELD
      self->state = START_FIELD;
      PD_FALLTHROUGH;
    }
    case START_FIELD:
      // expecting field
      if (IS_TERMINATOR(c)) {
        END_FIELD();
        END_LINE();
      } else if (IS_CARRIAGE(c)) {
        END_FIELD();
        self->state = EAT_CRNL;
      } else if (IS_QUOTE(c)) {
        // start quoted field
        self->state = IN_QUOTED_FIELD;
      } else if (IS_ESCAPE_CHAR(c)) {
        // possible escaped character
        self->state = ESCAPED_CHAR;
      } else if (IS_SKIPPABLE_SPACE(c)) {
        // ignore space at start of field
      } else if (IS_DELIMITER(c)) {
        if (self->delim_whitespace) {
          self->state = EAT_WHITESPACE;
        } else {
          // save empty field
          END_FIELD();
        }
      } else if (IS_COMMENT_CHAR(c)) {
        END_FIELD();
        self->state = EAT_COMMENT;
      } else {
        // begin new unquoted field
        PUSH_CHAR(c);
        self->state = IN_FIELD;
      }
      break;

    case ESCAPED_CHAR:
      PUSH_CHAR(c);
      self->state = IN_FIELD;
      break;

    case EAT_LINE_COMMENT:
      if (IS_TERMINATOR(c)) {
        self->file_lines++;
        self->state = START_RECORD;
      } else if (IS_CARRIAGE(c)) {
        self->file_lines++;
        self->state = EAT_CRNL_NOP;
      }
      break;

    case IN_FIELD:
      // in unquoted field
      if (IS_TERMINATOR(c)) {
        END_FIELD();
        END_LINE();
      } else if (IS_CARRIAGE(c)) {
        END_FIELD();
        self->state = EAT_CRNL;
      } else if (IS_ESCAPE_CHAR(c)) {
        // possible escaped character
        self->state = ESCAPED_CHAR;
      } else if (IS_DELIMITER(c)) {
        // end of field - end of line not reached yet
        END_FIELD();

        if (self->delim_whitespace) {
          self->state = EAT_WHITESPACE;
        } else {
          self->state = START_FIELD;
        }
      } else if (IS_COMMENT_CHAR(c)) {
        END_FIELD();
        self->state = EAT_COMMENT;
      } else {
        // normal character - save in field
        PUSH_CHAR(c);

        // SIMD bulk scan: process 16 bytes at a time, copying
        // normal characters directly without state-machine overhead.
#if defined(PANDAS_HAVE_NEON) || defined(PANDAS_HAVE_SSE2)
        if (!self->delim_whitespace) {
          size_t remaining = self->datalen - (i + 1);
          if (remaining >= 16) {
            size_t skip = fast_scan_simd(buf, remaining, vdelim, vterm, vcr,
                                         vquote, vescape, vcomment);
            if (skip > 0) {
              // in-bounds; see stream reservation at loop start
              memcpy(stream, buf, skip);
              stream += skip;
              slen += skip;
              buf += skip;
              i += skip;
            }
          }
        }
#endif
        // Scalar bulk scan fallback: copy remaining ordinary characters
        // directly, bypassing the per-char state machine overhead.
        while (i + 1 < self->datalen &&
               !(breaks_field_scan[(uint8_t)*buf] & 0x1)) {
          *stream++ = *buf++;
          slen++;
          i++;
        }
      }
      break;

    case IN_QUOTED_FIELD:
      // in quoted field
      if (IS_ESCAPE_CHAR(c)) {
        // possible escape character
        self->state = ESCAPE_IN_QUOTED_FIELD;
      } else if (IS_QUOTE(c)) {
        if (self->doublequote) {
          // double quote - " represented by ""
          self->state = QUOTE_IN_QUOTED_FIELD;
        } else {
          // end of quote part of field
          self->state = IN_FIELD;
        }
      } else {
        // normal character - save in field
        PUSH_CHAR(c);

        // SIMD bulk scan for quoted fields: only quote and escape
        // chars are special, so use a lighter scan.
#if defined(PANDAS_HAVE_NEON) || defined(PANDAS_HAVE_SSE2)
        {
          size_t remaining = self->datalen - (i + 1);
          if (remaining >= 16) {
            size_t skip =
                fast_scan_quoted_simd(buf, remaining, vquote, vescape);
            if (skip > 0) {
              memcpy(stream, buf, skip);
              stream += skip;
              slen += skip;
              buf += skip;
              i += skip;
            }
          }
        }
#endif
        // Scalar bulk scan fallback: copy remaining ordinary characters
        // directly, bypassing the per-char state machine overhead.
        while (i + 1 < self->datalen &&
               !(breaks_field_scan[(uint8_t)*buf] & 0x2)) {
          *stream++ = *buf++;
          slen++;
          i++;
        }
      }
      break;

    case ESCAPE_IN_QUOTED_FIELD:
      PUSH_CHAR(c);
      self->state = IN_QUOTED_FIELD;
      break;

    case QUOTE_IN_QUOTED_FIELD:
      // double quote - seen a quote in a quoted field
      if (IS_QUOTE(c)) {
        // save "" as "

        PUSH_CHAR(c);
        self->state = IN_QUOTED_FIELD;
      } else if (IS_DELIMITER(c)) {
        // end of field - end of line not reached yet
        END_FIELD();

        if (self->delim_whitespace) {
          self->state = EAT_WHITESPACE;
        } else {
          self->state = START_FIELD;
        }
      } else if (IS_TERMINATOR(c)) {
        END_FIELD();
        END_LINE();
      } else if (IS_CARRIAGE(c)) {
        END_FIELD();
        self->state = EAT_CRNL;
      } else {
        PUSH_CHAR(c);
        self->state = IN_FIELD;
      }
      break;

    case EAT_COMMENT:
      if (IS_TERMINATOR(c)) {
        END_LINE();
      } else if (IS_CARRIAGE(c)) {
        self->state = EAT_CRNL;
      }
      break;

    // only occurs with non-custom line terminator,
    // which is why we directly check for '\n'
    case EAT_CRNL:
      if (c == '\n') {
        END_LINE();
      } else if (IS_DELIMITER(c)) {
        if (self->delim_whitespace) {
          END_LINE_STATE(EAT_WHITESPACE);
        } else {
          // Handle \r-delimited files
          END_LINE_AND_FIELD_STATE(START_FIELD);
        }
      } else {
        if (self->delim_whitespace) {
          /* XXX
           * first character of a new record--need to back up and
           * reread
           * to handle properly...
           */
          i--;
          buf--; // back up one character (HACK!)
          END_LINE_STATE(START_RECORD);
        } else {
          // \r line terminator
          // UGH. we don't actually want
          // to consume the token. fix this later
          self->stream_len = slen;
          if (end_line(self) < 0) {
            goto parsingerror;
          }

          stream = self->stream + self->stream_len;
          slen = self->stream_len;
          self->state = START_RECORD;

          --i;
          buf--; // let's try this character again (HACK!)
          if (line_limit > 0 && self->lines == start_lines + line_limit) {
            goto linelimit;
          }
        }
      }
      break;

    // only occurs with non-custom line terminator,
    // which is why we directly check for '\n'
    case EAT_CRNL_NOP: // inside an ignored comment line
      self->state = START_RECORD;
      // \r line terminator -- parse this character again
      if (c != '\n' && !IS_DELIMITER(c)) {
        --i;
        --buf;
      }
      break;
    default:
      break;
    }
  }

  _TOKEN_CLEANUP();

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

static int parser_handle_eof(parser_t *self) {
  const size_t bufsize = 100;

  if (self->datalen != 0)
    return -1;

  switch (self->state) {
  case START_RECORD:
  case WHITESPACE_LINE:
  case EAT_CRNL_NOP:
  case EAT_LINE_COMMENT:
    return 0;

  case ESCAPE_IN_QUOTED_FIELD:
  case IN_QUOTED_FIELD:
    self->error_msg = (char *)malloc(bufsize);
    snprintf(self->error_msg, bufsize,
             "EOF inside string starting at row %" PRIu64, self->file_lines);
    return -1;

  case ESCAPED_CHAR:
    self->error_msg = (char *)malloc(bufsize);
    snprintf(self->error_msg, bufsize, "EOF following escape character");
    return -1;

  case IN_FIELD:
  case START_FIELD:
  case QUOTE_IN_QUOTED_FIELD:
    if (end_field(self) < 0)
      return -1;
    break;

  default:
    break;
  }

  if (end_line(self) < 0)
    return -1;
  else
    return 0;
}

int parser_consume_rows(parser_t *self, uint64_t nrows) {
  if (nrows > self->lines) {
    nrows = self->lines;
  }

  /* do nothing */
  if (nrows == 0)
    return 0;

  /* cannot guarantee that nrows + 1 has been observed */
  const int64_t word_deletions =
      self->line_start[nrows - 1] + self->line_fields[nrows - 1];

  uint64_t char_count;
  if (word_deletions < 1) {
    /* nothing deleted, so no data needs to be skipped */
    char_count = 0;
  } else if ((uint64_t)word_deletions < self->words_len) {
    /* start of the first surviving word, which equals the end (past the
     * trailing '\0') of the last deleted word */
    char_count = (uint64_t)self->word_starts[word_deletions];
  } else {
    /* every word is being deleted */
    char_count = self->stream_len;
  }

  /* move stream, only if something to move */
  if (char_count < self->stream_len) {
    memmove(self->stream, (self->stream + char_count),
            self->stream_len - char_count);
  }
  /* buffer counts */
  self->stream_len -= char_count;

  /* move token metadata */
  // Note: We should always have words_len < word_deletions, so this
  //  subtraction will remain appropriately-typed.
  int64_t offset;
  for (uint64_t i = 0; i < self->words_len - word_deletions; ++i) {
    offset = i + word_deletions;

    self->words[i] = self->words[offset] - char_count;
    self->word_starts[i] = self->word_starts[offset] - char_count;
  }
  self->words_len -= word_deletions;

  /* move current word pointer to stream */
  self->pword_start -= char_count;
  self->word_start -= char_count;

  /* move line metadata */
  // Note: We should always have self->lines - nrows + 1 >= 0, so this
  //  subtraction will remain appropriately-typed.
  for (uint64_t i = 0; i < self->lines - nrows + 1; ++i) {
    offset = i + nrows;
    self->line_start[i] = self->line_start[offset] - word_deletions;
    self->line_fields[i] = self->line_fields[offset];
  }
  self->lines -= nrows;

  return 0;
}

static size_t _next_pow2(size_t sz) {
  size_t result = 1;
  while (result < sz)
    result *= 2;
  return result;
}

int parser_trim_buffers(parser_t *self) {
  /*
    Free memory
   */

  /**
   * Before we free up space and trim, we should
   * save how many words we saw when parsing, if
   * it exceeds the maximum number we saw before.
   *
   * This is important for when we read in chunks,
   * so that we can inform subsequent chunk parsing
   * as to how many words we could possibly see.
   */
  if (self->words_cap > self->max_words_cap) {
    self->max_words_cap = self->words_cap;
  }

  /* trim words, word_starts */
  size_t new_cap = _next_pow2(self->words_len) + 1;
  if (new_cap < self->words_cap) {
    self->words = realloc(self->words, new_cap * sizeof(char *));
    if (self->words == NULL) {
      return PARSER_OUT_OF_MEMORY;
    }
    self->word_starts = realloc(self->word_starts, new_cap * sizeof(int64_t));
    if (self->word_starts == NULL) {
      return PARSER_OUT_OF_MEMORY;
    }
    self->words_cap = new_cap;
  }

  /* trim stream */
  new_cap = _next_pow2(self->stream_len) + 1;
  if (new_cap < self->stream_cap) {
    void *newptr = realloc(self->stream, new_cap);
    if (newptr == NULL) {
      return PARSER_OUT_OF_MEMORY;
    } else {
      // Update the pointers in the self->words array (char **) if
      // `realloc`
      //  moved the `self->stream` buffer. This block mirrors a similar
      //  block in
      //  `make_stream_space`.
      if (self->stream != newptr) {
        self->pword_start = (char *)newptr + self->word_start;

        for (uint64_t i = 0; i < self->words_len; ++i) {
          self->words[i] = (char *)newptr + self->word_starts[i];
        }
      }

      self->stream = newptr;
      self->stream_cap = new_cap;
    }
  }

  /* trim line_start, line_fields */
  new_cap = _next_pow2(self->lines) + 1;
  if (new_cap < self->lines_cap) {
    void *newptr = realloc(self->line_start, new_cap * sizeof(int64_t));
    if (newptr == NULL) {
      return PARSER_OUT_OF_MEMORY;
    } else {
      self->line_start = newptr;
    }
    newptr = realloc(self->line_fields, new_cap * sizeof(int64_t));
    if (newptr == NULL) {
      return PARSER_OUT_OF_MEMORY;
    } else {
      self->line_fields = newptr;
      self->lines_cap = new_cap;
    }
  }

  return 0;
}

/*
  nrows : number of rows to tokenize (or until reach EOF)
  all : tokenize all the data vs. certain number of rows
 */

static int _tokenize_helper(parser_t *self, uint64_t nrows, int all,
                            const char *encoding_errors) {
  int status = 0;
  const uint64_t start_lines = self->lines;

  if (self->state == FINISHED) {
    return 0;
  }

  while (1) {
    if (!all && self->lines - start_lines >= nrows)
      break;

    if (self->datapos == self->datalen) {
      if (self->source == NULL) {
        /* Pre-loaded buffer mode: the entire chunk has been consumed.
         * Signal EOF without invoking the Python I/O callback so the
         * GIL is never re-acquired during tokenisation. */
        self->datalen = 0;
        status = parser_handle_eof(self);
        self->state = FINISHED;
        break;
      }

      status = parser_buffer_bytes(self, self->chunksize, encoding_errors);

      if (status == REACHED_EOF) {
        // close out last line
        status = parser_handle_eof(self);
        self->state = FINISHED;
        break;
      } else if (status != 0) {
        return status;
      }
    }

    status = tokenize_bytes(self, nrows, start_lines);

    if (status < 0) {
      // XXX
      status = -1;
      break;
    }
  }
  return status;
}

int tokenize_nrows(parser_t *self, uint64_t nrows,
                   const char *encoding_errors) {
  return _tokenize_helper(self, nrows, 0, encoding_errors);
}

int tokenize_all_rows(parser_t *self, const char *encoding_errors) {
  return _tokenize_helper(self, 0, 1, encoding_errors);
}

/*
 * Function: to_boolean
 * --------------------
 *
 * Validate if item should be recognized as a boolean field.
 *
 * item: const char* representing parsed text
 * val : pointer to a uint8_t of boolean representation
 *
 * If item is determined to be boolean, this method will set
 * the appropriate value of val and return 0. A non-zero exit
 * status means that item was not inferred to be boolean, and
 * leaves the value of *val unmodified.
 */
int to_boolean(const char *item, uint8_t *val) {
  if (strcasecmp(item, "TRUE") == 0) {
    *val = 1;
    return 0;
  } else if (strcasecmp(item, "FALSE") == 0) {
    *val = 0;
    return 0;
  }

  return -1;
}

// ---------------------------------------------------------------------------
// Implementation of xstrtod

//
// strtod.c
//
// Convert string to double
//
// Copyright (C) 2002 Michael Ringgaard. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. Neither the name of the project nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.
//
// -----------------------------------------------------------------------

// Defined in fast_float_strtod.cpp — provides IEEE 754 correctly-rounded
// float parsing via the fast_float library.
int fast_float_strtod(const char *start, const char *end, double *value,
                      const char **endptr, char decimal);

double precise_xstrtod(const char *str, char **endptr, char decimal, char sci,
                       char tsep, int skip_trailing, int *error,
                       int *maybe_int) {
  // Use fast_float for standard format (no tsep, sci='e'/'E').
  // fast_float provides IEEE 754 correctly-rounded parsing.
  if (tsep == '\0' && (sci == 'e' || sci == 'E')) {
    const char *p = str;
    while (isspace_ascii(*p))
      p++;

    // Only try fast_float for numeric-looking input (digit, sign+digit,
    // or decimal point). This avoids fast_float parsing "nan"/"inf" which
    // the original precise_xstrtod did not handle.
    const char *q = p;
    if (*q == '-' || *q == '+')
      q++;
    if (!isdigit_ascii(*q) && !(*q == decimal && isdigit_ascii(*(q + 1))))
      goto fallback;

    // Find end of token (next whitespace or NUL).
    const char *end = p;
    while (*end && !isspace_ascii(*end))
      end++;

    double value;
    const char *parsed_end;
    if (fast_float_strtod(p, end, &value, &parsed_end, decimal) == 0) {
      // Determine maybe_int by checking if we saw decimal or 'e'/'E'.
      if (maybe_int != NULL) {
        *maybe_int = 1;
        for (const char *c = p; c < parsed_end; c++) {
          if (*c == decimal || *c == 'e' || *c == 'E') {
            *maybe_int = 0;
            break;
          }
        }
      }
      if (skip_trailing)
        while (isspace_ascii(*parsed_end))
          parsed_end++;
      if (endptr)
        *endptr = (char *)parsed_end;
      return value;
    }
  }

fallback:
    // Fallback for non-standard formats (custom tsep, or sci char).
    ;
  const char *p = str;
  const int max_digits = 17;

  if (maybe_int != NULL)
    *maybe_int = 1;
  // Cache powers of 10 in memory.
  static double e[] = {
      1.,    1e1,   1e2,   1e3,   1e4,   1e5,   1e6,   1e7,   1e8,   1e9,
      1e10,  1e11,  1e12,  1e13,  1e14,  1e15,  1e16,  1e17,  1e18,  1e19,
      1e20,  1e21,  1e22,  1e23,  1e24,  1e25,  1e26,  1e27,  1e28,  1e29,
      1e30,  1e31,  1e32,  1e33,  1e34,  1e35,  1e36,  1e37,  1e38,  1e39,
      1e40,  1e41,  1e42,  1e43,  1e44,  1e45,  1e46,  1e47,  1e48,  1e49,
      1e50,  1e51,  1e52,  1e53,  1e54,  1e55,  1e56,  1e57,  1e58,  1e59,
      1e60,  1e61,  1e62,  1e63,  1e64,  1e65,  1e66,  1e67,  1e68,  1e69,
      1e70,  1e71,  1e72,  1e73,  1e74,  1e75,  1e76,  1e77,  1e78,  1e79,
      1e80,  1e81,  1e82,  1e83,  1e84,  1e85,  1e86,  1e87,  1e88,  1e89,
      1e90,  1e91,  1e92,  1e93,  1e94,  1e95,  1e96,  1e97,  1e98,  1e99,
      1e100, 1e101, 1e102, 1e103, 1e104, 1e105, 1e106, 1e107, 1e108, 1e109,
      1e110, 1e111, 1e112, 1e113, 1e114, 1e115, 1e116, 1e117, 1e118, 1e119,
      1e120, 1e121, 1e122, 1e123, 1e124, 1e125, 1e126, 1e127, 1e128, 1e129,
      1e130, 1e131, 1e132, 1e133, 1e134, 1e135, 1e136, 1e137, 1e138, 1e139,
      1e140, 1e141, 1e142, 1e143, 1e144, 1e145, 1e146, 1e147, 1e148, 1e149,
      1e150, 1e151, 1e152, 1e153, 1e154, 1e155, 1e156, 1e157, 1e158, 1e159,
      1e160, 1e161, 1e162, 1e163, 1e164, 1e165, 1e166, 1e167, 1e168, 1e169,
      1e170, 1e171, 1e172, 1e173, 1e174, 1e175, 1e176, 1e177, 1e178, 1e179,
      1e180, 1e181, 1e182, 1e183, 1e184, 1e185, 1e186, 1e187, 1e188, 1e189,
      1e190, 1e191, 1e192, 1e193, 1e194, 1e195, 1e196, 1e197, 1e198, 1e199,
      1e200, 1e201, 1e202, 1e203, 1e204, 1e205, 1e206, 1e207, 1e208, 1e209,
      1e210, 1e211, 1e212, 1e213, 1e214, 1e215, 1e216, 1e217, 1e218, 1e219,
      1e220, 1e221, 1e222, 1e223, 1e224, 1e225, 1e226, 1e227, 1e228, 1e229,
      1e230, 1e231, 1e232, 1e233, 1e234, 1e235, 1e236, 1e237, 1e238, 1e239,
      1e240, 1e241, 1e242, 1e243, 1e244, 1e245, 1e246, 1e247, 1e248, 1e249,
      1e250, 1e251, 1e252, 1e253, 1e254, 1e255, 1e256, 1e257, 1e258, 1e259,
      1e260, 1e261, 1e262, 1e263, 1e264, 1e265, 1e266, 1e267, 1e268, 1e269,
      1e270, 1e271, 1e272, 1e273, 1e274, 1e275, 1e276, 1e277, 1e278, 1e279,
      1e280, 1e281, 1e282, 1e283, 1e284, 1e285, 1e286, 1e287, 1e288, 1e289,
      1e290, 1e291, 1e292, 1e293, 1e294, 1e295, 1e296, 1e297, 1e298, 1e299,
      1e300, 1e301, 1e302, 1e303, 1e304, 1e305, 1e306, 1e307, 1e308};

  // Skip leading whitespace.
  while (isspace_ascii(*p))
    p++;

  // Handle optional sign.
  int negative = 0;
  switch (*p) {
  case '-':
    negative = 1;
    PD_FALLTHROUGH; // Fall through to increment position.
  case '+':
    p++;
    break;
  }

  long int exponent = 0;
  long int num_digits = 0;
  long int num_decimals = 0;

  // Accumulate mantissa digits as an integer to avoid per-digit FP rounding.
  // max_digits=17 decimal digits fit safely in uint64_t (max ~9.9e17 < 2^64).
  uint64_t mantissa = 0;

  // Process string of digits.
  while (isdigit_ascii(*p)) {
    if (num_digits < max_digits) {
      mantissa = mantissa * 10 + (*p - '0');
      num_digits++;
    } else {
      ++exponent;
    }

    p++;
    p += (tsep != '\0' && *p == tsep);
  }

  // Process decimal part
  if (*p == decimal) {
    if (maybe_int != NULL)
      *maybe_int = 0;
    p++;

    while (num_digits < max_digits && isdigit_ascii(*p)) {
      mantissa = mantissa * 10 + (*p - '0');
      p++;
      num_digits++;
      num_decimals++;
    }

    if (num_digits >= max_digits) // Consume extra decimal digits.
      while (isdigit_ascii(*p))
        ++p;

    exponent -= num_decimals;
  }

  if (num_digits == 0) {
    *error = ERANGE;
    return 0.0;
  }

  // Single conversion from integer mantissa to double: at most one rounding,
  // compared to up to max_digits roundings in the old FP accumulation loop.
  double number = (double)mantissa;

  // Correct for sign.
  if (negative)
    number = -number;

  // Process an exponent string.
  if (toupper_ascii(*p) == toupper_ascii(sci)) {
    if (maybe_int != NULL)
      *maybe_int = 0;

    // move past scientific notation
    p++;

    char *tmp_ptr;
    long int n = strtol(p, &tmp_ptr, 10);

    if (errno == ERANGE || checked_add(exponent, n, &exponent)) {
      errno = 0;
      exponent = n;
    }

    // If no digits after the 'e'/'E', un-consume it.
    if (tmp_ptr == p)
      p--;
    else
      p = tmp_ptr;
  }

  if (exponent > 308) {
    number = number == 0 ? 0 : number < 0 ? -HUGE_VAL : HUGE_VAL;
  } else if (exponent > 0) {
    number *= e[exponent];
  } else if (exponent < -308) { // Subnormal
    if (exponent < -616) {      // Prevent invalid array access.
      number = 0.;
    } else {
      number /= e[-308 - exponent];
      number /= e[308];
    }

  } else {
    number /= e[-exponent];
  }

  if (skip_trailing) {
    // Skip trailing whitespace.
    while (isspace_ascii(*p))
      p++;
  }

  if (endptr)
    *endptr = (char *)p;
  return number;
}

// End of xstrtod code
// ---------------------------------------------------------------------------

void uint_state_init(uint_state *self) {
  self->seen_sint = 0;
  self->seen_uint = 0;
  self->seen_null = 0;
}

int uint64_conflict(uint_state *self) {
  return self->seen_uint && (self->seen_sint || self->seen_null);
}

/* Copy a numeric token without `tsep` into `output`.
 *
 * Returns the number of bytes written (excluding NUL), or -1 on overflow.
 * `*endptr` is set to the position in `str` where copying stopped (the first
 * non-digit / non-tsep char or the end of input) so the caller can detect
 * trailing garbage like "1 ," (GH#64631).
 */
static int copy_number_without_tsep(char output[PROCESSED_WORD_CAPACITY],
                                    const char *str, const char **endptr,
                                    size_t str_len, char tsep) {
  const char *p = str;
  const char *end = str + str_len;
  size_t bytes_written = 0;

  if (p < end && (*p == '+' || *p == '-')) {
    output[bytes_written++] = *p++;
  }

  while (p < end && (isdigit_ascii(*p) || (tsep != '\0' && *p == tsep))) {
    if (*p != tsep) {
      if (bytes_written + 1 >= PROCESSED_WORD_CAPACITY) {
        return -1;
      }
      output[bytes_written++] = *p;
    }
    p++;
  }

  output[bytes_written] = '\0';
  if (endptr != NULL) {
    *endptr = p;
  }
  return (int)bytes_written;
}

int64_t str_to_int64(const char *p_item, int64_t length, int *error,
                     char tsep) {
  const char *p = p_item;
  // Skip leading spaces.
  while (isspace_ascii(*p)) {
    ++p;
  }

  // Handle sign. std::from_chars accepts '-' but rejects '+', so strip '+'
  // after verifying the following char is a digit (not another sign).
  const bool has_sign = *p == '-' || *p == '+';
  const char *digit_start = has_sign ? p + 1 : p;
  if (!isdigit_ascii(*digit_start)) {
    *error = ERROR_NO_DIGITS;
    return 0;
  }
  if (*p == '+') {
    ++p;
  }

  char buffer[PROCESSED_WORD_CAPACITY];
  // length == strlen(p_item) supplied by the caller (-1 to compute here);
  // lets from_chars get its end pointer without a strlen scan.
  size_t str_len =
      length < 0 ? strlen(p) : (size_t)length - (size_t)(p - p_item);
  const char *number_end = NULL;
  if (tsep != '\0' && memchr(p, tsep, str_len) != NULL) {
    const int written =
        copy_number_without_tsep(buffer, p, &number_end, str_len, tsep);
    if (written < 0) {
      // Word is too big, probably will cause an overflow
      *error = ERROR_OVERFLOW;
      return 0;
    }
    p = buffer;
    str_len = (size_t)written;
  }

  int64_t number;
  const char *endptr;
  const pd_strtoi_status status = pd_strtoll(p, p + str_len, &number, &endptr);
  if (number_end != NULL) {
    // GH#64631: detect trailing junk in the original input that
    // copy_number_without_tsep stopped at (e.g. "1 ," with tsep=',').
    endptr = number_end;
  }
  if (status == PD_STRTOI_OVERFLOW) {
    // Overflow with trailing junk → INVALID_CHARS (so caller can fall through
    // to float parsing, e.g. "18446744073709551616.0"). Pure overflow (endptr
    // at NUL) → OVERFLOW (caller retries as uint64).
    *error = *endptr ? ERROR_INVALID_CHARS : ERROR_OVERFLOW;
    return 0;
  }
  if (status == PD_STRTOI_INVALID) {
    *error = ERROR_INVALID_CHARS;
    return 0;
  }

  // Skip trailing spaces.
  while (isspace_ascii(*endptr)) {
    ++endptr;
  }

  // Did we use up all the characters?
  if (*endptr) {
    *error = ERROR_INVALID_CHARS;
    return 0;
  }

  *error = 0;
  return number;
}

uint64_t str_to_uint64(uint_state *state, const char *p_item, int64_t length,
                       int *error, char tsep) {
  const char *p = p_item;
  // Skip leading spaces.
  while (isspace_ascii(*p)) {
    ++p;
  }

  // Handle sign.
  if (*p == '-') {
    state->seen_sint = 1;
    *error = 0;
    return 0;
  } else if (*p == '+') {
    p++;
  }

  // Check that there is a first digit.
  if (!isdigit_ascii(*p)) {
    *error = ERROR_NO_DIGITS;
    return 0;
  }

  char buffer[PROCESSED_WORD_CAPACITY];
  // length == strlen(p_item) supplied by the caller (-1 to compute here);
  // lets from_chars get its end pointer without a strlen scan.
  size_t str_len =
      length < 0 ? strlen(p) : (size_t)length - (size_t)(p - p_item);
  const char *number_end = NULL;
  if (tsep != '\0' && memchr(p, tsep, str_len) != NULL) {
    const int written =
        copy_number_without_tsep(buffer, p, &number_end, str_len, tsep);
    if (written < 0) {
      // Word is too big, probably will cause an overflow
      *error = ERROR_OVERFLOW;
      return 0;
    }
    p = buffer;
    str_len = (size_t)written;
  }

  uint64_t number;
  const char *endptr;
  const pd_strtoi_status status = pd_strtoull(p, p + str_len, &number, &endptr);
  if (number_end != NULL) {
    // GH#64631: detect trailing junk in the original input.
    endptr = number_end;
  }
  if (status == PD_STRTOI_OVERFLOW) {
    *error = *endptr ? ERROR_INVALID_CHARS : ERROR_OVERFLOW;
    return 0;
  }
  if (status == PD_STRTOI_INVALID) {
    *error = ERROR_INVALID_CHARS;
    return 0;
  }

  // Skip trailing spaces.
  while (isspace_ascii(*endptr)) {
    ++endptr;
  }

  // Did we use up all the characters?
  if (*endptr) {
    *error = ERROR_INVALID_CHARS;
    return 0;
  }

  if (number > (uint64_t)INT64_MAX) {
    state->seen_uint = 1;
  }

  *error = 0;
  return number;
}
