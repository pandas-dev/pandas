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


#include "tokenizer.h"

#include <ctype.h>
#include <math.h>
#include <float.h>


#define READ_ERROR_OUT_OF_MEMORY   1


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




static void *safe_realloc(void *buffer, size_t size) {
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


static  void free_if_not_null(void *ptr) {
     if (ptr != NULL) free(ptr);
 }



 /*

   Parser / tokenizer

 */


static void *grow_buffer(void *buffer, int length, int *capacity,
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
    self->strict = 0;

    self->expected_fields = -1;
    self->error_bad_lines = 0;
    self->warn_bad_lines = 0;

    self->commentchar = '#';
    self->thousands = '\0';

    self->skipset = NULL;
    self-> skip_first_N_rows = -1;
    self->skip_footer = 0;
}

int get_parser_memory_footprint(parser_t *self) {
    return 0;
}

parser_t* parser_new() {
    return (parser_t*) calloc(1, sizeof(parser_t));
}



int parser_clear_data_buffers(parser_t *self) {
    free_if_not_null(self->stream);
    free_if_not_null(self->words);
    free_if_not_null(self->word_starts);
    free_if_not_null(self->line_start);
    free_if_not_null(self->line_fields);

    return 0;
}

int parser_cleanup(parser_t *self) {
    if (self->cb_cleanup == NULL) {
        return 0;
    }

    if (self->cb_cleanup(self->source) < 0) {
        return -1;
    }

    if (parser_clear_data_buffers(self) < 0) {
        return -1;
    }

    // XXX where to put this
    free_if_not_null(self->error_msg);
    free_if_not_null(self->warn_msg);

    if (self->skipset != NULL)
        kh_destroy_int64((kh_int64_t*) self->skipset);

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
    self->warn_msg = NULL;

    self->commentchar = '\0';

    return 0;
}


void parser_free(parser_t *self) {
    // opposite of parser_init
    parser_cleanup(self);
    free(self);
}

static int make_stream_space(parser_t *self, size_t nbytes) {
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
    /*
    printf("Line_start: ");

    for (j = 0; j < self->lines + 1; ++j) {
         printf("%d ", self->line_fields[j]);
     }
    printf("\n");

    printf("lines_cap: %d\n", self->lines_cap);
    */
    cap = self->lines_cap;
    self->line_start = (int*) grow_buffer((void *) self->line_start,
                                          self->lines + 1,
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


static int push_char(parser_t *self, char c) {
    /* TRACE(("pushing %c \n", c)) */
    self->stream[self->stream_len++] = c;
    return 0;
}

static int P_INLINE end_field(parser_t *self) {
    // XXX cruft
    self->numeric_field = 0;

    // null terminate token
    push_char(self, '\0');

    // set pointer and metadata
    self->words[self->words_len] = self->pword_start;

    TRACE(("Char diff: %d\n", self->pword_start - self->words[0]));

    TRACE(("Saw word %s at: %d. Total: %d\n",
           self->pword_start, self->word_start, self->words_len + 1))

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
    int ex_length;
    int length = strlen(msg);

    if (self->warn_msg == NULL) {
        self->warn_msg = (char*) malloc(length + 1);
        strcpy(self->warn_msg, msg);
    } else {
        ex_length = strlen(self->warn_msg);
        self->warn_msg = (char*) safe_realloc(self->warn_msg,
                                              ex_length + length + 1);
        strcpy(self->warn_msg + ex_length, msg);
    }
}

static int end_line(parser_t *self) {
    int fields;
    khiter_t k;  /* for hash set detection */
    int ex_fields = self->expected_fields;
    char *msg;

    fields = self->line_fields[self->lines];

    TRACE(("Line end, nfields: %d\n", fields));

    if (self->lines > 0) {
        if (self->expected_fields >= 0) {
            ex_fields = self->expected_fields;
        } else {
            ex_fields = self->line_fields[self->lines - 1];
        }
    }

    if (self->state == SKIP_LINE) {
        TRACE(("Skipping row %d\n", self->file_lines));
        // increment file line count
        self->file_lines++;
        
        // skip the tokens from this bad line
        self->line_start[self->lines] += fields;

        // reset field count
        self->line_fields[self->lines] = 0;
        return 0;
    }

    /* printf("Line: %d, Fields: %d, Ex-fields: %d\n", self->lines, fields, ex_fields); */

    if (!(self->lines <= self->header_end + 1)
        && (self->expected_fields < 0 && fields > ex_fields)) {
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
                // pass up error message
                msg = (char*) malloc(100);
                sprintf(msg, "Skipping line %d: expected %d fields, saw %d\n",
                        self->file_lines, ex_fields, fields);
                append_warning(self, msg);
                free(msg);
            }
        }
    }
    else {
        /* missing trailing delimiters */
        if ((self->lines >= self->header_end + 1) && fields < ex_fields) {

            /* Might overrun the buffer when closing fields */
            if (make_stream_space(self, ex_fields - fields) < 0) {
                self->error_msg = "out of memory";
                return -1;
            }

            while (fields < ex_fields){
                end_field(self);
                /* printf("Prior word: %s\n", self->words[self->words_len - 2]); */
                fields++;
            }
        }

        // increment both line counts
        self->file_lines++;

        self->lines++;

        /* coliter_t it; */
        /* coliter_setup(&it, self, 5, self->lines - 1); */
        /* printf("word at column 5: %s\n", COLITER_NEXT(it)); */

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

int parser_set_skipfirstnrows(parser_t *self, int64_t nrows) {
    // self->file_lines is zero based so subtract 1 from nrows
    if (nrows > 0) {
        self->skip_first_N_rows = nrows - 1;
    }

    return 0;
}

static int parser_buffer_bytes(parser_t *self, size_t nbytes) {
    int status;
    size_t bytes_read;
    void *src = self->source;

    status = 0;
    self->datapos = 0;
    self->data = self->cb_io(self->source, nbytes, &bytes_read, &status);
    self->datalen = bytes_read;

    if (status != REACHED_EOF && self->data == NULL) {
        self->error_msg = (char*) malloc(200);

        if (status == CALLING_READ_FAILED) {
            sprintf(self->error_msg, ("Calling read(nbytes) on source failed. "
                                      "Try engine='python'."));
        } else {
            sprintf(self->error_msg, "Unknown error in IO callback");
        }
        return -1;
    }

    TRACE(("datalen: %d\n", self->datalen));

    return status;
}


/*

  Tokenization macros and state machine code

*/

//    printf("pushing %c\n", c);

#if defined(VERBOSE)
#define PUSH_CHAR(c)                                \
    printf("Pushing %c, slen now: %d\n", c, slen);  \
    *stream++ = c;                                  \
    slen++;
#else
#define PUSH_CHAR(c)                           \
    *stream++ = c;                             \
    slen++;
#endif



// This is a little bit of a hack but works for now

#define END_FIELD()                                \
    self->stream_len = slen;                   \
    if (end_field(self) < 0) {                 \
        goto parsingerror;                     \
    }                                          \
    stream = self->stream + self->stream_len;  \
    slen = self->stream_len;

#define END_LINE_STATE(STATE)                                           \
    self->stream_len = slen;                                            \
    if (end_line(self) < 0) {                                           \
        goto parsingerror;                                              \
    }                                                                   \
    stream = self->stream + self->stream_len;                           \
    slen = self->stream_len;                                            \
    self->state = STATE;                                                \
    if (line_limit > 0 && self->lines == start_lines + line_limit) {    \
        goto linelimit;                                                 \
                                                                        \
    }

#define END_LINE_AND_FIELD_STATE(STATE)                                 \
    self->stream_len = slen;                                            \
    if (end_line(self) < 0) {                                           \
        goto parsingerror;                                              \
    }                                                                   \
    if (end_field(self) < 0) {                                          \
        goto parsingerror;                                              \
    }                                                                   \
    stream = self->stream + self->stream_len;                           \
    slen = self->stream_len;                                            \
    self->state = STATE;                                                \
    if (line_limit > 0 && self->lines == start_lines + line_limit) {    \
        goto linelimit;                                                 \
                                                                        \
    }

#define END_LINE() END_LINE_STATE(START_RECORD)

#define IS_WHITESPACE(c) ((c == ' ' || c == '\t'))

typedef int (*parser_op)(parser_t *self, size_t line_limit);

#define _TOKEN_CLEANUP()                                                \
    self->stream_len = slen;                                            \
    self->datapos = i;                                                  \
    TRACE(("datapos: %d, datalen: %d\n", self->datapos, self->datalen));


int skip_this_line(parser_t *self, int64_t rownum) {
    if (self->skipset != NULL) {
        return ( kh_get_int64((kh_int64_t*) self->skipset, self->file_lines) != 
                 ((kh_int64_t*)self->skipset)->n_buckets );
    }
    else {
        return ( rownum <= self->skip_first_N_rows );
    }
}

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

        case SKIP_LINE:
//            TRACE(("tokenize_delimited SKIP_LINE %c, state %d\n", c, self->state));
            if (c == '\n') {
                END_LINE();
            }
            break;

        case START_RECORD:
            // start of record
            if (skip_this_line(self, self->file_lines)) {
                if (c == '\n') {
                    END_LINE()
                }
                else {
                    self->state = SKIP_LINE;
                }
                break;
            }
            else if (c == '\n') {
                // \n\r possible?
                if (self->skip_empty_lines)
                {
                    self->file_lines++;
                }
                else
                {
                    END_LINE();
                }
                break;
            }
            else if (c == '\r') {
                if (self->skip_empty_lines)
                {
                    self->file_lines++;
                    self->state = EAT_CRNL_NOP;
                }
                else
                    self->state = EAT_CRNL;
                break;
            } 
            else if (c == self->commentchar) {
                self->state = EAT_LINE_COMMENT;
                break;
            }
            else if (IS_WHITESPACE(c) && c != self->delimiter && self->skip_empty_lines) {
                self->state = WHITESPACE_LINE;
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
            else if (c == self->commentchar) {
                END_FIELD();
                self->state = EAT_COMMENT;
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

        case WHITESPACE_LINE: // check if line is whitespace-only
            if (c == '\n') {
                self->file_lines++;
                self->state = START_RECORD; // ignore empty line
            }
            else if (c == '\r') {
                self->file_lines++;
                self->state = EAT_CRNL_NOP;
            }
            else if (IS_WHITESPACE(c) && c != self->delimiter)
                ;
            else { // backtrack
                /* We have to use i + 1 because buf has been incremented but not i */
                while (i + 1 > self->datapos && *buf != '\n') {
                    --buf;
                    --i;
                }
                if (i + 1 > self->datapos) // reached a newline rather than the beginning
                {
                    ++buf; // move pointer to first char after newline
                    ++i;
                }
                self->state = START_FIELD;
            }
            break;

        case ESCAPED_CHAR:
            /* if (c == '\0') */
            /*  c = '\n'; */

            PUSH_CHAR(c);
            self->state = IN_FIELD;
            break;

        case EAT_LINE_COMMENT:
            if (c == '\n') {
                self->file_lines++;
                self->state = START_RECORD;
            } else if (c == '\r') {
                self->file_lines++;
                self->state = EAT_CRNL_NOP;
            }
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
            else if (c == self->commentchar) {
                END_FIELD();
                self->state = EAT_COMMENT;
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

        case EAT_COMMENT:
            if (c == '\n') {
                END_LINE();
            } else if (c == '\r') {
                self->state = EAT_CRNL;
            }
            break;

        case EAT_CRNL:
            if (c == '\n') {
                END_LINE();
                /* self->state = START_RECORD; */
            } else if (c == self->delimiter){
                // Handle \r-delimited files
                END_LINE_AND_FIELD_STATE(START_FIELD);
            } else {
                /* \r line terminator */

                /* UGH. we don't actually want to consume the token. fix this later */
                self->stream_len = slen;
                if (end_line(self) < 0) {
                    goto parsingerror;
                }
                stream = self->stream + self->stream_len;
                slen = self->stream_len;
                self->state = START_RECORD;

                /* HACK, let's try this one again */
                --i; buf--;
                if (line_limit > 0 && self->lines == start_lines + line_limit) {
                    goto linelimit;
                }

            }
            break;

        case EAT_CRNL_NOP: /* inside an ignored comment line */
            self->state = START_RECORD;
            /* \r line terminator -- parse this character again */
            if (c != '\n' && c != self->delimiter) {
                --i;
                --buf;
            }
            break;
        default:
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

/* custom line terminator */
int tokenize_delim_customterm(parser_t *self, size_t line_limit)
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

        case SKIP_LINE:
//            TRACE(("tokenize_delim_customterm SKIP_LINE %c, state %d\n", c, self->state));
            if (c == self->lineterminator) {
                END_LINE();
            }
            break;

        case START_RECORD:
            // start of record
            if (skip_this_line(self, self->file_lines)) {
                if (c == self->lineterminator) {
                    END_LINE()
                }
                else {
                    self->state = SKIP_LINE;
                }
                break;
            }
            else if (c == self->lineterminator) {
                // \n\r possible?
                if (self->skip_empty_lines)
                {
                    self->file_lines++;
                }
                else
                {
                    END_LINE();
                }
                break;
            }
            else if (c == self->commentchar) {
                self->state = EAT_LINE_COMMENT;
                break;
            }
            else if (IS_WHITESPACE(c) && c != self->delimiter && self->skip_empty_lines)
            {
                self->state = WHITESPACE_LINE;
                break;
            }
            /* normal character - handle as START_FIELD */
            self->state = START_FIELD;
            /* fallthru */
        case START_FIELD:
            /* expecting field */
            if (c == self->lineterminator) {
                END_FIELD();
                END_LINE();
                /* self->state = START_RECORD; */
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
            else if (c == self->commentchar) {
                END_FIELD();
                self->state = EAT_COMMENT;
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

        case WHITESPACE_LINE: // check if line is whitespace-only
            if (c == self->lineterminator) {
                self->file_lines++;
                self->state = START_RECORD; // ignore empty line
            }
            else if (IS_WHITESPACE(c) && c != self->delimiter)
                ;
            else { // backtrack
                /* We have to use i + 1 because buf has been incremented but not i */
                while (i + 1 > self->datapos && *buf != self->lineterminator) {
                    --buf;
                    --i;
                }
                if (i + 1 > self->datapos) // reached a newline rather than the beginning
                {
                    ++buf; // move pointer to first char after newline
                    ++i;
                }
                self->state = START_FIELD;
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
            if (c == self->lineterminator) {
                END_FIELD();
                END_LINE();
                /* self->state = START_RECORD; */
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
            else if (c == self->commentchar) {
                END_FIELD();
                self->state = EAT_COMMENT;
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
            else if (c == self->lineterminator) {
                END_FIELD();
                END_LINE();
                /* self->state = START_RECORD; */
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

        case EAT_LINE_COMMENT:
            if (c == self->lineterminator) {
                self->file_lines++;
                self->state = START_RECORD;
            }
            break;

        case EAT_COMMENT:
            if (c == self->lineterminator) {
                END_LINE();
            }
            break;

        default:
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

        TRACE(("Iter: %d Char: %c Line %d field_count %d, state %d\n",
               i, c, self->file_lines + 1, self->line_fields[self->lines],
               self->state));

        switch(self->state) {

        case SKIP_LINE:
//            TRACE(("tokenize_whitespace SKIP_LINE %c, state %d\n", c, self->state));
            if (c == '\n') {
                END_LINE();
            }
            break;

        case WHITESPACE_LINE:
            if (c == '\n') {
                self->file_lines++;
                self->state = START_RECORD;
                break;
            }
            else if (c == '\r') {
                self->file_lines++;
                self->state = EAT_CRNL_NOP;
                break;
            }
            // fall through

        case EAT_WHITESPACE:
            if (c == '\n') {
                END_LINE();
                self->state = START_RECORD;
                break;
            } else if (c == '\r') {
                self->state = EAT_CRNL;
                break;
            } else if (!IS_WHITESPACE(c)) {
                // END_FIELD();
                self->state = START_FIELD;
                // Fall through to subsequent state
            } else {
                // if whitespace char, keep slurping
                break;
            }

        case START_RECORD:
            // start of record
            if (skip_this_line(self, self->file_lines)) {
                if (c == '\n') {
                    END_LINE()
                }
                else {
                    self->state = SKIP_LINE;
                }
                break;
            } else  if (c == '\n') {
                if (self->skip_empty_lines)
                // \n\r possible?
                {
                    self->file_lines++;
                }
                else
                {
                    END_LINE();
                }
                break;
            } else if (c == '\r') {
                if (self->skip_empty_lines)
                {
                    self->file_lines++;
                    self->state = EAT_CRNL_NOP;
                }
                else
                    self->state = EAT_CRNL;
                break;
            } else if (IS_WHITESPACE(c)) {
                /*if (self->skip_empty_lines)
                    self->state = WHITESPACE_LINE;
                    else*/
                    self->state = EAT_WHITESPACE;
                break;
            } else if (c == self->commentchar) {
                self->state = EAT_LINE_COMMENT;
                break;
            } else {
                /* normal character - handle as START_FIELD */
                self->state = START_FIELD;
            }
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
            else if (c == self->commentchar) {
                END_FIELD();
                self->state = EAT_COMMENT;
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

        case EAT_LINE_COMMENT:
            if (c == '\n') {
                self->file_lines++;
                self->state = START_RECORD;
            } else if (c == '\r') {
                self->file_lines++;
                self->state = EAT_CRNL_NOP;
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
            else if (c == self->commentchar) {
                END_FIELD();
                self->state = EAT_COMMENT;
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
            } else if (IS_WHITESPACE(c)){
                // Handle \r-delimited files
                END_LINE_STATE(EAT_WHITESPACE);
            } else {
                /* XXX
                 * first character of a new record--need to back up and reread
                 * to handle properly...
                 */
                i--; buf--; /* back up one character (HACK!) */
                END_LINE_STATE(START_RECORD);
            }
            break;

        case EAT_CRNL_NOP: // inside an ignored comment line
            self->state = START_RECORD;
            /* \r line terminator -- parse this character again */
            if (c != '\n' && c != self->delimiter) {
                --i;
                --buf;
            }
            break;

        case EAT_COMMENT:
            if (c == '\n') {
                END_LINE();
            } else if (c == '\r') {
                self->state = EAT_CRNL;
            }
            break;

        default:
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


static int parser_handle_eof(parser_t *self) {
    TRACE(("handling eof, datalen: %d, pstate: %d\n", self->datalen, self->state))
    if (self->datalen == 0 && (self->state != START_RECORD)) {
        // test cases needed here
        // TODO: empty field at end of line
        TRACE(("handling eof\n"));

        if (self->state == IN_FIELD || self->state == START_FIELD) {
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
    else if (self->datalen == 0 && (self->state == START_RECORD)) {
        return 0;
    }

    return -1;
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
    /*
    printf("Line_start: ");
    for (i = 0; i < self->lines + 1; ++i) {
         printf("%d ", self->line_fields[i]);
     }
    printf("\n");
    */
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
    new_cap = _next_pow2(self->stream_len) + 1;
    if (new_cap < self->stream_cap) {
        self->stream = safe_realloc((void*) self->stream, new_cap);
        self->stream_cap = new_cap;
    }

    /* trim words, word_starts */
    new_cap = _next_pow2(self->words_len) + 1;
    if (new_cap < self->words_cap) {
        self->words = (char**) safe_realloc((void*) self->words,
                                            new_cap * sizeof(char*));
        self->word_starts = (int*) safe_realloc((void*) self->word_starts,
                                                new_cap * sizeof(int));
        self->words_cap = new_cap;
    }

    /* trim line_start, line_fields */
    new_cap = _next_pow2(self->lines) + 1;
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
    int j, line;
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
    } else if (self->lineterminator == '\0') {
        tokenize_bytes = tokenize_delimited;
    } else {
        tokenize_bytes = tokenize_delim_customterm;
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
            } else if (status != 0) {
                return status;
            }
        }

        TRACE(("Trying to process %d bytes\n", self->datalen - self->datapos));
        /* TRACE(("sourcetype: %c, status: %d\n", self->sourcetype, status)); */

        status = tokenize_bytes(self, nrows);

        /* debug_print_parser(self); */

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



P_INLINE void lowercase(char *p) {
    for ( ; *p; ++p) *p = tolower(*p);
}

P_INLINE void uppercase(char *p) {
    for ( ; *p; ++p) *p = toupper(*p);
}


int P_INLINE to_complex(char *item, double *p_real, double *p_imag, char sci, char decimal)
{
    char *p_end;

    *p_real = xstrtod(item, &p_end, decimal, sci, '\0', FALSE);
    if (*p_end == '\0') {
        *p_imag = 0.0;
        return errno == 0;
    }
    if (*p_end == 'i' || *p_end == 'j') {
        *p_imag = *p_real;
        *p_real = 0.0;
        ++p_end;
    }
    else {
        if (*p_end == '+') {
            ++p_end;
        }
        *p_imag = xstrtod(p_end, &p_end, decimal, sci, '\0', FALSE);
        if (errno || ((*p_end != 'i') && (*p_end != 'j'))) {
            return FALSE;
        }
        ++p_end;
    }
    while(*p_end == ' ') {
        ++p_end;
    }
    return *p_end == '\0';
}


int P_INLINE to_longlong(char *item, long long *p_value)
{
    char *p_end;

    // Try integer conversion.  We explicitly give the base to be 10. If
    // we used 0, strtoll() would convert '012' to 10, because the leading 0 in
    // '012' signals an octal number in C.  For a general purpose reader, that
    // would be a bug, not a feature.
    *p_value = strtoll(item, &p_end, 10);

    // Allow trailing spaces.
    while (isspace(*p_end)) ++p_end;

    return (errno == 0) && (!*p_end);
}

int P_INLINE to_longlong_thousands(char *item, long long *p_value, char tsep)
{
    int i, pos, status, n = strlen(item), count = 0;
    char *tmp;
    char *p_end;

    for (i = 0; i < n; ++i)
    {
        if (*(item + i) == tsep) {
            count++;
        }
    }

    if (count == 0) {
        return to_longlong(item, p_value);
    }

    tmp = (char*) malloc((n - count + 1) * sizeof(char));
    if (tmp == NULL) {
        return 0;
    }

    pos = 0;
    for (i = 0; i < n; ++i)
    {
        if (item[i] != tsep)
            tmp[pos++] = item[i];
    }

    tmp[pos] = '\0';

    status = to_longlong(tmp, p_value);
    free(tmp);

    return status;
}

int to_boolean(char *item, uint8_t *val) {
    char *tmp;
    int i, status = 0;

    static const char *tstrs[1] = {"TRUE"};
    static const char *fstrs[1] = {"FALSE"};

    tmp = malloc(sizeof(char) * (strlen(item) + 1));
    strcpy(tmp, item);
    uppercase(tmp);

    for (i = 0; i < 1; ++i)
    {
        if (strcmp(tmp, tstrs[i]) == 0) {
            *val = 1;
            goto done;
        }
    }

    for (i = 0; i < 1; ++i)
    {
        if (strcmp(tmp, fstrs[i]) == 0) {
            *val = 0;
            goto done;
        }
    }

    status = -1;

done:
    free(tmp);
    return status;
}

// #define TEST

#ifdef TEST

int main(int argc, char *argv[])
{
    double x, y;
    long long xi;
    int status;
    char *s;

    //s = "0.10e-3-+5.5e2i";
    // s = "1-0j";
    // status = to_complex(s, &x, &y, 'e', '.');
    s = "123,789";
    status = to_longlong_thousands(s, &xi, ',');
    printf("s = '%s'\n", s);
    printf("status = %d\n", status);
    printf("x = %d\n", (int) xi);

    // printf("x = %lg,  y = %lg\n", x, y);

    return 0;
}
#endif

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
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.
//
// -----------------------------------------------------------------------
// Modifications by Warren Weckesser, March 2011:
// * Rename strtod() to xstrtod().
// * Added decimal and sci arguments.
// * Skip trailing spaces.
// * Commented out the other functions.
// Modifications by Richard T Guy, August 2013:
// * Add tsep argument for thousands separator
//

double xstrtod(const char *str, char **endptr, char decimal,
                      char sci, char tsep, int skip_trailing)
{
  double number;
  int exponent;
  int negative;
  char *p = (char *) str;
  double p10;
  int n;
  int num_digits;
  int num_decimals;

  errno = 0;

  // Skip leading whitespace
  while (isspace(*p)) p++;

  // Handle optional sign
  negative = 0;
  switch (*p)
  {
    case '-': negative = 1; // Fall through to increment position
    case '+': p++;
  }

  number = 0.;
  exponent = 0;
  num_digits = 0;
  num_decimals = 0;

  // Process string of digits
  while (isdigit(*p))
  {
    number = number * 10. + (*p - '0');
    p++;
    num_digits++;

    p += (tsep != '\0' & *p == tsep);
  }

  // Process decimal part
  if (*p == decimal)
  {
    p++;

    while (isdigit(*p))
    {
      number = number * 10. + (*p - '0');
      p++;
      num_digits++;
      num_decimals++;
    }

    exponent -= num_decimals;
  }

  if (num_digits == 0)
  {
    errno = ERANGE;
    return 0.0;
  }

  // Correct for sign
  if (negative) number = -number;

  // Process an exponent string
  if (toupper(*p) == toupper(sci))
  {
    // Handle optional sign
    negative = 0;
    switch (*++p)
    {
      case '-': negative = 1;   // Fall through to increment pos
      case '+': p++;
    }

    // Process string of digits
    n = 0;
    while (isdigit(*p))
    {
      n = n * 10 + (*p - '0');
      p++;
    }

    if (negative)
      exponent -= n;
    else
      exponent += n;
  }


  if (exponent < DBL_MIN_EXP  || exponent > DBL_MAX_EXP)
  {

    errno = ERANGE;
    return HUGE_VAL;
  }

  // Scale the result
  p10 = 10.;
  n = exponent;
  if (n < 0) n = -n;
  while (n)
  {
    if (n & 1)
    {
      if (exponent < 0)
        number /= p10;
      else
        number *= p10;
    }
    n >>= 1;
    p10 *= p10;
  }


  if (number == HUGE_VAL) {
      errno = ERANGE;
  }

  if (skip_trailing) {
      // Skip trailing whitespace
      while (isspace(*p)) p++;
  }

  if (endptr) *endptr = p;


  return number;
}

double precise_xstrtod(const char *str, char **endptr, char decimal,
               char sci, char tsep, int skip_trailing)
{
    double number;
    int exponent;
    int negative;
    char *p = (char *) str;
    int num_digits;
    int num_decimals;
    int max_digits = 17;
    int n;
    // Cache powers of 10 in memory
    static double e[] = {1., 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10,
                         1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20,
                         1e21, 1e22, 1e23, 1e24, 1e25, 1e26, 1e27, 1e28, 1e29, 1e30,
                         1e31, 1e32, 1e33, 1e34, 1e35, 1e36, 1e37, 1e38, 1e39, 1e40,
                         1e41, 1e42, 1e43, 1e44, 1e45, 1e46, 1e47, 1e48, 1e49, 1e50,
                         1e51, 1e52, 1e53, 1e54, 1e55, 1e56, 1e57, 1e58, 1e59, 1e60,
                         1e61, 1e62, 1e63, 1e64, 1e65, 1e66, 1e67, 1e68, 1e69, 1e70,
                         1e71, 1e72, 1e73, 1e74, 1e75, 1e76, 1e77, 1e78, 1e79, 1e80,
                         1e81, 1e82, 1e83, 1e84, 1e85, 1e86, 1e87, 1e88, 1e89, 1e90,
                         1e91, 1e92, 1e93, 1e94, 1e95, 1e96, 1e97, 1e98, 1e99, 1e100,
                         1e101, 1e102, 1e103, 1e104, 1e105, 1e106, 1e107, 1e108, 1e109, 1e110,
                         1e111, 1e112, 1e113, 1e114, 1e115, 1e116, 1e117, 1e118, 1e119, 1e120,
                         1e121, 1e122, 1e123, 1e124, 1e125, 1e126, 1e127, 1e128, 1e129, 1e130,
                         1e131, 1e132, 1e133, 1e134, 1e135, 1e136, 1e137, 1e138, 1e139, 1e140,
                         1e141, 1e142, 1e143, 1e144, 1e145, 1e146, 1e147, 1e148, 1e149, 1e150,
                         1e151, 1e152, 1e153, 1e154, 1e155, 1e156, 1e157, 1e158, 1e159, 1e160,
                         1e161, 1e162, 1e163, 1e164, 1e165, 1e166, 1e167, 1e168, 1e169, 1e170,
                         1e171, 1e172, 1e173, 1e174, 1e175, 1e176, 1e177, 1e178, 1e179, 1e180,
                         1e181, 1e182, 1e183, 1e184, 1e185, 1e186, 1e187, 1e188, 1e189, 1e190,
                         1e191, 1e192, 1e193, 1e194, 1e195, 1e196, 1e197, 1e198, 1e199, 1e200,
                         1e201, 1e202, 1e203, 1e204, 1e205, 1e206, 1e207, 1e208, 1e209, 1e210,
                         1e211, 1e212, 1e213, 1e214, 1e215, 1e216, 1e217, 1e218, 1e219, 1e220,
                         1e221, 1e222, 1e223, 1e224, 1e225, 1e226, 1e227, 1e228, 1e229, 1e230,
                         1e231, 1e232, 1e233, 1e234, 1e235, 1e236, 1e237, 1e238, 1e239, 1e240,
                         1e241, 1e242, 1e243, 1e244, 1e245, 1e246, 1e247, 1e248, 1e249, 1e250,
                         1e251, 1e252, 1e253, 1e254, 1e255, 1e256, 1e257, 1e258, 1e259, 1e260,
                         1e261, 1e262, 1e263, 1e264, 1e265, 1e266, 1e267, 1e268, 1e269, 1e270,
                         1e271, 1e272, 1e273, 1e274, 1e275, 1e276, 1e277, 1e278, 1e279, 1e280,
                         1e281, 1e282, 1e283, 1e284, 1e285, 1e286, 1e287, 1e288, 1e289, 1e290,
                         1e291, 1e292, 1e293, 1e294, 1e295, 1e296, 1e297, 1e298, 1e299, 1e300,
                         1e301, 1e302, 1e303, 1e304, 1e305, 1e306, 1e307, 1e308};
    errno = 0;

    // Skip leading whitespace
    while (isspace(*p)) p++;

    // Handle optional sign
    negative = 0;
    switch (*p)
    {
    case '-': negative = 1; // Fall through to increment position
    case '+': p++;
    }

    number = 0.;
    exponent = 0;
    num_digits = 0;
    num_decimals = 0;

    // Process string of digits
    while (isdigit(*p))
    {
        if (num_digits < max_digits)
        {
            number = number * 10. + (*p - '0');
            num_digits++;
        }
        else
            ++exponent;

        p++;
        p += (tsep != '\0' & *p == tsep);
    }

    // Process decimal part
    if (*p == decimal)
    {
        p++;

        while (num_digits < max_digits && isdigit(*p))
        {
            number = number * 10. + (*p - '0');
            p++;
            num_digits++;
            num_decimals++;
        }
        
        if (num_digits >= max_digits) // consume extra decimal digits
            while (isdigit(*p))
                ++p;

        exponent -= num_decimals;
    }

    if (num_digits == 0)
    {
        errno = ERANGE;
        return 0.0;
    }

    // Correct for sign
    if (negative) number = -number;

    // Process an exponent string
    if (toupper(*p) == toupper(sci))
    {
        // Handle optional sign
        negative = 0;
        switch (*++p)
        {
        case '-': negative = 1;   // Fall through to increment pos
        case '+': p++;
        }

        // Process string of digits
        n = 0;
        while (isdigit(*p))
        {
            n = n * 10 + (*p - '0');
            p++;
        }

        if (negative)
            exponent -= n;
        else
            exponent += n;
    }

    if (exponent > 308)
    {
        errno = ERANGE;
        return HUGE_VAL;
    }
    else if (exponent > 0)
        number *= e[exponent];
    else if (exponent < -308) // subnormal
    {
        if (exponent < -616) // prevent invalid array access
            number = 0.;
        number /= e[-308 - exponent];
        number /= e[308];
    }
    else
        number /= e[-exponent];

    if (number == HUGE_VAL || number == -HUGE_VAL)
        errno = ERANGE;

    if (skip_trailing) {
        // Skip trailing whitespace
        while (isspace(*p)) p++;
    }

    if (endptr) *endptr = p;
    return number;
}

double round_trip(const char *p, char **q, char decimal, char sci,
                  char tsep, int skip_trailing)
{
#if PY_VERSION_HEX >= 0x02070000
    return PyOS_string_to_double(p, q, 0);
#else
    return strtod(p, q);
#endif
}

/*
float strtof(const char *str, char **endptr)
{
  return (float) strtod(str, endptr);
}


long double strtold(const char *str, char **endptr)
{
  return strtod(str, endptr);
}

double atof(const char *str)
{
  return strtod(str, NULL);
}
*/

// End of xstrtod code
// ---------------------------------------------------------------------------

int64_t str_to_int64(const char *p_item, int64_t int_min, int64_t int_max,
                     int *error, char tsep)
{
    const char *p = (const char *) p_item;
    int isneg = 0;
    int64_t number = 0;
    int d;

    // Skip leading spaces.
    while (isspace(*p)) {
        ++p;
    }

    // Handle sign.
    if (*p == '-') {
        isneg = 1;
        ++p;
    }
    else if (*p == '+') {
        p++;
    }

    // Check that there is a first digit.
    if (!isdigit(*p)) {
        // Error...
        *error = ERROR_NO_DIGITS;
        return 0;
    }

    if (isneg) {
        // If number is greater than pre_min, at least one more digit
        // can be processed without overflowing.
        int dig_pre_min = -(int_min % 10);
        int64_t pre_min = int_min / 10;

        // Process the digits.
        d = *p;
        if (tsep != '\0') {
            while (1) {
                if (d == tsep) {
                    d = *++p;
                    continue;
                } else if (!isdigit(d)) {
                    break;
                }
                if ((number > pre_min) ||
                    ((number == pre_min) && (d - '0' <= dig_pre_min))) {

                    number = number * 10 - (d - '0');
                    d = *++p;
                }
                else {
                    *error = ERROR_OVERFLOW;
                    return 0;
                }
            }
        } else {
            while (isdigit(d)) {
                if ((number > pre_min) ||
                    ((number == pre_min) && (d - '0' <= dig_pre_min))) {

                    number = number * 10 - (d - '0');
                    d = *++p;
                }
                else {
                    *error = ERROR_OVERFLOW;
                    return 0;
                }
            }
        }
    }
    else {
        // If number is less than pre_max, at least one more digit
        // can be processed without overflowing.
        int64_t pre_max = int_max / 10;
        int dig_pre_max = int_max % 10;

        //printf("pre_max = %lld  dig_pre_max = %d\n", pre_max, dig_pre_max);

        // Process the digits.
        d = *p;
        if (tsep != '\0') {
            while (1) {
                if (d == tsep) {
                    d = *++p;
                    continue;
                } else if (!isdigit(d)) {
                    break;
                }
                if ((number < pre_max) ||
                    ((number == pre_max) && (d - '0' <= dig_pre_max))) {

                    number = number * 10 + (d - '0');
                    d = *++p;

                }
                else {
                    *error = ERROR_OVERFLOW;
                    return 0;
                }
            }
        } else {
            while (isdigit(d)) {
                if ((number < pre_max) ||
                    ((number == pre_max) && (d - '0' <= dig_pre_max))) {

                    number = number * 10 + (d - '0');
                    d = *++p;

                }
                else {
                    *error = ERROR_OVERFLOW;
                    return 0;
                }
            }
        }
    }

    // Skip trailing spaces.
    while (isspace(*p)) {
        ++p;
    }

    // Did we use up all the characters?
    if (*p) {
        *error = ERROR_INVALID_CHARS;
        return 0;
    }

    *error = 0;
    return number;
}


uint64_t str_to_uint64(const char *p_item, uint64_t uint_max, int *error)
{
    int d, dig_pre_max;
    uint64_t pre_max;
    const char *p = (const char *) p_item;
    uint64_t number = 0;

    // Skip leading spaces.
    while (isspace(*p)) {
        ++p;
    }

    // Handle sign.
    if (*p == '-') {
        *error = ERROR_MINUS_SIGN;
        return 0;
    }
    if (*p == '+') {
        p++;
    }

    // Check that there is a first digit.
    if (!isdigit(*p)) {
        // Error...
        *error = ERROR_NO_DIGITS;
        return 0;
    }

    // If number is less than pre_max, at least one more digit
    // can be processed without overflowing.
    pre_max = uint_max / 10;
    dig_pre_max = uint_max % 10;

    // Process the digits.
    d = *p;
    while (isdigit(d)) {
        if ((number < pre_max) || ((number == pre_max) && (d - '0' <= dig_pre_max))) {
            number = number * 10 + (d - '0');
            d = *++p;
        }
        else {
            *error = ERROR_OVERFLOW;
            return 0;
        }
    }

    // Skip trailing spaces.
    while (isspace(*p)) {
        ++p;
    }

    // Did we use up all the characters?
    if (*p) {
        *error = ERROR_INVALID_CHARS;
        return 0;
    }

    *error = 0;
    return number;
}
