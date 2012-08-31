from numpy cimport *
import numpy as np

import_array()

from khash cimport *

from cpython cimport PyString_FromString, Py_INCREF, PyString_AsString

cdef extern from "stdint.h":
    enum: INT64_MAX
    enum: INT64_MIN
    enum: INT32_MAX
    enum: INT32_MIN


cdef extern from "Python.h":
    ctypedef struct FILE
    FILE* PyFile_AsFile(object)

cdef extern from "parser/conversions.h":
    inline int to_double(char *item, double *p_value,
                         char sci, char decimal)
    inline int to_complex(char *item, double *p_real,
                          double *p_imag, char sci, char decimal)
    inline int to_longlong(char *item, long long *p_value)


cdef extern from "parser/common.h":

    ctypedef enum ParserState:
        START_RECORD
        START_FIELD
        ESCAPED_CHAR
        IN_FIELD
        IN_QUOTED_FIELD
        ESCAPE_IN_QUOTED_FIELD
        QUOTE_IN_QUOTED_FIELD
        EAT_CRNL
        FINISHED

    ctypedef struct table_chunk:
        void **columns
        int ncols

    ctypedef struct parser_t:
        void *source
        char sourcetype   # 'M' for mmap, 'F' for FILE, 'A' for array

        int chunksize  # Number of bytes to prepare for each chunk
        char *data     # pointer to data to be processed
        int datalen    # amount of data available

        # where to write out tokenized data
        char *stream
        int stream_len
        int stream_cap

        # Store words in (potentially ragged) matrix for now, hmm
        char **words
        int *word_starts # where we are in the stream
        int words_len
        int words_cap

        char *pword_start    # pointer to stream start of current field
        int word_start       # position start of current field

        int *line_start      # position in words for start of line
        int *line_fields     # Number of fields in each line
        int lines            # Number of lines observed
        int lines_cap        # Vector capacity

        # Tokenizing stuff
        ParserState state
        int doublequote            # is " represented by ""? */
        char delimiter             # field separator */
        char quotechar             # quote character */
        char escapechar            # escape character */
        int skipinitialspace       # ignore spaces following delimiter? */
        int quoting                # style of quoting to write */

        # hmm =/
        int numeric_field

        char commentchar
        int allow_embedded_newline
        int strict                 # raise exception on bad CSV */

        int error_bad_lines
        int warn_bad_lines

        int infer_types

        # floating point options
        char decimal
        char sci

        # thousands separator (comma, period)
        char thousands

        int header # Boolean: 1: has header, 0: no header

        int skiprows
        int skip_footer

        table_chunk *chunks
        int nchunks

        void **columns
        int ncols

        # PyObject *converters

        #  error handling
        char *error_msg

    ctypedef struct coliter_t:
        char **words
        int *line_start
        int col
        int line

    void coliter_setup(coliter_t *it, parser_t *parser, int i)
    char* COLITER_NEXT(coliter_t it)

    parser_t* parser_new()

    int parser_init(parser_t *self)
    void parser_free(parser_t *self)

    void set_parser_default_options(parser_t *self)

    int parser_file_source_init(parser_t *self, FILE* fp)
    int parser_mmap_init(parser_t *self, FILE* fp)
    int parser_array_source_init(parser_t *self, char *bytes, size_t length)

    void debug_print_parser(parser_t *self)

    int tokenize_all_rows(parser_t *self)
    int tokenize_nrows(parser_t *self, size_t nrows)

    int64_t str_to_int64(char *p_item, int64_t int_min,
                         int64_t int_max, int *error)
    uint64_t str_to_uint64(char *p_item, uint64_t uint_max, int *error)



DEFAULT_CHUNKSIZE = 256 * 1024

cdef class TextReader:
    '''

    # source: StringIO or file object

    '''

    cdef:
        parser_t *parser
        object file_handle, should_close
        bint factorize

    cdef public:
        object delimiter, na_values, converters, thousands, delim_whitespace
        object memory_map

    def __cinit__(self, source, delimiter=',', header=0,
                  memory_map=False,
                  chunksize=DEFAULT_CHUNKSIZE,
                  delim_whitespace=False,
                  na_values=None,
                  converters=None,
                  thousands=None,
                  factorize=True):
        self.parser = parser_new()
        self.parser.chunksize = chunksize

        self._setup_parser_source(source)
        set_parser_default_options(self.parser)

        parser_init(self.parser)

        if delim_whitespace:
            raise NotImplementedError
        else:
            if len(delimiter) > 1:
                raise ValueError('only length-1 separators excluded right now')
            self.parser.delimiter = (<char*> delimiter)[0]

        self.factorize = factorize

        # TODO: no header vs. header is not the first row
        self.parser.header = header

        self.should_close = False

        self.delimiter = delimiter
        self.delim_whitespace = delim_whitespace

        self.memory_map = memory_map
        self.na_values = na_values
        self.converters = converters
        self.thousands = thousands

    def __dealloc__(self):
        parser_free(self.parser)

    def __del__(self):
        if self.should_close:
            self.file_handle.close()

    cdef _setup_parser_source(self, source):
        cdef:
            int status

        if isinstance(source, (basestring, file)):
            if isinstance(source, basestring):
                source = open(source, 'rb')
                self.should_close = True

            self.file_handle = source

            if self.memory_map:
                status = parser_mmap_init(self.parser,
                                          PyFile_AsFile(source))
            else:
                status = parser_file_source_init(self.parser,
                                                 PyFile_AsFile(source))

            if status != 0:
                raise Exception('Initializing from file failed')

        elif hasattr(source, 'read'):
            # e.g., StringIO

            bytes = source.read()

            # TODO: unicode
            if isinstance(bytes, unicode):
                raise ValueError('Only ascii/bytes supported at the moment')

            status = parser_array_source_init(self.parser,
                                              PyString_AsString(bytes),
                                              len(bytes))
            if status != 0:
                raise Exception('Initializing parser from file-like '
                                'object failed')

    def _parse_table_header(self):
        pass

    def read(self, rows=None):
        """
        rows=None --> read all rows
        """
        cdef:
            int prior_lines
            int status

        if rows is not None:
            raise NotImplementedError
        else:
            status = tokenize_all_rows(self.parser)

        if status < 0:
            raise_parser_error('Error tokenizing data', self.parser)

        result = self._convert_column_data()

        # debug_print_parser(self.parser)
        return result

    def _convert_column_data(self):
        cdef:
            Py_ssize_t i, ncols
            cast_func func

        ncols = self.parser.line_fields[0]

        results = {}
        for i in range(ncols):
            col_res = None
            for func in cast_func_order:
                col_res = func(self.parser, i, 0, self.parser.lines)
                if col_res is not None:
                    results[i] = col_res
                    break

            if col_res is None:
                raise Exception('Unable to parse column %d' % i)

            # col_res = _try_double(self.parser, i, 0, self.parser.lines)

            # if col_res is None:
            #     col_res = _string_box_factorize(self.parser, i,
            #                                     0, self.parser.lines)

            results[i] = col_res

        return results

class CParserError(Exception):
    pass


# ----------------------------------------------------------------------
# Type conversions / inference support code

ctypedef object (*cast_func)(parser_t *parser, int col,
                             int line_start, int line_end)

cdef _string_box_factorize(parser_t *parser, int col,
                           int line_start, int line_end):
    cdef:
        int error
        Py_ssize_t i
        size_t lines
        coliter_t it
        char *word
        ndarray[object] result

        int ret = 0
        kh_strbox_t *table
        kh_iter_t

        object pyval


    table = kh_init_strbox()

    lines = line_end - line_start
    result = np.empty(lines, dtype=np.object_)

    coliter_setup(&it, parser, col)
    for i in range(lines):
        word = COLITER_NEXT(it)

        k = kh_get_strbox(table, word)

        # in the hash table
        if k != table.n_buckets:
            # this increments the refcount, but need to test
            pyval = <object> table.vals[k]
        else:
            # box it. new ref?
            pyval = PyString_FromString(word)

            k = kh_put_strbox(table, word, &ret)
            table.vals[k] = <PyObject*> pyval

        result[i] = pyval

    return result



cdef _try_double(parser_t *parser, int col, int line_start, int line_end):
    cdef:
        int error
        size_t i, lines
        coliter_t it
        char *word
        double *data
        ndarray result

    lines = line_end - line_start

    result = np.empty(lines, dtype=np.float64)

    data = <double *> result.data

    coliter_setup(&it, parser, col)
    for i in range(lines):
        word = COLITER_NEXT(it)
        error = to_double(word, data, parser.sci, parser.decimal)

        if error != 1:
            return None

        data += 1

    return result

cdef _try_int64(parser_t *parser, int col, int line_start, int line_end):
    cdef:
        int error
        size_t i, lines
        coliter_t it
        char *word
        int64_t *data
        ndarray result

    lines = line_end - line_start

    result = np.empty(lines, dtype=np.int64)

    data = <int64_t *> result.data

    coliter_setup(&it, parser, col)
    for i in range(lines):
        word = COLITER_NEXT(it)
        data[i] = str_to_int64(word, INT64_MIN, INT64_MAX, &error);

        if error != 0:
            return None

    return result


cdef cast_func cast_func_order[3]
cast_func_order[0] = _try_int64
cast_func_order[1] = _try_double
cast_func_order[2] = _string_box_factorize

cdef raise_parser_error(object base, parser_t *parser):
    message = '%s. C error: ' % base
    if parser.error_msg != NULL:
        message += parser.error_msg
    else:
        message += 'no error message set'

    raise CParserError(message)
