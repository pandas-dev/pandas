# Copyright (c) 2012, Lambda Foundry, Inc.
# See LICENSE for the license

from libc.stdlib cimport malloc, free

from cpython cimport (PyObject, PyString_FromString,
                      PyString_AsString, PyString_Check)

cdef extern from "stdlib.h":
    void memcpy(void *dst, void *src, size_t n)

cimport numpy as cnp

from numpy cimport ndarray, uint8_t, uint64_t

import numpy as np

# cdef extern from "Python.h":
#     void Py_INCREF(PyObject*)
#     void Py_XDECREF(PyObject*)

cimport util

import pandas.lib as lib

import time

cnp.import_array()

from khash cimport *

import sys

cdef bint PY3 = (sys.version_info[0] >= 3)

cdef extern from "stdint.h":
    enum: UINT8_MAX
    enum: UINT16_MAX
    enum: UINT32_MAX
    enum: INT8_MIN
    enum: INT8_MAX
    enum: INT16_MIN
    enum: INT16_MAX
    enum: INT32_MAX
    enum: INT32_MIN
    enum: INT64_MAX
    enum: INT64_MIN

try:
    basestring
except NameError:
    basestring = str

cdef extern from "Python.h":
    ctypedef struct FILE
    FILE* PyFile_AsFile(object)

cdef extern from "parser/conversions.h":
    inline int to_double(char *item, double *p_value,
                         char sci, char decimal)
    inline int to_complex(char *item, double *p_real,
                          double *p_imag, char sci, char decimal)
    inline int to_longlong(char *item, long long *p_value)
    inline int to_longlong_thousands(char *item, long long *p_value,
                                     char tsep)
    inline int to_boolean(char *item, uint8_t *val)


cdef extern from "parser/parser.h":

    ctypedef enum ParserState:
        START_RECORD
        START_FIELD
        ESCAPED_CHAR
        IN_FIELD
        IN_QUOTED_FIELD
        ESCAPE_IN_QUOTED_FIELD
        QUOTE_IN_QUOTED_FIELD
        EAT_CRNL
        EAT_WHITESPACE
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
        int datapos

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
        int file_lines       # Number of file lines observed (with bad/skipped)
        int lines_cap        # Vector capacity

        # Tokenizing stuff
        ParserState state
        int doublequote            # is " represented by ""? */
        char delimiter             # field separator */
        int delim_whitespace       # consume tabs / spaces instead
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

        void *skipset
        int skip_footer

        table_chunk *chunks
        int nchunks

        void **columns
        int ncols

        #  error handling
        char *error_msg

    ctypedef struct coliter_t:
        char **words
        int *line_start
        int col

    void coliter_setup(coliter_t *it, parser_t *parser, int i, int start)
    char* COLITER_NEXT(coliter_t it)

    parser_t* parser_new()

    int parser_init(parser_t *self) nogil
    void parser_free(parser_t *self) nogil
    int parser_add_skiprow(parser_t *self, int64_t row)

    void parser_set_default_options(parser_t *self)

    int parser_file_source_init(parser_t *self, FILE* fp)
    int parser_mmap_init(parser_t *self, FILE* fp)
    int parser_array_source_init(parser_t *self, char *bytes, size_t length)

    int parser_consume_rows(parser_t *self, size_t nrows)

    int parser_trim_buffers(parser_t *self)

    void debug_print_parser(parser_t *self)

    int tokenize_all_rows(parser_t *self) nogil
    int tokenize_nrows(parser_t *self, size_t nrows) nogil

    int64_t str_to_int64(char *p_item, int64_t int_min,
                         int64_t int_max, int *error, char tsep)
    uint64_t str_to_uint64(char *p_item, uint64_t uint_max, int *error)



DEFAULT_CHUNKSIZE = 1024 * 1024

# common NA values
# no longer excluding inf representations
# '1.#INF','-1.#INF', '1.#INF000000',
_NA_VALUES = ['-1.#IND', '1.#QNAN', '1.#IND', '-1.#QNAN',
              '#N/A N/A', 'NA', '#NA', 'NULL', 'NaN',
              'nan', '']


cdef class TextReader:
    '''

    # source: StringIO or file object

    '''

    cdef:
        parser_t *parser
        object file_handle, should_close
        bint factorize, na_filter, verbose
        int parser_start
        list clocks

    cdef public:
        int leading_cols, table_width, skip_footer, buffer_lines
        object delimiter, na_values, converters, delim_whitespace
        object memory_map
        object as_recarray
        object header, names
        object low_memory
        object skiprows
        object compact_ints, use_unsigned
        set noconvert

    def __cinit__(self, source,
                  delimiter=b',',

                  header=0,
                  names=None,

                  memory_map=False,
                  tokenize_chunksize=DEFAULT_CHUNKSIZE,
                  delim_whitespace=False,

                  converters=None,

                  factorize=True,
                  as_recarray=False,

                  skipinitialspace=False,
                  escapechar=None,
                  doublequote=True,
                  quotechar=b'"',
                  quoting=0,

                  comment=None,
                  decimal=b'.',
                  thousands=None,

                  error_bad_lines=True,
                  warn_bad_lines=True,

                  na_filter=True,
                  na_values=None,
                  compact_ints=False,
                  use_unsigned=False,
                  low_memory=False,
                  buffer_lines=2**16,
                  skiprows=None,
                  skip_footer=0,
                  verbose=False):

        self.parser = parser_new()
        self.parser.chunksize = tokenize_chunksize

        # For timekeeping
        self.clocks = []

        self._setup_parser_source(source)
        parser_set_default_options(self.parser)

        parser_init(self.parser)

        if delim_whitespace:
            self.parser.delim_whitespace = delim_whitespace
        else:
            if len(delimiter) > 1:
                raise ValueError('only length-1 separators excluded right now')
            self.parser.delimiter = ord(delimiter)

        self.factorize = factorize

        #----------------------------------------
        # parser options

        self.parser.doublequote = doublequote
        self.parser.skipinitialspace = skipinitialspace

        if len(decimal) != 1:
            raise ValueError('Only length-1 decimal markers supported')
        self.parser.decimal = ord(decimal)

        if thousands is not None:
            if len(thousands) != 1:
                raise ValueError('Only length-1 decimal markers supported')
            self.parser.thousands = ord(thousands)

        if escapechar is not None:
            if len(escapechar) != 1:
                raise ValueError('Only length-1 escapes  supported')
            self.parser.escapechar = ord(escapechar)

        self.parser.quotechar = ord(quotechar)
        self.parser.quoting = quoting

        # error handling of bad lines
        self.parser.error_bad_lines = int(error_bad_lines)
        self.parser.warn_bad_lines = int(warn_bad_lines)

        self.skiprows = skiprows
        if skiprows is not None:
            self._make_skiprow_set()

        self.skip_footer = skip_footer

        # KLUDGE
        if skip_footer > 0:
            self.parser.error_bad_lines = 0
            self.parser.warn_bad_lines = 0

        self.should_close = False

        self.delimiter = delimiter
        self.delim_whitespace = delim_whitespace

        self.memory_map = memory_map
        self.na_values = na_values
        self.converters = converters

        self.na_filter = na_filter
        self.as_recarray = as_recarray

        self.compact_ints = compact_ints
        self.use_unsigned = use_unsigned

        self.verbose = verbose
        self.low_memory = low_memory
        self.buffer_lines = buffer_lines

        # XXX
        self.noconvert = set()

        #----------------------------------------
        # header stuff

        self.names = names
        self.leading_cols = 0
        if names is not None:
            header = None

        # TODO: no header vs. header is not the first row
        if header is None:
            # sentinel value
            self.parser.header = -1
            self.parser_start = 0
        else:
            self.parser.header = header
            self.parser_start = header + 1

        self.header, self.table_width = self._get_header()

        # print 'Header: %s, width: %d' % (str(self.header), self.table_width)

    def __init__(self, *args, **kwards):
        pass

    def __dealloc__(self):
        parser_free(self.parser)

    def __del__(self):
        if self.should_close:
            self.file_handle.close()

    def set_error_bad_lines(self, int status):
        self.parser.error_bad_lines = status

    cdef _make_skiprow_set(self):
        if isinstance(self.skiprows, (int, np.integer)):
            self.skiprows = range(self.skiprows)

        for i in self.skiprows:
            parser_add_skiprow(self.parser, i)

    cdef _setup_parser_source(self, source):
        cdef:
            int status

        if _is_file_like(source):
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

    cdef _get_header(self):
        cdef:
            size_t i, start, data_line, field_count, passed_count
            char *word
            object name
            int status

        header = []

        if self.parser.header >= 0:
            # Header is in the file

            if self.parser.lines < self.parser.header + 1:
                self._tokenize_rows(self.parser.header + 2)

            # e.g., if header=3 and file only has 2 lines
            if self.parser.lines < self.parser.header + 1:
                raise CParserError('Passed header=%d but only %d lines in file'
                                   % (self.parser.header, self.parser.lines))

            field_count = self.parser.line_fields[self.parser.header]
            start = self.parser.line_start[self.parser.header]

            # TODO: Py3 vs. Py2
            counts = {}
            for i in range(field_count):
                word = self.parser.words[start + i]
                name = PyString_FromString(word)

                if name == '':
                    name = 'Unnamed: %d' % i

                count = counts.get(name, 0)
                if count > 0:
                    header.append('%s.%d' % (name, count))
                else:
                    header.append(name)
                counts[name] = count + 1

            data_line = self.parser.header + 1

        elif self.names is not None:
            # Names passed
            if self.parser.lines < 1:
                self._tokenize_rows(1)

            header = self.names
            data_line = 0
        else:
            # No header passed nor to be found in the file
            if self.parser.lines < 1:
                self._tokenize_rows(1)

            return None, self.parser.line_fields[0]

        # Corner case, not enough lines in the file
        if self.parser.lines < data_line + 1:
            field_count = len(header)
        else:
            field_count = self.parser.line_fields[data_line]
            passed_count = len(header)

            if passed_count > field_count:
                raise CParserError('Column names have %d fields, data has %d'
                                   ' fields' % (passed_count, field_count))

            self.leading_cols = field_count - passed_count

        return header, field_count

    cdef _implicit_index_count(self):
        pass

    def read(self, rows=None):
        """
        rows=None --> read all rows
        """
        cdef:
            int status

        if self.low_memory:
            # Conserve intermediate space
            columns = self._read_low_memory(rows)
        else:
            # Don't care about memory usage
            columns = self._read_high_memory(rows, 1)

        if self.as_recarray:
            self._start_clock()
            result = _to_structured_array(columns, self.header)
            self._end_clock('Conversion to structured array')

            return result
        else:
            return columns

    cdef _read_low_memory(self, rows):
        cdef:
            size_t rows_read = 0
            chunks = []

        if rows is None:
            while True:
                try:
                    chunk = self._read_high_memory(self.buffer_lines, 0)
                    if len(chunk) == 0:
                        break
                except StopIteration:
                    break
                else:
                    chunks.append(chunk)
        else:
            while rows_read < rows:
                try:
                    crows = min(self.buffer_lines, rows - rows_read)

                    chunk = self._read_high_memory(crows, 0)
                    if len(chunk) == 0:
                        break

                    rows_read += len(chunk.values()[0])
                except StopIteration:
                    break
                else:
                    chunks.append(chunk)

        parser_trim_buffers(self.parser)

        if len(chunks) == 0:
            raise StopIteration

        # destructive to chunks
        return _concatenate_chunks(chunks)

    cdef _tokenize_rows(self, size_t nrows):
        cdef int status
        with nogil:
            status = tokenize_nrows(self.parser, nrows)
        if status < 0:
            raise_parser_error('Error tokenizing data', self.parser)

    cdef _read_high_memory(self, rows, bint trim):
        cdef:
            int buffered_lines
            int irows, footer = 0

        self._start_clock()

        if rows is not None:
            irows = rows
            buffered_lines = self.parser.lines - self.parser_start
            if buffered_lines < irows:
                self._tokenize_rows(irows - buffered_lines)

            if self.skip_footer > 0:
                raise ValueError('skip_footer can only be used to read '
                                 'the whole file')
        else:
            with nogil:
                status = tokenize_all_rows(self.parser)
            if status < 0:
                raise_parser_error('Error tokenizing data', self.parser)
            footer = self.skip_footer

        if self.parser_start == self.parser.lines:
            raise StopIteration
        self._end_clock('Tokenization')


        self._start_clock()
        columns = self._convert_column_data(rows=rows,
                                            footer=footer,
                                            upcast_na=not self.as_recarray)
        self._end_clock('Type conversion')

        self._start_clock()
        if len(columns) > 0:
            rows_read = len(columns.values()[0])
            # trim
            parser_consume_rows(self.parser, rows_read)
            if trim:
                parser_trim_buffers(self.parser)
            self.parser_start -= rows_read

        self._end_clock('Parser memory cleanup')

        return columns

    def debug_print(self):
        debug_print_parser(self.parser)

    cdef _start_clock(self):
        self.clocks.append(time.time())

    cdef _end_clock(self, what):
        if self.verbose:
            elapsed = time.time() - self.clocks.pop(-1)
            print '%s took: %.2f ms' % (what, elapsed * 1000)

    def set_noconvert(self, i):
        self.noconvert.add(i)

    def remove_noconvert(self, i):
        self.noconvert.remove(i)

    def _convert_column_data(self, rows=None, upcast_na=False, footer=0):
        cdef:
            Py_ssize_t i, ncols
            cast_func func
            kh_str_t *na_hashset = NULL
            int start, end
            bint na_filter = 0

        start = self.parser_start

        if rows is None:
            end = self.parser.lines
        else:
            end = min(start + rows, self.parser.lines)

        # # skip footer
        # if footer > 0:
        #     end -= footer

        results = {}
        for i in range(self.table_width):
            name = self._get_column_name(i)
            conv = self._get_converter(i, name)

            # XXX
            if self.na_filter:
                na_list = self._get_na_list(i, name)
                if na_list is None:
                    na_filter = 0
                else:
                    na_filter = 1
                    na_hashset = kset_from_list(na_list)
            else:
                na_filter = 0

            if conv:
                results[i] = _apply_converter(conv, self.parser, i, start, end)
                continue

            if i in self.noconvert:
                func = _string_box_factorize
                col_res, na_count = func(self.parser, i, start, end,
                                         na_filter, na_hashset)
                results[i] = col_res
            else:
                col_res = None
                for func in cast_func_order:
                    col_res, na_count = func(self.parser, i, start, end,
                                             na_filter, na_hashset)
                    if col_res is not None:
                        results[i] = col_res
                        break

            if na_filter:
                self._free_na_set(na_hashset)

            if upcast_na and na_count > 0:
                col_res = _maybe_upcast(col_res)

            if issubclass(col_res.dtype.type, np.integer) and self.compact_ints:
                col_res = downcast_int64(col_res, self.use_unsigned)

            if col_res is None:
                raise Exception('Unable to parse column %d' % i)

            results[i] = col_res

        self.parser_start += end - start

        return results

    def _get_converter(self, i, name):
        if self.converters is None:
            return None

        if name is not None and name in self.converters:
            return self.converters[name]

        # Converter for position, if any
        return self.converters.get(i)

    cdef _get_na_list(self, i, name):
        if self.na_values is None:
            return None

        if isinstance(self.na_values, dict):
            values = None
            if name is not None and name in self.na_values:
                values = self.na_values[name]
                if values is not None and not isinstance(values, list):
                    values = list(values)
            else:
                if i in self.na_values:
                    return self.na_values[i]
                else:
                    return _NA_VALUES

            return values
        else:
            if not isinstance(self.na_values, list):
                self.na_values = list(self.na_values)

            return self.na_values

    cdef _free_na_set(self, kh_str_t *table):
        kh_destroy_str(table)

    cdef _get_column_name(self, i):
        if self.header is not None:
            return self.header[i - self.leading_cols]
        else:
            return None

class CParserError(Exception):
    pass

def _is_file_like(obj):
    if PY3:
        import io
        if isinstance(obj, io.TextIOWrapper):
            raise CParserError('Cannot handle open unicode files (yet)')

        # BufferedReader is a byte reader for Python 3
        file = io.BufferedReader
    else:
        import __builtin__
        file = __builtin__.file

    return isinstance(obj, (basestring, file))


def _maybe_upcast(arr):
    """

    """
    if issubclass(arr.dtype.type, np.integer):
        na_value = na_values[arr.dtype]
        arr = arr.astype(float)
        np.putmask(arr, arr == na_value, np.nan)

    return arr

# ----------------------------------------------------------------------
# Type conversions / inference support code

ctypedef object (*cast_func)(parser_t *parser, int col,
                             int line_start, int line_end,
                             bint na_filter, kh_str_t *na_hashset)

cdef _string_box_factorize(parser_t *parser, int col,
                           int line_start, int line_end,
                           bint na_filter, kh_str_t *na_hashset):
    cdef:
        int error, na_count = 0
        Py_ssize_t i
        size_t lines
        coliter_t it
        char *word
        ndarray[object] result

        int ret = 0
        kh_strbox_t *table

        object pyval

        object NA = na_values[np.object_]
        khiter_t k

    table = kh_init_strbox()
    lines = line_end - line_start
    result = np.empty(lines, dtype=np.object_)
    coliter_setup(&it, parser, col, line_start)

    for i in range(lines):
        word = COLITER_NEXT(it)

        if na_filter:
            k = kh_get_str(na_hashset, word)
            # in the hash table
            if k != na_hashset.n_buckets:
                na_count += 1
                result[i] = NA
                continue

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

    kh_destroy_strbox(table)

    return result, na_count



cdef _try_double(parser_t *parser, int col, int line_start, int line_end,
                 bint na_filter, kh_str_t *na_hashset):
    cdef:
        int error, na_count = 0
        size_t i, lines
        coliter_t it
        char *word
        double *data
        double NA = na_values[np.float64]
        ndarray result
        khiter_t k

    lines = line_end - line_start
    result = np.empty(lines, dtype=np.float64)
    data = <double *> result.data
    coliter_setup(&it, parser, col, line_start)

    if na_filter:
        for i in range(lines):
            word = COLITER_NEXT(it)

            k = kh_get_str(na_hashset, word)
            # in the hash table
            if k != na_hashset.n_buckets:
                na_count += 1
                data[0] = NA
            else:
                error = to_double(word, data, parser.sci, parser.decimal)
                if error != 1:
                    return None, None
            data += 1
    else:
        for i in range(lines):
            word = COLITER_NEXT(it)
            error = to_double(word, data, parser.sci, parser.decimal)
            if error != 1:
                return None, None
            data += 1

    return result, na_count

cdef _try_int64(parser_t *parser, int col, int line_start, int line_end,
                bint na_filter, kh_str_t *na_hashset):
    cdef:
        int error, na_count = 0
        size_t i, lines
        coliter_t it
        char *word
        int64_t *data
        ndarray result

        int64_t NA = na_values[np.int64]
        khiter_t k

    lines = line_end - line_start
    result = np.empty(lines, dtype=np.int64)
    data = <int64_t *> result.data
    coliter_setup(&it, parser, col, line_start)

    if na_filter:
        for i in range(lines):
            word = COLITER_NEXT(it)

            k = kh_get_str(na_hashset, word)
            # in the hash table
            if k != na_hashset.n_buckets:
                na_count += 1
                data[i] = NA
                continue

            data[i] = str_to_int64(word, INT64_MIN, INT64_MAX,
                                   &error, parser.thousands);
            if error != 0:
                return None, None
    else:
        for i in range(lines):
            word = COLITER_NEXT(it)
            data[i] = str_to_int64(word, INT64_MIN, INT64_MAX,
                                   &error, parser.thousands);
            if error != 0:
                return None, None

    return result, na_count

cdef _try_bool(parser_t *parser, int col, int line_start, int line_end,
               bint na_filter, kh_str_t *na_hashset):
    cdef:
        int error, na_count = 0
        size_t i, lines
        coliter_t it
        char *word
        uint8_t *data
        ndarray result

        uint8_t NA = na_values[np.bool_]
        khiter_t k

    lines = line_end - line_start
    result = np.empty(lines, dtype=np.uint8)
    data = <uint8_t *> result.data
    coliter_setup(&it, parser, col, line_start)

    if na_filter:
        for i in range(lines):
            word = COLITER_NEXT(it)

            k = kh_get_str(na_hashset, word)
            # in the hash table
            if k != na_hashset.n_buckets:
                na_count += 1
                data[i] = NA
                continue

            error = to_boolean(word, data)
            if error != 0:
                return None, None
            data += 1
    else:
        for i in range(lines):
            word = COLITER_NEXT(it)

            error = to_boolean(word, data)
            if error != 0:
                return None, None
            data += 1

    if na_count > 0:
        return result, na_count
    else:
        return result.view(np.bool_), na_count

cdef _get_na_mask(parser_t *parser, int col, int line_start, int line_end,
                  kh_str_t *na_hashset):
    cdef:
        int error
        Py_ssize_t i
        size_t lines
        coliter_t it
        char *word
        ndarray[uint8_t, cast=True] result
        khiter_t k

    lines = line_end - line_start
    result = np.empty(lines, dtype=np.bool_)

    coliter_setup(&it, parser, col, line_start)
    for i in range(lines):
        word = COLITER_NEXT(it)

        k = kh_get_str(na_hashset, word)
        # in the hash table
        if k != na_hashset.n_buckets:
            result[i] = 1
        else:
            result[i] = 0

    return result

cdef kh_str_t* kset_from_list(list values) except NULL:
    # caller takes responsibility for freeing the hash table
    cdef:
        Py_ssize_t i
        khiter_t k
        kh_str_t *table
        int ret = 0

        object val

    table = kh_init_str()

    for i in range(len(values)):
        val = values[i]

        # None creeps in sometimes, which isn't possible here
        if not PyString_Check(val):
            raise TypeError('must be string, was %s' % type(val))

        k = kh_put_str(table, PyString_AsString(val), &ret)

    return table


# if at first you don't succeed...

cdef cast_func cast_func_order[4]
cast_func_order[0] = _try_int64
cast_func_order[1] = _try_double
cast_func_order[2] = _try_bool
cast_func_order[3] = _string_box_factorize

cdef raise_parser_error(object base, parser_t *parser):
    message = '%s. C error: ' % base
    if parser.error_msg != NULL:
        message += parser.error_msg
    else:
        message += 'no error message set'

    raise CParserError(message)

def downcast_int64(ndarray[int64_t] arr, bint use_unsigned=0):
    cdef:
        Py_ssize_t i, n = len(arr)
        int64_t mx = INT64_MIN + 1, mn = INT64_MAX
        int64_t NA = na_values[np.int64]
        int64_t val
        ndarray[uint8_t] mask
        int na_count = 0

    _mask = np.empty(n, dtype=bool)
    mask = _mask.view(np.uint8)

    for i in range(n):
        val = arr[i]

        if val == NA:
            mask[i] = 1
            na_count += 1
            continue

        # not NA
        mask[i] = 0

        if val > mx:
            mx = val

        if val < mn:
            mn = val

    if mn >= 0 and use_unsigned:
        if mx <= UINT8_MAX - 1:
            result = arr.astype(np.uint8)
            if na_count:
                np.putmask(result, _mask, na_values[np.uint8])
            return result

        if mx <= UINT16_MAX - 1:
            result = arr.astype(np.uint16)
            if na_count:
                np.putmask(result, _mask, na_values[np.uint16])
            return result

        if mx <= UINT32_MAX - 1:
            result = arr.astype(np.uint32)
            if na_count:
                np.putmask(result, _mask, na_values[np.uint32])
            return result

    else:
        if mn >= INT8_MIN + 1 and mx <= INT8_MAX:
            result = arr.astype(np.int8)
            if na_count:
                np.putmask(result, _mask, na_values[np.int8])
            return result

        if mn >= INT16_MIN + 1 and mx <= INT16_MAX:
            result = arr.astype(np.int16)
            if na_count:
                np.putmask(result, _mask, na_values[np.int16])
            return result

        if mn >= INT32_MIN + 1 and mx <= INT32_MAX:
            result = arr.astype(np.int32)
            if na_count:
                np.putmask(result, _mask, na_values[np.int32])
            return result

    return arr


def _concatenate_chunks(list chunks):
    cdef:
        Py_ssize_t i, j, ncols = len(chunks[0])

    result = {}
    for i in range(ncols):
        arrs = [chunk.pop(i) for chunk in chunks]
        result[i] = np.concatenate(arrs)
    return result

#----------------------------------------------------------------------

# NA values

na_values = {
    np.float64 : np.nan,
    np.int64 : INT64_MIN,
    np.int32 : INT32_MIN,
    np.bool_ : UINT8_MAX,
    np.object_ : np.nan    # oof
}

for k in list(na_values):
    na_values[np.dtype(k)] = na_values[k]


cdef _apply_converter(object f, parser_t *parser, int col,
                       int line_start, int line_end):
    cdef:
        int error
        Py_ssize_t i
        size_t lines
        coliter_t it
        char *word
        ndarray[object] result
        object val

    lines = line_end - line_start
    result = np.empty(lines, dtype=np.object_)

    coliter_setup(&it, parser, col, line_start)
    for i in range(lines):
        word = COLITER_NEXT(it)
        val = PyString_FromString(word)
        result[i] = f(val)

    values = lib.maybe_convert_objects(result)

    if issubclass(values.dtype.type, (np.number, np.bool_)):
        return values

    # XXX
    na_values = set([''])
    try:
        return lib.maybe_convert_numeric(values, na_values, False)
    except Exception:
        na_count = lib.sanitize_objects(values, na_values, False)
        return result

def _to_structured_array(dict columns, object names):
    cdef:
        ndarray recs, column
        cnp.dtype dt
        dict fields

        object name, fnames, field_type
        Py_ssize_t i, offset, nfields, length
        int stride, elsize
        char *buf

    if names is None:
        names = ['%d' % i for i in range(len(columns))]

    dt = np.dtype([(str(name), columns[i].dtype)
                   for i, name in enumerate(names)])
    fnames = dt.names
    fields = dt.fields

    nfields = len(fields)

    length = len(columns.values()[0])
    stride = dt.itemsize

    # start = time.time()

    # we own the data
    buf = <char*> malloc(length * stride)

    recs = util.sarr_from_data(dt, length, buf)
    assert(recs.flags.owndata)

    # buf = <char*> recs.data
    # end = time.time()
    # print 'took %.4f' % (end - start)

    for i in range(nfields):
        # start = time.clock()
        # name = names[i]

        # XXX
        field_type = fields[fnames[i]]

        # (dtype, stride) tuple
        offset = field_type[1]
        elsize = field_type[0].itemsize
        column = columns[i]

        _fill_structured_column(buf + offset, <char*> column.data,
                                elsize, stride, length,
                                field_type[0] == np.object_)

        # print 'Transfer of %s took %.4f' % (str(field_type),
        #                                     time.clock() - start)

    return recs

cdef _fill_structured_column(char *dst, char* src, int elsize,
                             int stride, int length, bint incref):
    cdef:
        size_t i

    if incref:
        util.transfer_object_column(dst, src, stride, length)
    else:
        for i in range(length):
            memcpy(dst, src, elsize)
            dst += stride
            src += elsize

