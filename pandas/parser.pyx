# Copyright (c) 2012, Lambda Foundry, Inc.
# See LICENSE for the license

from libc.stdio cimport fopen, fclose
from libc.stdlib cimport malloc, free
from libc.string cimport strncpy, strlen, strcmp, strcasecmp
cimport libc.stdio as stdio
import warnings

from csv import QUOTE_MINIMAL, QUOTE_NONNUMERIC, QUOTE_NONE
from cpython cimport (PyObject, PyBytes_FromString,
                      PyBytes_AsString, PyBytes_Check,
                      PyUnicode_Check, PyUnicode_AsUTF8String,
                      PyErr_Occurred, PyErr_Fetch)
from cpython.ref cimport PyObject, Py_XDECREF
from io.common import CParserError, DtypeWarning, EmptyDataError


cdef extern from "Python.h":
    object PyUnicode_FromString(char *v)

    object PyUnicode_Decode(char *v, Py_ssize_t size, char *encoding,
                            char *errors)

cdef extern from "stdlib.h":
    void memcpy(void *dst, void *src, size_t n)

cimport cython
cimport numpy as cnp

from numpy cimport ndarray, uint8_t, uint64_t

import numpy as np
cimport util

import pandas.lib as lib
import pandas.compat as compat
from pandas.types.common import (is_categorical_dtype, CategoricalDtype,
                                 is_integer_dtype, is_float_dtype,
                                 is_bool_dtype, is_object_dtype,
                                 is_string_dtype, is_datetime64_dtype,
                                 pandas_dtype)
from pandas.core.categorical import Categorical
from pandas.core.algorithms import take_1d
from pandas.types.concat import union_categoricals
from pandas import Index

import time
import os

cnp.import_array()

from khash cimport *

import sys

cdef bint PY3 = (sys.version_info[0] >= 3)

cdef double INF = <double> np.inf
cdef double NEGINF = -INF

cdef extern from "headers/stdint.h":
    enum: UINT8_MAX
    enum: UINT16_MAX
    enum: UINT32_MAX
    enum: UINT64_MAX
    enum: INT8_MIN
    enum: INT8_MAX
    enum: INT16_MIN
    enum: INT16_MAX
    enum: INT32_MAX
    enum: INT32_MIN
    enum: INT64_MAX
    enum: INT64_MIN

cdef extern from "headers/portable.h":
    pass

cdef extern from "errno.h":
    int errno

try:
    basestring
except NameError:
    basestring = str

cdef extern from "parser/tokenizer.h":

    ctypedef enum ParserState:
        START_RECORD
        START_FIELD
        ESCAPED_CHAR
        IN_FIELD
        IN_QUOTED_FIELD
        ESCAPE_IN_QUOTED_FIELD
        QUOTE_IN_QUOTED_FIELD
        EAT_CRNL
        EAT_CRNL_NOP
        EAT_WHITESPACE
        EAT_COMMENT
        EAT_LINE_COMMENT
        WHITESPACE_LINE
        SKIP_LINE
        FINISHED

    enum: ERROR_OVERFLOW

    ctypedef void* (*io_callback)(void *src, size_t nbytes, size_t *bytes_read,
                                  int *status)
    ctypedef int (*io_cleanup)(void *src)

    ctypedef struct parser_t:
        void *source
        io_callback cb_io
        io_cleanup cb_cleanup

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
        char lineterminator
        int skipinitialspace       # ignore spaces following delimiter? */
        int quoting                # style of quoting to write */

        # hmm =/
#        int numeric_field

        char commentchar
        int allow_embedded_newline
        int strict                 # raise exception on bad CSV */

        int usecols

        int expected_fields
        int error_bad_lines
        int warn_bad_lines

        # floating point options
        char decimal
        char sci

        # thousands separator (comma, period)
        char thousands

        int header # Boolean: 1: has header, 0: no header
        int header_start # header row start
        int header_end # header row end

        void *skipset
        int64_t skip_first_N_rows
        int skipfooter
        double (*converter)(const char *, char **, char, char, char, int) nogil

        #  error handling
        char *warn_msg
        char *error_msg

        int skip_empty_lines

    ctypedef struct coliter_t:
        char **words
        int *line_start
        int col

    void coliter_setup(coliter_t *it, parser_t *parser, int i, int start) nogil
    void COLITER_NEXT(coliter_t, const char *) nogil

    parser_t* parser_new()

    int parser_init(parser_t *self) nogil
    void parser_free(parser_t *self) nogil
    int parser_add_skiprow(parser_t *self, int64_t row)

    int parser_set_skipfirstnrows(parser_t *self, int64_t nrows)

    void parser_set_default_options(parser_t *self)

    int parser_consume_rows(parser_t *self, size_t nrows)

    int parser_trim_buffers(parser_t *self)

    void debug_print_parser(parser_t *self)

    int tokenize_all_rows(parser_t *self) nogil
    int tokenize_nrows(parser_t *self, size_t nrows) nogil

    int64_t str_to_int64(char *p_item, int64_t int_min,
                         int64_t int_max, int *error, char tsep) nogil
#    uint64_t str_to_uint64(char *p_item, uint64_t uint_max, int *error)

    double xstrtod(const char *p, char **q, char decimal, char sci,
                   char tsep, int skip_trailing) nogil
    double precise_xstrtod(const char *p, char **q, char decimal, char sci,
                           char tsep, int skip_trailing) nogil
    double round_trip(const char *p, char **q, char decimal, char sci,
                      char tsep, int skip_trailing) nogil

#    inline int to_complex(char *item, double *p_real,
#                          double *p_imag, char sci, char decimal)
    inline int to_longlong(char *item, long long *p_value) nogil
#    inline int to_longlong_thousands(char *item, long long *p_value,
#                                     char tsep)
    int to_boolean(const char *item, uint8_t *val) nogil


cdef extern from "parser/io.h":
    void *new_mmap(char *fname)
    int del_mmap(void *src)
    void* buffer_mmap_bytes(void *source, size_t nbytes,
                            size_t *bytes_read, int *status)

    void *new_file_source(char *fname, size_t buffer_size)

    void *new_rd_source(object obj)

    int del_file_source(void *src)
    int del_rd_source(void *src)

    void* buffer_file_bytes(void *source, size_t nbytes,
                            size_t *bytes_read, int *status)

    void* buffer_rd_bytes(void *source, size_t nbytes,
                          size_t *bytes_read, int *status)


DEFAULT_CHUNKSIZE = 256 * 1024

# common NA values
# no longer excluding inf representations
# '1.#INF','-1.#INF', '1.#INF000000',
_NA_VALUES = [b'-1.#IND', b'1.#QNAN', b'1.#IND', b'-1.#QNAN',
              b'#N/A N/A', b'NA', b'#NA', b'NULL', b'NaN',
              b'nan', b'']


cdef class TextReader:
    """

    # source: StringIO or file object

    """

    cdef:
        parser_t *parser
        object file_handle, na_fvalues
        object true_values, false_values
        object handle
        bint na_filter, verbose, has_usecols, has_mi_columns
        int parser_start
        list clocks
        char *c_encoding
        kh_str_t *false_set
        kh_str_t *true_set

    cdef public:
        int leading_cols, table_width, skipfooter, buffer_lines
        object allow_leading_cols
        object delimiter, converters, delim_whitespace
        object na_values
        object memory_map
        object as_recarray
        object header, orig_header, names, header_start, header_end
        object index_col
        object low_memory
        object skiprows
        object compact_ints, use_unsigned
        object dtype
        object encoding
        object compression
        object mangle_dupe_cols
        object tupleize_cols
        list dtype_cast_order
        set noconvert, usecols

    def __cinit__(self, source,
                  delimiter=b',',

                  header=0,
                  header_start=0,
                  header_end=0,
                  index_col=None,
                  names=None,

                  memory_map=False,
                  tokenize_chunksize=DEFAULT_CHUNKSIZE,
                  delim_whitespace=False,

                  compression=None,

                  converters=None,

                  as_recarray=False,

                  skipinitialspace=False,
                  escapechar=None,
                  doublequote=True,
                  quotechar=b'"',
                  quoting=0,
                  lineterminator=None,

                  encoding=None,

                  comment=None,
                  decimal=b'.',
                  thousands=None,

                  dtype=None,
                  usecols=None,
                  error_bad_lines=True,
                  warn_bad_lines=True,

                  na_filter=True,
                  na_values=None,
                  na_fvalues=None,
                  true_values=None,
                  false_values=None,

                  compact_ints=False,
                  allow_leading_cols=True,
                  use_unsigned=False,
                  low_memory=False,
                  buffer_lines=None,
                  skiprows=None,
                  skipfooter=0,
                  verbose=False,
                  mangle_dupe_cols=True,
                  tupleize_cols=False,
                  float_precision=None,
                  skip_blank_lines=True):

        self.parser = parser_new()
        self.parser.chunksize = tokenize_chunksize

        self.mangle_dupe_cols=mangle_dupe_cols
        self.tupleize_cols=tupleize_cols

        # For timekeeping
        self.clocks = []

        self.compression = compression
        self.memory_map = memory_map

        self.parser.usecols = (usecols is not None)

        self._setup_parser_source(source)
        parser_set_default_options(self.parser)

        parser_init(self.parser)

        if delim_whitespace:
            self.parser.delim_whitespace = delim_whitespace
        else:
            if len(delimiter) > 1:
                raise ValueError('only length-1 separators excluded right now')
            self.parser.delimiter = ord(delimiter)

        #----------------------------------------
        # parser options

        self.parser.doublequote = doublequote
        self.parser.skipinitialspace = skipinitialspace
        self.parser.skip_empty_lines = skip_blank_lines

        if lineterminator is not None:
            if len(lineterminator) != 1:
                raise ValueError('Only length-1 line terminators supported')
            self.parser.lineterminator = ord(lineterminator)

        if len(decimal) != 1:
            raise ValueError('Only length-1 decimal markers supported')
        self.parser.decimal = ord(decimal)

        if thousands is not None:
            if len(thousands) != 1:
                raise ValueError('Only length-1 thousands markers supported')
            self.parser.thousands = ord(thousands)

        if escapechar is not None:
            if len(escapechar) != 1:
                raise ValueError('Only length-1 escapes  supported')
            self.parser.escapechar = ord(escapechar)

        self._set_quoting(quotechar, quoting)

        dtype_order = ['int64', 'float64', 'bool', 'object']
        if quoting == QUOTE_NONNUMERIC:
            # consistent with csv module semantics, cast all to float
            dtype_order = dtype_order[1:]
        self.dtype_cast_order = [np.dtype(x) for x in dtype_order]

        if comment is not None:
            if len(comment) > 1:
                raise ValueError('Only length-1 comment characters supported')
            self.parser.commentchar = ord(comment)

        # error handling of bad lines
        self.parser.error_bad_lines = int(error_bad_lines)
        self.parser.warn_bad_lines = int(warn_bad_lines)

        self.skiprows = skiprows
        if skiprows is not None:
            self._make_skiprow_set()

        self.skipfooter = skipfooter

        # suboptimal
        if usecols is not None:
            self.has_usecols = 1
            self.usecols = set(usecols)

        # XXX
        if skipfooter > 0:
            self.parser.error_bad_lines = 0
            self.parser.warn_bad_lines = 0

        self.delimiter = delimiter
        self.delim_whitespace = delim_whitespace

        self.na_values = na_values
        if na_fvalues is None:
            na_fvalues = set()
        self.na_fvalues = na_fvalues

        self.true_values = _maybe_encode(true_values) + _true_values
        self.false_values = _maybe_encode(false_values) + _false_values

        self.true_set = kset_from_list(self.true_values)
        self.false_set = kset_from_list(self.false_values)

        self.converters = converters

        self.na_filter = na_filter
        self.as_recarray = as_recarray

        self.compact_ints = compact_ints
        self.use_unsigned = use_unsigned

        self.verbose = verbose
        self.low_memory = low_memory
        self.parser.converter = xstrtod
        if float_precision == 'high':
            self.parser.converter = precise_xstrtod
        elif float_precision == 'round_trip':
            self.parser.converter = round_trip

        # encoding
        if encoding is not None:
            if not isinstance(encoding, bytes):
                encoding = encoding.encode('utf-8')
            encoding = encoding.lower()
            self.c_encoding = <char*> encoding
        else:
            self.c_encoding = NULL

        self.encoding = encoding

        if isinstance(dtype, dict):
            dtype = {k: pandas_dtype(dtype[k])
                     for k in dtype}
        elif dtype is not None:
            dtype = pandas_dtype(dtype)

        self.dtype = dtype

        # XXX
        self.noconvert = set()

        self.index_col = index_col

        #----------------------------------------
        # header stuff

        self.allow_leading_cols = allow_leading_cols
        self.leading_cols = 0

        # TODO: no header vs. header is not the first row
        self.has_mi_columns = 0
        self.orig_header = header
        if header is None:
            # sentinel value
            self.parser.header_start = -1
            self.parser.header_end = -1
            self.parser.header = -1
            self.parser_start = 0
            self.header = []
        else:
            if isinstance(header, list) and len(header):
                # need to artifically skip the final line
                # which is still a header line
                header = list(header)
                header.append(header[-1] + 1)

                self.parser.header_start = header[0]
                self.parser.header_end = header[-1]
                self.parser.header = header[0]
                self.parser_start = header[-1] + 1
                self.has_mi_columns = 1
                self.header = header
            else:
                self.parser.header_start = header
                self.parser.header_end = header
                self.parser.header = header
                self.parser_start = header + 1
                self.header = [ header ]

        self.names = names
        self.header, self.table_width = self._get_header()

        if not self.table_width:
            raise EmptyDataError("No columns to parse from file")

        # compute buffer_lines as function of table width
        heuristic = 2**20 // self.table_width
        self.buffer_lines = 1
        while self.buffer_lines * 2< heuristic:
            self.buffer_lines *= 2

    def __init__(self, *args, **kwards):
        pass

    def __dealloc__(self):
        parser_free(self.parser)
        kh_destroy_str(self.true_set)
        kh_destroy_str(self.false_set)

    def close(self):
        # we need to properly close an open derived
        # filehandle here, e.g. and UTFRecoder
        if self.handle is not None:
            try:
                self.handle.close()
            except:
                pass

    def set_error_bad_lines(self, int status):
        self.parser.error_bad_lines = status

    def _set_quoting(self, quote_char, quoting):
        if not isinstance(quoting, int):
            raise TypeError('"quoting" must be an integer')

        if not QUOTE_MINIMAL <= quoting <= QUOTE_NONE:
            raise TypeError('bad "quoting" value')

        if not isinstance(quote_char, (str, compat.text_type,
                                       bytes)) and quote_char is not None:
            dtype = type(quote_char).__name__
            raise TypeError('"quotechar" must be string, '
                            'not {dtype}'.format(dtype=dtype))

        if quote_char is None or quote_char == '':
            if quoting != QUOTE_NONE:
                raise TypeError("quotechar must be set if quoting enabled")
            self.parser.quoting = quoting
            self.parser.quotechar = -1
        elif len(quote_char) > 1:  # 0-len case handled earlier
            raise TypeError('"quotechar" must be a 1-character string')
        else:
            self.parser.quoting = quoting
            self.parser.quotechar = ord(quote_char)

    cdef _make_skiprow_set(self):
        if isinstance(self.skiprows, (int, np.integer)):
            parser_set_skipfirstnrows(self.parser, self.skiprows)
        else:
            for i in self.skiprows:
                parser_add_skiprow(self.parser, i)

    cdef _setup_parser_source(self, source):
        cdef:
            int status
            void *ptr

        self.parser.cb_io = NULL
        self.parser.cb_cleanup = NULL

        if self.compression:
            if self.compression == 'gzip':
                import gzip
                if isinstance(source, basestring):
                    source = gzip.GzipFile(source, 'rb')
                else:
                    source = gzip.GzipFile(fileobj=source)
            elif self.compression == 'bz2':
                import bz2
                if isinstance(source, basestring) or PY3:
                    source = bz2.BZ2File(source, 'rb')
                else:
                    raise ValueError('Python 2 cannot read bz2 from open file '
                                     'handle')
            elif self.compression == 'zip':
                import zipfile
                zip_file = zipfile.ZipFile(source)
                zip_names = zip_file.namelist()

                if len(zip_names) == 1:
                    file_name = zip_names.pop()
                    source = zip_file.open(file_name)

                elif len(zip_names) == 0:
                    raise ValueError('Zero files found in compressed '
                                     'zip file %s', source)
                else:
                    raise ValueError('Multiple files found in compressed '
                                     'zip file %s', str(zip_names))
            elif self.compression == 'xz':
                lzma = compat.import_lzma()

                if isinstance(source, basestring):
                    source = lzma.LZMAFile(source, 'rb')
                else:
                    source = lzma.LZMAFile(filename=source)
            else:
                raise ValueError('Unrecognized compression type: %s' %
                                 self.compression)
            self.handle = source

        if isinstance(source, basestring):
            if not isinstance(source, bytes):
                source = source.encode(sys.getfilesystemencoding() or 'utf-8')

            if self.memory_map:
                ptr = new_mmap(source)
                if ptr == NULL:
                    # fall back
                    ptr = new_file_source(source, self.parser.chunksize)
                    self.parser.cb_io = &buffer_file_bytes
                    self.parser.cb_cleanup = &del_file_source
                else:
                    self.parser.cb_io = &buffer_mmap_bytes
                    self.parser.cb_cleanup = &del_mmap
            else:
                ptr = new_file_source(source, self.parser.chunksize)
                self.parser.cb_io = &buffer_file_bytes
                self.parser.cb_cleanup = &del_file_source

            if ptr == NULL:
                if not os.path.exists(source):
                    raise compat.FileNotFoundError(
                        'File %s does not exist' % source)
                raise IOError('Initializing from file failed')

            self.parser.source = ptr

        elif hasattr(source, 'read'):
            # e.g., StringIO

            ptr = new_rd_source(source)
            if ptr == NULL:
                raise IOError('Initializing parser from file-like '
                              'object failed')

            self.parser.source = ptr
            self.parser.cb_io = &buffer_rd_bytes
            self.parser.cb_cleanup = &del_rd_source
        else:
            raise IOError('Expected file path name or file-like object,'
                          ' got %s type' % type(source))

    cdef _get_header(self):
        # header is now a list of lists, so field_count should use header[0]

        cdef:
            size_t i, start, data_line, field_count, passed_count, hr, unnamed_count  # noqa
            char *word
            object name
            int status
            Py_ssize_t size
            char *errors = "strict"
            cdef StringPath path = _string_path(self.c_encoding)

        header = []

        if self.parser.header_start >= 0:

            # Header is in the file
            for level, hr in enumerate(self.header):

                this_header = []

                if self.parser.lines < hr + 1:
                    self._tokenize_rows(hr + 2)

                if self.parser.lines == 0:
                    field_count = 0
                    start = self.parser.line_start[0]

                # e.g., if header=3 and file only has 2 lines
                elif (self.parser.lines < hr + 1
                      and not isinstance(self.orig_header, list)) or (
                          self.parser.lines < hr):
                    msg = self.orig_header
                    if isinstance(msg, list):
                        msg = "[%s], len of %d," % (
                            ','.join([ str(m) for m in msg ]), len(msg))
                    raise CParserError(
                        'Passed header=%s but only %d lines in file'
                        % (msg, self.parser.lines))

                else:
                    field_count = self.parser.line_fields[hr]
                    start = self.parser.line_start[hr]

                counts = {}
                unnamed_count = 0
                for i in range(field_count):
                    word = self.parser.words[start + i]

                    if path == CSTRING:
                        name = PyBytes_FromString(word)
                    elif path == UTF8:
                        name = PyUnicode_FromString(word)
                    elif path == ENCODED:
                        name = PyUnicode_Decode(word, strlen(word),
                                                self.c_encoding, errors)

                    if name == '':
                        if self.has_mi_columns:
                            name = 'Unnamed: %d_level_%d' % (i, level)
                        else:
                            name = 'Unnamed: %d' % i
                        unnamed_count += 1

                    count = counts.get(name, 0)
                    if (count > 0 and self.mangle_dupe_cols
                        and not self.has_mi_columns):
                        this_header.append('%s.%d' % (name, count))
                    else:
                        this_header.append(name)
                    counts[name] = count + 1

                if self.has_mi_columns:

                    # If we have grabbed an extra line, but it's not in our
                    # format, save in the buffer, and create an blank extra
                    # line for the rest of the parsing code.
                    if hr == self.header[-1]:
                        lc = len(this_header)
                        ic = (len(self.index_col) if self.index_col
                                                     is not None else 0)
                        if lc != unnamed_count and lc - ic > unnamed_count:
                            hr -= 1
                            self.parser_start -= 1
                            this_header = [None] * lc

                data_line = hr + 1
                header.append(this_header)

            if self.names is not None:
                header = [ self.names ]

        elif self.names is not None:
            # Enforce this unless usecols
            if not self.has_usecols:
                self.parser.expected_fields = len(self.names)

            # Names passed
            if self.parser.lines < 1:
                self._tokenize_rows(1)

            header = [ self.names ]
            data_line = 0

            if self.parser.lines < 1:
                field_count = len(header[0])
            else:
                field_count = self.parser.line_fields[data_line]
        else:
            # No header passed nor to be found in the file
            if self.parser.lines < 1:
                self._tokenize_rows(1)

            return None, self.parser.line_fields[0]

        # Corner case, not enough lines in the file
        if self.parser.lines < data_line + 1:
            field_count = len(header[0])
        else: # not self.has_usecols:

            field_count = self.parser.line_fields[data_line]

            # #2981
            if self.names is not None:
                field_count = max(field_count, len(self.names))

            passed_count = len(header[0])

            # if passed_count > field_count:
            #     raise CParserError('Column names have %d fields, '
            #                        'data has %d fields'
            #                        % (passed_count, field_count))

            if self.has_usecols and self.allow_leading_cols:
                nuse = len(self.usecols)
                if nuse == passed_count:
                    self.leading_cols = 0
                elif self.names is None and nuse < passed_count:
                    self.leading_cols = field_count - passed_count
                elif passed_count != field_count:
                    raise ValueError('Passed header names '
                                     'mismatches usecols')
            # oh boy, #2442, #2981
            elif self.allow_leading_cols and passed_count < field_count:
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
            columns = self._read_rows(rows, 1)

        if self.as_recarray:
            self._start_clock()
            result = _to_structured_array(columns, self.header, self.usecols)
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
                    chunk = self._read_rows(self.buffer_lines, 0)
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

                    chunk = self._read_rows(crows, 0)
                    if len(chunk) == 0:
                        break

                    rows_read += len(list(chunk.values())[0])
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

        if self.parser.warn_msg != NULL:
            print >> sys.stderr, self.parser.warn_msg
            free(self.parser.warn_msg)
            self.parser.warn_msg = NULL

        if status < 0:
            raise_parser_error('Error tokenizing data', self.parser)

    cdef _read_rows(self, rows, bint trim):
        cdef:
            int buffered_lines
            int irows, footer = 0

        self._start_clock()

        if rows is not None:
            irows = rows
            buffered_lines = self.parser.lines - self.parser_start
            if buffered_lines < irows:
                self._tokenize_rows(irows - buffered_lines)

            if self.skipfooter > 0:
                raise ValueError('skipfooter can only be used to read '
                                 'the whole file')
        else:
            with nogil:
                status = tokenize_all_rows(self.parser)

            if self.parser.warn_msg != NULL:
                print >> sys.stderr, self.parser.warn_msg
                free(self.parser.warn_msg)
                self.parser.warn_msg = NULL

            if status < 0:
                raise_parser_error('Error tokenizing data', self.parser)
            footer = self.skipfooter

        if self.parser_start >= self.parser.lines:
            raise StopIteration
        self._end_clock('Tokenization')

        self._start_clock()
        columns = self._convert_column_data(rows=rows,
                                            footer=footer,
                                            upcast_na=not self.as_recarray)
        self._end_clock('Type conversion')

        self._start_clock()
        if len(columns) > 0:
            rows_read = len(list(columns.values())[0])
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
            Py_ssize_t i, nused
            kh_str_t *na_hashset = NULL
            int start, end
            object name, na_flist
            bint na_filter = 0
            Py_ssize_t num_cols

        start = self.parser_start

        if rows is None:
            end = self.parser.lines
        else:
            end = min(start + rows, self.parser.lines)

        # # skip footer
        # if footer > 0:
        #     end -= footer

        num_cols = -1
        for i in range(self.parser.lines):
            num_cols = (num_cols < self.parser.line_fields[i]) * \
                self.parser.line_fields[i] + \
                (num_cols >= self.parser.line_fields[i]) * num_cols

        if self.table_width - self.leading_cols > num_cols:
            raise CParserError(
                "Too many columns specified: expected %s and found %s" %
                (self.table_width - self.leading_cols, num_cols))

        results = {}
        nused = 0
        for i in range(self.table_width):
            if i < self.leading_cols:
                # Pass through leading columns always
                name = i
            elif self.usecols and nused == len(self.usecols):
                # Once we've gathered all requested columns, stop. GH5766
                break
            else:
                name = self._get_column_name(i, nused)
                if self.has_usecols and not (i in self.usecols or
                                             name in self.usecols):
                    continue
                nused += 1

            conv = self._get_converter(i, name)

            # XXX
            na_flist = set()
            if self.na_filter:
                na_list, na_flist = self._get_na_list(i, name)
                if na_list is None:
                    na_filter = 0
                else:
                    na_filter = 1
                    na_hashset = kset_from_list(na_list)
            else:
                na_filter = 0

            if conv:
                results[i] = _apply_converter(conv, self.parser, i, start, end,
                                              self.c_encoding)
                continue

            # Should return as the desired dtype (inferred or specified)
            col_res, na_count = self._convert_tokens(
                i, start, end, name, na_filter, na_hashset, na_flist)

            if na_filter:
                self._free_na_set(na_hashset)

            if upcast_na and na_count > 0:
                col_res = _maybe_upcast(col_res)

            if issubclass(col_res.dtype.type,
                          np.integer) and self.compact_ints:
                col_res = lib.downcast_int64(col_res, na_values,
                                             self.use_unsigned)

            if col_res is None:
                raise CParserError('Unable to parse column %d' % i)

            results[i] = col_res

        self.parser_start += end - start

        return results

    cdef inline _convert_tokens(self, Py_ssize_t i, int start, int end,
                                object name, bint na_filter,
                                kh_str_t *na_hashset,
                                object na_flist):
        cdef:
            object col_dtype = None

        if self.dtype is not None:
            if isinstance(self.dtype, dict):
                if name in self.dtype:
                    col_dtype = self.dtype[name]
                elif i in self.dtype:
                    col_dtype = self.dtype[i]
            else:
                if self.dtype.names:
                    # structured array
                    col_dtype = np.dtype(self.dtype.descr[i][1])
                else:
                    col_dtype = self.dtype

            if col_dtype is not None:
                col_res, na_count = self._convert_with_dtype(
                    col_dtype, i, start, end, na_filter,
                    1, na_hashset, na_flist)

                # Fallback on the parse (e.g. we requested int dtype,
                # but its actually a float).
                if col_res is not None:
                    return col_res, na_count

        if i in self.noconvert:
            return self._string_convert(i, start, end, na_filter, na_hashset)
        else:
            col_res = None
            for dt in self.dtype_cast_order:
                try:
                    col_res, na_count = self._convert_with_dtype(
                        dt, i, start, end, na_filter, 0, na_hashset, na_flist)
                except OverflowError:
                    col_res, na_count = self._convert_with_dtype(
                        np.dtype('object'), i, start, end, na_filter,
                        0, na_hashset, na_flist)

                if col_res is not None:
                    break

        # we had a fallback parse on the dtype, so now try to cast
        # only allow safe casts, eg. with a nan you cannot safely cast to int
        if col_res is not None and col_dtype is not None:
            try:
                col_res = col_res.astype(col_dtype, casting='safe')
            except TypeError:

                # float -> int conversions can fail the above
                # even with no nans
                col_res_orig = col_res
                col_res = col_res.astype(col_dtype)
                if (col_res != col_res_orig).any():
                    raise ValueError(
                        "cannot safely convert passed user dtype of "
                        "{col_dtype} for {col_res} dtyped data in "
                        "column {column}".format(
                            col_dtype=col_dtype,
                            col_res=col_res_orig.dtype.name,
                            column=i))

        return col_res, na_count

    cdef _convert_with_dtype(self, object dtype, Py_ssize_t i,
                             int start, int end,
                             bint na_filter,
                             bint user_dtype,
                             kh_str_t *na_hashset,
                             object na_flist):
        if is_integer_dtype(dtype):
            result, na_count = _try_int64(self.parser, i, start,
                                          end, na_filter, na_hashset)
            if user_dtype and na_count is not None:
                if na_count > 0:
                    raise ValueError("Integer column has NA values in "
                                     "column {column}".format(column=i))

            if result is not None and dtype != 'int64':
                result = result.astype(dtype)

            return result, na_count

        elif is_float_dtype(dtype):
            result, na_count = _try_double(self.parser, i, start, end,
                                           na_filter, na_hashset, na_flist)

            if result is not None and dtype != 'float64':
                result = result.astype(dtype)
            return result, na_count

        elif is_bool_dtype(dtype):
            result, na_count = _try_bool_flex(self.parser, i, start, end,
                                              na_filter, na_hashset,
                                              self.true_set, self.false_set)
            return result, na_count
        elif dtype.kind == 'S':
            # TODO: na handling
            width = dtype.itemsize
            if width > 0:
                result = _to_fw_string(self.parser, i, start, end, width)
                return result, 0

            # treat as a regular string parsing
            return self._string_convert(i, start, end, na_filter,
                                        na_hashset)
        elif dtype.kind == 'U':
            width = dtype.itemsize
            if width > 0:
                raise TypeError("the dtype %s is not "
                                "supported for parsing" % dtype)

            # unicode variable width
            return self._string_convert(i, start, end, na_filter,
                                        na_hashset)
        elif is_categorical_dtype(dtype):
            codes, cats, na_count = _categorical_convert(
                self.parser, i, start, end, na_filter,
                na_hashset, self.c_encoding)
            # sort categories and recode if necessary
            cats = Index(cats)
            if not cats.is_monotonic_increasing:
                unsorted = cats.copy()
                cats = cats.sort_values()
                indexer = cats.get_indexer(unsorted)
                codes = take_1d(indexer, codes, fill_value=-1)

            return Categorical(codes, categories=cats, ordered=False,
                               fastpath=True), na_count
        elif is_object_dtype(dtype):
            return self._string_convert(i, start, end, na_filter,
                                        na_hashset)
        elif is_datetime64_dtype(dtype):
            raise TypeError("the dtype %s is not supported "
                            "for parsing, pass this column "
                            "using parse_dates instead" % dtype)
        else:
            raise TypeError("the dtype %s is not "
                            "supported for parsing" % dtype)

    cdef _string_convert(self, Py_ssize_t i, int start, int end,
                         bint na_filter, kh_str_t *na_hashset):

        cdef StringPath path = _string_path(self.c_encoding)

        if path == UTF8:
            return _string_box_utf8(self.parser, i, start, end, na_filter,
                                    na_hashset)
        elif path == ENCODED:
            return _string_box_decode(self.parser, i, start, end,
                                      na_filter, na_hashset, self.c_encoding)
        elif path == CSTRING:
            return _string_box_factorize(self.parser, i, start, end,
                                         na_filter, na_hashset)

    def _get_converter(self, i, name):
        if self.converters is None:
            return None

        if name is not None and name in self.converters:
            return self.converters[name]

        # Converter for position, if any
        return self.converters.get(i)

    cdef _get_na_list(self, i, name):
        if self.na_values is None:
            return None, set()

        if isinstance(self.na_values, dict):
            key = None
            values = None

            if name is not None and name in self.na_values:
                key = name
            elif i in self.na_values:
                key = i
            else:  # No na_values provided for this column.
                return _NA_VALUES, set()

            values = self.na_values[key]
            if values is not None and not isinstance(values, list):
                values = list(values)

            fvalues = self.na_fvalues[key]
            if fvalues is not None and not isinstance(fvalues, set):
                fvalues = set(fvalues)

            return _ensure_encoded(values), fvalues
        else:
            if not isinstance(self.na_values, list):
                self.na_values = list(self.na_values)
            if not isinstance(self.na_fvalues, set):
                self.na_fvalues = set(self.na_fvalues)

            return _ensure_encoded(self.na_values), self.na_fvalues

    cdef _free_na_set(self, kh_str_t *table):
        kh_destroy_str(table)

    cdef _get_column_name(self, Py_ssize_t i, Py_ssize_t nused):
        if self.has_usecols and self.names is not None:
            if len(self.names) == len(self.usecols):
                return self.names[nused]
            else:
                return self.names[i - self.leading_cols]
        else:
            if self.header is not None:
                j = i - self.leading_cols
                # hack for #2442
                if j == len(self.header[0]):
                    return j
                else:
                    return self.header[0][j]
            else:
                return None

cdef object _true_values = [b'True', b'TRUE', b'true']
cdef object _false_values = [b'False', b'FALSE', b'false']


def _ensure_encoded(list lst):
    cdef list result = []
    for x in lst:
        if PyUnicode_Check(x):
            x = PyUnicode_AsUTF8String(x)
        elif not PyBytes_Check(x):
            x = asbytes(x)

        result.append(x)
    return result

cdef asbytes(object o):
    if PY3:
        return str(o).encode('utf-8')
    else:
        return str(o)


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
    elif arr.dtype == np.bool_:
        mask = arr.view(np.uint8) == na_values[np.uint8]
        arr = arr.astype(object)
        np.putmask(arr, mask, np.nan)

    return arr

cdef enum StringPath:
    CSTRING
    UTF8
    ENCODED

# factored out logic to pick string converter
cdef inline StringPath _string_path(char *encoding):
    if encoding != NULL and encoding != b"utf-8":
        return ENCODED
    elif PY3 or encoding != NULL:
        return UTF8
    else:
        return CSTRING
# ----------------------------------------------------------------------
# Type conversions / inference support code

cdef _string_box_factorize(parser_t *parser, int col,
                           int line_start, int line_end,
                           bint na_filter, kh_str_t *na_hashset):
    cdef:
        int error, na_count = 0
        Py_ssize_t i
        size_t lines
        coliter_t it
        const char *word = NULL
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
        COLITER_NEXT(it, word)

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
            pyval = PyBytes_FromString(word)

            k = kh_put_strbox(table, word, &ret)
            table.vals[k] = <PyObject*> pyval

        result[i] = pyval

    kh_destroy_strbox(table)

    return result, na_count

cdef _string_box_utf8(parser_t *parser, int col,
                      int line_start, int line_end,
                      bint na_filter, kh_str_t *na_hashset):
    cdef:
        int error, na_count = 0
        Py_ssize_t i
        size_t lines
        coliter_t it
        const char *word = NULL
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
        COLITER_NEXT(it, word)

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
            pyval = PyUnicode_FromString(word)

            k = kh_put_strbox(table, word, &ret)
            table.vals[k] = <PyObject *> pyval

        result[i] = pyval

    kh_destroy_strbox(table)

    return result, na_count

cdef _string_box_decode(parser_t *parser, int col,
                        int line_start, int line_end,
                        bint na_filter, kh_str_t *na_hashset,
                        char *encoding):
    cdef:
        int error, na_count = 0
        Py_ssize_t i, size
        size_t lines
        coliter_t it
        const char *word = NULL
        ndarray[object] result

        int ret = 0
        kh_strbox_t *table

        char *errors = "strict"

        object pyval

        object NA = na_values[np.object_]
        khiter_t k

    table = kh_init_strbox()
    lines = line_end - line_start
    result = np.empty(lines, dtype=np.object_)
    coliter_setup(&it, parser, col, line_start)

    for i in range(lines):
        COLITER_NEXT(it, word)

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
            size = strlen(word)
            pyval = PyUnicode_Decode(word, size, encoding, errors)

            k = kh_put_strbox(table, word, &ret)
            table.vals[k] = <PyObject *> pyval

        result[i] = pyval

    kh_destroy_strbox(table)

    return result, na_count


@cython.boundscheck(False)
cdef _categorical_convert(parser_t *parser, int col,
                          int line_start, int line_end,
                          bint na_filter, kh_str_t *na_hashset,
                          char *encoding):
    "Convert column data into codes, categories"
    cdef:
        int error, na_count = 0
        Py_ssize_t i, size
        size_t lines
        coliter_t it
        const char *word = NULL

        int64_t NA = -1
        int64_t[:] codes
        int64_t current_category = 0

        char *errors = "strict"
        cdef StringPath path = _string_path(encoding)

        int ret = 0
        kh_str_t *table
        khiter_t k

    lines = line_end - line_start
    codes = np.empty(lines, dtype=np.int64)

    # factorize parsed values, creating a hash table
    # bytes -> category code
    with nogil:
        table = kh_init_str()
        coliter_setup(&it, parser, col, line_start)

        for i in range(lines):
            COLITER_NEXT(it, word)

            if na_filter:
                k = kh_get_str(na_hashset, word)
                # is in NA values
                if k != na_hashset.n_buckets:
                    na_count += 1
                    codes[i] = NA
                    continue

            k = kh_get_str(table, word)
            # not in the hash table
            if k == table.n_buckets:
                k = kh_put_str(table, word, &ret)
                table.vals[k] = current_category
                current_category += 1

            codes[i] = table.vals[k]

    # parse and box categories to python strings
    result = np.empty(table.n_occupied, dtype=np.object_)
    if path == ENCODED:
        for k in range(table.n_buckets):
            if kh_exist_str(table, k):
                size = strlen(table.keys[k])
                result[table.vals[k]] = PyUnicode_Decode(
                    table.keys[k], size, encoding, errors)
    elif path == UTF8:
        for k in range(table.n_buckets):
            if kh_exist_str(table, k):
                result[table.vals[k]] = PyUnicode_FromString(table.keys[k])
    elif path == CSTRING:
        for k in range(table.n_buckets):
            if kh_exist_str(table, k):
                result[table.vals[k]] = PyBytes_FromString(table.keys[k])

    kh_destroy_str(table)
    return np.asarray(codes), result, na_count

cdef _to_fw_string(parser_t *parser, int col, int line_start,
                   int line_end, size_t width):
    cdef:
        Py_ssize_t i
        coliter_t it
        const char *word = NULL
        char *data
        ndarray result

    result = np.empty(line_end - line_start, dtype='|S%d' % width)
    data = <char*> result.data

    with nogil:
        _to_fw_string_nogil(parser, col, line_start, line_end, width, data)

    return result

cdef inline void _to_fw_string_nogil(parser_t *parser, int col,
                                     int line_start, int line_end,
                                     size_t width, char *data) nogil:
    cdef:
        Py_ssize_t i
        coliter_t it
        const char *word = NULL

    coliter_setup(&it, parser, col, line_start)

    for i in range(line_end - line_start):
        COLITER_NEXT(it, word)
        strncpy(data, word, width)
        data += width

cdef char* cinf = b'inf'
cdef char* cposinf = b'+inf'
cdef char* cneginf = b'-inf'

cdef _try_double(parser_t *parser, int col, int line_start, int line_end,
                 bint na_filter, kh_str_t *na_hashset, object na_flist):
    cdef:
        int error, na_count = 0
        size_t i, lines
        coliter_t it
        const char *word = NULL
        char *p_end
        double *data
        double NA = na_values[np.float64]
        kh_float64_t *na_fset
        ndarray result
        khiter_t k
        bint use_na_flist = len(na_flist) > 0

    lines = line_end - line_start
    result = np.empty(lines, dtype=np.float64)
    data = <double *> result.data
    na_fset = kset_float64_from_list(na_flist)
    with nogil:
        error = _try_double_nogil(parser, col, line_start, line_end,
                                  na_filter, na_hashset, use_na_flist,
                                  na_fset, NA, data, &na_count)
    kh_destroy_float64(na_fset)
    if error != 0:
        return None, None
    return result, na_count

cdef inline int _try_double_nogil(parser_t *parser, int col,
                                  int line_start, int line_end,
                                  bint na_filter, kh_str_t *na_hashset,
                                  bint use_na_flist,
                                  const kh_float64_t *na_flist,
                                  double NA, double *data,
                                  int *na_count) nogil:
    cdef:
        int error,
        size_t i
        size_t lines = line_end - line_start
        coliter_t it
        const char *word = NULL
        char *p_end
        khiter_t k, k64

    global errno

    na_count[0] = 0
    coliter_setup(&it, parser, col, line_start)

    if na_filter:
        for i in range(lines):
            COLITER_NEXT(it, word)

            k = kh_get_str(na_hashset, word)
            # in the hash table
            if k != na_hashset.n_buckets:
                na_count[0] += 1
                data[0] = NA
            else:
                data[0] = parser.converter(word, &p_end, parser.decimal,
                                           parser.sci, parser.thousands, 1)
                if errno != 0 or p_end[0] or p_end == word:
                    if (strcasecmp(word, cinf) == 0 or
                                strcasecmp(word, cposinf) == 0):
                        data[0] = INF
                    elif strcasecmp(word, cneginf) == 0:
                        data[0] = NEGINF
                    else:
                        # Just return a non-zero value since
                        # the errno is never consumed.
                        return 1
                if use_na_flist:
                    k64 = kh_get_float64(na_flist, data[0])
                    if k64 != na_flist.n_buckets:
                        na_count[0] += 1
                        data[0] = NA
            data += 1
    else:
        for i in range(lines):
            COLITER_NEXT(it, word)
            data[0] = parser.converter(word, &p_end, parser.decimal,
                                       parser.sci, parser.thousands, 1)
            if errno != 0 or p_end[0] or p_end == word:
                if (strcasecmp(word, cinf) == 0 or
                            strcasecmp(word, cposinf) == 0):
                    data[0] = INF
                elif strcasecmp(word, cneginf) == 0:
                    data[0] = NEGINF
                else:
                    # Just return a non-zero value since
                    # the errno is never consumed.
                    return 1
            data += 1

    return 0

cdef _try_int64(parser_t *parser, int col, int line_start, int line_end,
                bint na_filter, kh_str_t *na_hashset):
    cdef:
        int error, na_count = 0
        size_t i, lines
        coliter_t it
        int64_t *data
        ndarray result

        int64_t NA = na_values[np.int64]
        khiter_t k

    lines = line_end - line_start
    result = np.empty(lines, dtype=np.int64)
    data = <int64_t *> result.data
    coliter_setup(&it, parser, col, line_start)
    with nogil:
        error = _try_int64_nogil(parser, col, line_start, line_end,
                                 na_filter, na_hashset, NA, data, &na_count)
    if error != 0:
        if error == ERROR_OVERFLOW:
            # Can't get the word variable
            raise OverflowError('Overflow')
        return None, None

    return result, na_count

cdef inline int _try_int64_nogil(parser_t *parser, int col, int line_start,
                                 int line_end, bint na_filter,
                                 const kh_str_t *na_hashset, int64_t NA,
                                 int64_t *data, int *na_count) nogil:
    cdef:
        int error
        size_t i
        size_t lines = line_end - line_start
        coliter_t it
        const char *word = NULL
        khiter_t k

    na_count[0] = 0
    coliter_setup(&it, parser, col, line_start)

    if na_filter:
        for i in range(lines):
            COLITER_NEXT(it, word)
            k = kh_get_str(na_hashset, word)
            # in the hash table
            if k != na_hashset.n_buckets:
                na_count[0] += 1
                data[i] = NA
                continue

            data[i] = str_to_int64(word, INT64_MIN, INT64_MAX,
                                   &error, parser.thousands)
            if error != 0:
                return error
    else:
        for i in range(lines):
            COLITER_NEXT(it, word)
            data[i] = str_to_int64(word, INT64_MIN, INT64_MAX,
                                   &error, parser.thousands)
            if error != 0:
                return error

    return 0

cdef _try_bool(parser_t *parser, int col, int line_start, int line_end,
               bint na_filter, kh_str_t *na_hashset):
    cdef:
        int na_count
        size_t lines = line_end - line_start
        uint8_t *data
        cnp.ndarray[cnp.uint8_t, ndim=1] result

        uint8_t NA = na_values[np.bool_]

    result = np.empty(lines)
    data = <uint8_t *> result.data

    with nogil:
        error = _try_bool_nogil(parser, col, line_start,
                                line_end, na_filter,
                                na_hashset, NA, data,
                                &na_count)
    if error != 0:
        return None, None
    return result.view(np.bool_), na_count

cdef inline int _try_bool_nogil(parser_t *parser, int col, int line_start,
                                int line_end, bint na_filter,
                                const kh_str_t *na_hashset, uint8_t NA,
                                uint8_t *data, int *na_count) nogil:
    cdef:
        int error
        size_t lines = line_end - line_start
        size_t i
        coliter_t it
        const char *word = NULL
        khiter_t k
    na_count[0] = 0

    coliter_setup(&it, parser, col, line_start)

    if na_filter:
        for i in range(lines):
            COLITER_NEXT(it, word)

            k = kh_get_str(na_hashset, word)
            # in the hash table
            if k != na_hashset.n_buckets:
                na_count[0] += 1
                data[0] = NA
                data += 1
                continue

            error = to_boolean(word, data)
            if error != 0:
                return error
            data += 1
    else:
        for i in range(lines):
            COLITER_NEXT(it, word)

            error = to_boolean(word, data)
            if error != 0:
                return error
            data += 1
    return 0

cdef _try_bool_flex(parser_t *parser, int col, int line_start, int line_end,
                    bint na_filter, const kh_str_t *na_hashset,
                    const kh_str_t *true_hashset,
                    const kh_str_t *false_hashset):
    cdef:
        int error, na_count = 0
        size_t i, lines
        coliter_t it
        const char *word = NULL
        uint8_t *data
        ndarray result

        uint8_t NA = na_values[np.bool_]
        khiter_t k

    lines = line_end - line_start
    result = np.empty(lines, dtype=np.uint8)
    data = <uint8_t *> result.data
    with nogil:
        error = _try_bool_flex_nogil(parser, col, line_start, line_end,
                                     na_filter, na_hashset, true_hashset,
                                     false_hashset, NA, data, &na_count)
    if error != 0:
        return None, None
    return result.view(np.bool_), na_count

cdef inline int _try_bool_flex_nogil(parser_t *parser, int col, int line_start,
                                     int line_end, bint na_filter,
                                     const kh_str_t *na_hashset,
                                     const kh_str_t *true_hashset,
                                     const kh_str_t *false_hashset,
                                     uint8_t NA, uint8_t *data,
                                     int *na_count) nogil:
    cdef:
        int error = 0
        size_t i
        size_t lines = line_end - line_start
        coliter_t it
        const char *word = NULL
        khiter_t k

    na_count[0] = 0
    coliter_setup(&it, parser, col, line_start)

    if na_filter:
        for i in range(lines):
            COLITER_NEXT(it, word)

            k = kh_get_str(na_hashset, word)
            # in the hash table
            if k != na_hashset.n_buckets:
                na_count[0] += 1
                data[0] = NA
                data += 1
                continue

            k = kh_get_str(true_hashset, word)
            if k != true_hashset.n_buckets:
                data[0] = 1
                data += 1
                continue
            k = kh_get_str(false_hashset, word)
            if k != false_hashset.n_buckets:
                data[0] = 0
                data += 1
                continue

            error = to_boolean(word, data)
            if error != 0:
                return error
            data += 1
    else:
        for i in range(lines):
            COLITER_NEXT(it, word)

            k = kh_get_str(true_hashset, word)
            if k != true_hashset.n_buckets:
                data[0] = 1
                data += 1
                continue

            k = kh_get_str(false_hashset, word)
            if k != false_hashset.n_buckets:
                data[0] = 0
                data += 1
                continue

            error = to_boolean(word, data)
            if error != 0:
                return error
            data += 1

    return 0


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
        if not PyBytes_Check(val):
            raise ValueError('Must be all encoded bytes')

        k = kh_put_str(table, PyBytes_AsString(val), &ret)

    return table

cdef kh_float64_t* kset_float64_from_list(values) except NULL:
    # caller takes responsibility for freeing the hash table
    cdef:
        Py_ssize_t i
        khiter_t k
        kh_float64_t *table
        int ret = 0
        cnp.float64_t val
        object value

    table = kh_init_float64()

    for value in values:
        val = float(value)

        k = kh_put_float64(table, val, &ret)

    return table


cdef raise_parser_error(object base, parser_t *parser):
    cdef:
        object old_exc
        PyObject *type
        PyObject *value
        PyObject *traceback

    if PyErr_Occurred():
        PyErr_Fetch(&type, &value, &traceback);
        Py_XDECREF(type)
        Py_XDECREF(traceback)
        if value != NULL:
            old_exc = <object> value
            Py_XDECREF(value)
            raise old_exc
    message = '%s. C error: ' % base
    if parser.error_msg != NULL:
        if PY3:
            message += parser.error_msg.decode('utf-8')
        else:
            message += parser.error_msg
    else:
        message += 'no error message set'

    raise CParserError(message)


def _concatenate_chunks(list chunks):
    cdef:
        list names = list(chunks[0].keys())
        object name
        list warning_columns
        object warning_names
        object common_type

    result = {}
    warning_columns = list()
    for name in names:
        arrs = [chunk.pop(name) for chunk in chunks]
        # Check each arr for consistent types.
        dtypes = set([a.dtype for a in arrs])
        if len(dtypes) > 1:
            common_type = np.find_common_type(dtypes, [])
            if common_type == np.object:
                warning_columns.append(str(name))

        if is_categorical_dtype(dtypes.pop()):
            result[name] = union_categoricals(arrs, sort_categories=True)
        else:
            result[name] = np.concatenate(arrs)

    if warning_columns:
        warning_names = ','.join(warning_columns)
        warning_message = " ".join([
            "Columns (%s) have mixed types." % warning_names,
            "Specify dtype option on import or set low_memory=False."
          ])
        warnings.warn(warning_message, DtypeWarning, stacklevel=8)
    return result


# ----------------------------------------------------------------------
# NA values
def _compute_na_values():
    int64info = np.iinfo(np.int64)
    int32info = np.iinfo(np.int32)
    int16info = np.iinfo(np.int16)
    int8info = np.iinfo(np.int8)
    uint64info = np.iinfo(np.uint64)
    uint32info = np.iinfo(np.uint32)
    uint16info = np.iinfo(np.uint16)
    uint8info = np.iinfo(np.uint8)
    na_values = {
        np.float64: np.nan,
        np.int64: int64info.min,
        np.int32: int32info.min,
        np.int16: int16info.min,
        np.int8: int8info.min,
        np.uint64: uint64info.max,
        np.uint32: uint32info.max,
        np.uint16: uint16info.max,
        np.uint8: uint8info.max,
        np.bool_: uint8info.max,
        np.object_: np.nan   # oof
    }
    return na_values

na_values = _compute_na_values()

for k in list(na_values):
    na_values[np.dtype(k)] = na_values[k]


cdef _apply_converter(object f, parser_t *parser, int col,
                      int line_start, int line_end,
                      char* c_encoding):
    cdef:
        int error
        Py_ssize_t i
        size_t lines
        coliter_t it
        const char *word = NULL
        char *errors = "strict"
        ndarray[object] result
        object val

    lines = line_end - line_start
    result = np.empty(lines, dtype=np.object_)

    coliter_setup(&it, parser, col, line_start)

    if not PY3 and c_encoding == NULL:
        for i in range(lines):
            COLITER_NEXT(it, word)
            val = PyBytes_FromString(word)
            result[i] = f(val)
    elif ((PY3 and c_encoding == NULL) or c_encoding == b'utf-8'):
        for i in range(lines):
            COLITER_NEXT(it, word)
            val = PyUnicode_FromString(word)
            result[i] = f(val)
    else:
        for i in range(lines):
            COLITER_NEXT(it, word)
            val = PyUnicode_Decode(word, strlen(word),
                                   c_encoding, errors)
            result[i] = f(val)

    return lib.maybe_convert_objects(result)


def _to_structured_array(dict columns, object names, object usecols):
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
    else:
        # single line header
        names = names[0]

    if usecols is not None:
        names = [n for i, n in enumerate(names)
                 if i in usecols or n in usecols]

    dt = np.dtype([(str(name), columns[i].dtype)
                   for i, name in enumerate(names)])
    fnames = dt.names
    fields = dt.fields

    nfields = len(fields)

    if PY3:
        length = len(list(columns.values())[0])
    else:
        length = len(columns.values()[0])

    stride = dt.itemsize

    # We own the data.
    buf = <char*> malloc(length * stride)

    recs = util.sarr_from_data(dt, length, buf)
    assert(recs.flags.owndata)

    for i in range(nfields):
        # XXX
        field_type = fields[fnames[i]]

        # (dtype, stride) tuple
        offset = field_type[1]
        elsize = field_type[0].itemsize
        column = columns[i]

        _fill_structured_column(buf + offset, <char*> column.data,
                                elsize, stride, length,
                                field_type[0] == np.object_)

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


def _maybe_encode(values):
    if values is None:
        return []
    return [x.encode('utf-8') if isinstance(x, unicode) else x for x in values]
