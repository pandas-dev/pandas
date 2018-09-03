# Copyright (c) 2012, Lambda Foundry, Inc.
# See LICENSE for the license
import os
import sys
import time
import warnings

from csv import QUOTE_MINIMAL, QUOTE_NONNUMERIC, QUOTE_NONE

from libc.stdlib cimport free
from libc.string cimport strncpy, strlen, strcasecmp

cimport cython
from cython cimport Py_ssize_t

from cpython cimport (PyObject, PyBytes_FromString,
                      PyBytes_AsString, PyBytes_Check,
                      PyUnicode_Check, PyUnicode_AsUTF8String,
                      PyErr_Occurred, PyErr_Fetch)
from cpython.ref cimport Py_XDECREF


cdef extern from "Python.h":
    object PyUnicode_FromString(char *v)

    object PyUnicode_Decode(char *v, Py_ssize_t size, char *encoding,
                            char *errors)


import numpy as np
cimport numpy as cnp
from numpy cimport ndarray, uint8_t, uint64_t, int64_t, float64_t
cnp.import_array()

from util cimport UINT64_MAX, INT64_MAX, INT64_MIN
import lib

from khash cimport (
    khiter_t,
    kh_str_t, kh_init_str, kh_put_str, kh_exist_str,
    kh_get_str, kh_destroy_str,
    kh_float64_t, kh_get_float64, kh_destroy_float64,
    kh_put_float64, kh_init_float64,
    kh_strbox_t, kh_put_strbox, kh_get_strbox, kh_init_strbox,
    kh_destroy_strbox)

import pandas.compat as compat
from pandas.core.dtypes.common import (
    is_categorical_dtype,
    is_integer_dtype, is_float_dtype,
    is_bool_dtype, is_object_dtype,
    is_datetime64_dtype,
    pandas_dtype)
from pandas.core.arrays import Categorical
from pandas.core.dtypes.concat import union_categoricals
import pandas.io.common as icom

from pandas.errors import (ParserError, DtypeWarning,
                           EmptyDataError, ParserWarning)

# Import CParserError as alias of ParserError for backwards compatibility.
# Ultimately, we want to remove this import. See gh-12665 and gh-14479.
CParserError = ParserError


cdef bint PY3 = (sys.version_info[0] >= 3)

cdef double INF = <double> np.inf
cdef double NEGINF = -INF


cdef extern from "errno.h":
    int errno

cdef extern from "headers/portable.h":
    # I *think* this is here so that strcasecmp is defined on Windows
    # so we don't get
    # `parsers.obj : error LNK2001: unresolved external symbol strcasecmp`
    # in Appveyor.
    # In a sane world, the `from libc.string cimport` above would fail
    # loudly.
    pass

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

        int64_t chunksize  # Number of bytes to prepare for each chunk
        char *data         # pointer to data to be processed
        int64_t datalen    # amount of data available
        int64_t datapos

        # where to write out tokenized data
        char *stream
        int64_t stream_len
        int64_t stream_cap

        # Store words in (potentially ragged) matrix for now, hmm
        char **words
        int64_t *word_starts  # where we are in the stream
        int64_t words_len
        int64_t words_cap

        char *pword_start        # pointer to stream start of current field
        int64_t word_start       # position start of current field

        int64_t *line_start      # position in words for start of line
        int64_t *line_fields     # Number of fields in each line
        int64_t lines            # Number of lines observed
        int64_t file_lines       # Number of lines observed (with bad/skipped)
        int64_t lines_cap        # Vector capacity

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
        # int numeric_field

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

        int header                  # Boolean: 1: has header, 0: no header
        int64_t header_start        # header row start
        int64_t header_end          # header row end

        void *skipset
        PyObject *skipfunc
        int64_t skip_first_N_rows
        int64_t skipfooter
        # pick one, depending on whether the converter requires GIL
        double (*double_converter_nogil)(const char *, char **,
                                         char, char, char, int) nogil
        double (*double_converter_withgil)(const char *, char **,
                                           char, char, char, int)

        #  error handling
        char *warn_msg
        char *error_msg

        int64_t skip_empty_lines

    ctypedef struct coliter_t:
        char **words
        int64_t *line_start
        int64_t col

    ctypedef struct uint_state:
        int seen_sint
        int seen_uint
        int seen_null

    void uint_state_init(uint_state *self)
    int uint64_conflict(uint_state *self)

    void coliter_setup(coliter_t *it, parser_t *parser,
                       int64_t i, int64_t start) nogil
    void COLITER_NEXT(coliter_t, const char *) nogil

    parser_t* parser_new()

    int parser_init(parser_t *self) nogil
    void parser_free(parser_t *self) nogil
    void parser_del(parser_t *self) nogil
    int parser_add_skiprow(parser_t *self, int64_t row)

    int parser_set_skipfirstnrows(parser_t *self, int64_t nrows)

    void parser_set_default_options(parser_t *self)

    int parser_consume_rows(parser_t *self, size_t nrows)

    int parser_trim_buffers(parser_t *self)

    int tokenize_all_rows(parser_t *self) nogil
    int tokenize_nrows(parser_t *self, size_t nrows) nogil

    int64_t str_to_int64(char *p_item, int64_t int_min,
                         int64_t int_max, int *error, char tsep) nogil
    uint64_t str_to_uint64(uint_state *state, char *p_item, int64_t int_max,
                           uint64_t uint_max, int *error, char tsep) nogil

    double xstrtod(const char *p, char **q, char decimal, char sci,
                   char tsep, int skip_trailing) nogil
    double precise_xstrtod(const char *p, char **q, char decimal, char sci,
                           char tsep, int skip_trailing) nogil
    double round_trip(const char *p, char **q, char decimal, char sci,
                      char tsep, int skip_trailing) nogil

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


cdef class TextReader:
    """

    # source: StringIO or file object

    """

    cdef:
        parser_t *parser
        object file_handle, na_fvalues
        object true_values, false_values
        object handle
        bint na_filter, keep_default_na, verbose, has_usecols, has_mi_columns
        int64_t parser_start
        list clocks
        char *c_encoding
        kh_str_t *false_set
        kh_str_t *true_set

    cdef public:
        int64_t leading_cols, table_width, skipfooter, buffer_lines
        object allow_leading_cols
        object delimiter, converters, delim_whitespace
        object na_values
        object memory_map
        object header, orig_header, names, header_start, header_end
        object index_col
        object low_memory
        object skip_rows
        object dtype
        object encoding
        object compression
        object mangle_dupe_cols
        object tupleize_cols
        object usecols
        list dtype_cast_order
        set noconvert

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
                  keep_default_na=True,

                  true_values=None,
                  false_values=None,
                  allow_leading_cols=True,
                  low_memory=False,
                  skip_rows=None,
                  skipfooter=0,
                  verbose=False,
                  mangle_dupe_cols=True,
                  tupleize_cols=False,
                  float_precision=None,
                  skip_blank_lines=True):

        # set encoding for native Python and C library
        if encoding is not None:
            if not isinstance(encoding, bytes):
                encoding = encoding.encode('utf-8')
            encoding = encoding.lower()
            self.c_encoding = <char*> encoding
        else:
            self.c_encoding = NULL

        self.encoding = encoding

        self.parser = parser_new()
        self.parser.chunksize = tokenize_chunksize

        self.mangle_dupe_cols = mangle_dupe_cols
        self.tupleize_cols = tupleize_cols

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

        # ----------------------------------------
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
                raise ValueError('Only length-1 escapes supported')
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

        self.skip_rows = skip_rows
        if skip_rows is not None:
            self._make_skiprow_set()

        self.skipfooter = skipfooter

        # suboptimal
        if usecols is not None:
            self.has_usecols = 1
            # GH-20558, validate usecols at higher level and only pass clean
            # usecols into TextReader.
            self.usecols = usecols

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

        self.keep_default_na = keep_default_na
        self.converters = converters
        self.na_filter = na_filter

        self.verbose = verbose
        self.low_memory = low_memory
        self.parser.double_converter_nogil = xstrtod
        self.parser.double_converter_withgil = NULL
        if float_precision == 'high':
            self.parser.double_converter_nogil = precise_xstrtod
            self.parser.double_converter_withgil = NULL
        elif float_precision == 'round_trip':  # avoid gh-15140
            self.parser.double_converter_nogil = NULL
            self.parser.double_converter_withgil = round_trip

        if isinstance(dtype, dict):
            dtype = {k: pandas_dtype(dtype[k])
                     for k in dtype}
        elif dtype is not None:
            dtype = pandas_dtype(dtype)

        self.dtype = dtype

        # XXX
        self.noconvert = set()

        self.index_col = index_col

        # ----------------------------------------
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
            if isinstance(header, list):
                if len(header) > 1:
                    # need to artificially skip the final line
                    # which is still a header line
                    header = list(header)
                    header.append(header[-1] + 1)
                    self.parser.header_end = header[-1]
                    self.has_mi_columns = 1
                else:
                    self.parser.header_end = header[0]

                self.parser_start = header[-1] + 1
                self.parser.header_start = header[0]
                self.parser.header = header[0]
                self.header = header
            else:
                self.parser.header_start = header
                self.parser.header_end = header
                self.parser_start = header + 1
                self.parser.header = header
                self.header = [ header ]

        self.names = names
        self.header, self.table_width = self._get_header()

        if not self.table_width:
            raise EmptyDataError("No columns to parse from file")

        # Compute buffer_lines as function of table width.
        heuristic = 2**20 // self.table_width
        self.buffer_lines = 1
        while self.buffer_lines * 2 < heuristic:
            self.buffer_lines *= 2

    def __init__(self, *args, **kwargs):
        pass

    def __dealloc__(self):
        parser_free(self.parser)
        if self.true_set:
            kh_destroy_str(self.true_set)
            self.true_set = NULL
        if self.false_set:
            kh_destroy_str(self.false_set)
            self.false_set = NULL
        parser_del(self.parser)

    def close(self):
        # we need to properly close an open derived
        # filehandle here, e.g. and UTFRecoder
        if self.handle is not None:
            try:
                self.handle.close()
            except:
                pass
        # also preemptively free all allocated memory
        parser_free(self.parser)
        if self.true_set:
            kh_destroy_str(self.true_set)
            self.true_set = NULL
        if self.false_set:
            kh_destroy_str(self.false_set)
            self.false_set = NULL

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
        if isinstance(self.skip_rows, (int, np.integer)):
            parser_set_skipfirstnrows(self.parser, self.skip_rows)
        elif not callable(self.skip_rows):
            for i in self.skip_rows:
                parser_add_skiprow(self.parser, i)
        else:
            self.parser.skipfunc = <PyObject *> self.skip_rows

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
                    content = source.read()
                    source.close()
                    source = compat.StringIO(bz2.decompress(content))
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

            if b'utf-16' in (self.encoding or b''):
                # we need to read utf-16 through UTF8Recoder.
                # if source is utf-16, convert source to utf-8 by UTF8Recoder.
                source = icom.UTF8Recoder(source,
                                          self.encoding.decode('utf-8'))
                self.encoding = b'utf-8'
                self.c_encoding = <char*> self.encoding

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
                        'File {source} does not exist'.format(source=source))
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
            Py_ssize_t i, start, field_count, passed_count, unnamed_count  # noqa
            char *word
            object name
            int status
            int64_t hr, data_line
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
                            ','.join(str(m) for m in msg), len(msg))
                    raise ParserError(
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
                            name = ('Unnamed: {i}_level_{lvl}'
                                    .format(i=i, lvl=level))
                        else:
                            name = 'Unnamed: {i}'.format(i=i)
                        unnamed_count += 1

                    count = counts.get(name, 0)

                    if not self.has_mi_columns and self.mangle_dupe_cols:
                        while count > 0:
                            counts[name] = count + 1
                            name = '%s.%d' % (name, count)
                            count = counts.get(name, 0)

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
        else:  # not self.has_usecols:

            field_count = self.parser.line_fields[data_line]

            # #2981
            if self.names is not None:
                field_count = max(field_count, len(self.names))

            passed_count = len(header[0])

            # if passed_count > field_count:
            #     raise ParserError('Column names have %d fields, '
            #                        'data has %d fields'
            #                        % (passed_count, field_count))

            if (self.has_usecols and self.allow_leading_cols and
                    not callable(self.usecols)):
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
            int64_t buffered_lines
            int64_t irows, footer = 0

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
                                            upcast_na=True)
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

    cdef _start_clock(self):
        self.clocks.append(time.time())

    cdef _end_clock(self, what):
        if self.verbose:
            elapsed = time.time() - self.clocks.pop(-1)
            print('%s took: %.2f ms' % (what, elapsed * 1000))

    def set_noconvert(self, i):
        self.noconvert.add(i)

    def remove_noconvert(self, i):
        self.noconvert.remove(i)

    def _convert_column_data(self, rows=None, upcast_na=False, footer=0):
        cdef:
            int64_t i
            int nused
            kh_str_t *na_hashset = NULL
            int64_t start, end
            object name, na_flist, col_dtype = None
            bint na_filter = 0
            int64_t num_cols

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
            raise ParserError(
                "Too many columns specified: expected {expected} and "
                "found {found}"
                .format(expected=self.table_width - self.leading_cols,
                        found=num_cols))

        results = {}
        nused = 0
        for i in range(self.table_width):
            if i < self.leading_cols:
                # Pass through leading columns always
                name = i
            elif (self.usecols and not callable(self.usecols) and
                    nused == len(self.usecols)):
                # Once we've gathered all requested columns, stop. GH5766
                break
            else:
                name = self._get_column_name(i, nused)
                usecols = set()
                if callable(self.usecols):
                    if self.usecols(name):
                        usecols = {i}
                else:
                    usecols = self.usecols
                if self.has_usecols and not (i in usecols or
                                             name in usecols):
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

            col_dtype = None
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

            if conv:
                if col_dtype is not None:
                    warnings.warn(("Both a converter and dtype were specified "
                                   "for column {0} - only the converter will "
                                   "be used").format(name), ParserWarning,
                                  stacklevel=5)
                results[i] = _apply_converter(conv, self.parser, i, start, end,
                                              self.c_encoding)
                continue

            # Should return as the desired dtype (inferred or specified)
            col_res, na_count = self._convert_tokens(
                i, start, end, name, na_filter, na_hashset,
                na_flist, col_dtype)

            if na_filter:
                self._free_na_set(na_hashset)

            if upcast_na and na_count > 0:
                col_res = _maybe_upcast(col_res)

            if col_res is None:
                raise ParserError('Unable to parse column {i}'.format(i=i))

            results[i] = col_res

        self.parser_start += end - start

        return results

    cdef inline _convert_tokens(self, Py_ssize_t i, int start, int end,
                                object name, bint na_filter,
                                kh_str_t *na_hashset,
                                object na_flist, object col_dtype):

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
                except ValueError:
                    # This error is raised from trying to convert to uint64,
                    # and we discover that we cannot convert to any numerical
                    # dtype successfully. As a result, we leave the data
                    # column AS IS with object dtype.
                    col_res, na_count = self._convert_with_dtype(
                        np.dtype('object'), i, start, end, 0,
                        0, na_hashset, na_flist)
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
                             int64_t start, int64_t end,
                             bint na_filter,
                             bint user_dtype,
                             kh_str_t *na_hashset,
                             object na_flist):
        if is_integer_dtype(dtype):
            try:
                result, na_count = _try_int64(self.parser, i, start,
                                              end, na_filter, na_hashset)
                if user_dtype and na_count is not None:
                    if na_count > 0:
                        raise ValueError("Integer column has NA values in "
                                         "column {column}".format(column=i))
            except OverflowError:
                result = _try_uint64(self.parser, i, start, end,
                                     na_filter, na_hashset)
                na_count = 0

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
                raise TypeError("the dtype {dtype} is not "
                                "supported for parsing".format(dtype=dtype))

            # unicode variable width
            return self._string_convert(i, start, end, na_filter,
                                        na_hashset)
        elif is_categorical_dtype(dtype):
            # TODO: I suspect that _categorical_convert could be
            # optimized when dtype is an instance of CategoricalDtype
            codes, cats, na_count = _categorical_convert(
                self.parser, i, start, end, na_filter,
                na_hashset, self.c_encoding)
            cat = Categorical._from_inferred_categories(cats, codes, dtype)
            return cat, na_count

        elif is_object_dtype(dtype):
            return self._string_convert(i, start, end, na_filter,
                                        na_hashset)
        elif is_datetime64_dtype(dtype):
            raise TypeError("the dtype {dtype} is not supported "
                            "for parsing, pass this column "
                            "using parse_dates instead".format(dtype=dtype))
        else:
            raise TypeError("the dtype {dtype} is not "
                            "supported for parsing".format(dtype=dtype))

    cdef _string_convert(self, Py_ssize_t i, int64_t start, int64_t end,
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
                if self.keep_default_na:
                    return _NA_VALUES, set()

                return list(), set()

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
        cdef int64_t j
        if self.has_usecols and self.names is not None:
            if (not callable(self.usecols) and
                    len(self.names) == len(self.usecols)):
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


# common NA values
# no longer excluding inf representations
# '1.#INF','-1.#INF', '1.#INF000000',
_NA_VALUES = _ensure_encoded(list(icom._NA_VALUES))


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


cdef _string_box_factorize(parser_t *parser, int64_t col,
                           int64_t line_start, int64_t line_end,
                           bint na_filter, kh_str_t *na_hashset):
    cdef:
        int error, na_count = 0
        Py_ssize_t i, lines
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


cdef _string_box_utf8(parser_t *parser, int64_t col,
                      int64_t line_start, int64_t line_end,
                      bint na_filter, kh_str_t *na_hashset):
    cdef:
        int error, na_count = 0
        Py_ssize_t i, lines
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


cdef _string_box_decode(parser_t *parser, int64_t col,
                        int64_t line_start, int64_t line_end,
                        bint na_filter, kh_str_t *na_hashset,
                        char *encoding):
    cdef:
        int error, na_count = 0
        Py_ssize_t i, size, lines
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
cdef _categorical_convert(parser_t *parser, int64_t col,
                          int64_t line_start, int64_t line_end,
                          bint na_filter, kh_str_t *na_hashset,
                          char *encoding):
    "Convert column data into codes, categories"
    cdef:
        int error, na_count = 0
        Py_ssize_t i, size, lines
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


cdef _to_fw_string(parser_t *parser, int64_t col, int64_t line_start,
                   int64_t line_end, int64_t width):
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


cdef inline void _to_fw_string_nogil(parser_t *parser, int64_t col,
                                     int64_t line_start, int64_t line_end,
                                     size_t width, char *data) nogil:
    cdef:
        int64_t i
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


cdef _try_double(parser_t *parser, int64_t col,
                 int64_t line_start, int64_t line_end,
                 bint na_filter, kh_str_t *na_hashset, object na_flist):
    cdef:
        int error, na_count = 0
        Py_ssize_t i, lines
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
    if parser.double_converter_nogil != NULL:  # if it can run without the GIL
        with nogil:
            error = _try_double_nogil(parser, parser.double_converter_nogil,
                                      col, line_start, line_end,
                                      na_filter, na_hashset, use_na_flist,
                                      na_fset, NA, data, &na_count)
    else:
        assert parser.double_converter_withgil != NULL
        error = _try_double_nogil(parser,
                                  <double (*)(const char *, char **,
                                              char, char, char, int)
                                  nogil>parser.double_converter_withgil,
                                  col, line_start, line_end,
                                  na_filter, na_hashset, use_na_flist,
                                  na_fset, NA, data, &na_count)
    kh_destroy_float64(na_fset)
    if error != 0:
        return None, None
    return result, na_count


cdef inline int _try_double_nogil(parser_t *parser,
                                  double (*double_converter)(
                                      const char *, char **, char,
                                      char, char, int) nogil,
                                  int col, int line_start, int line_end,
                                  bint na_filter, kh_str_t *na_hashset,
                                  bint use_na_flist,
                                  const kh_float64_t *na_flist,
                                  double NA, double *data,
                                  int *na_count) nogil:
    cdef:
        int error,
        Py_ssize_t i, lines = line_end - line_start
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
                data[0] = double_converter(word, &p_end, parser.decimal,
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
            data[0] = double_converter(word, &p_end, parser.decimal,
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


cdef _try_uint64(parser_t *parser, int64_t col,
                 int64_t line_start, int64_t line_end,
                 bint na_filter, kh_str_t *na_hashset):
    cdef:
        int error
        Py_ssize_t i, lines
        coliter_t it
        uint64_t *data
        ndarray result
        khiter_t k
        uint_state state

    lines = line_end - line_start
    result = np.empty(lines, dtype=np.uint64)
    data = <uint64_t *> result.data

    uint_state_init(&state)
    coliter_setup(&it, parser, col, line_start)
    with nogil:
        error = _try_uint64_nogil(parser, col, line_start, line_end,
                                  na_filter, na_hashset, data, &state)
    if error != 0:
        if error == ERROR_OVERFLOW:
            # Can't get the word variable
            raise OverflowError('Overflow')
        return None

    if uint64_conflict(&state):
        raise ValueError('Cannot convert to numerical dtype')

    if state.seen_sint:
        raise OverflowError('Overflow')

    return result


cdef inline int _try_uint64_nogil(parser_t *parser, int64_t col,
                                  int64_t line_start,
                                  int64_t line_end, bint na_filter,
                                  const kh_str_t *na_hashset,
                                  uint64_t *data, uint_state *state) nogil:
    cdef:
        int error
        Py_ssize_t i, lines = line_end - line_start
        coliter_t it
        const char *word = NULL
        khiter_t k

    coliter_setup(&it, parser, col, line_start)

    if na_filter:
        for i in range(lines):
            COLITER_NEXT(it, word)
            k = kh_get_str(na_hashset, word)
            # in the hash table
            if k != na_hashset.n_buckets:
                state.seen_null = 1
                data[i] = 0
                continue

            data[i] = str_to_uint64(state, word, INT64_MAX, UINT64_MAX,
                                    &error, parser.thousands)
            if error != 0:
                return error
    else:
        for i in range(lines):
            COLITER_NEXT(it, word)
            data[i] = str_to_uint64(state, word, INT64_MAX, UINT64_MAX,
                                    &error, parser.thousands)
            if error != 0:
                return error

    return 0


cdef _try_int64(parser_t *parser, int64_t col,
                int64_t line_start, int64_t line_end,
                bint na_filter, kh_str_t *na_hashset):
    cdef:
        int error, na_count = 0
        Py_ssize_t i, lines
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


cdef inline int _try_int64_nogil(parser_t *parser, int64_t col,
                                 int64_t line_start,
                                 int64_t line_end, bint na_filter,
                                 const kh_str_t *na_hashset, int64_t NA,
                                 int64_t *data, int *na_count) nogil:
    cdef:
        int error
        Py_ssize_t i, lines = line_end - line_start
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


cdef _try_bool_flex(parser_t *parser, int64_t col,
                    int64_t line_start, int64_t line_end,
                    bint na_filter, const kh_str_t *na_hashset,
                    const kh_str_t *true_hashset,
                    const kh_str_t *false_hashset):
    cdef:
        int error, na_count = 0
        Py_ssize_t i, lines
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


cdef inline int _try_bool_flex_nogil(parser_t *parser, int64_t col,
                                     int64_t line_start,
                                     int64_t line_end, bint na_filter,
                                     const kh_str_t *na_hashset,
                                     const kh_str_t *true_hashset,
                                     const kh_str_t *false_hashset,
                                     uint8_t NA, uint8_t *data,
                                     int *na_count) nogil:
    cdef:
        int error = 0
        Py_ssize_t i, lines = line_end - line_start
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
        float64_t val
        object value

    table = kh_init_float64()

    for value in values:
        val = float(value)

        k = kh_put_float64(table, val, &ret)

    return table


cdef raise_parser_error(object base, parser_t *parser):
    cdef:
        object old_exc
        object exc_type
        PyObject *type
        PyObject *value
        PyObject *traceback

    if PyErr_Occurred():
        PyErr_Fetch(&type, &value, &traceback)
        Py_XDECREF(traceback)

        if value != NULL:
            old_exc = <object> value
            Py_XDECREF(value)

            # PyErr_Fetch only returned the error message in *value,
            # so the Exception class must be extracted from *type.
            if isinstance(old_exc, compat.string_types):
                if type != NULL:
                    exc_type = <object> type
                else:
                    exc_type = ParserError

                Py_XDECREF(type)
                raise exc_type(old_exc)
            else:
                Py_XDECREF(type)
                raise old_exc

    message = '{base}. C error: '.format(base=base)
    if parser.error_msg != NULL:
        if PY3:
            message += parser.error_msg.decode('utf-8')
        else:
            message += parser.error_msg
    else:
        message += 'no error message set'

    raise ParserError(message)


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
        dtypes = {a.dtype for a in arrs}
        numpy_dtypes = {x for x in dtypes if not is_categorical_dtype(x)}
        if len(numpy_dtypes) > 1:
            common_type = np.find_common_type(numpy_dtypes, [])
            if common_type == np.object:
                warning_columns.append(str(name))

        dtype = dtypes.pop()
        if is_categorical_dtype(dtype):
            sort_categories = isinstance(dtype, str)
            result[name] = union_categoricals(arrs,
                                              sort_categories=sort_categories)
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


cdef _apply_converter(object f, parser_t *parser, int64_t col,
                      int64_t line_start, int64_t line_end,
                      char* c_encoding):
    cdef:
        int error
        Py_ssize_t i, lines
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


def _maybe_encode(values):
    if values is None:
        return []
    return [x.encode('utf-8') if isinstance(x, unicode) else x for x in values]


def sanitize_objects(ndarray[object] values, set na_values,
                     convert_empty=True):
    """
    Convert specified values, including the given set na_values and empty
    strings if convert_empty is True, to np.nan.

    Parameters
    ----------
    values : ndarray[object]
    na_values : set
    convert_empty : bool (default True)
    """
    cdef:
        Py_ssize_t i, n
        object val, onan
        Py_ssize_t na_count = 0
        dict memo = {}

    n = len(values)
    onan = np.nan

    for i in range(n):
        val = values[i]
        if (convert_empty and val == '') or (val in na_values):
            values[i] = onan
            na_count += 1
        elif val in memo:
            values[i] = memo[val]
        else:
            memo[val] = val

    return na_count
