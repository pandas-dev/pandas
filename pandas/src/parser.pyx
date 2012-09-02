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

cdef extern from "stdint.h":
    enum: UINT8_MAX
    enum: INT8_MIN
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

    int parser_init(parser_t *self) nogil
    void parser_free(parser_t *self) nogil

    void parser_set_default_options(parser_t *self)

    int parser_file_source_init(parser_t *self, FILE* fp)
    int parser_mmap_init(parser_t *self, FILE* fp)
    int parser_array_source_init(parser_t *self, char *bytes, size_t length)

    void debug_print_parser(parser_t *self)

    int tokenize_all_rows(parser_t *self) nogil
    int tokenize_nrows(parser_t *self, size_t nrows) nogil

    int64_t str_to_int64(char *p_item, int64_t int_min,
                         int64_t int_max, int *error, char tsep)
    uint64_t str_to_uint64(char *p_item, uint64_t uint_max, int *error)



DEFAULT_CHUNKSIZE = 256 * 1024

cdef class TextReader:
    '''

    # source: StringIO or file object

    '''

    cdef:
        parser_t *parser
        object file_handle, should_close
        bint factorize, na_filter

    cdef public:
        object delimiter, na_values, converters, delim_whitespace
        object memory_map
        object as_recarray

    def __cinit__(self, source, delimiter=b',', header=0,
                  memory_map=False,
                  tokenize_chunksize=DEFAULT_CHUNKSIZE,
                  delim_whitespace=False,
                  na_values=None,
                  converters=None,
                  thousands=None,
                  factorize=True,
                  as_recarray=False,
                  skipinitialspace=False,
                  decimal=b'.',
                  na_filter=True):
        self.parser = parser_new()
        self.parser.chunksize = tokenize_chunksize

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

        # TODO: no header vs. header is not the first row
        self.parser.header = header
        self.parser.skipinitialspace = skipinitialspace

        if len(decimal) != 1:
            raise ValueError('Only length-1 decimal markers supported')
        self.parser.decimal = ord(decimal)

        if thousands is not None:
            if len(thousands) != 1:
                raise ValueError('Only length-1 decimal markers supported')
            self.parser.thousands = ord(thousands)

        self.should_close = False

        self.delimiter = delimiter
        self.delim_whitespace = delim_whitespace

        self.memory_map = memory_map
        self.na_values = na_values
        self.converters = converters

        self.na_filter = na_filter
        self.as_recarray = as_recarray

    def __init__(self, *args, **kwards):
        pass

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

        # start = time.clock()

        if rows is not None:
            raise NotImplementedError
        else:
            with nogil:
                status = tokenize_all_rows(self.parser)

        # end = time.clock()
        # print 'Tokenization took %.4f sec' % (end - start)

        if status < 0:
            raise_parser_error('Error tokenizing data', self.parser)

        # start = time.clock()
        columns, names = self._convert_column_data()
        # end = time.clock()
        # print 'Type conversion took %.4f sec' % (end - start)

        # debug_print_parser(self.parser)

        if self.as_recarray:
            # start = time.clock()
            result = _to_structured_array(columns, names)
            # end = time.clock()
            # print 'to_structured_array took %.4f sec' % (end - start)

            return result
        else:
            return columns

    def _convert_column_data(self):
        cdef:
            Py_ssize_t i, ncols
            cast_func func
            kh_str_t *table
            int start, end

        ncols = self.parser.line_fields[0]

        na_values = ['NA', 'nan', 'NaN']
        table = kset_from_list(na_values)

        start = 0
        end = self.parser.lines

        names = []

        results = {}
        for i in range(ncols):
            # XXX
            if self.na_filter:
                na_mask = _get_na_mask(self.parser, i, start, end, table)
            else:
                na_mask = None

            conv = self._get_converter(i)

            if conv:
                col_res = _apply_converter(conv, self.parser, i, start, end)
                results[i] = col_res
                continue

            col_res = None
            for func in cast_func_order:
                col_res, na_count = func(self.parser, i, start, end,
                                         na_mask, self.na_filter)
                if col_res is not None:
                    results[i] = col_res
                    break

            if col_res is None:
                raise Exception('Unable to parse column %d' % i)

            names.append(i)
            results[i] = col_res

        # XXX: needs to be done elsewhere
        kh_destroy_str(table)

        return results, names

    def _get_converter(self, col):
        if self.converters is None:
            return None

        return self.converters.get(col)

    def _get_col_name(self, col):
        pass

class CParserError(Exception):
    pass


# ----------------------------------------------------------------------
# Type conversions / inference support code

ctypedef object (*cast_func)(parser_t *parser, int col,
                             int line_start, int line_end,
                             object _na_mask, bint na_filter)

cdef _string_box_factorize(parser_t *parser, int col,
                           int line_start, int line_end,
                           object _na_mask, bint na_filter):
    cdef:
        int error, na_count = 0
        Py_ssize_t i
        size_t lines
        coliter_t it
        char *word
        ndarray[object] result
        ndarray[uint8_t, cast=True] na_mask

        int ret = 0
        kh_strbox_t *table
        khiter_t

        object pyval

        object NA = na_values[np.object_]

    if na_filter:
        na_mask = _na_mask

    table = kh_init_strbox()

    lines = line_end - line_start
    result = np.empty(lines, dtype=np.object_)

    coliter_setup(&it, parser, col)

    for i in range(lines):
        word = COLITER_NEXT(it)

        if na_filter and na_mask[i]:
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
                 object _na_mask, bint na_filter):
    cdef:
        int error, na_count = 0
        size_t i, lines
        coliter_t it
        char *word
        double *data
        double NA = na_values[np.float64]
        ndarray result
        ndarray[uint8_t, cast=True] na_mask

    if na_filter:
        na_mask = _na_mask

    lines = line_end - line_start
    result = np.empty(lines, dtype=np.float64)
    data = <double *> result.data
    coliter_setup(&it, parser, col)

    if na_filter:
        for i in range(lines):
            word = COLITER_NEXT(it)
            if na_mask[i]:
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
                object _na_mask, bint na_filter):
    cdef:
        int error, na_count = 0
        size_t i, lines
        coliter_t it
        char *word
        int64_t *data
        ndarray result

        ndarray[uint8_t, cast=True] na_mask
        int64_t NA = na_values[np.int64]

    if na_filter:
        na_mask = _na_mask

    lines = line_end - line_start

    result = np.empty(lines, dtype=np.int64)

    data = <int64_t *> result.data

    coliter_setup(&it, parser, col)

    if na_filter:
        for i in range(lines):
            word = COLITER_NEXT(it)

            if na_mask[i]:
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
               object _na_mask, bint na_filter):
    cdef:
        int error, na_count = 0
        size_t i, lines
        coliter_t it
        char *word
        uint8_t *data
        ndarray result

        ndarray[uint8_t, cast=True] na_mask
        uint8_t NA = na_values[np.bool_]

    if na_filter:
        na_mask = _na_mask

    lines = line_end - line_start
    result = np.empty(lines, dtype=np.uint8)
    data = <uint8_t *> result.data
    coliter_setup(&it, parser, col)

    if na_filter:
        for i in range(lines):
            word = COLITER_NEXT(it)

            if na_mask[i]:
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
                  kh_str_t *na_table):
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

    coliter_setup(&it, parser, col)
    for i in range(lines):
        word = COLITER_NEXT(it)

        # length 0
        if word[0] == '\x00':
            result[i] = 1
            continue

        k = kh_get_str(na_table, word)
        # in the hash table
        if k != na_table.n_buckets:
            result[i] = 1
        else:
            result[i] = 0

    return result

cdef kh_str_t* kset_from_list(list values):
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

#----------------------------------------------------------------------
# NA values

na_values = {
    np.float64 : np.nan,
    np.int64 : INT64_MIN,
    np.int32 : INT32_MIN,
    np.bool_ : UINT8_MAX,
    np.object_ : np.nan    # oof
}


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

    coliter_setup(&it, parser, col)
    for i in range(lines):
        word = COLITER_NEXT(it)
        val = PyString_FromString(word)
        result[i] = f(val)

    return lib.maybe_convert_objects(result)

def _to_structured_array(dict columns, object colnames):
    cdef:
        ndarray recs, column
        cnp.dtype dt
        dict fields

        object name, fnames, field_type
        Py_ssize_t i, offset, nfields, length
        int stride, elsize
        char *buf

    dt = np.dtype([(str(name), columns[name].dtype) for name in colnames])
    fnames = dt.names
    fields = dt.fields

    nfields = len(fields)

    length = len(columns.values()[0])
    stride = dt.itemsize

    # start = time.clock()

    # we own the data
    buf = <char*> malloc(length * stride)

    recs = util.sarr_from_data(dt, length, buf)
    assert(recs.flags.owndata)

    # buf = <char*> recs.data
    # end = time.clock()
    # print 'took %.4f' % (end - start)

    for i in range(nfields):
        # start = time.clock()
        name = colnames[i]

        # XXX
        field_type = fields[fnames[i]]

        # (dtype, stride) tuple
        offset = field_type[1]
        elsize = field_type[0].itemsize
        column = columns[name]

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
