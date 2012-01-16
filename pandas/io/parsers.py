"""
Module contains tools for processing files into DataFrames or other objects
"""
from StringIO import StringIO
import re

import numpy as np

from pandas.core.index import Index, MultiIndex
from pandas.core.frame import DataFrame
import pandas._tseries as lib

from pandas.util.decorators import Appender

_parser_params = """Also supports optionally iterating or breaking of the file
into chunks.

Parameters
----------
filepath_or_buffer : string or file handle / StringIO
%s
header : int, default 0
    Row to use for the column labels of the parsed DataFrame
skiprows : list-like
    Row numbers to skip (0-indexed)
index_col : int or sequence, default None
    Column to use as the row labels of the DataFrame. If a sequence is
    given, a MultiIndex is used.
names : array-like
    List of column names
na_values : list-like, default None
    List of additional strings to recognize as NA/NaN
parse_dates : boolean, default False
    Attempt to parse dates in the index column(s)
date_parser : function
    Function to use for converting dates to strings. Defaults to
    dateutil.parser
nrows : int, default None
    Number of rows of file to read. Useful for reading pieces of large files
iterator : boolean, default False
    Return TextParser object
chunksize : int, default None
    Return TextParser object for iteration
skip_footer : int, default 0
    Number of line at bottom of file to skip
converters : dict. optional
    Dict of functions for converting values in certain columns. Keys can either
    be integers or column labels
verbose : boolean, default False
    Indicate number of NA values placed in non-numeric columns

Returns
-------
result : DataFrame or TextParser
"""

_csv_sep = """sep : string, default None
    Delimiter to use. By default will try to automatically determine
    this"""

_table_sep = """sep : string, default \\t (tab-stop)
    Delimiter to use"""

_read_csv_doc = """
Read CSV (comma-separated) file into DataFrame

%s

Returns
-------
parsed : DataFrame
""" % (_parser_params % _csv_sep)

_read_csv_doc = """
Read CSV (comma-separated) file into DataFrame

%s

Returns
-------
parsed : DataFrame
""" % (_parser_params % _csv_sep)

_read_table_doc = """
Read delimited file into DataFrame

%s

Returns
-------
parsed : DataFrame
""" % (_parser_params % _table_sep)

@Appender(_read_csv_doc)
def read_csv(filepath_or_buffer, sep=None, header=0, index_col=None, names=None,
             skiprows=None, na_values=None, parse_dates=False,
             date_parser=None, nrows=None, iterator=False, chunksize=None,
             skip_footer=0, converters=None, verbose=False):
    if hasattr(filepath_or_buffer, 'read'):
        f = filepath_or_buffer
    else:
        try:
            # universal newline mode
            f = open(filepath_or_buffer, 'U')
        except Exception: # pragma: no cover
            f = open(filepath_or_buffer, 'r')

    if date_parser is not None:
        parse_dates = True

    parser = TextParser(f, header=header, index_col=index_col,
                        names=names, na_values=na_values,
                        parse_dates=parse_dates,
                        date_parser=date_parser,
                        skiprows=skiprows,
                        delimiter=sep,
                        chunksize=chunksize,
                        skip_footer=skip_footer,
                        converters=converters,
                        verbose=verbose)

    if nrows is not None:
        return parser.get_chunk(nrows)
    elif chunksize or iterator:
        return parser

    return parser.get_chunk()

@Appender(_read_table_doc)
def read_table(filepath_or_buffer, sep='\t', header=0, index_col=None,
               names=None, skiprows=None, na_values=None, parse_dates=False,
               date_parser=None, nrows=None, iterator=False, chunksize=None,
               skip_footer=0, converters=None, verbose=False):
    return read_csv(filepath_or_buffer, sep=sep, header=header,
                    skiprows=skiprows, index_col=index_col,
                    na_values=na_values, date_parser=date_parser,
                    names=names, parse_dates=parse_dates,
                    nrows=nrows, iterator=iterator, chunksize=chunksize,
                    skip_footer=skip_footer, converters=converters,
                    verbose=verbose)

def read_clipboard(**kwargs):  # pragma: no cover
    """
    Read text from clipboard and pass to read_table. See read_table for the full
    argument list

    Returns
    -------
    parsed : DataFrame
    """
    from pandas.util.clipboard import clipboard_get
    text = clipboard_get()
    return read_table(StringIO(text), **kwargs)


class BufferedReader(object):
    """
    For handling different kinds of files, e.g. zip files where reading out a
    chunk of lines is faster than reading out one line at a time.
    """

    def __init__(self, fh, delimiter=','):
        pass # pragma: no coverage

class BufferedCSVReader(BufferedReader):
    pass

class TextParser(object):
    """
    Converts lists of lists/tuples into DataFrames with proper type inference
    and optional (e.g. string to datetime) conversion. Also enables iterating
    lazily over chunks of large files

    Parameters
    ----------
    data : file-like object or list
    names : sequence, default
    header : int, default 0
        Row to use to parse column labels. Defaults to the first row. Prior
        rows will be discarded
    index_col : int or list, default None
        Column or columns to use as the (possibly hierarchical) index
    na_values : iterable, defualt None
        Custom NA values
    parse_dates : boolean, default False
    date_parser : function, default None
    skiprows : list of integers
        Row numbers to skip
    skip_footer : int
        Number of line at bottom of file to skip
    """

    # common NA values
    # no longer excluding inf representations
    # '1.#INF','-1.#INF', '1.#INF000000',
    NA_VALUES = set(['-1.#IND', '1.#QNAN', '1.#IND', '-1.#QNAN',
                     '#N/A N/A', 'NA', '#NA', 'NULL', 'NaN',
                     'nan', ''])

    def __init__(self, f, delimiter=None, names=None, header=0,
                 index_col=None, na_values=None, parse_dates=False,
                 date_parser=None, chunksize=None, skiprows=None,
                 skip_footer=0, converters=None, verbose=False):
        """
        Workhorse function for processing nested list into DataFrame

        Should be replaced by np.genfromtxt eventually?
        """
        self.data = None
        self.buf = []
        self.pos = 0
        self.names = list(names) if names is not None else names
        self.header = header
        self.index_col = index_col
        self.parse_dates = parse_dates
        self.date_parser = date_parser
        self.chunksize = chunksize
        self.passed_names = names is not None
        self.skiprows = set() if skiprows is None else set(skiprows)
        self.skip_footer = skip_footer
        self.delimiter = delimiter
        self.verbose = verbose

        if converters is not None:
            assert(isinstance(converters, dict))
            self.converters = converters
        else:
            self.converters = {}

        assert(self.skip_footer >= 0)

        if na_values is None:
            self.na_values = self.NA_VALUES
        else:
            self.na_values = set(list(na_values)) | self.NA_VALUES

        if hasattr(f, 'readline'):
            self._make_reader(f)
        else:
            self.data = f
        self.columns = self._infer_columns()
        self.index_name = self._get_index_name()
        self._first_chunk = True

    def _make_reader(self, f):
        import csv

        sep = self.delimiter

        if sep is None or len(sep) == 1:
            sniff_sep = True
            # default dialect
            dia = csv.excel
            if sep is not None:
                sniff_sep = False
                dia.delimiter = sep
            # attempt to sniff the delimiter
            if sniff_sep:
                line = f.readline()
                while self.pos in self.skiprows:
                    self.pos += 1
                    line = f.readline()

                self.pos += 1
                sniffed = csv.Sniffer().sniff(line)
                dia.delimiter = sniffed.delimiter
                self.buf.extend(list(csv.reader(StringIO(line), dialect=dia)))
            reader = csv.reader(f, dialect=dia)
        else:
            reader = (re.split(sep, line.strip()) for line in f)

        self.data = reader

    def _infer_columns(self):
        names = self.names
        passed_names = self.names is not None
        if passed_names:
            self.header = None

        if self.header is not None:
            if len(self.buf) > 0:
                line = self.buf[0]
            else:
                line = self._next_line()

            while self.pos <= self.header:
                line = self._next_line()

            columns = []
            for i, c in enumerate(line):
                if c == '':
                    columns.append('Unnamed: %d' % i)
                else:
                    columns.append(c)

            counts = {}
            for i, col in enumerate(columns):
                cur_count = counts.get(col, 0)
                if cur_count > 0:
                    columns[i] = '%s.%d' % (col, cur_count)
                counts[col] = cur_count + 1
            self._clear_buffer()
        else:
            line = self._next_line()

            ncols = len(line)
            if not names:
                columns = ['X.%d' % (i + 1) for i in range(ncols)]
            else:
                columns = names

        return columns

    def _next_line(self):
        if isinstance(self.data, list):
            while self.pos in self.skiprows:
                self.pos += 1

            line = self.data[self.pos]
        else:
            while self.pos in self.skiprows:
                self.data.next()
                self.pos += 1
            line = self.data.next()
        self.pos += 1
        self.buf.append(line)

        return line

    def _clear_buffer(self):
        self.buf = []

    def __iter__(self):
        try:
            while True:
                yield self.get_chunk(self.chunksize)
        except StopIteration:
            pass

    def _get_index_name(self):
        columns = self.columns

        try:
            line = self._next_line()
        except StopIteration:
            line = None

        try:
            next_line = self._next_line()
        except StopIteration:
            next_line = None

        index_name = None

        # implicitly index_col=0 b/c 1 fewer column names
        implicit_first_cols = 0
        if line is not None:
            implicit_first_cols = len(line) - len(columns)
            if next_line is not None:
                if len(next_line) == len(line) + len(columns):
                    implicit_first_cols = 0
                    self.index_col = range(len(line))
                    self.buf = self.buf[1:]
                    return line

        if implicit_first_cols > 0:
            if self.index_col is None:
                if implicit_first_cols == 1:
                    self.index_col = 0
                else:
                    self.index_col = range(implicit_first_cols)
            index_name = None
        elif np.isscalar(self.index_col):
            index_name = columns.pop(self.index_col)
            if 'Unnamed' in index_name:
                index_name = None
        elif self.index_col is not None:
            cp_cols = list(columns)
            index_name = []
            for i in self.index_col:
                name = cp_cols[i]
                columns.remove(name)
                index_name.append(name)

        return index_name

    def get_chunk(self, rows=None):
        if rows is not None and self.skip_footer:
            raise ValueError('skip_footer not supported for iteration')

        try:
            content = self._get_lines(rows)
        except StopIteration:
            if self._first_chunk:
                content = []
            else:
                raise

        # done with first read, next time raise StopIteration
        self._first_chunk = False

        if len(content) == 0: # pragma: no cover
            if self.index_col is not None:
                if np.isscalar(self.index_col):
                    index = Index([], name=self.index_name)
                else:
                    index = MultiIndex.from_arrays([[]] * len(self.index_col),
                                                   names=self.index_name)
            else:
                index = Index([])

            return DataFrame(index=index, columns=self.columns)

        zipped_content = list(lib.to_object_array(content).T)

        # no index column specified, so infer that's what is wanted
        if self.index_col is not None:
            if np.isscalar(self.index_col):
                index = zipped_content.pop(self.index_col)
            else: # given a list of index
                index = []
                for idx in self.index_col:
                    index.append(zipped_content[idx])
                # remove index items from content and columns, don't pop in loop
                for i in reversed(sorted(self.index_col)):
                    zipped_content.pop(i)

            if np.isscalar(self.index_col):
                if self.parse_dates:
                    index = lib.try_parse_dates(index, parser=self.date_parser)
                index, na_count = _convert_types(index, self.na_values)
                index = Index(index, name=self.index_name)
                if self.verbose and na_count:
                    print 'Found %d NA values in the index' % na_count
            else:
                arrays = []
                for arr in index:
                    if self.parse_dates:
                        arr = lib.try_parse_dates(arr, parser=self.date_parser)
                    arr, _ = _convert_types(arr, self.na_values)
                    arrays.append(arr)
                index = MultiIndex.from_arrays(arrays, names=self.index_name)
        else:
            index = Index(np.arange(len(content)))

        if not index._verify_integrity():
            dups = index.get_duplicates()
            raise Exception('Index has duplicates: %s' % str(dups))

        if len(self.columns) != len(zipped_content):
            raise Exception('wrong number of columns')

        data = dict((k, v) for k, v in zip(self.columns, zipped_content))

        # apply converters
        for col, f in self.converters.iteritems():
            if isinstance(col, int) and col not in self.columns:
                col = self.columns[col]
            result = np.vectorize(f)(data[col])
            if issubclass(result.dtype.type, (basestring, unicode)):
                result = result.astype('O')
            data[col] = result

        data = _convert_to_ndarrays(data, self.na_values, self.verbose)

        return DataFrame(data=data, columns=self.columns, index=index)

    def _get_lines(self, rows=None):
        source = self.data
        lines = self.buf

        # already fetched some number
        if rows is not None:
            rows -= len(self.buf)

        if isinstance(source, list):
            if self.pos > len(source):
                raise StopIteration
            if rows is None:
                lines.extend(source[self.pos:])
                self.pos = len(source)
            else:
                lines.extend(source[self.pos:self.pos+rows])
                self.pos += rows
        else:
            try:
                if rows is not None:
                    for _ in xrange(rows):
                        lines.append(source.next())
                else:
                    while True:
                        lines.append(source.next())
            except StopIteration:
                if len(lines) == 0:
                    raise
            self.pos += len(lines)

        self.buf = []

        if self.skip_footer:
            lines = lines[:-self.skip_footer]

        return lines

def _convert_to_ndarrays(dct, na_values, verbose=False):
    result = {}
    for c, values in dct.iteritems():
        cvals, na_count = _convert_types(values, na_values)
        result[c] = cvals
        if verbose and na_count:
            print 'Filled %d NA values in column %s' % (na_count, str(c))
    return result

def _convert_types(values, na_values):
    na_count = 0
    if issubclass(values.dtype.type, (np.number, np.bool_)):
        return values, na_count

    try:
        result = lib.maybe_convert_numeric(values, na_values)
    except Exception:
        na_count = lib.sanitize_objects(values, na_values)
        result = values

    if result.dtype == np.object_:
        result = lib.maybe_convert_bool(values)

    return result, na_count

#-------------------------------------------------------------------------------
# ExcelFile class


class ExcelFile(object):
    """
    Class for parsing tabular .xls sheets into DataFrame objects, uses xlrd. See
    ExcelFile.parse for more documentation

    Parameters
    ----------
    path : string
        Path to xls file
    """
    def __init__(self, path):
        import xlrd
        self.path = path
        self.book = xlrd.open_workbook(path)

    def __repr__(self):
        return object.__repr__(self)

    def parse(self, sheetname, header=0, skiprows=None, index_col=None,
              parse_dates=False, date_parser=None, na_values=None,
              chunksize=None):
        """
        Read Excel table into DataFrame

        Parameters
        ----------
        sheetname : string
            Name of Excel sheet
        header : int, default 0
            Row to use for the column labels of the parsed DataFrame
        skiprows : list-like
            Row numbers to skip (0-indexed)
        index_col : int, default 0
            Column to use as the row labels of the DataFrame. Pass None if there
            is no such column
        na_values : list-like, default None
            List of additional strings to recognize as NA/NaN

        Returns
        -------
        parsed : DataFrame
        """
        from datetime import MINYEAR, time, datetime
        from xlrd import xldate_as_tuple, XL_CELL_DATE

        datemode = self.book.datemode
        sheet = self.book.sheet_by_name(sheetname)

        data = []
        for i in range(sheet.nrows):
            row = []
            for value, typ in zip(sheet.row_values(i), sheet.row_types(i)):
                if typ == XL_CELL_DATE:
                    dt = xldate_as_tuple(value, datemode)
                    # how to produce this first case?
                    if dt[0] < MINYEAR: # pragma: no cover
                        value = time(*dt[3:])
                    else:
                        value = datetime(*dt)
                row.append(value)
            data.append(row)

        parser = TextParser(data, header=header, index_col=index_col,
                            na_values=na_values,
                            parse_dates=parse_dates,
                            date_parser=date_parser,
                            skiprows=skiprows,
                            chunksize=chunksize)

        return parser.get_chunk()
