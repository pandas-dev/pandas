"""
Module contains tools for processing files into DataFrames or other objects
"""
from StringIO import StringIO
import re
from itertools import izip
import urlparse
import csv

import numpy as np

import sys
import struct
from pandas.core.index import Index, MultiIndex
from pandas.core.frame import DataFrame
from pandas.core.series import Series
from pandas.core.categorical import Categorical
import datetime
import pandas.core.common as com
from pandas.util import py3compat
from pandas.io.date_converters import generic_parser
from pandas import isnull

from pandas.util.decorators import Appender

import pandas.lib as lib
import pandas.tslib as tslib
import pandas._parser as _parser
from pandas.tseries.period import Period
import json


class DateConversionError(Exception):
    pass

_parser_params = """Also supports optionally iterating or breaking of the file
into chunks.

Parameters
----------
filepath_or_buffer : string or file handle / StringIO. The string could be
    a URL. Valid URL schemes include http, ftp, s3, and file. For file URLs, a host
    is expected. For instance, a local file could be
    file ://localhost/path/to/table.csv
%s
lineterminator : string (length 1), default None
    Character to break file into lines. Only valid with C parser
quotechar : string
quoting : string
skipinitialspace : boolean, default False
    Skip spaces after delimiter
escapechar : string
dtype : Type name or dict of column -> type
    Data type for data or columns. E.g. {'a': np.float64, 'b': np.int32}
compression : {'gzip', 'bz2', None}, default None
    For on-the-fly decompression of on-disk data
dialect : string or csv.Dialect instance, default None
    If None defaults to Excel dialect. Ignored if sep longer than 1 char
    See csv.Dialect documentation for more details
header : int, default 0 if names parameter not specified, otherwise None
    Row to use for the column labels of the parsed DataFrame. Specify None if
    there is no header row.
skiprows : list-like or integer
    Row numbers to skip (0-indexed) or number of rows to skip (int)
    at the start of the file
index_col : int or sequence or False, default None
    Column to use as the row labels of the DataFrame. If a sequence is given, a
    MultiIndex is used. If you have a malformed file with delimiters at the end
    of each line, you might consider index_col=False to force pandas to _not_
    use the first column as the index (row names)
names : array-like
    List of column names to use. If file contains no header row, then you
    should explicitly pass header=None
prefix : string or None (default)
    Prefix to add to column numbers when no header, e.g 'X' for X0, X1, ...
na_values : list-like or dict, default None
    Additional strings to recognize as NA/NaN. If dict passed, specific
    per-column NA values
true_values : list
    Values to consider as True
false_values : list
    Values to consider as False
keep_default_na : bool, default True
    If na_values are specified and keep_default_na is False the default NaN
    values are overridden, otherwise they're appended to
parse_dates : boolean, list of ints or names, list of lists, or dict
    If True -> try parsing the index.
    If [1, 2, 3] -> try parsing columns 1, 2, 3 each as a separate date column.
    If [[1, 3]] -> combine columns 1 and 3 and parse as a single date column.
    {'foo' : [1, 3]} -> parse columns 1, 3 as date and call result 'foo'
keep_date_col : boolean, default False
    If True and parse_dates specifies combining multiple columns then
    keep the original columns.
date_parser : function
    Function to use for converting a sequence of string columns to an
    array of datetime instances. The default uses dateutil.parser.parser
    to do the conversion.
dayfirst : boolean, default False
    DD/MM format dates, international and European format
thousands : str, default None
    Thousands separator
comment : str, default None
    Indicates remainder of line should not be parsed
    Does not support line commenting (will return empty line)
decimal : str, default '.'
    Character to recognize as decimal point. E.g. use ',' for European data
nrows : int, default None
    Number of rows of file to read. Useful for reading pieces of large files
iterator : boolean, default False
    Return TextParser object
chunksize : int, default None
    Return TextParser object for iteration
skipfooter : int, default 0
    Number of line at bottom of file to skip
converters : dict. optional
    Dict of functions for converting values in certain columns. Keys can either
    be integers or column labels
verbose : boolean, default False
    Indicate number of NA values placed in non-numeric columns
delimiter : string, default None
    Alternative argument name for sep. Regular expressions are accepted.
encoding : string, default None
    Encoding to use for UTF when reading/writing (ex. 'utf-8')
squeeze : boolean, default False
    If the parsed data only contains one column then return a Series
na_filter: boolean, default True
    Detect missing value markers (empty strings and the value of na_values). In
    data without any NAs, passing na_filter=False can improve the performance
    of reading a large file

Returns
-------
result : DataFrame or TextParser
"""

_csv_sep = """sep : string, default ','
    Delimiter to use. If sep is None, will try to automatically determine
    this. Regular expressions are accepted.
"""

_table_sep = """sep : string, default \\t (tab-stop)
    Delimiter to use. Regular expressions are accepted."""

_read_csv_doc = """
Read CSV (comma-separated) file into DataFrame

%s
""" % (_parser_params % _csv_sep)

_read_table_doc = """
Read general delimited file into DataFrame

%s
""" % (_parser_params % _table_sep)

_read_stata_doc = """
Read Stata file into DataFrame

%s
""" % (_parser_params)

_fwf_widths = """\
colspecs : a list of pairs (tuples), giving the extents
    of the fixed-width fields of each line as half-open internals
    (i.e.,  [from, to[  ).
widths : a list of field widths, which can be used instead of
    'colspecs' if the intervals are contiguous.
"""

_read_fwf_doc = """
Read a table of fixed-width formatted lines into DataFrame

%s

Also, 'delimiter' is used to specify the filler character of the
fields if it is not spaces (e.g., '~').
""" % (_parser_params % _fwf_widths)


_VALID_URLS = set(urlparse.uses_relative + urlparse.uses_netloc +
                  urlparse.uses_params)
_VALID_URLS.discard('')


def _is_url(url):
    """Check to see if a URL has a valid protocol.

    Parameters
    ----------
    url : str or unicode

    Returns
    -------
    isurl : bool
        If `url` has a valid protocol return True otherwise False.
    """
    try:
        return urlparse.urlparse(url).scheme in _VALID_URLS
    except:
        return False

def _is_s3_url(url):
    """ Check for an s3 url """
    try:
        return urlparse.urlparse(url).scheme == 's3'
    except:
        return False

def _read(filepath_or_buffer, kwds):
    "Generic reader of line files."
    encoding = kwds.get('encoding', None)
    skipfooter = kwds.pop('skipfooter', None)
    if skipfooter is not None:
        kwds['skip_footer'] = skipfooter

    if isinstance(filepath_or_buffer, basestring):
        if _is_url(filepath_or_buffer):
            from urllib2 import urlopen
            filepath_or_buffer = urlopen(filepath_or_buffer)
            if py3compat.PY3:  # pragma: no cover
                if encoding:
                    errors = 'strict'
                else:
                    errors = 'replace'
                    encoding = 'utf-8'
                bytes = filepath_or_buffer.read()
                filepath_or_buffer = StringIO(bytes.decode(encoding, errors))

        if _is_s3_url(filepath_or_buffer):
            try:
                import boto
            except:
                raise ImportError("boto is required to handle s3 files")
            # Assuming AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY
            # are environment variables
            parsed_url = urlparse.urlparse(filepath_or_buffer)
            conn = boto.connect_s3()
            b = conn.get_bucket(parsed_url.netloc)
            k = boto.s3.key.Key(b)
            k.key = parsed_url.path
            filepath_or_buffer = StringIO(k.get_contents_as_string())

    if kwds.get('date_parser', None) is not None:
        if isinstance(kwds['parse_dates'], bool):
            kwds['parse_dates'] = True

    # Extract some of the arguments (pass chunksize on).
    iterator = kwds.pop('iterator', False)
    nrows = kwds.pop('nrows', None)
    chunksize = kwds.get('chunksize', None)

    # Create the parser.
    parser = TextFileReader(filepath_or_buffer, **kwds)

    if nrows is not None:
        return parser.read(nrows)
    elif chunksize or iterator:
        return parser

    return parser.read()

_parser_defaults = {
    'delimiter': None,

    'doublequote': True,
    'escapechar': None,
    'quotechar': '"',
    'quoting': csv.QUOTE_MINIMAL,
    'skipinitialspace': False,
    'lineterminator': None,

    'header': 'infer',
    'index_col': None,
    'names': None,
    'prefix': None,
    'skiprows': None,
    'na_values': None,
    'true_values': None,
    'false_values': None,
    'skip_footer': 0,
    'converters': None,

    'keep_default_na': True,
    'thousands': None,
    'comment': None,

    # 'engine': 'c',
    'parse_dates': False,
    'keep_date_col': False,
    'dayfirst': False,
    'date_parser': None,

    'usecols': None,

    # 'nrows': None,
    # 'iterator': False,
    'chunksize': None,
    'verbose': False,
    'encoding': None,
    'squeeze': False,
    'compression': None,
    'mangle_dupe_cols': True,
}


_c_parser_defaults = {
    'delim_whitespace': False,
    'as_recarray': False,
    'na_filter': True,
    'compact_ints': False,
    'use_unsigned': False,
    'low_memory': True,
    'memory_map': False,
    'buffer_lines': None,
    'error_bad_lines': True,
    'warn_bad_lines': True,
    'factorize': True,
    'dtype': None,
    'decimal': b'.'
}

_fwf_defaults = {
    'colspecs': None,
    'widths': None
}

_c_unsupported = set(['skip_footer'])
_python_unsupported = set(_c_parser_defaults.keys())


def _make_parser_function(name, sep=','):

    def parser_f(filepath_or_buffer,
                 sep=sep,
                 dialect=None,
                 compression=None,

                 doublequote=True,
                 escapechar=None,
                 quotechar='"',
                 quoting=csv.QUOTE_MINIMAL,
                 skipinitialspace=False,
                 lineterminator=None,

                 header='infer',
                 index_col=None,
                 names=None,
                 prefix=None,
                 skiprows=None,
                 skipfooter=None,
                 skip_footer=0,
                 na_values=None,
                 true_values=None,
                 false_values=None,
                 delimiter=None,
                 converters=None,
                 dtype=None,
                 usecols=None,

                 engine='c',
                 delim_whitespace=False,
                 as_recarray=False,
                 na_filter=True,
                 compact_ints=False,
                 use_unsigned=False,
                 low_memory=_c_parser_defaults['low_memory'],
                 buffer_lines=None,
                 warn_bad_lines=True,
                 error_bad_lines=True,

                 keep_default_na=True,
                 thousands=None,
                 comment=None,
                 decimal=b'.',

                 parse_dates=False,
                 keep_date_col=False,
                 dayfirst=False,
                 date_parser=None,

                 memory_map=False,
                 nrows=None,
                 iterator=False,
                 chunksize=None,

                 verbose=False,
                 encoding=None,
                 squeeze=False,
                 mangle_dupe_cols=True
                 ):

        # Alias sep -> delimiter.
        if delimiter is None:
            delimiter = sep

        kwds = dict(delimiter=delimiter,
                    engine=engine,
                    dialect=dialect,
                    compression=compression,

                    doublequote=doublequote,
                    escapechar=escapechar,
                    quotechar=quotechar,
                    quoting=quoting,
                    skipinitialspace=skipinitialspace,
                    lineterminator=lineterminator,

                    header=header,
                    index_col=index_col,
                    names=names,
                    prefix=prefix,
                    skiprows=skiprows,
                    na_values=na_values,
                    true_values=true_values,
                    false_values=false_values,
                    keep_default_na=keep_default_na,
                    thousands=thousands,
                    comment=comment,
                    decimal=decimal,

                    parse_dates=parse_dates,
                    keep_date_col=keep_date_col,
                    dayfirst=dayfirst,
                    date_parser=date_parser,

                    nrows=nrows,
                    iterator=iterator,
                    chunksize=chunksize,
                    skipfooter=skipfooter or skip_footer,
                    converters=converters,
                    dtype=dtype,
                    usecols=usecols,
                    verbose=verbose,
                    encoding=encoding,
                    squeeze=squeeze,
                    memory_map=memory_map,

                    na_filter=na_filter,
                    compact_ints=compact_ints,
                    use_unsigned=use_unsigned,
                    delim_whitespace=delim_whitespace,
                    as_recarray=as_recarray,
                    warn_bad_lines=warn_bad_lines,
                    error_bad_lines=error_bad_lines,
                    low_memory=low_memory,
                    buffer_lines=buffer_lines,
                    mangle_dupe_cols=mangle_dupe_cols
            )

        return _read(filepath_or_buffer, kwds)

    parser_f.__name__ = name

    return parser_f

read_csv = _make_parser_function('read_csv', sep=',')
read_csv = Appender(_read_csv_doc)(read_csv)

read_table = _make_parser_function('read_table', sep='\t')
read_table = Appender(_read_table_doc)(read_table)


@Appender(_read_stata_doc)
def read_stata(filepath_or_buffer, convert_dates=True, convert_categoricals=True, encoding=None, index=None):
    reader = StataReader(filepath_or_buffer, encoding)

    return reader.data(convert_dates, convert_categoricals, index)


@Appender(_read_fwf_doc)
def read_fwf(filepath_or_buffer, colspecs=None, widths=None, **kwds):
    # Check input arguments.
    if bool(colspecs is None) == bool(widths is None):
        raise ValueError("You must specify only one of 'widths' and "
                         "'colspecs'")

    # Compute 'colspec' from 'widths', if specified.
    if widths is not None:
        colspecs, col = [], 0
        for w in widths:
            colspecs.append((col, col + w))
            col += w

    kwds['colspecs'] = colspecs
    kwds['engine'] = 'python-fwf'
    return _read(filepath_or_buffer, kwds)


def read_clipboard(**kwargs):  # pragma: no cover
    """
    Read text from clipboard and pass to read_table. See read_table for the
    full argument list

    Returns
    -------
    parsed : DataFrame
    """
    from pandas.util.clipboard import clipboard_get
    text = clipboard_get()
    return read_table(StringIO(text), **kwargs)


def to_clipboard(obj):  # pragma: no cover
    """
    Attempt to write text representation of object to the system clipboard

    Notes
    -----
    Requirements for your platform
      - Linux: xsel command line tool
      - Windows: Python win32 extensions
      - OS X:
    """
    from pandas.util.clipboard import clipboard_set
    clipboard_set(str(obj))


# common NA values
# no longer excluding inf representations
# '1.#INF','-1.#INF', '1.#INF000000',
_NA_VALUES = set(['-1.#IND', '1.#QNAN', '1.#IND', '-1.#QNAN',
                 '#N/A N/A', 'NA', '#NA', 'NULL', 'NaN',
                 'nan', ''])


class TextFileReader(object):
    """

    Passed dialect overrides any of the related parser options

    """

    def __init__(self, f, engine='python', **kwds):

        self.f = f

        if kwds.get('dialect') is not None:
            dialect = kwds['dialect']
            kwds['delimiter'] = dialect.delimiter
            kwds['doublequote'] = dialect.doublequote
            kwds['escapechar'] = dialect.escapechar
            kwds['skipinitialspace'] = dialect.skipinitialspace
            kwds['quotechar'] = dialect.quotechar
            kwds['quoting'] = dialect.quoting

        if kwds.get('header', 'infer') == 'infer':
            kwds['header'] = 0 if kwds.get('names') is None else None

        self.orig_options = kwds

        # miscellanea
        self.engine = engine
        self._engine = None

        options = self._get_options_with_defaults(engine)

        self.chunksize = options.pop('chunksize', None)
        self.squeeze = options.pop('squeeze', False)

        # might mutate self.engine
        self.options, self.engine = self._clean_options(options, engine)
        if 'has_index_names' in kwds:
            self.options['has_index_names'] = kwds['has_index_names']

        self._make_engine(self.engine)

    def _get_options_with_defaults(self, engine):
        kwds = self.orig_options

        options = {}
        for argname, default in _parser_defaults.iteritems():
            if argname in kwds:
                value = kwds[argname]
            else:
                value = default

            options[argname] = value

        for argname, default in _c_parser_defaults.iteritems():
            if argname in kwds:
                value = kwds[argname]
                if engine != 'c' and value != default:
                    raise ValueError('%s is not supported with %s parser' %
                                     (argname, engine))
            options[argname] = value

        if engine == 'python-fwf':
            for argname, default in _fwf_defaults.iteritems():
                if argname in kwds:
                    value = kwds[argname]
                options[argname] = value

        return options

    def _clean_options(self, options, engine):
        result = options.copy()

        sep = options['delimiter']
        if (sep is None and not options['delim_whitespace']):
            if engine == 'c':
                print 'Using Python parser to sniff delimiter'
                engine = 'python'
        elif sep is not None and len(sep) > 1:
            # wait until regex engine integrated
            engine = 'python'

        # C engine not supported yet
        if engine == 'c':
            if options['skip_footer'] > 0:
                engine = 'python'

        if engine == 'c':
            for arg in _c_unsupported:
                del result[arg]

        if 'python' in engine:
            for arg in _python_unsupported:
                del result[arg]

        index_col = options['index_col']
        names = options['names']
        converters = options['converters']
        na_values = options['na_values']
        skiprows = options['skiprows']

        # really delete this one
        keep_default_na = result.pop('keep_default_na')

        if _is_index_col(index_col):
            if not isinstance(index_col, (list, tuple, np.ndarray)):
                index_col = [index_col]
        result['index_col'] = index_col

        names = list(names) if names is not None else names

        # type conversion-related
        if converters is not None:
            if not (isinstance(converters, dict)):
                raise AssertionError()
        else:
            converters = {}

        # Converting values to NA
        na_values = _clean_na_values(na_values, keep_default_na)

        if com.is_integer(skiprows):
            skiprows = range(skiprows)
        skiprows = set() if skiprows is None else set(skiprows)

        # put stuff back
        result['index_col'] = index_col
        result['names'] = names
        result['converters'] = converters
        result['na_values'] = na_values
        result['skiprows'] = skiprows

        return result, engine

    def __iter__(self):
        try:
            while True:
                yield self.read(self.chunksize)
        except StopIteration:
            pass

    def _make_engine(self, engine='c'):
        if engine == 'c':
            self._engine = CParserWrapper(self.f, **self.options)
        else:
            if engine == 'python':
                klass = PythonParser
            elif engine == 'python-fwf':
                klass = FixedWidthFieldParser
            self._engine = klass(self.f, **self.options)

    def _failover_to_python(self):
        raise NotImplementedError

    def read(self, nrows=None):
        suppressed_warnings = False
        if nrows is not None:
            if self.options.get('skip_footer'):
                raise ValueError('skip_footer not supported for iteration')

            # # XXX hack
            # if isinstance(self._engine, CParserWrapper):
            #     suppressed_warnings = True
            #     self._engine.set_error_bad_lines(False)

        ret = self._engine.read(nrows)

        if self.options.get('as_recarray'):
            return ret

        index, columns, col_dict = ret

        # May alter columns / col_dict
        # index, columns, col_dict = self._create_index(col_dict, columns)

        df = DataFrame(col_dict, columns=columns, index=index)

        if self.squeeze and len(df.columns) == 1:
            return df[df.columns[0]]
        return df

    def _create_index(self, col_dict, columns):
        pass

    def get_chunk(self, size=None):
        if size is None:
            size = self.chunksize
        return self.read(nrows=size)

def _is_index_col(col):
    return col is not None and col is not False


class ParserBase(object):

    def __init__(self, kwds):
        self.names = kwds.get('names')
        self.orig_names = None
        self.prefix = kwds.pop('prefix', None)

        self.index_col = kwds.pop('index_col', None)
        self.index_names = None

        self.parse_dates = kwds.pop('parse_dates', False)
        self.date_parser = kwds.pop('date_parser', None)
        self.dayfirst = kwds.pop('dayfirst', False)
        self.keep_date_col = kwds.pop('keep_date_col', False)

        self.na_values = kwds.get('na_values')
        self.true_values = kwds.get('true_values')
        self.false_values = kwds.get('false_values')

        self._date_conv = _make_date_converter(date_parser=self.date_parser,
                                               dayfirst=self.dayfirst)

        self._name_processed = False

    @property
    def _has_complex_date_col(self):
        return (isinstance(self.parse_dates, dict) or
                (isinstance(self.parse_dates, list) and
                 len(self.parse_dates) > 0 and
                 isinstance(self.parse_dates[0], list)))

    def _should_parse_dates(self, i):
        if isinstance(self.parse_dates, bool):
            return self.parse_dates
        else:
            name = self.index_names[i]
            j = self.index_col[i]

            if np.isscalar(self.parse_dates):
                return (j == self.parse_dates) or (name == self.parse_dates)
            else:
                return (j in self.parse_dates) or (name in self.parse_dates)

    def _make_index(self, data, alldata, columns):
        if not _is_index_col(self.index_col) or len(self.index_col) == 0:
            index = None

        elif not self._has_complex_date_col:
            index = self._get_simple_index(alldata, columns)
            index = self._agg_index(index)

        elif self._has_complex_date_col:
            if not self._name_processed:
                (self.index_names, _,
                 self.index_col) = _clean_index_names(list(columns),
                                                      self.index_col)
                self._name_processed = True
            index = self._get_complex_date_index(data, columns)
            index = self._agg_index(index, try_parse_dates=False)

        return index

    _implicit_index = False

    def _get_simple_index(self, data, columns):
        def ix(col):
            if not isinstance(col, basestring):
                return col
            raise ValueError('Index %s invalid' % col)
        index = None

        to_remove = []
        index = []
        for idx in self.index_col:
            i = ix(idx)
            to_remove.append(i)
            index.append(data[i])

        # remove index items from content and columns, don't pop in
        # loop
        for i in reversed(sorted(to_remove)):
            data.pop(i)
            if not self._implicit_index:
                columns.pop(i)

        return index

    def _get_complex_date_index(self, data, col_names):
        def _get_name(icol):
            if isinstance(icol, basestring):
                return icol

            if col_names is None:
                raise ValueError(('Must supply column order to use %s as '
                                  'index') % str(icol))

            for i, c in enumerate(col_names):
                if i == icol:
                    return c

        index = None

        to_remove = []
        index = []
        for idx in self.index_col:
            name = _get_name(idx)
            to_remove.append(name)
            index.append(data[name])

        # remove index items from content and columns, don't pop in
        # loop
        for c in reversed(sorted(to_remove)):
            data.pop(c)
            col_names.remove(c)

        return index

    def _agg_index(self, index, try_parse_dates=True):
        arrays = []
        for i, arr in enumerate(index):

            if (try_parse_dates and self._should_parse_dates(i)):
                arr = self._date_conv(arr)

            col_na_values = self.na_values

            if isinstance(self.na_values, dict):
                col_name = self.index_names[i]
                if col_name is not None:
                    col_na_values = _get_na_values(col_name,
                                                   self.na_values)

            arr, _ = self._convert_types(arr, col_na_values)
            arrays.append(arr)

        index = MultiIndex.from_arrays(arrays, names=self.index_names)

        return index

    def _convert_to_ndarrays(self, dct, na_values, verbose=False,
                             converters=None):
        result = {}
        for c, values in dct.iteritems():
            conv_f = None if converters is None else converters.get(c, None)
            col_na_values = _get_na_values(c, na_values)
            coerce_type = True
            if conv_f is not None:
                values = lib.map_infer(values, conv_f)
                coerce_type = False
            cvals, na_count = self._convert_types(values, col_na_values,
                                                  coerce_type)
            result[c] = cvals
            if verbose and na_count:
                print 'Filled %d NA values in column %s' % (na_count, str(c))
        return result

    def _convert_types(self, values, na_values, try_num_bool=True):
        na_count = 0
        if issubclass(values.dtype.type, (np.number, np.bool_)):
            mask = lib.ismember(values, na_values)
            na_count = mask.sum()
            if na_count > 0:
                if com.is_integer_dtype(values):
                    values = values.astype(np.float64)
                np.putmask(values, mask, np.nan)
            return values, na_count

        if try_num_bool:
            try:
                result = lib.maybe_convert_numeric(values, na_values, False)
            except Exception:
                result = values
                if values.dtype == np.object_:
                    na_count = lib.sanitize_objects(result, na_values, False)
        else:
            result = values
            if values.dtype == np.object_:
                na_count = lib.sanitize_objects(values, na_values, False)

        if result.dtype == np.object_ and try_num_bool:
            result = lib.maybe_convert_bool(values,
                                            true_values=self.true_values,
                                            false_values=self.false_values)

        return result, na_count

    def _do_date_conversions(self, names, data):
        # returns data, columns
        if self.parse_dates is not None:
            data, names = _process_date_conversion(
                data, self._date_conv, self.parse_dates, self.index_col,
                self.index_names, names, keep_date_col=self.keep_date_col)

        return names, data

    def _exclude_implicit_index(self, alldata):

        if self._implicit_index:
            excl_indices = self.index_col

            data = {}
            offset = 0
            for i, col in enumerate(self.orig_names):
                while i + offset in excl_indices:
                    offset += 1
                data[col] = alldata[i + offset]
        else:
            data = dict((k, v) for k, v in izip(self.orig_names, alldata))

        return data


class CParserWrapper(ParserBase):
    """

    """

    def __init__(self, src, **kwds):
        self.kwds = kwds
        kwds = kwds.copy()

        self.as_recarray = kwds.get('as_recarray', False)
        ParserBase.__init__(self, kwds)

        if 'utf-16' in (kwds.get('encoding') or ''):
            if isinstance(src, basestring):
                src = open(src, 'rb')
            src = com.UTF8Recoder(src, kwds['encoding'])
            kwds['encoding'] = 'utf-8'

        # #2442
        kwds['allow_leading_cols'] = self.index_col is not False
        self._reader = _parser.TextReader(src, **kwds)

        # XXX
        self.usecols = self._reader.usecols

        passed_names = self.names is None

        if self._reader.header is None:
            self.names = None
        else:
            self.names = list(self._reader.header)

        if self.names is None:
            if self.prefix:
                self.names = ['X%d' % i
                              for i in range(self._reader.table_width)]
            else:
                self.names = range(self._reader.table_width)

        # XXX
        self._set_noconvert_columns()

        self.orig_names = self.names

        if not self._has_complex_date_col:
            if (self._reader.leading_cols == 0 and
                    _is_index_col(self.index_col)):

                self._name_processed = True
                (self.index_names, self.names,
                 self.index_col) = _clean_index_names(self.names,
                                                      self.index_col)

            if self._reader.header is None and not passed_names:
                self.index_names = [None] * len(self.index_names)

        self._implicit_index = self._reader.leading_cols > 0

    def _set_noconvert_columns(self):
        names = self.names

        def _set(x):
            if com.is_integer(x):
                self._reader.set_noconvert(x)
            else:
                self._reader.set_noconvert(names.index(x))

        if isinstance(self.parse_dates, list):
            for val in self.parse_dates:
                if isinstance(val, list):
                    for k in val:
                        _set(k)
                else:
                    _set(val)

    def set_error_bad_lines(self, status):
        self._reader.set_error_bad_lines(int(status))

    def read(self, nrows=None):
        if self.as_recarray:
            # what to do if there are leading columns?
            return self._reader.read(nrows)

        try:
            data = self._reader.read(nrows)
        except StopIteration:
            if nrows is None:
                return None, self.names, {}
            else:
                raise

        names = self.names

        if self._reader.leading_cols:
            if self._has_complex_date_col:
                raise NotImplementedError('file structure not yet supported')

            # implicit index, no index names
            arrays = []

            for i in range(self._reader.leading_cols):
                if self.index_col is None:
                    values = data.pop(i)
                else:
                    values = data.pop(self.index_col[i])

                values = self._maybe_parse_dates(values, i,
                                                 try_parse_dates=True)
                arrays.append(values)

            index = MultiIndex.from_arrays(arrays)

            if self.usecols is not None:
                names = self._filter_usecols(names)

            # rename dict keys
            data = sorted(data.items())
            data = dict((k, v) for k, (i, v) in zip(names, data))

            names, data = self._do_date_conversions(names, data)

        else:
            # rename dict keys
            data = sorted(data.items())

            # ugh, mutation
            names = list(self.orig_names)

            if self.usecols is not None:
                names = self._filter_usecols(names)

            # columns as list
            alldata = [x[1] for x in data]

            data = dict((k, v) for k, (i, v) in zip(names, data))

            names, data = self._do_date_conversions(names, data)
            index = self._make_index(data, alldata, names)

        return index, names, data

    def _filter_usecols(self, names):
        # hackish
        if self.usecols is not None and len(names) != len(self.usecols):
            names = [name for i, name in enumerate(names)
                     if i in self.usecols or name in self.usecols]
        return names

    def _get_index_names(self):
        names = list(self._reader.header)
        idx_names = None

        if self._reader.leading_cols == 0 and self.index_col is not None:
            (idx_names, names,
             self.index_col) = _clean_index_names(names, self.index_col)

        return names, idx_names

    def _maybe_parse_dates(self, values, index, try_parse_dates=True):
        if try_parse_dates and self._should_parse_dates(index):
            values = self._date_conv(values)
        return values


def TextParser(*args, **kwds):
    """
    Converts lists of lists/tuples into DataFrames with proper type inference
    and optional (e.g. string to datetime) conversion. Also enables iterating
    lazily over chunks of large files

    Parameters
    ----------
    data : file-like object or list
    delimiter : separator character to use
    dialect : str or csv.Dialect instance, default None
        Ignored if delimiter is longer than 1 character
    names : sequence, default
    header : int, default 0
        Row to use to parse column labels. Defaults to the first row. Prior
        rows will be discarded
    index_col : int or list, default None
        Column or columns to use as the (possibly hierarchical) index
    has_index_names: boolean, default False
        True if the cols defined in index_col have an index name and are
        not in the header
    na_values : iterable, default None
        Custom NA values
    keep_default_na : bool, default True
    thousands : str, default None
        Thousands separator
    comment : str, default None
        Comment out remainder of line
    parse_dates : boolean, default False
    keep_date_col : boolean, default False
    date_parser : function, default None
    skiprows : list of integers
        Row numbers to skip
    skip_footer : int
        Number of line at bottom of file to skip
    encoding : string, default None
        Encoding to use for UTF when reading/writing (ex. 'utf-8')
    squeeze : boolean, default False
        returns Series if only one column
    """
    kwds['engine'] = 'python'
    return TextFileReader(*args, **kwds)

# delimiter=None, dialect=None, names=None, header=0,
# index_col=None,
# na_values=None,
# na_filter=True,
# thousands=None,
# quotechar='"',
# escapechar=None,
# doublequote=True,
# skipinitialspace=False,
# quoting=csv.QUOTE_MINIMAL,
# comment=None, parse_dates=False, keep_date_col=False,
# date_parser=None, dayfirst=False,
# chunksize=None, skiprows=None, skip_footer=0, converters=None,
# verbose=False, encoding=None, squeeze=False):


def count_empty_vals(vals):
    return sum([1 for v in vals if v == '' or v is None])


def _wrap_compressed(f, compression):
    compression = compression.lower()
    if compression == 'gzip':
        import gzip
        return gzip.GzipFile(fileobj=f)
    elif compression == 'bz2':
        raise ValueError('Python cannot read bz2 data from file handle')
    else:
        raise ValueError('do not recognize compression method %s'
                         % compression)


class PythonParser(ParserBase):

    def __init__(self, f, **kwds):
        """
        Workhorse function for processing nested list into DataFrame

        Should be replaced by np.genfromtxt eventually?
        """
        ParserBase.__init__(self, kwds)

        self.data = None
        self.buf = []
        self.pos = 0

        if kwds['usecols'] is not None:
            raise Exception("usecols not supported with engine='python'"
                            " or multicharacter separators (yet).")

        self.header = kwds['header']
        self.encoding = kwds['encoding']
        self.compression = kwds['compression']
        self.skiprows = kwds['skiprows']

        self.skip_footer = kwds['skip_footer']
        self.delimiter = kwds['delimiter']

        self.quotechar = kwds['quotechar']
        self.escapechar = kwds['escapechar']
        self.doublequote = kwds['doublequote']
        self.skipinitialspace = kwds['skipinitialspace']
        self.lineterminator = kwds['lineterminator']
        self.quoting = kwds['quoting']
        self.mangle_dupe_cols = kwds.get('mangle_dupe_cols',True)

        self.has_index_names = False
        if 'has_index_names' in kwds:
            self.has_index_names = kwds['has_index_names']

        self.verbose = kwds['verbose']
        self.converters = kwds['converters']

        self.thousands = kwds['thousands']
        self.comment = kwds['comment']
        self._comment_lines = []

        if isinstance(f, basestring):
            f = com._get_handle(f, 'r', encoding=self.encoding,
                                compression=self.compression)
        elif self.compression:
            f = _wrap_compressed(f, self.compression)

        if hasattr(f, 'readline'):
            self._make_reader(f)
        else:
            self.data = f
        self.columns = self._infer_columns()

        # get popped off for index
        self.orig_names = list(self.columns)

        # needs to be cleaned/refactored
        # multiple date column thing turning into a real spaghetti factory

        if not self._has_complex_date_col:
            (self.index_names,
             self.orig_names, _) = self._get_index_name(self.columns)
            self._name_processed = True
        self._first_chunk = True

    def _make_reader(self, f):
        sep = self.delimiter

        if sep is None or len(sep) == 1:
            if self.lineterminator:
                raise ValueError('Custom line terminators not supported in '
                                 'python parser (yet)')

            class MyDialect(csv.Dialect):
                delimiter = self.delimiter
                quotechar = self.quotechar
                escapechar = self.escapechar
                doublequote = self.doublequote
                skipinitialspace = self.skipinitialspace
                quoting = self.quoting
                lineterminator = '\n'

            dia = MyDialect

            sniff_sep = True

            if sep is not None:
                sniff_sep = False
                dia.delimiter = sep
            # attempt to sniff the delimiter
            if sniff_sep:
                line = f.readline()
                while self.pos in self.skiprows:
                    self.pos += 1
                    line = f.readline()

                line = self._check_comments([line])[0]

                self.pos += 1
                sniffed = csv.Sniffer().sniff(line)
                dia.delimiter = sniffed.delimiter
                if self.encoding is not None:
                    self.buf.extend(list(
                        com.UnicodeReader(StringIO(line),
                                          dialect=dia,
                                          encoding=self.encoding)))
                else:
                    self.buf.extend(list(csv.reader(StringIO(line),
                                                    dialect=dia)))

            if self.encoding is not None:
                reader = com.UnicodeReader(f, dialect=dia,
                                           encoding=self.encoding,
                                           strict=True)
            else:
                reader = csv.reader(f, dialect=dia,
                                    strict=True)

        else:
            def _read():
                line = next(f)
                pat = re.compile(sep)
                if (py3compat.PY3 and isinstance(line, bytes)):
                    yield pat.split(line.decode('utf-8').strip())
                    for line in f:
                        yield pat.split(line.decode('utf-8').strip())
                else:
                    yield pat.split(line.strip())
                    for line in f:
                        yield pat.split(line.strip())
            reader = _read()

        self.data = reader

    def read(self, rows=None):
        try:
            content = self._get_lines(rows)
        except StopIteration:
            if self._first_chunk:
                content = []
            else:
                raise

        # done with first read, next time raise StopIteration
        self._first_chunk = False

        columns = list(self.orig_names)
        if len(content) == 0:  # pragma: no cover
            # DataFrame with the right metadata, even though it's length 0
            return _get_empty_meta(self.orig_names,
                                   self.index_col,
                                   self.index_names)

        # handle new style for names in index
        count_empty_content_vals = count_empty_vals(content[0])
        indexnamerow = None
        if self.has_index_names and count_empty_content_vals == len(columns):
            indexnamerow = content[0]
            content = content[1:]

        alldata = self._rows_to_cols(content)
        data = self._exclude_implicit_index(alldata)

        columns, data = self._do_date_conversions(self.columns, data)

        data = self._convert_data(data)
        index = self._make_index(data, alldata, columns)
        if indexnamerow:
            coffset = len(indexnamerow) - len(columns)
            index.names = indexnamerow[:coffset]

        return index, columns, data

    # legacy
    def get_chunk(self, size=None):
        if size is None:
            size = self.chunksize
        return self.read(nrows=size)

    def _convert_data(self, data):
        # apply converters
        clean_conv = {}

        for col, f in self.converters.iteritems():
            if isinstance(col, int) and col not in self.orig_names:
                col = self.orig_names[col]
            clean_conv[col] = f

        return self._convert_to_ndarrays(data, self.na_values, self.verbose,
                                         clean_conv)

    def _infer_columns(self):
        names = self.names

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

            if self.mangle_dupe_cols:
                counts = {}
                for i, col in enumerate(columns):
                    cur_count = counts.get(col, 0)
                    if cur_count > 0:
                        columns[i] = '%s.%d' % (col, cur_count)
                    counts[col] = cur_count + 1

            self._clear_buffer()

            if names is not None:
                if len(names) != len(columns):
                    raise Exception('Number of passed names did not match '
                                    'number of header fields in the file')
                columns = names
        else:
            if len(self.buf) > 0:
                line = self.buf[0]
            else:
                line = self._next_line()

            ncols = len(line)
            if not names:
                if self.prefix:
                    columns = ['X%d' % i for i in range(ncols)]
                else:
                    columns = range(ncols)
            else:
                columns = names

        return columns

    def _next_line(self):
        if isinstance(self.data, list):
            while self.pos in self.skiprows:
                self.pos += 1

            try:
                line = self.data[self.pos]
            except IndexError:
                raise StopIteration
        else:
            while self.pos in self.skiprows:
                next(self.data)
                self.pos += 1

            line = next(self.data)

        line = self._check_comments([line])[0]
        line = self._check_thousands([line])[0]

        self.pos += 1
        self.buf.append(line)

        return line

    def _check_comments(self, lines):
        if self.comment is None:
            return lines
        ret = []
        for l in lines:
            rl = []
            for x in l:
                if (not isinstance(x, basestring) or
                        self.comment not in x):
                    rl.append(x)
                else:
                    x = x[:x.find(self.comment)]
                    if len(x) > 0:
                        rl.append(x)
                    break
            ret.append(rl)
        return ret

    def _check_thousands(self, lines):
        if self.thousands is None:
            return lines
        nonnum = re.compile('[^-^0-9^%s^.]+' % self.thousands)
        ret = []
        for l in lines:
            rl = []
            for x in l:
                if (not isinstance(x, basestring) or
                    self.thousands not in x or
                        nonnum.search(x.strip())):
                    rl.append(x)
                else:
                    rl.append(x.replace(',', ''))
            ret.append(rl)
        return ret

    def _clear_buffer(self):
        self.buf = []

    _implicit_index = False

    def _get_index_name(self, columns):
        orig_names = list(columns)
        columns = list(columns)

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
            # leave it 0, #2442
            if self.index_col is not False:
                implicit_first_cols = len(line) - len(columns)

            if next_line is not None:
                if len(next_line) == len(line) + len(columns):
                    # column and index names on diff rows
                    implicit_first_cols = 0

                    self.index_col = range(len(line))
                    self.buf = self.buf[1:]

                    for c in reversed(line):
                        columns.insert(0, c)

                    return line, columns, orig_names

        if implicit_first_cols > 0:
            self._implicit_index = True
            if self.index_col is None:
                self.index_col = range(implicit_first_cols)
            index_name = None

        else:
            (index_name, columns,
             self.index_col) = _clean_index_names(columns, self.index_col)

        return index_name, orig_names, columns

    def _rows_to_cols(self, content):
        zipped_content = list(lib.to_object_array(content).T)

        col_len = len(self.orig_names)
        zip_len = len(zipped_content)

        if self._implicit_index:
            col_len += len(self.index_col)

        if not ((self.skip_footer >= 0)):
            raise AssertionError()

        if col_len != zip_len and self.index_col is not False:
            row_num = -1
            i = 0
            for (i, l) in enumerate(content):
                if len(l) != col_len:
                    break

            footers = 0
            if self.skip_footer:
                footers = self.skip_footer

            row_num = self.pos - (len(content) - i + footers)

            msg = ('Expected %d fields in line %d, saw %d' %
                   (col_len, row_num + 1, zip_len))
            raise ValueError(msg)

        return zipped_content

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
                lines.extend(source[self.pos:self.pos + rows])
                self.pos += rows
        else:
            new_rows = []
            try:
                if rows is not None:
                    for _ in xrange(rows):
                        new_rows.append(next(source))
                    lines.extend(new_rows)
                else:
                    rows = 0
                    while True:
                        try:
                            new_rows.append(next(source))
                            rows += 1
                        except csv.Error, inst:
                            if 'newline inside string' in str(inst):
                                row_num = str(self.pos + rows)
                                msg = ('EOF inside string starting with line '
                                       + row_num)
                                raise Exception(msg)
                            raise
            except StopIteration:
                lines.extend(new_rows)
                if len(lines) == 0:
                    raise
            self.pos += len(new_rows)

        self.buf = []

        if self.skip_footer:
            lines = lines[:-self.skip_footer]

        lines = self._check_comments(lines)
        return self._check_thousands(lines)


def _make_date_converter(date_parser=None, dayfirst=False):
    def converter(*date_cols):
        if date_parser is None:
            strs = _concat_date_cols(date_cols)
            try:
                return tslib.array_to_datetime(com._ensure_object(strs),
                                               utc=None, dayfirst=dayfirst)
            except:
                return lib.try_parse_dates(strs, dayfirst=dayfirst)
        else:
            try:
                result = date_parser(*date_cols)
                if isinstance(result, datetime.datetime):
                    raise Exception('scalar parser')
                return result
            except Exception:
                try:
                    return lib.try_parse_dates(_concat_date_cols(date_cols),
                                               parser=date_parser,
                                               dayfirst=dayfirst)
                except Exception:
                    return generic_parser(date_parser, *date_cols)

    return converter


def _process_date_conversion(data_dict, converter, parse_spec,
                             index_col, index_names, columns,
                             keep_date_col=False):
    def _isindex(colspec):
        return ((isinstance(index_col, list) and
                 colspec in index_col)
                or (isinstance(index_names, list) and
                    colspec in index_names))

    new_cols = []
    new_data = {}

    orig_names = columns
    columns = list(columns)

    date_cols = set()

    if parse_spec is None or isinstance(parse_spec, bool):
        return data_dict, columns

    if isinstance(parse_spec, list):
        # list of column lists
        for colspec in parse_spec:
            if np.isscalar(colspec):
                if isinstance(colspec, int) and colspec not in data_dict:
                    colspec = orig_names[colspec]
                if _isindex(colspec):
                    continue
                data_dict[colspec] = converter(data_dict[colspec])
            else:
                new_name, col, old_names = _try_convert_dates(
                    converter, colspec, data_dict, orig_names)
                if new_name in data_dict:
                    raise ValueError('New date column already in dict %s' %
                                     new_name)
                new_data[new_name] = col
                new_cols.append(new_name)
                date_cols.update(old_names)

    elif isinstance(parse_spec, dict):
        # dict of new name to column list
        for new_name, colspec in parse_spec.iteritems():
            if new_name in data_dict:
                raise ValueError('Date column %s already in dict' %
                                 new_name)

            _, col, old_names = _try_convert_dates(converter, colspec,
                                                   data_dict, orig_names)

            new_data[new_name] = col
            new_cols.append(new_name)
            date_cols.update(old_names)

    data_dict.update(new_data)
    new_cols.extend(columns)

    if not keep_date_col:
        for c in list(date_cols):
            data_dict.pop(c)
            new_cols.remove(c)

    return data_dict, new_cols


def _try_convert_dates(parser, colspec, data_dict, columns):
    colset = set(columns)
    colnames = []

    for c in colspec:
        if c in colset:
            colnames.append(c)
        elif isinstance(c, int) and c not in columns:
            colnames.append(str(columns[c]))
        else:
            colnames.append(c)

    new_name = '_'.join([str(x) for x in colnames])
    to_parse = [data_dict[c] for c in colnames if c in data_dict]

    try:
        new_col = parser(*to_parse)
    except DateConversionError:
        new_col = parser(_concat_date_cols(to_parse))
    return new_name, new_col, colnames


def _clean_na_values(na_values, keep_default_na=True):
    if na_values is None and keep_default_na:
        na_values = _NA_VALUES
    elif isinstance(na_values, dict):
        if keep_default_na:
            for k, v in na_values.iteritems():
                v = set(list(v)) | _NA_VALUES
                na_values[k] = v
    else:
        if not com.is_list_like(na_values):
            na_values = [na_values]
        na_values = set(list(na_values))
        if keep_default_na:
            na_values = na_values | _NA_VALUES

    return na_values


def _clean_index_names(columns, index_col):
    if not _is_index_col(index_col):
        return None, columns, index_col

    columns = list(columns)

    cp_cols = list(columns)
    index_names = []

    # don't mutate
    index_col = list(index_col)

    for i, c in enumerate(index_col):
        if isinstance(c, basestring):
            index_names.append(c)
            for j, name in enumerate(cp_cols):
                if name == c:
                    index_col[i] = j
                    columns.remove(name)
                    break
        else:
            name = cp_cols[c]
            columns.remove(name)
            index_names.append(name)

    # hack
    if isinstance(index_names[0], basestring) and 'Unnamed' in index_names[0]:
        index_names[0] = None

    return index_names, columns, index_col


def _get_empty_meta(columns, index_col, index_names):
    columns = list(columns)

    if index_col is not None:
        index = MultiIndex.from_arrays([[]] * len(index_col),
                                       names=index_names)
        for n in index_col:
            columns.pop(n)
    else:
        index = Index([])

    return index, columns, {}


def _get_na_values(col, na_values):
    if isinstance(na_values, dict):
        if col in na_values:
            return set(list(na_values[col]))
        else:
            return _NA_VALUES
    else:
        return na_values


def _get_col_names(colspec, columns):
    colset = set(columns)
    colnames = []
    for c in colspec:
        if c in colset:
            colnames.append(c)
        elif isinstance(c, int):
            colnames.append(columns[c])
    return colnames


def _concat_date_cols(date_cols):
    if len(date_cols) == 1:
        if py3compat.PY3:
            return np.array([unicode(x) for x in date_cols[0]], dtype=object)
        else:
            return np.array([str(x) if not isinstance(x, basestring) else x
                             for x in date_cols[0]], dtype=object)

    # stripped = [map(str.strip, x) for x in date_cols]
    rs = np.array([' '.join([unicode(y) for y in x])
                   for x in zip(*date_cols)], dtype=object)
    return rs


class FixedWidthReader(object):
    """
    A reader of fixed-width lines.
    """
    def __init__(self, f, colspecs, filler, thousands=None):
        self.f = f
        self.colspecs = colspecs
        self.filler = filler  # Empty characters between fields.
        self.thousands = thousands

        if not ( isinstance(colspecs, (tuple, list))):
            raise AssertionError()

        for colspec in colspecs:
            if not ( isinstance(colspec, (tuple, list)) and
                       len(colspec) == 2 and
                       isinstance(colspec[0], int) and
                       isinstance(colspec[1], int) ):
                raise AssertionError()

    def next(self):
        line = next(self.f)
        # Note: 'colspecs' is a sequence of half-open intervals.
        return [line[fromm:to].strip(self.filler or ' ')
                for (fromm, to) in self.colspecs]

    # Iterator protocol in Python 3 uses __next__()
    __next__ = next


class FixedWidthFieldParser(PythonParser):
    """
    Specialization that Converts fixed-width fields into DataFrames.
    See PythonParser for details.
    """
    def __init__(self, f, **kwds):
        # Support iterators, convert to a list.
        self.colspecs = list(kwds.pop('colspecs'))

        PythonParser.__init__(self, f, **kwds)

    def _make_reader(self, f):
        self.data = FixedWidthReader(f, self.colspecs, self.delimiter)


#----------------------------------------------------------------------
# ExcelFile class

class ExcelFile(object):
    """
    Class for parsing tabular excel sheets into DataFrame objects.
    Uses xlrd. See ExcelFile.parse for more documentation

    Parameters
    ----------
    path : string or file-like object
        Path to xls or xlsx file
    """
    def __init__(self, path_or_buf, kind=None, **kwds):
        self.kind = kind

        import xlrd # throw an ImportError if we need to
        ver = tuple(map(int,xlrd.__VERSION__.split(".")[:2]))
        if ver < (0, 9):
            raise ImportError("pandas requires xlrd >= 0.9.0 for excel support")

        self.path_or_buf = path_or_buf
        self.tmpfile = None

        if isinstance(path_or_buf, basestring):
            self.book = xlrd.open_workbook(path_or_buf)
        else:
            data = path_or_buf.read()
            self.book = xlrd.open_workbook(file_contents=data)

    def __repr__(self):
        return object.__repr__(self)

    def parse(self, sheetname, header=0, skiprows=None, skip_footer=0,
              index_col=None, parse_cols=None, parse_dates=False,
              date_parser=None, na_values=None, thousands=None, chunksize=None,
              **kwds):
        """
        Read Excel table into DataFrame

        Parameters
        ----------
        sheetname : string
            Name of Excel sheet
        header : int, default 0
            Row to use for the column labels of the parsed DataFrame
        skiprows : list-like
            Rows to skip at the beginning (0-indexed)
        skip_footer : int, default 0
            Rows at the end to skip (0-indexed)
        index_col : int, default None
            Column to use as the row labels of the DataFrame. Pass None if
            there is no such column
        parse_cols : int or list, default None
            If None then parse all columns,
            If int then indicates last column to be parsed
            If list of ints then indicates list of column numbers to be parsed
            If string then indicates comma separated list of column names and
                column ranges (e.g. "A:E" or "A,C,E:F")
        na_values : list-like, default None
            List of additional strings to recognize as NA/NaN

        Returns
        -------
        parsed : DataFrame
        """

        # has_index_names: boolean, default False
        #     True if the cols defined in index_col have an index name and are
        #     not in the header
        has_index_names = False  # removed as new argument of API function

        skipfooter = kwds.pop('skipfooter', None)
        if skipfooter is not None:
            skip_footer = skipfooter

        return  self._parse_excel(sheetname, header=header,
                                     skiprows=skiprows, index_col=index_col,
                                     has_index_names=has_index_names,
                                     parse_cols=parse_cols,
                                     parse_dates=parse_dates,
                                     date_parser=date_parser,
                                     na_values=na_values,
                                     thousands=thousands,
                                     chunksize=chunksize,
                                     skip_footer=skip_footer)

    def _should_parse(self, i, parse_cols):

        def _range2cols(areas):
            """
            Convert comma separated list of column names and column ranges to a
            list of 0-based column indexes.

            >>> _range2cols('A:E')
            [0, 1, 2, 3, 4]
            >>> _range2cols('A,C,Z:AB')
            [0, 2, 25, 26, 27]
            """
            def _excel2num(x):
                "Convert Excel column name like 'AB' to 0-based column index"
                return reduce(lambda s, a: s * 26 + ord(a) - ord('A') + 1, x.upper().strip(), 0) - 1

            cols = []
            for rng in areas.split(','):
                if ':' in rng:
                    rng = rng.split(':')
                    cols += range(_excel2num(rng[0]), _excel2num(rng[1]) + 1)
                else:
                    cols.append(_excel2num(rng))
            return cols

        if isinstance(parse_cols, int):
            return i <= parse_cols
        elif isinstance(parse_cols, basestring):
            return i in _range2cols(parse_cols)
        else:
            return i in parse_cols

    def _parse_excel(self, sheetname, header=0, skiprows=None,
                   skip_footer=0, index_col=None, has_index_names=None,
                   parse_cols=None, parse_dates=False, date_parser=None,
                   na_values=None, thousands=None, chunksize=None):
        from xlrd import (xldate_as_tuple, XL_CELL_DATE,
                          XL_CELL_ERROR, XL_CELL_BOOLEAN)

        datemode = self.book.datemode
        sheet = self.book.sheet_by_name(sheetname)

        data = []
        should_parse = {}
        for i in range(sheet.nrows):
            row = []
            for j, (value, typ) in enumerate(izip(sheet.row_values(i),
                                                  sheet.row_types(i))):
                if parse_cols is not None and j not in should_parse:
                    should_parse[j] = self._should_parse(j, parse_cols)

                if parse_cols is None or should_parse[j]:
                    if typ == XL_CELL_DATE:
                        dt = xldate_as_tuple(value, datemode)
                        # how to produce this first case?
                        if dt[0] < datetime.MINYEAR:  # pragma: no cover
                            value = datetime.time(*dt[3:])
                        else:
                            value = datetime.datetime(*dt)
                    elif typ == XL_CELL_ERROR:
                        value = np.nan
                    elif typ == XL_CELL_BOOLEAN:
                        value = bool(value)
                    row.append(value)

            data.append(row)

        if header is not None:
            data[header] = _trim_excel_header(data[header])

        parser = TextParser(data, header=header, index_col=index_col,
                            has_index_names=has_index_names,
                            na_values=na_values,
                            thousands=thousands,
                            parse_dates=parse_dates,
                            date_parser=date_parser,
                            skiprows=skiprows,
                            skip_footer=skip_footer,
                            chunksize=chunksize)

        return parser.read()

    @property
    def sheet_names(self):
            return self.book.sheet_names()


def _trim_excel_header(row):
    # trim header row so auto-index inference works
    # xlrd uses '' , openpyxl None
    while len(row) > 0 and (row[0] == '' or row[0] is None):
        row = row[1:]
    return row


class CellStyleConverter(object):
    """
    Utility Class which converts a style dict to xlrd or openpyxl style
    """

    @staticmethod
    def to_xls(style_dict, num_format_str=None):
        """
        converts a style_dict to an xlwt style object
        Parameters
        ----------
        style_dict: style dictionary to convert
        """
        import xlwt

        def style_to_xlwt(item, firstlevel=True, field_sep=',', line_sep=';'):
            """helper wich recursively generate an xlwt easy style string
            for example:

              hstyle = {"font": {"bold": True},
              "border": {"top": "thin",
                        "right": "thin",
                        "bottom": "thin",
                        "left": "thin"},
              "align": {"horiz": "center"}}
              will be converted to
              font: bold on; \
                      border: top thin, right thin, bottom thin, left thin; \
                      align: horiz center;
            """
            if hasattr(item, 'items'):
                if firstlevel:
                    it = ["%s: %s" % (key, style_to_xlwt(value, False))
                          for key, value in item.items()]
                    out = "%s " % (line_sep).join(it)
                    return out
                else:
                    it = ["%s %s" % (key, style_to_xlwt(value, False))
                          for key, value in item.items()]
                    out = "%s " % (field_sep).join(it)
                    return out
            else:
                item = "%s" % item
                item = item.replace("True", "on")
                item = item.replace("False", "off")
                return item

        if style_dict:
            xlwt_stylestr = style_to_xlwt(style_dict)
            style = xlwt.easyxf(xlwt_stylestr, field_sep=',', line_sep=';')
        else:
            style = xlwt.XFStyle()
        if num_format_str is not None:
            style.num_format_str = num_format_str

        return style

    @staticmethod
    def to_xlsx(style_dict):
        """
        converts a style_dict to an openpyxl style object
        Parameters
        ----------
        style_dict: style dictionary to convert
        """

        from openpyxl.style import Style
        xls_style = Style()
        for key, value in style_dict.items():
            for nk, nv in value.items():
                if key == "borders":
                    (xls_style.borders.__getattribute__(nk)
                     .__setattr__('border_style', nv))
                else:
                    xls_style.__getattribute__(key).__setattr__(nk, nv)

        return xls_style


def _conv_value(val):
    # convert value for excel dump
    if isinstance(val, np.int64):
        val = int(val)
    elif isinstance(val, np.bool8):
        val = bool(val)
    elif isinstance(val, Period):
        val = "%s" % val

    return val


class ExcelWriter(object):
    """
    Class for writing DataFrame objects into excel sheets, uses xlwt for xls,
    openpyxl for xlsx.  See DataFrame.to_excel for typical usage.

    Parameters
    ----------
    path : string
        Path to xls file
    """
    def __init__(self, path):
        self.use_xlsx = True
        if path.endswith('.xls'):
            self.use_xlsx = False
            import xlwt
            self.book = xlwt.Workbook()
            self.fm_datetime = xlwt.easyxf(
                num_format_str='YYYY-MM-DD HH:MM:SS')
            self.fm_date = xlwt.easyxf(num_format_str='YYYY-MM-DD')
        else:
            from openpyxl.workbook import Workbook
            self.book = Workbook()  # optimized_write=True)
            # open pyxl 1.6.1 adds a dummy sheet remove it
            if self.book.worksheets:
                self.book.remove_sheet(self.book.worksheets[0])
        self.path = path
        self.sheets = {}
        self.cur_sheet = None

    def save(self):
        """
        Save workbook to disk
        """
        self.book.save(self.path)

    def write_cells(self, cells, sheet_name=None, startrow=0, startcol=0):
        """
        Write given formated cells into Excel an excel sheet

        Parameters
        ----------
        cells : generator
            cell of formated data to save to Excel sheet
        sheet_name : string, default None
            Name of Excel sheet, if None, then use self.cur_sheet
        startrow: upper left cell row to dump data frame
        startcol: upper left cell column to dump data frame
        """
        if sheet_name is None:
            sheet_name = self.cur_sheet
        if sheet_name is None:  # pragma: no cover
            raise Exception('Must pass explicit sheet_name or set '
                            'cur_sheet property')
        if self.use_xlsx:
            self._writecells_xlsx(cells, sheet_name, startrow, startcol)
        else:
            self._writecells_xls(cells, sheet_name, startrow, startcol)

    def _writecells_xlsx(self, cells, sheet_name, startrow, startcol):

        from openpyxl.cell import get_column_letter

        if sheet_name in self.sheets:
            wks = self.sheets[sheet_name]
        else:
            wks = self.book.create_sheet()
            wks.title = sheet_name
            self.sheets[sheet_name] = wks

        for cell in cells:
            colletter = get_column_letter(startcol + cell.col + 1)
            xcell = wks.cell("%s%s" % (colletter, startrow + cell.row + 1))
            xcell.value = _conv_value(cell.val)
            if cell.style:
                style = CellStyleConverter.to_xlsx(cell.style)
                for field in style.__fields__:
                    xcell.style.__setattr__(field,
                                            style.__getattribute__(field))

            if isinstance(cell.val, datetime.datetime):
                xcell.style.number_format.format_code = "YYYY-MM-DD HH:MM:SS"
            elif isinstance(cell.val, datetime.date):
                xcell.style.number_format.format_code = "YYYY-MM-DD"

            # merging requires openpyxl latest (works on 1.6.1)
            # todo add version check
            if cell.mergestart is not None and cell.mergeend is not None:
                cletterstart = get_column_letter(startcol + cell.col + 1)
                cletterend = get_column_letter(startcol + cell.mergeend + 1)

                wks.merge_cells('%s%s:%s%s' % (cletterstart,
                                               startrow + cell.row + 1,
                                               cletterend,
                                               startrow + cell.mergestart + 1))

    def _writecells_xls(self, cells, sheet_name, startrow, startcol):
        if sheet_name in self.sheets:
            wks = self.sheets[sheet_name]
        else:
            wks = self.book.add_sheet(sheet_name)
            self.sheets[sheet_name] = wks

        style_dict = {}

        for cell in cells:
            val = _conv_value(cell.val)

            num_format_str = None
            if isinstance(cell.val, datetime.datetime):
                num_format_str = "YYYY-MM-DD HH:MM:SS"
            if isinstance(cell.val, datetime.date):
                num_format_str = "YYYY-MM-DD"

            stylekey = json.dumps(cell.style)
            if num_format_str:
                stylekey += num_format_str

            if stylekey in style_dict:
                style = style_dict[stylekey]
            else:
                style = CellStyleConverter.to_xls(cell.style, num_format_str)
                style_dict[stylekey] = style

            if cell.mergestart is not None and cell.mergeend is not None:
                wks.write_merge(startrow + cell.row,
                                startrow + cell.mergestart,
                                startcol + cell.col,
                                startcol + cell.mergeend,
                                val, style)
            else:
                wks.write(startrow + cell.row,
                          startcol + cell.col,
                          val, style)


"""
The StataReader below was originally written by Joe Presbrey as part of PyDTA.
It has been extended and improved by Skipper Seabold from the Statsmodels project
who also developed the StataWriter and was finally added to pandas in an once again
improved version.

You can find more information on http://presbrey.mit.edu/PyDTA and
http://statsmodels.sourceforge.net/devel/
"""


def is_py3():
    if sys.version_info[0] == 3:
        return True
    return False
PY3 = is_py3()


_date_formats = ["%tc", "%tC", "%td", "%tw", "%tm", "%tq", "%th", "%ty"]


def _stata_elapsed_date_to_datetime(date, fmt):
    """
    Convert from SIF to datetime. http://www.stata.com/help.cgi?datetime

    Parameters
    ----------
    date : int
        The Stata Internal Format date to convert to datetime according to fmt
    fmt : str
        The format to convert to. Can be, tc, td, tw, tm, tq, th, ty

    Examples
    --------
    >>> _stata_elapsed_date_to_datetime(52, "%tw")                                datetime.datetime(1961, 1, 1, 0, 0)

    Notes
    -----
    datetime/c - tc
        milliseconds since 01jan1960 00:00:00.000, assuming 86,400 s/day
    datetime/C - tC - NOT IMPLEMENTED
        milliseconds since 01jan1960 00:00:00.000, adjusted for leap seconds
    date - td
        days since 01jan1960 (01jan1960 = 0)
    weekly date - tw
        weeks since 1960w1
        This assumes 52 weeks in a year, then adds 7 * remainder of the weeks.
        The datetime value is the start of the week in terms of days in the
        year, not ISO calendar weeks.
    monthly date - tm
        months since 1960m1
    quarterly date - tq
        quarters since 1960q1
    half-yearly date - th
        half-years since 1960h1 yearly
    date - ty
        years since 0000

    If you don't have pandas with datetime support, then you can't do
    milliseconds accurately.
    """
    #NOTE: we could run into overflow / loss of precision situations here
    # casting to int, but I'm not sure what to do. datetime won't deal with
    # numpy types and numpy datetime isn't mature enough / we can't rely on
    # pandas version > 0.7.1
    #TODO: IIRC relative delta doesn't play well with np.datetime?
    if np.isnan(date):
        return np.datetime64('nat')

    date = int(date)
    stata_epoch = datetime.datetime(1960, 1, 1)
    if fmt in ["%tc", "tc"]:
        from dateutil.relativedelta import relativedelta
        return stata_epoch + relativedelta(microseconds=date * 1000)
    elif fmt in ["%tC", "tC"]:
        from warnings import warn
        warn("Encountered %tC format. Leaving in Stata Internal Format.")
        return date
    elif fmt in ["%td", "td"]:
        return stata_epoch + datetime.timedelta(int(date))
    elif fmt in ["%tw", "tw"]:  # does not count leap days - 7 days is a week
        year = datetime.datetime(stata_epoch.year + date // 52, 1, 1)
        day_delta = (date % 52) * 7
        return year + datetime.timedelta(int(day_delta))
    elif fmt in ["%tm", "tm"]:
        year = stata_epoch.year + date // 12
        month_delta = (date % 12) + 1
        return datetime.datetime(year, month_delta, 1)
    elif fmt in ["%tq", "tq"]:
        year = stata_epoch.year + date // 4
        month_delta = (date % 4) * 3 + 1
        return datetime.datetime(year, month_delta, 1)
    elif fmt in ["%th", "th"]:
        year = stata_epoch.year + date // 2
        month_delta = (date % 2) * 6 + 1
        return datetime.datetime(year, month_delta, 1)
    elif fmt in ["%ty", "ty"]:
        if date > 0:
            return datetime.datetime(date, 1, 1)
        else:  # don't do negative years bc can't mix dtypes in column
            raise ValueError("Year 0 and before not implemented")
    else:
        raise ValueError("Date fmt %s not understood" % fmt)


def _datetime_to_stata_elapsed(date, fmt):
    """
    Convert from datetime to SIF. http://www.stata.com/help.cgi?datetime

    Parameters
    ----------
    date : datetime.datetime
        The date to convert to the Stata Internal Format given by fmt
    fmt : str
        The format to convert to. Can be, tc, td, tw, tm, tq, th, ty
    """
    if not isinstance(date, datetime.datetime):
        raise ValueError("date should be datetime.datetime format")
    stata_epoch = datetime.datetime(1960, 1, 1)
    if fmt in ["%tc", "tc"]:
        delta = date - stata_epoch
        return (delta.days * 86400000 + delta.seconds*1000 +
                delta.microseconds/1000)
    elif fmt in ["%tC", "tC"]:
        from warnings import warn
        warn("Stata Internal Format tC not supported.")
        return date
    elif fmt in ["%td", "td"]:
        return (date - stata_epoch).days
    elif fmt in ["%tw", "tw"]:
        return (52*(date.year-stata_epoch.year) +
                (date - datetime.datetime(date.year, 1, 1)).days / 7)
    elif fmt in ["%tm", "tm"]:
        return (12 * (date.year - stata_epoch.year) + date.month - 1)
    elif fmt in ["%tq", "tq"]:
        return 4*(date.year-stata_epoch.year) + int((date.month - 1)/3)
    elif fmt in ["%th", "th"]:
        return 2 * (date.year - stata_epoch.year) + int(date.month > 6)
    elif fmt in ["%ty", "ty"]:
        return date.year
    else:
        raise ValueError("fmt %s not understood" % fmt)


class StataMissingValue(object):
    """
    An observation's missing value.

    Parameters
    -----------
    offset
    value

    Attributes
    ----------
    string
    value

    Notes
    -----
    More information: <http://www.stata.com/help.cgi?missing>
    """

    def __init__(self, offset, value):
        self._value = value
        if type(value) is int or type(value) is long:
            self._str = value - offset is 1 and \
                '.' or ('.' + chr(value - offset + 96))
        else:
            self._str = '.'
    string = property(lambda self: self._str, doc="The Stata representation of the missing value: '.', '.a'..'.z'")
    value = property(lambda self: self._value, doc='The binary representation of the missing value.')

    def __str__(self):
        return self._str

    __str__.__doc__ = string.__doc__


class StataParser(object):
    def __init__(self):
        #type          code
        #--------------------
        #str1        1 = 0x01
        #str2        2 = 0x02
        #...
        #str244    244 = 0xf4
        #byte      251 = 0xfb  (sic)
        #int       252 = 0xfc
        #long      253 = 0xfd
        #float     254 = 0xfe
        #double    255 = 0xff
        #--------------------
        #NOTE: the byte type seems to be reserved for categorical variables
        # with a label, but the underlying variable is -127 to 100
        # we're going to drop the label and cast to int
        self.DTYPE_MAP = \
            dict(
                zip(range(1, 245), ['a' + str(i) for i in range(1, 245)]) +
                [
                    (251, np.int16),
                    (252, np.int32),
                    (253, np.int64),
                    (254, np.float32),
                    (255, np.float64)
                ]
            )
        self.TYPE_MAP = range(251) + list('bhlfd')
        #NOTE: technically, some of these are wrong. there are more numbers
        # that can be represented. it's the 27 ABOVE and BELOW the max listed
        # numeric data type in [U] 12.2.2 of the 11.2 manual
        self.MISSING_VALUES = \
            {
                'b': (-127, 100),
                'h': (-32767, 32740),
                'l': (-2147483647, 2147483620),
                'f': (-1.701e+38, +1.701e+38),
                'd': (-1.798e+308, +8.988e+307)
            }

        self.OLD_TYPE_MAPPING = \
            {
                'i': 252,
                'f': 254,
                'b': 251
            }


class StataReader(StataParser):
    """
    Class for working with a Stata dataset. There are two possibilities for usage:

     * The from_dta() method on the DataFrame class.
       This will return a DataFrame with the Stata dataset. Note that when using the
       from_dta() method, you will not have access to meta-information like variable
       labels or the data label.

     * Work with this object directly. Upon instantiation, the header of the Stata data
       file is read, giving you access to attributes like variable_labels(), data_label(),
       nobs(), ... A DataFrame with the data is returned by the read() method; this will
       also fill up the value_labels. Note that calling the value_labels() method will
       result in an error if the read() method has not been called yet. This is because
       the value labels are stored at the end of a Stata dataset, after the data.

    Parameters
    ----------
    path_or_buf : string or file-like object
        Path to .dta file or object implementing a binary read() functions
    encoding : string, None or encoding
        Encoding used to parse the files. Note that Stata doesn't
        support unicode. None defaults to cp1252.
    """
    def __init__(self, path_or_buf, encoding=None):
        super(StataReader, self).__init__()
        self.col_sizes = ()
        self._has_string_data = False
        self._missing_values = False
        self._data_read = False
        self._value_labels_read = False
        if(encoding is None):
            self.encoding = 'cp1252'
        else:
            self.encoding = encoding
        if isinstance(path_or_buf, str) and _is_url(path_or_buf):
            from urllib.request import urlopen
            path_or_buf = urlopen(path_or_buf)
            if py3compat.PY3:  # pragma: no cover
                if self.encoding:
                    errors = 'strict'
                else:
                    errors = 'replace'
                    self.encoding = 'cp1252'
                bytes = path_or_buf.read()
                self.path_or_buf = StringIO(bytes.decode(self.encoding, errors))
        elif type(path_or_buf) is str:
            self.path_or_buf = open(path_or_buf, 'rb')
        else:
            self.path_or_buf = path_or_buf

        self._read_header()

    def _read_header(self):
        # header
        self.format_version = struct.unpack('b', self.path_or_buf.read(1))[0]
        if self.format_version not in [104, 105, 108, 113, 114, 115]:
            raise ValueError("Version of given Stata file is not 104, 105, 108, 113 (Stata 8/9), 114 (Stata 10/11) or 115 (Stata 12)")
        self.byteorder = self.path_or_buf.read(1) == 0x1 and '>' or '<'
        self.filetype = struct.unpack('b', self.path_or_buf.read(1))[0]
        self.path_or_buf.read(1)  # unused

        self.nvar = struct.unpack(self.byteorder + 'H', self.path_or_buf.read(2))[0]
        self.nobs = struct.unpack(self.byteorder + 'I', self.path_or_buf.read(4))[0]
        if self.format_version > 105:
            self.data_label = self.path_or_buf.read(81)
        else:
            self.data_label = self.path_or_buf.read(32)
        if self.format_version > 104:
            self.time_stamp = self.path_or_buf.read(18)

        # descriptors
        if self.format_version > 108:
            typlist = [ord(self.path_or_buf.read(1)) for i in range(self.nvar)]
        else:
            typlist = [self.OLD_TYPE_MAPPING[self.path_or_buf.read(1).decode(self.encoding)] for i in range(self.nvar)]
        self.typlist = [self.TYPE_MAP[typ] for typ in typlist]
        self.dtyplist = [self.DTYPE_MAP[typ] for typ in typlist]
        if self.format_version > 108:
            self.varlist = [self._null_terminate(self.path_or_buf.read(33)) for i in range(self.nvar)]
        else:
            self.varlist = [self._null_terminate(self.path_or_buf.read(9)) for i in range(self.nvar)]
        self.srtlist = struct.unpack(self.byteorder + ('h' * (self.nvar + 1)), self.path_or_buf.read(2 * (self.nvar + 1)))[:-1]
        if self.format_version > 113:
            self.fmtlist = [self._null_terminate(self.path_or_buf.read(49)) for i in range(self.nvar)]
        elif self.format_version > 104:
            self.fmtlist = [self._null_terminate(self.path_or_buf.read(12)) for i in range(self.nvar)]
        else:
            self.fmtlist = [self._null_terminate(self.path_or_buf.read(7)) for i in range(self.nvar)]
        if self.format_version > 108:
            self.lbllist = [self._null_terminate(self.path_or_buf.read(33)) for i in range(self.nvar)]
        else:
            self.lbllist = [self._null_terminate(self.path_or_buf.read(9)) for i in range(self.nvar)]
        if self.format_version > 105:
            self.vlblist = [self._null_terminate(self.path_or_buf.read(81)) for i in range(self.nvar)]
        else:
            self.vlblist = [self._null_terminate(self.path_or_buf.read(32)) for i in range(self.nvar)]

        # ignore expansion fields (Format 105 and later)
        # When reading, read five bytes; the last four bytes now tell you the
        # size of the next read, which you discard.  You then continue like
        # this until you read 5 bytes of zeros.

        if self.format_version > 104:
            while True:
                data_type = struct.unpack(self.byteorder + 'b', self.path_or_buf.read(1))[0]
                if self.format_version > 108:
                    data_len = struct.unpack(self.byteorder + 'i', self.path_or_buf.read(4))[0]
                else:
                    data_len = struct.unpack(self.byteorder + 'h', self.path_or_buf.read(2))[0]
                if data_type == 0:
                    break
                self.path_or_buf.read(data_len)

        # necessary data to continue parsing
        self.data_location = self.path_or_buf.tell()
        self.has_string_data = len([x for x in self.typlist if type(x) is int]) > 0
        self._col_size()

    def _calcsize(self, fmt):
        return type(fmt) is int and fmt or struct.calcsize(self.byteorder + fmt)

    def _col_size(self, k=None):
        """Calculate size of a data record."""
        if len(self.col_sizes) == 0:
            self.col_sizes = map(lambda x: self._calcsize(x), self.typlist)
        if k is None:
            return self.col_sizes
        else:
            return self.col_sizes[k]

    def _unpack(self, fmt, byt):
        d = struct.unpack(self.byteorder + fmt, byt)[0]
        if fmt[-1] in self.MISSING_VALUES:
            nmin, nmax = self.MISSING_VALUES[fmt[-1]]
            if d < nmin or d > nmax:
                if self._missing_values:
                    return StataMissingValue(nmax, d)
                else:
                    return None
        return d

    def _null_terminate(self, s):
        if PY3:  # have bytes not strings, so must decode
            null_byte = 0x00
            try:
                s = s[:s.index(null_byte)]
            except:
                pass
            return s.decode(self.encoding)
        else:
            null_byte = 0x00
            try:
                return s.lstrip(null_byte)[:s.index(null_byte)]
            except:
                return s

    def _next(self):
        typlist = self.typlist
        if self._has_string_data:
            data = [None] * self.nvar
            for i in range(len(data)):
                if type(typlist[i]) is int:
                    data[i] = self._null_terminate(self.path_or_buf.read(typlist[i]))
                else:
                    data[i] = self._unpack(typlist[i], self.path_or_buf.read(self._col_size(i)))
            return data
        else:
            return map(lambda i: self._unpack(typlist[i],
                                              self.path_or_buf.read(self._col_size(i))),
                       range(self.nvar))

    def _dataset(self):
        """
        Returns a Python generator object for iterating over the dataset.


        Parameters
        ----------

        Returns
        -------
        Generator object for iterating over the dataset.  Yields each row of
        observations as a list by default.

        Notes
        -----
        If missing_values is True during instantiation of StataReader then
        observations with _StataMissingValue(s) are not filtered and should
        be handled by your applcation.
        """

        try:
            self._file.seek(self._data_location)
        except Exception:
            pass

        for i in range(self.nobs):
            yield self._next()

    def _read_value_labels(self):
        if not self._data_read:
            raise Exception("Data has not been read. Because of the layout of Stata files, this is necessary before reading value labels.")
        if self._value_labels_read:
            raise Exception("Value labels have already been read.")

        self.value_label_dict = dict()

        if self.format_version <= 108:
            return  # Value labels are not supported in version 108 and earlier.

        while True:
            slength = self.path_or_buf.read(4)
            if not slength:
                break  # end of variable lable table
            labname = self._null_terminate(self.path_or_buf.read(33))
            self.path_or_buf.read(3)  # padding

            n = struct.unpack(self.byteorder + 'I', self.path_or_buf.read(4))[0]
            txtlen = struct.unpack(self.byteorder + 'I', self.path_or_buf.read(4))[0]
            off = []
            for i in range(n):
                off.append(struct.unpack(self.byteorder + 'I', self.path_or_buf.read(4))[0])
            val = []
            for i in range(n):
                val.append(struct.unpack(self.byteorder + 'I', self.path_or_buf.read(4))[0])
            txt = self.path_or_buf.read(txtlen)
            self.value_label_dict[labname] = dict()
            for i in range(n):
                self.value_label_dict[labname][val[i]] = self._null_terminate(txt[off[i]:])
        self._value_labels_read = True

    def data(self, convert_dates=True, convert_categoricals=True, index=None):
        """
        Reads observations from Stata file, converting them into a dataframe

        Parameters
        ----------
        convert_dates : boolean, defaults to True
            Convert date variables to DataFrame time values
        convert_categoricals : boolean, defaults to True
            Read value labels and convert columns to Categorical/Factor variables
        index : identifier of index column
            identifier of column that should be used as index of the DataFrame

        Returns
        -------
        y : DataFrame instance
        """
        if self._data_read:
            raise Exception("Data has already been read.")
        self._data_read = True

        stata_dta = self._dataset()

        data = []
        for rownum, line in enumerate(stata_dta):
            # doesn't handle missing value objects, just casts
            # None will only work without missing value object.
            for i, val in enumerate(line):
                #NOTE: This will only be scalar types because missing strings
                # are empty not None in Stata
                if val is None:
                    line[i] = np.nan
            data.append(tuple(line))

        if convert_categoricals:
            self._read_value_labels()

        data = DataFrame(data, columns=self.varlist, index=index)

        cols_ = np.where(self.dtyplist)[0]
        for i in cols_:
            if self.dtyplist[i] is not None:
                col = data.columns[i]
                data[col] = Series(data[col], data[col].index, self.dtyplist[i])

        if convert_dates:
            cols = np.where(map(lambda x: x in _date_formats, self.fmtlist))[0]
            for i in cols:
                col = data.columns[i]
                data[col] = data[col].apply(_stata_elapsed_date_to_datetime, args=(self.fmtlist[i],))

        if convert_categoricals:
            cols = np.where(map(lambda x: x in self.value_label_dict.iterkeys(), self.lbllist))[0]
            for i in cols:
                col = data.columns[i]
                labeled_data = np.copy(data[col])
                labeled_data = labeled_data.astype(object)
                for k, v in self.value_label_dict[self.lbllist[i]].iteritems():
                    labeled_data[data[col] == k] = v
                data[col] = Categorical.from_array(labeled_data)

        return data

    def data_label(self):
        """Returns data label of Stata file"""
        return self.data_label

    def variable_labels(self):
        """Returns variable labels as a dict, associating each variable name with corresponding label"""
        return dict(zip(self.varlist, self.vlblist))

    def value_labels(self):
        """Returns a dict, associating each variable name a dict, associating each value its corresponding label"""
        if not self._value_labels_read:
            self._read_value_labels()

        return self.value_label_dict


def _open_file_binary_write(fname, encoding):
    if hasattr(fname, 'write'):
        #if 'b' not in fname.mode:
        return fname
    return open(fname, "wb")


def _set_endianness(endianness):
    if endianness.lower() in ["<", "little"]:
        return "<"
    elif endianness.lower() in [">", "big"]:
        return ">"
    else:  # pragma : no cover
        raise ValueError("Endianness %s not understood" % endianness)


def _pad_bytes(name, length):
    """
    Takes a char string and pads it wih null bytes until it's length chars
    """
    return name + "\x00" * (length - len(name))


def _default_names(nvar):
    """
    Returns default Stata names v1, v2, ... vnvar
    """
    return ["v%d" % i for i in range(1, nvar+1)]


def _convert_datetime_to_stata_type(fmt):
    """
    Converts from one of the stata date formats to a type in TYPE_MAP
    """
    if fmt in ["tc", "%tc", "td", "%td", "tw", "%tw", "tm", "%tm", "tq",
               "%tq", "th", "%th", "ty", "%ty"]:
        return np.float64  # Stata expects doubles for SIFs
    else:
        raise ValueError("fmt %s not understood" % fmt)


def _maybe_convert_to_int_keys(convert_dates, varlist):
    new_dict = {}
    for key in convert_dates:
        if not convert_dates[key].startswith("%"):  # make sure proper fmts
            convert_dates[key] = "%" + convert_dates[key]
        if key in varlist:
            new_dict.update({varlist.index(key): convert_dates[key]})
        else:
            if not isinstance(key, int):
                raise ValueError("convery_dates key is not in varlist and is not an int")
            new_dict.update({key: convert_dates[key]})
    return new_dict


def _dtype_to_stata_type(dtype):
    """
    Converts dtype types to stata types. Returns the byte of the given ordinal.
    See TYPE_MAP and comments for an explanation. This is also explained in
    the dta spec.
    1 - 244 are strings of this length
    251 - chr(251) - for int8 and int16, byte
    252 - chr(252) - for int32, int
    253 - chr(253) - for int64, long
    254 - chr(254) - for float32, float
    255 - chr(255) - double, double

    If there are dates to convert, then dtype will already have the correct
    type inserted.
    """
    #TODO: expand to handle datetime to integer conversion
    if dtype.type == np.string_:
        return chr(dtype.itemsize)
    elif dtype.type == np.object_:  # try to coerce it to the biggest string
                                    # not memory efficient, what else could we do?
        return chr(244)
    elif dtype == np.float64:
        return chr(255)
    elif dtype == np.float32:
        return chr(254)
    elif dtype == np.int64:
        return chr(253)
    elif dtype == np.int32:
        return chr(252)
    elif dtype == np.int8 or dtype == np.int16:
        return chr(251)
    else:  # pragma : no cover
        raise ValueError("Data type %s not currently understood. "
                         "Please report an error to the developers." % dtype)


def _dtype_to_default_stata_fmt(dtype):
    """
    Maps numpy dtype to stata's default format for this type. Not terribly
    important since users can change this in Stata. Semantics are

    string  -> "%DDs" where DD is the length of the string
    float64 -> "%10.0g"
    float32 -> "%9.0g"
    int64   -> "%9.0g"
    int32   -> "%12.0g"
    int16   -> "%8.0g"
    int8    -> "%8.0g"
    """
    #TODO: expand this to handle a default datetime format?
    if dtype.type == np.string_:
        return "%" + str(dtype.itemsize) + "s"
    elif dtype.type == np.object_:
        return "%244s"
    elif dtype == np.float64:
        return "%10.0g"
    elif dtype == np.float32:
        return "%9.0g"
    elif dtype == np.int64:
        return "%9.0g"
    elif dtype == np.int32:
        return "%12.0g"
    elif dtype == np.int8 or dtype == np.int16:
        return "%8.0g"
    else:  # pragma : no cover
        raise ValueError("Data type %s not currently understood. "
                         "Please report an error to the developers." % dtype)


class StataWriter(StataParser):
    """
    A class for writing Stata binary dta files from array-like objects

    Parameters
    ----------
    fname : file path or buffer
        Where to save the dta file.
    data : array-like
        Array-like input to save. Pandas objects are also accepted.
    convert_dates : dict
        Dictionary mapping column of datetime types to the stata internal
        format that you want to use for the dates. Options are
        'tc', 'td', 'tm', 'tw', 'th', 'tq', 'ty'. Column can be either a
        number or a name.
    encoding : str
        Default is latin-1. Note that Stata does not support unicode.
    byteorder : str
        Can be ">", "<", "little", or "big". The default is None which uses
        `sys.byteorder`

    Returns
    -------
    writer : StataWriter instance
        The StataWriter instance has a write_file method, which will
        write the file to the given `fname`.

    Examples
    --------
    >>> writer = StataWriter('./data_file.dta', data)
    >>> writer.write_file()

    Or with dates

    >>> writer = StataWriter('./date_data_file.dta', date, {2 : 'tw'})
    >>> writer.write_file()
    """
    def __init__(self, fname, data, convert_dates=None, write_index=True, encoding="latin-1",
                 byteorder=None):
        super(StataWriter, self).__init__()
        self._convert_dates = convert_dates
        self._write_index = write_index
        # attach nobs, nvars, data, varlist, typlist
        self._prepare_pandas(data)

        if byteorder is None:
            byteorder = sys.byteorder
        self._byteorder = _set_endianness(byteorder)
        self._encoding = encoding
        self._file = _open_file_binary_write(fname, encoding)
        self.type_converters = {253: np.long, 252: int}

    def _write(self, to_write):
        """
        Helper to call asbytes before writing to file for Python 3 compat.
        """
        self._file.write(to_write.encode(self._encoding))

    def _prepare_pandas(self, data):
        #NOTE: we might need a different API / class for pandas objects so
        # we can set different semantics - handle this with a PR to pandas.io
        class DataFrameRowIter(object):
            def __init__(self, data):
                self.data = data

            def __iter__(self):
                for i, row in data.iterrows():
                    yield row

        if self._write_index:
            data = data.reset_index()
        self.datarows = DataFrameRowIter(data)
        self.nobs, self.nvar = data.shape
        self.data = data
        self.varlist = data.columns.tolist()
        dtypes = data.dtypes
        if self._convert_dates is not None:
            self._convert_dates = _maybe_convert_to_int_keys(self._convert_dates, self.varlist)
            for key in self._convert_dates:
                new_type = _convert_datetime_to_stata_type(self._convert_dates[key])
                dtypes[key] = np.dtype(new_type)
        self.typlist = [_dtype_to_stata_type(dt) for dt in dtypes]
        self.fmtlist = [_dtype_to_default_stata_fmt(dt) for dt in dtypes]
        # set the given format for the datetime cols
        if self._convert_dates is not None:
            for key in self._convert_dates:
                self.fmtlist[key] = self._convert_dates[key]

    def write_file(self):
        self._write_header()
        self._write_descriptors()
        self._write_variable_labels()
        # write 5 zeros for expansion fields
        self._write(_pad_bytes("", 5))
        if self._convert_dates is None:
            self._write_data_nodates()
        else:
            self._write_data_dates()
        #self._write_value_labels()
        self._file.close()

    def _write_header(self, data_label=None, time_stamp=None):
        byteorder = self._byteorder
        # ds_format - just use 114
        self._file.write(struct.pack("b", 114))
        # byteorder
        self._write(byteorder == ">" and "\x01" or "\x02")
        # filetype
        self._write("\x01")
        # unused
        self._write("\x00")
        # number of vars, 2 bytes
        self._file.write(struct.pack(byteorder+"h", self.nvar)[:2])
        # number of obs, 4 bytes
        self._file.write(struct.pack(byteorder+"i", self.nobs)[:4])
        # data label 81 bytes, char, null terminated
        if data_label is None:
            self._file.write(self._null_terminate(_pad_bytes("", 80)))
        else:
            self._file.write(self._null_terminate(_pad_bytes(data_label[:80], 80)))
        # time stamp, 18 bytes, char, null terminated
        # format dd Mon yyyy hh:mm
        if time_stamp is None:
            time_stamp = datetime.datetime.now()
        elif not isinstance(time_stamp, datetime):
            raise ValueError("time_stamp should be datetime type")
        self._file.write(self._null_terminate(time_stamp.strftime("%d %b %Y %H:%M")))

    def _write_descriptors(self, typlist=None, varlist=None, srtlist=None,
                           fmtlist=None, lbllist=None):
        nvar = self.nvar
        # typlist, length nvar, format byte array
        for typ in self.typlist:
            self._write(typ)

        # varlist, length 33*nvar, char array, null terminated
        for name in self.varlist:
            name = self._null_terminate(name)
            name = _pad_bytes(name[:32].decode(self._encoding), 33)
            self._write(name)

        # srtlist, 2*(nvar+1), int array, encoded by byteorder
        srtlist = _pad_bytes("", (2*(nvar+1)))
        self._write(srtlist)

        # fmtlist, 49*nvar, char array
        for fmt in self.fmtlist:
            self._write(_pad_bytes(fmt, 49))

        # lbllist, 33*nvar, char array
        #NOTE: this is where you could get fancy with pandas categorical type
        for i in range(nvar):
            self._write(_pad_bytes("", 33))

    def _write_variable_labels(self, labels=None):
        nvar = self.nvar
        if labels is None:
            for i in range(nvar):
                self._write(_pad_bytes("", 81))

    def _write_data_nodates(self):
        data = self.datarows
        byteorder = self._byteorder
        TYPE_MAP = self.TYPE_MAP
        typlist = self.typlist
        for row in data:
            #row = row.squeeze().tolist() # needed for structured arrays
            for i, var in enumerate(row):
                typ = ord(typlist[i])
                if typ <= 244:  # we've got a string
                    if len(var) < typ:
                        var = _pad_bytes(var.decode(self._encoding), len(var) + 1)
                    self._write(var)
                else:
                    try:
                        self._file.write(struct.pack(byteorder + TYPE_MAP[typ], var))
                    except struct.error:
                        # have to be strict about type pack won't do any
                        # kind of casting
                        self._file.write(struct.pack(byteorder+TYPE_MAP[typ],
                                         self.type_converters[typ](var)))

    def _write_data_dates(self):
        convert_dates = self._convert_dates
        data = self.datarows
        byteorder = self._byteorder
        TYPE_MAP = self.TYPE_MAP
        MISSING_VALUES = self.MISSING_VALUES
        typlist = self.typlist
        for row in data:
            #row = row.squeeze().tolist() # needed for structured arrays
            for i, var in enumerate(row):
                typ = ord(typlist[i])
                #NOTE: If anyone finds this terribly slow, there is
                # a vectorized way to convert dates, see genfromdta for going
                # from int to datetime and reverse it. will copy data though
                if i in convert_dates:
                    var = _datetime_to_stata_elapsed(var, self.fmtlist[i])
                if typ <= 244:  # we've got a string
                    if isnull(var):
                        var = ""  # missing string
                    if len(var) < typ:
                        var = _pad_bytes(var, len(var) + 1)
                    self._write(var)
                else:
                    if isnull(var):  # this only matters for floats
                        var = MISSING_VALUES[typ]
                    self._write(struct.pack(byteorder+TYPE_MAP[typ], var))

    def _null_terminate(self, s):
        null_byte = '\x00'
        if PY3:
            s += null_byte
            return s.encode(self._encoding)
        else:
            s += null_byte
            return s
