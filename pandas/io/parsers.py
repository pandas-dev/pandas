"""
Module contains tools for processing files into DataFrames or other objects
"""
from StringIO import StringIO
import re
from itertools import izip
import csv
from warnings import warn

import numpy as np

from pandas.core.index import Index, MultiIndex
from pandas.core.frame import DataFrame
import datetime
import pandas.core.common as com
from pandas.util import py3compat
from pandas.io.date_converters import generic_parser
from pandas.io.common import get_filepath_or_buffer

from pandas.util.decorators import Appender

import pandas.lib as lib
import pandas.tslib as tslib
import pandas.parser as _parser
from pandas.tseries.period import Period

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
    The character to used to denote the start and end of a quoted item. Quoted items can include the delimiter and it will be ignored.
quoting : int
    Controls whether quotes should be recognized. Values are taken from
    `csv.QUOTE_*` values. Acceptable values are 0, 1, 2, and 3 for
    QUOTE_MINIMAL, QUOTE_ALL, QUOTE_NONE, and QUOTE_NONNUMERIC, respectively.
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
header : int, default 0 if names parameter not specified,
    Row to use for the column labels of the parsed DataFrame. Specify None if
    there is no header row. Can be a list of integers that specify row
    locations for a multi-index on the columns E.g. [0,1,3]. Interveaning
    rows that are not specified (E.g. 2 in this example are skipped)
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
    Return TextFileReader object
chunksize : int, default None
    Return TextFileReader object for iteration
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
usecols : array-like
    Return a subset of the columns.
    Results in much faster parsing time and lower memory usage.
mangle_dupe_cols: boolean, default True
    Duplicate columns will be specified as 'X.0'...'X.N', rather than 'X'...'X'
tupleize_cols: boolean, default False
    Leave a list of tuples on columns as is (default is to convert to
    a Multi Index on the columns)

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


def _read(filepath_or_buffer, kwds):
    "Generic reader of line files."
    encoding = kwds.get('encoding', None)
    skipfooter = kwds.pop('skipfooter', None)
    if skipfooter is not None:
        kwds['skip_footer'] = skipfooter

    filepath_or_buffer, _ = get_filepath_or_buffer(filepath_or_buffer)

    if kwds.get('date_parser', None) is not None:
        if isinstance(kwds['parse_dates'], bool):
            kwds['parse_dates'] = True

    # Extract some of the arguments (pass chunksize on).
    iterator = kwds.get('iterator', False)
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
    'tupleize_cols':True,
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
                 na_fvalues=None,
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
                 mangle_dupe_cols=True,
                 tupleize_cols=True,
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
                    na_fvalues=na_fvalues,
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
                    mangle_dupe_cols=mangle_dupe_cols,
                    tupleize_cols=tupleize_cols,
            )

        return _read(filepath_or_buffer, kwds)

    parser_f.__name__ = name

    return parser_f

read_csv = _make_parser_function('read_csv', sep=',')
read_csv = Appender(_read_csv_doc)(read_csv)

read_table = _make_parser_function('read_table', sep='\t')
read_table = Appender(_read_table_doc)(read_table)


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
                print ('Using Python parser to sniff delimiter')
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
        na_values, na_fvalues = _clean_na_values(na_values, keep_default_na)

        if com.is_integer(skiprows):
            skiprows = range(skiprows)
        skiprows = set() if skiprows is None else set(skiprows)

        # put stuff back
        result['index_col'] = index_col
        result['names'] = names
        result['converters'] = converters
        result['na_values'] = na_values
        result['na_fvalues'] = na_fvalues
        result['skiprows'] = skiprows

        return result, engine

    def __iter__(self):
        try:
            if self.chunksize:
                while True:
                    yield self.read(self.chunksize)
            else:
                yield self.read()
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

        # May alter columns / col_dict
        index, columns, col_dict = self._create_index(ret)

        df = DataFrame(col_dict, columns=columns, index=index)

        if self.squeeze and len(df.columns) == 1:
            return df[df.columns[0]]
        return df

    def _create_index(self, ret):
        index, columns, col_dict = ret
        return index, columns, col_dict

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
        self.col_names = None

        self.parse_dates = kwds.pop('parse_dates', False)
        self.date_parser = kwds.pop('date_parser', None)
        self.dayfirst = kwds.pop('dayfirst', False)
        self.keep_date_col = kwds.pop('keep_date_col', False)

        self.na_values = kwds.get('na_values')
        self.na_fvalues = kwds.get('na_fvalues')
        self.true_values = kwds.get('true_values')
        self.false_values = kwds.get('false_values')
        self.tupleize_cols = kwds.get('tupleize_cols',True)

        self._date_conv = _make_date_converter(date_parser=self.date_parser,
                                               dayfirst=self.dayfirst)

        # validate header options for mi
        self.header = kwds.get('header')
        if isinstance(self.header,(list,tuple,np.ndarray)):
            if kwds.get('as_recarray'):
                raise Exception("cannot specify as_recarray when "
                                "specifying a multi-index header")
            if kwds.get('usecols'):
                raise Exception("cannot specify usecols when "
                                "specifying a multi-index header")
            if kwds.get('names'):
                raise Exception("cannot specify names when "
                                "specifying a multi-index header")

            # validate index_col that only contains integers
            if self.index_col is not None:
                if not (isinstance(self.index_col,(list,tuple,np.ndarray)) and all(
                        [ com.is_integer(i) for i in self.index_col ]) or com.is_integer(self.index_col)):
                    raise Exception("index_col must only contain row numbers "
                                    "when specifying a multi-index header")

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


    def _extract_multi_indexer_columns(self, header, index_names, col_names, passed_names=False):
        """ extract and return the names, index_names, col_names
            header is a list-of-lists returned from the parsers """
        if len(header) < 2:
            return header[0], index_names, col_names, passed_names

        # the names are the tuples of the header that are not the index cols
        # 0 is the name of the index, assuming index_col is a list of column
        # numbers
        ic = self.index_col
        if ic is None:
            ic = []

        if not isinstance(ic, (list,tuple,np.ndarray)):
            ic = [ ic ]
        sic = set(ic)

        orig_header = list(header)

        # clean the index_names
        index_names = header.pop(-1)
        (index_names, names,
         index_col) = _clean_index_names(index_names, self.index_col)

        # extract the columns
        field_count = len(header[0])
        def extract(r):
            return tuple([ r[i] for i in range(field_count) if i not in sic ])
        columns = zip(*[ extract(r) for r in header ])
        names = ic + columns

        # if we find 'Unnamed' all of a single level, then our header was too long
        for n in range(len(columns[0])):
            if all([ 'Unnamed' in c[n] for c in columns ]):
                raise _parser.CParserError("Passed header=[%s] are too many rows for this "
                                           "multi_index of columns" % ','.join([ str(x) for x in self.header ]))

        # clean the column names (if we have an index_col)
        if len(ic):
            col_names = [ r[0] if len(r[0]) and 'Unnamed' not in r[0] else None for r in header ]
        else:
            col_names = [ None ] * len(header)

        passed_names = True

        return names, index_names, col_names, passed_names

    def _maybe_make_multi_index_columns(self, columns, col_names=None):
        # possibly create a column mi here
        if not self.tupleize_cols and len(columns) and not isinstance(
            columns, MultiIndex) and all([ isinstance(c,tuple) for c in columns]):
            columns = MultiIndex.from_tuples(columns,names=col_names)
        return columns

    def _make_index(self, data, alldata, columns, indexnamerow=False):
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

        # add names for the index
        if indexnamerow:
            coffset = len(indexnamerow) - len(columns)
            index.names = indexnamerow[:coffset]

        # maybe create a mi on the columns
        columns = self._maybe_make_multi_index_columns(columns, self.col_names)

        return index, columns

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
            col_na_fvalues = self.na_fvalues

            if isinstance(self.na_values, dict):
                col_name = self.index_names[i]
                if col_name is not None:
                    col_na_values, col_na_fvalues = _get_na_values(col_name,
                                                                   self.na_values,
                                                                   self.na_fvalues)
                    
            arr, _ = self._convert_types(arr, col_na_values | col_na_fvalues)
            arrays.append(arr)

        index = MultiIndex.from_arrays(arrays, names=self.index_names)

        return index

    def _convert_to_ndarrays(self, dct, na_values, na_fvalues, verbose=False,
                             converters=None):
        result = {}
        for c, values in dct.iteritems():
            conv_f = None if converters is None else converters.get(c, None)
            col_na_values, col_na_fvalues = _get_na_values(c, na_values, na_fvalues)
            coerce_type = True
            if conv_f is not None:
                values = lib.map_infer(values, conv_f)
                coerce_type = False
            cvals, na_count = self._convert_types(values,
                                                  set(col_na_values) | col_na_fvalues,
                                                  coerce_type)
            result[c] = cvals
            if verbose and na_count:
                print ('Filled %d NA values in column %s' % (na_count, str(c)))
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
            if len(self._reader.header) > 1:
                # we have a multi index in the columns
                self.names, self.index_names, self.col_names, passed_names = self._extract_multi_indexer_columns(
                    self._reader.header, self.index_names, self.col_names, passed_names)
            else:
                self.names = list(self._reader.header[0])

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
                (index_names, self.names,
                 self.index_col) = _clean_index_names(self.names, self.index_col)

                if self.index_names is None:
                    self.index_names = index_names

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
            index, names = self._make_index(data, alldata, names)

        # maybe create a mi on the columns
        names = self._maybe_make_multi_index_columns(names, self.col_names)

        return index, names, data

    def _filter_usecols(self, names):
        # hackish
        if self.usecols is not None and len(names) != len(self.usecols):
            names = [name for i, name in enumerate(names)
                     if i in self.usecols or name in self.usecols]
        return names

    def _get_index_names(self):
        names = list(self._reader.header[0])
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

        # we are processing a multi index column
        if len(self.columns) > 1:
            self.columns, self.index_names, self.col_names, _ = self._extract_multi_indexer_columns(
                self.columns, self.index_names, self.col_names)
        else:
            self.columns = self.columns[0]

        # get popped off for index
        self.orig_names = list(self.columns)

        # needs to be cleaned/refactored
        # multiple date column thing turning into a real spaghetti factory

        if not self._has_complex_date_col:
            (index_names,
             self.orig_names, _) = self._get_index_name(self.columns)
            self._name_processed = True
            if self.index_names is None:
                self.index_names = index_names
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
        index, columns = self._make_index(data, alldata, columns, indexnamerow)

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

        return self._convert_to_ndarrays(data, self.na_values, self.na_fvalues, self.verbose,
                                         clean_conv)

    def _infer_columns(self):
        names = self.names

        if self.header is not None:
            header = self.header

            # we have a mi columns, so read and extra line
            if isinstance(header,(list,tuple,np.ndarray)):
                have_mi_columns = True
                header = list(header) + [header[-1]+1]
            else:
                have_mi_columns = False
                header = [ header ]

            columns = []
            for level, hr in enumerate(header):

                if len(self.buf) > 0:
                    line = self.buf[0]
                else:
                    line = self._next_line()

                while self.pos <= hr:
                    line = self._next_line()

                this_columns = []
                for i, c in enumerate(line):
                    if c == '':
                        if have_mi_columns:
                            this_columns.append('Unnamed: %d_level_%d' % (i,level))
                        else:
                            this_columns.append('Unnamed: %d' % i)
                    else:
                        this_columns.append(c)

                if not have_mi_columns:
                    if self.mangle_dupe_cols:
                        counts = {}
                        for i, col in enumerate(this_columns):
                            cur_count = counts.get(col, 0)
                            if cur_count > 0:
                                this_columns[i] = '%s.%d' % (col, cur_count)
                            counts[col] = cur_count + 1

                columns.append(this_columns)

            self._clear_buffer()

            if names is not None:
                if len(names) != len(columns[0]):
                    raise Exception('Number of passed names did not match '
                                    'number of header fields in the file')
                if len(columns) > 1:
                    raise Exception('Cannot pass names with multi-index columns')
                columns = [ names ]

        else:
            if len(self.buf) > 0:
                line = self.buf[0]
            else:
                line = self._next_line()

            ncols = len(line)
            if not names:
                if self.prefix:
                    columns = [ ['X%d' % i for i in range(ncols)] ]
                else:
                    columns = [ range(ncols) ]
            else:
                columns = [ names ]

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
        new_rows = None

        # already fetched some number
        if rows is not None:

            # we already have the lines in the buffer
            if len(self.buf) >= rows:
                new_rows, self.buf = self.buf[:rows], self.buf[rows:]

            # need some lines
            else:
                rows -= len(self.buf)

        if new_rows is None:
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
        else:
            lines = new_rows

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

    new_col = parser(*to_parse)
    return new_name, new_col, colnames


def _clean_na_values(na_values, keep_default_na=True):

    if na_values is None and keep_default_na:
        na_values = _NA_VALUES
        na_fvalues = set()
    elif isinstance(na_values, dict):
        if keep_default_na:
            for k, v in na_values.iteritems():
                v = set(list(v)) | _NA_VALUES
                na_values[k] = v
        na_fvalues = dict([ (k, _floatify_na_values(v)) for k, v in na_values.items() ])
    else:
        if not com.is_list_like(na_values):
            na_values = [na_values]
        na_values = _stringify_na_values(na_values)
        if keep_default_na:
            na_values = na_values | _NA_VALUES

        na_fvalues = _floatify_na_values(na_values)

    return na_values, na_fvalues

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


def _floatify_na_values(na_values):
    # create float versions of the na_values
    result = set()
    for v in na_values:
        try:
            v = float(v)
            if not np.isnan(v):
                result.add(v)
        except:
            pass
    return result

def _stringify_na_values(na_values):
    """ return a stringified and numeric for these values """
    result = []
    for x in na_values:
        result.append(str(x))
        result.append(x)
        try:
            v = float(x)

            # we are like 999 here
            if v == int(v):
                v = int(v)
                result.append("%s.0" % v)
                result.append(str(v))

            result.append(v)
        except:
            pass
        try:
            result.append(int(x))
        except:
            pass
    return set(result)

def _get_na_values(col, na_values, na_fvalues):
    if isinstance(na_values, dict):
        if col in na_values:
            values = na_values[col]
            fvalues = na_fvalues[col]
            return na_values[col], na_fvalues[col]
        else:
            return _NA_VALUES, set()
    else:
        return na_values, na_fvalues


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


##### deprecations in 0.12 #####
##### remove in 0.12         #####

from pandas.io import clipboard
def read_clipboard(**kwargs):
    warn("read_clipboard is now a top-level accessible via pandas.read_clipboard", FutureWarning)
    clipboard.read_clipboard(**kwargs)

def to_clipboard(obj):
    warn("to_clipboard is now an object level method accessible via obj.to_clipboard()", FutureWarning)
    clipboard.to_clipboard(obj)

from pandas.io import excel
class ExcelWriter(excel.ExcelWriter):
    def __init__(self, path):
        warn("ExcelWriter can now be imported from: pandas.io.excel", FutureWarning)
        super(ExcelWriter, self).__init__(path)

class ExcelFile(excel.ExcelFile):
    def __init__(self, path_or_buf, kind=None, **kwds):
        warn("ExcelFile can now be imported from: pandas.io.excel", FutureWarning)
        super(ExcelFile, self).__init__(path_or_buf, kind=kind, **kwds)
