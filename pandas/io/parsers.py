"""
Module contains tools for processing files into DataFrames or other objects
"""

import numpy as np

from pandas.core.index import Index, MultiIndex
from pandas.core.frame import DataFrame
import pandas._tseries as lib


def read_csv(filepath_or_buffer, sep=None, header=0, index_col=None, names=None,
             skiprows=None, na_values=None, parse_dates=False,
             date_parser=None):
    import csv

    if hasattr(filepath_or_buffer, 'read'):
        f = filepath_or_buffer
    else:
        try:
            # universal newline mode
            f = open(filepath_or_buffer, 'U')
        except Exception: # pragma: no cover
            f = open(filepath_or_buffer, 'r')

    sniff_sep = True
    # default dialect
    dia = csv.excel
    if sep is not None:
        sniff_sep = False
        dia.delimiter = sep
    # attempt to sniff the delimiter
    if sniff_sep:
        sample = f.readline()
        sniffed = csv.Sniffer().sniff(sample)
        dia.delimiter = sniffed.delimiter
        f.seek(0)

    reader = csv.reader(f, dialect=dia)

    if skiprows is not None:
        skiprows = set(skiprows)
        lines = [l for i, l in enumerate(reader) if i not in skiprows]
    else:
        lines = [l for l in reader]
    f.close()

    if date_parser is not None:
        parse_dates = True

    return _simple_parser(lines,
                          header=header,
                          index_col=index_col,
                          names=names,
                          na_values=na_values,
                          parse_dates=parse_dates,
                          date_parser=date_parser)


def read_table(filepath_or_buffer, sep='\t', header=0, index_col=None,
               names=None, skiprows=None, na_values=None, parse_dates=False,
               date_parser=None):
    return read_csv(filepath_or_buffer, sep=sep, header=header,
                    skiprows=skiprows, index_col=index_col,
                    na_values=na_values, date_parser=date_parser,
                    names=names, parse_dates=parse_dates)

_parser_params = """Parameters
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
na_values : list-like, default None
    List of additional strings to recognize as NA/NaN
parse_dates : boolean, default False
    Attempt to parse dates in the index column(s)
date_parser : function
    Function to use for converting dates to strings. Defaults to
    dateutil.parser
names : array-like
    List of column names"""

_csv_sep = """sep : string, default None
    Delimiter to use. By default will try to automatically determine
    this"""

_table_sep = """sep : string, default \\t (tab-stop)
    Delimiter to use"""

read_csv.__doc__ = """
Read CSV (comma-separated) file into DataFrame

%s

Returns
-------
parsed : DataFrame
""" % (_parser_params % _csv_sep)

read_table.__doc__ = """
Read delimited file into DataFrame

%s

Returns
-------
parsed : DataFrame
""" % (_parser_params % _table_sep)


def _simple_parser(lines, names=None, header=0, index_col=0,
                   na_values=None, date_parser=None, parse_dates=True):
    """
    Workhorse function for processing nested list into DataFrame

    Should be replaced by np.genfromtxt eventually?
    """
    passed_names = names is not None
    if passed_names:
        names = list(names)
        header = None

    if header is not None:
        columns = []
        for i, c in enumerate(lines[header]):
            if c == '':
                columns.append('Unnamed: %d' % i)
            else:
                columns.append(c)

        content = lines[header+1:]

        counts = {}
        for i, col in enumerate(columns):
            cur_count = counts.get(col, 0)
            if cur_count > 0:
                columns[i] = '%s.%d' % (col, cur_count)
            counts[col] = cur_count + 1
    else:
        ncols = len(lines[0])
        if not names:
            columns = ['X.%d' % (i + 1) for i in range(ncols)]
        else:
            columns = names
        content = lines

    # spaghetti

    # implicitly index_col=0 b/c 1 fewer column names
    index_name = None
    implicit_first_col = (len(content) > 0 and
                          len(content[0]) == len(columns) + 1)

    if implicit_first_col:
        if index_col is None:
            index_col = 0
        index_name = None
    elif np.isscalar(index_col):
        if passed_names:
            index_name = None
        else:
            index_name = columns.pop(index_col)
    elif index_col is not None:
        if not passed_names:
            cp_cols = list(columns)
            index_name = []
            for i in index_col:
                name = cp_cols[i]
                columns.remove(name)
                index_name.append(name)
        else:
            index_name=None

    if len(content) == 0: # pragma: no cover
        if index_col is not None:
            if np.isscalar(index_col):
                index = Index([], name=index_name)
            else:
                index = MultiIndex.fromarrays([[]] * len(index_col),
                                              names=index_name)
        else:
            index = Index([])

        return DataFrame(index=index, columns=columns)

    # common NA values
    # no longer excluding inf representations
    # '1.#INF','-1.#INF', '1.#INF000000',
    NA_VALUES = set(['-1.#IND', '1.#QNAN', '1.#IND', '-1.#QNAN',
                     '#N/A N/A', 'NA', '#NA', 'NULL', 'NaN',
                     'nan', ''])
    if na_values is None:
        na_values = NA_VALUES
    else:
        na_values = set(list(na_values)) | NA_VALUES

    zipped_content = list(lib.to_object_array(content).T)

    # no index column specified, so infer that's what is wanted
    if index_col is not None:
        if np.isscalar(index_col):
            index = zipped_content.pop(index_col)
        else: # given a list of index
            index = []
            for idx in index_col:
                index.append(zipped_content[idx])
            #remove index items from content and columns, don't pop in loop
            for i in range(len(index_col)):
                zipped_content.remove(index[i])

        if np.isscalar(index_col):
            if parse_dates:
                index = lib.try_parse_dates(index, parser=date_parser)
            index = Index(_convert_types(index, na_values),
                          name=index_name)
        else:
            arrays = _maybe_convert_int_mindex(index, parse_dates,
                                               date_parser)
            index = MultiIndex.from_arrays(arrays, names=index_name)
    else:
        index = Index(np.arange(len(content)))

    if not index._verify_integrity():
        dups = index._get_duplicates()
        raise Exception('Index has duplicates: %s' % str(dups))

    if len(columns) != len(zipped_content):
        raise Exception('wrong number of columns')

    data = dict((k, v) for k, v in zip(columns, zipped_content))
    data = _convert_to_ndarrays(data, na_values)
    return DataFrame(data=data, columns=columns, index=index)


def _maybe_convert_int_mindex(index, parse_dates, date_parser):
    for i in range(len(index)):
        try:
            int(index[i][0])
            index[i] = map(int, index[i])
        except ValueError:
            if parse_dates:
                index[i] = lib.try_parse_dates(index[i], date_parser)

    return index

def _convert_to_ndarrays(dct, na_values):
    result = {}
    for c, values in dct.iteritems():
        result[c] = _convert_types(values, na_values)
    return result

def _convert_types(values, na_values):
    try:
        values = lib.maybe_convert_numeric(values, na_values)
    except Exception:
        lib.sanitize_objects(values)

    if values.dtype == np.object_:
        return lib.maybe_convert_bool(values)
    return values

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
              parse_dates=False, date_parser=None, na_values=None):
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

        if skiprows is None:
            skiprows = set()
        else:
            skiprows = set(skiprows)

        data = []
        for i in range(sheet.nrows):
            if i in skiprows:
                continue

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
        return _simple_parser(data, header=header, index_col=index_col,
                              parse_dates=parse_dates, date_parser=date_parser,
                              na_values=na_values)
