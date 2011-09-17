"""
Module contains tools for processing files into DataFrames or other objects
"""

from datetime import datetime
from itertools import izip
import re
import string

import numpy as np

from pandas.core.index import Index, MultiIndex
from pandas.core.frame import DataFrame

def read_csv(filepath_or_buffer, sep=None, header=0, skiprows=None, index_col=0,
             na_values=None, date_parser=None, names=None):
    """
    Read CSV file into DataFrame

    Parameters
    ----------
    filepath_or_buffer : string or file handle / StringIO
    sep : string, default None
        Delimiter to use. By default will try to automatically determine
        this
    header : int, default 0
        Row to use for the column labels of the parsed DataFrame
    skiprows : list-like
        Row numbers to skip (0-indexed)
    index_col : int or sequence., default 0
        Column to use as the row labels of the DataFrame. Pass None if there is
        no such column. If a sequence is given, a MultiIndex is used.
    na_values : list-like, default None
        List of additional strings to recognize as NA/NaN
    date_parser : function
        Function to use for converting dates to strings. Defaults to
        dateutil.parser
    names : array-like
        List of column names

    Returns
    -------
    parsed : DataFrame
    """
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
    return _simple_parser(lines, header=header, indexCol=index_col,
                          colNames=names, na_values=na_values,
                          date_parser=date_parser)

def read_table(filepath_or_buffer, sep='\t', header=0, skiprows=None,
               index_col=0, na_values=None, date_parser=None, names=None):
    """
    Read delimited file into DataFrame

    Parameters
    ----------
    filepath_or_buffer : string or file handle
    sep : string, default '\t'
        Delimiter to use
    header : int, default 0
        Row to use for the column labels of the parsed DataFrame
    skiprows : list-like
        Row numbers to skip (0-indexed)
    index_col : int or sequence, default 0
        Column to use as the row labels of the DataFrame. Pass None if there is
        no such column. If a sequence is given, a MultiIndex is used.
    na_values : list-like, default None
        List of additional strings to recognize as NA/NaN
    date_parser : function
        Function to use for converting dates to strings. Defaults to
        dateutil.parser
    names : array-like
        List of column names

    Returns
    -------
    parsed : DataFrame
    """
    return read_csv(filepath_or_buffer, sep, header, skiprows,
                    index_col, na_values, date_parser, names)

def _simple_parser(lines, colNames=None, header=0, indexCol=0,
                   na_values=None, date_parser=None, parse_dates=True):
    """
    Workhorse function for processing nested list into DataFrame

    Should be replaced by np.genfromtxt eventually?
    """
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
        if not colNames:
            columns = ['X.%d' % (i + 1) for i in range(ncols)]
        else:
            assert(len(colNames) == ncols)
            columns = colNames
        content = lines

    zipped_content = zip(*content)

    if len(content) == 0: # pragma: no cover
        raise Exception('No content to parse')

    # no index column specified, so infer that's what is wanted
    if indexCol is not None:
        if np.isscalar(indexCol):
            if indexCol == 0 and len(content[0]) == len(columns) + 1:
                index = zipped_content[0]
                zipped_content = zipped_content[1:]
            else:
                index = zipped_content.pop(indexCol)
                columns.pop(indexCol)
        else: # given a list of index
            idx_names = []
            index = []
            for idx in indexCol:
                idx_names.append(columns[idx])
                index.append(zipped_content[idx])
            #remove index items from content and columns, don't pop in loop
            for i in range(len(indexCol)):
                columns.remove(idx_names[i])
                zipped_content.remove(index[i])


        if np.isscalar(indexCol):
            if parse_dates:
                index = _try_parse_dates(index, parser=date_parser)
            index = Index(_maybe_convert_int(np.array(index, dtype=object)))
        else:
            index = MultiIndex.from_arrays(_maybe_convert_int_mindex(index,
                                                 parse_dates, date_parser),
                                                 names=idx_names)
    else:
        index = Index(np.arange(len(content)))

    if len(columns) != len(zipped_content):
        raise Exception('wrong number of columns')

    data = dict(izip(columns, zipped_content))
    data = _floatify(data, na_values=na_values)
    data = _convert_to_ndarrays(data)
    return DataFrame(data=data, columns=columns, index=index)

def _floatify(data_dict, na_values=None):
    """

    """
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

    def _convert_float(val):
        if val in na_values:
            return np.nan
        else:
            try:
                return np.float64(val)
            except Exception:
                return val

    result = {}
    for col, values in data_dict.iteritems():
        result[col] = [_convert_float(val) for val in values]

    return result

def _maybe_convert_int(arr):
    if len(arr) == 0: # pragma: no cover
        return arr

    try:
        if arr.dtype == np.object_:
            return arr.astype(int)

        if abs(arr[0] - int(arr[0])) < 1e-10:
            casted = arr.astype(int)
            if (np.abs(casted - arr) < 1e-10).all():
                return casted
    except (TypeError, ValueError):
        pass

    return arr

def _maybe_convert_int_mindex(index, parse_dates, date_parser):
    if len(index) == 0:
        return index

    for i in range(len(index)):
        try:
            int(index[i][0])
            index[i] = map(int, index[i])
        except ValueError:
            if parse_dates:
                index[i] = _try_parse_dates(index[i], date_parser)

    return index

def _convert_to_ndarrays(dct):
    result = {}
    for c, values in dct.iteritems():
        try:
            values = np.array(values, dtype=float)
        except Exception:
            values = np.array(values, dtype=object)
        result[c] = _maybe_convert_int(values)
    return result

def _try_parse_dates(values, parser=None):
    if parser is None:
        try:
            from dateutil import parser
            parse_date = parser.parse
        except ImportError: # pragma: no cover
            def parse_date(s):
                try:
                    return datetime.strptime(s, '%m/%d/%Y')
                except Exception:
                    return s
    else:
        parse_date = parser

    # EAFP
    try:
        return [parse_date(val) for val in values]
    except Exception:
        # failed
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

    def parse(self, sheetname, header=0, skiprows=None, index_col=0,
              na_values=None):
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
        return _simple_parser(data, header=header, indexCol=index_col,
                              na_values=na_values)

#-------------------------------------------------------------------------------
# Deprecated stuff

import warnings

def parseCSV(filepath, header=0, skiprows=None, indexCol=0,
             na_values=None): # pragma: no cover
    """
    Parse CSV file into a DataFrame object. Try to parse dates if possible.
    """
    warnings.warn("parseCSV is deprecated. Use read_csv instead", FutureWarning)
    return read_csv(filepath, header=header, skiprows=skiprows,
                    index_col=indexCol, na_values=na_values)

def parseText(filepath, sep='\t', header=0,
              indexCol=0, colNames=None): # pragma: no cover
    """
    Parse whitespace separated file into a DataFrame object.
    Try to parse dates if possible.
    """
    warnings.warn("parseText is deprecated. Use read_table instead",
                  FutureWarning)
    return read_table(filepath, sep=sep, header=header, index_col=indexCol,
                      names=colNames)


def parseExcel(filepath, header=None, indexCol=0,
               sheetname=None, **kwds): # pragma: no cover
    """

    """
    warnings.warn("parseExcel is deprecated. Use the ExcelFile class instead",
                  FutureWarning)
    excel_file = ExcelFile(filepath)
    return excel_file.parse(sheetname, header=header, index_col=indexCol)


