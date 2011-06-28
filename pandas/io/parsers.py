"""
Module contains tools for processing files into DataFrames or other objects
"""

from datetime import datetime
from itertools import izip
import re
import string

import numpy as np

from pandas.core.index import Index
from pandas.core.frame import DataFrame

def read_csv(filepath_or_buffer, header=0, skiprows=None, index_col=0,
             na_values=None, date_parser=None):
    """
    Read CSV file into DataFrame

    Parameters
    ----------
    filepath_or_buffer : string or file handle / StringIO

    header : int, default 0
        Row to use for the column labels of the parsed DataFrame
    skiprows : list-like
        Row numbers to skip (0-indexed)
    index_col : int, default 0
        Column to use as the row labels of the DataFrame. Pass None if there is
        no such column
    na_values : list-like, default None
        List of additional strings to recognize as NA/NaN
    date_parser : function
        Function to use for converting dates to strings. Defaults to
        dateutil.parser
    """
    import csv

    if hasattr(filepath_or_buffer, 'read'):
        f = filepath_or_buffer
    else:
        try:
            # universal newline mode
            f = open(filepath_or_buffer, 'U')
        except Exception:
            f = open(filepath_or_buffer, 'r')

    reader = csv.reader(f, dialect='excel')

    if skiprows is not None:
        skiprows = set(skiprows)
        lines = [l for i, l in enumerate(reader) if i not in skiprows]
    else:
        lines = [l for l in reader]
    f.close()
    return _simple_parser(lines, header=header, indexCol=index_col,
                          na_values=na_values, date_parser=date_parser)

def read_table(filepath_or_buffer, sep='\t', header=0, skiprows=None, index_col=0,
               na_values=None, names=None, date_parser=None):
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
    index_col : int, default 0
        Column to use as the row labels of the DataFrame. Pass None if there is
        no such column
    na_values : list-like, default None
        List of additional strings to recognize as NA/NaN
    date_parser : function
        Function to use for converting dates to strings. Defaults to
        dateutil.parser
    """
    reader = open(filepath_or_buffer,'rb')

    if skiprows is not None:
        skiprows = set(skiprows)
        lines = [l for i, l in enumerate(reader) if i not in skiprows]
    else:
        lines = [l for l in reader]

    lines = [re.split(sep, l.rstrip()) for l in lines]
    return _simple_parser(lines, header=header, indexCol=index_col,
                          colNames=names, na_values=na_values,
                          date_parser=date_parser)

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

        colCounts = dict(((col, 0) for col in columns))
        for i, col in enumerate(columns):
            if columns.count(col) > 1:
                columns[i] = col + str(colCounts[col])
                colCounts[col] += 1
    else:
        if not colNames:
            # columns = list(string.ascii_uppercase[:len(lines[0])])
            columns = ['X.%d' % (i + 1) for i in range(len(lines[0]))]
        else:
            columns = colNames
        content = lines

    zipped_content = zip(*content)

    if len(content) == 0:
        raise Exception('No content to parse')

    # no index column specified, so infer that's what is wanted
    if indexCol is not None:
        if indexCol == 0 and len(content[0]) == len(columns) + 1:
            index = zipped_content[0]
            zipped_content = zipped_content[1:]
        else:
            index = zipped_content.pop(indexCol)
            columns.pop(indexCol)

        if parse_dates:
            index = _try_parse_dates(index, parser=date_parser)

    else:
        index = np.arange(len(content))

    data = dict(izip(columns, zipped_content))
    data = _floatify(data, na_values=na_values)
    data = _convert_to_ndarrays(data)
    return DataFrame(data=data, columns=columns, index=Index(index))

def _floatify(data_dict, na_values=None):
    # common NA values
    NA_VALUES = set(['-1.#IND', '1.#QNAN', '1.#IND',
                     '-1.#QNAN','1.#INF','-1.#INF', '1.#INF000000',
                     'NA', '#NA', 'NULL', 'NaN', 'nan', ''])
    if na_values is None:
        na_values = NA_VALUES
    else:
        na_values = set(list(na_values)) | NA_VALUES

    def _convert_float(val):
        if val in na_values:
            return np.nan
        else:
            try:
                parsed = np.float64(val)
                if np.isinf(parsed):
                    return val
                return parsed
            except Exception:
                return val

    result = {}
    for col, values in data_dict.iteritems():
        result[col] = [_convert_float(val) for val in values]

    return result

def _convert_to_ndarrays(dct):
    result = {}
    for c, values in dct.iteritems():
        try:
            result[c] = np.array(values, dtype=float)
        except Exception:
            result[c] = np.array(values, dtype=object)

    return result

def _try_parse_dates(values, parser=None):
    if parser is None:
        try:
            from dateutil import parser
            parse_date = parser.parse
        except ImportError:
            def parse_date(s):
                try:
                    return datetime.strptime(s, '%m/%d/%Y')
                except Exception:
                    return s
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

    def old_parse(self, sheetname, header=None, index_col=0, date_col=0):
        from pandas.core.datetools import ole2datetime
        sheet = self.book.sheet_by_name(sheetname)

        data = [sheet.row_values(i) for i in range(sheet.nrows)]
        if date_col is not None:
            for row in data:
                try:
                    row[date_col] = ole2datetime(row[date_col])
                except Exception:
                    pass
        return _simple_parser(data, header=header, indexCol=index_col)

    def parse(self, sheetname, header=None, skiprows=None, index_col=0,
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
                    if dt[0] < MINYEAR:
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
             na_values=None):
    """
    Parse CSV file into a DataFrame object. Try to parse dates if possible.
    """
    warnings.warn("parseCSV is deprecated. Use read_csv instead", FutureWarning)
    return read_csv(filepath, header=header, skiprows=skiprows,
                    index_col=indexCol, na_values=na_values)

def parseText(filepath, sep='\t', header=0, indexCol=0, colNames=None):
    """
    Parse whitespace separated file into a DataFrame object.
    Try to parse dates if possible.
    """
    warnings.warn("parseText is deprecated. Use read_table instead",
                  FutureWarning)
    return read_table(filepath, sep=sep, header=header, index_col=indexCol,
                      names=colNames)


def parseExcel(filepath, header=None, indexCol=0, sheetname=None, **kwds):
    """

    """
    warnings.warn("parseExcel is deprecated. Use the ExcelFile class instead",
                  FutureWarning)
    excel_file = ExcelFile(filepath)
    return excel_file.parse(sheetname, header=header, index_col=indexCol)


