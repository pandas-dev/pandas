"""
Module contains tools for processing files into DataFrames or other objects
"""

from datetime import datetime, timedelta
from itertools import izip
import re
import string

import numpy as np

from pandas.core.index import Index
from pandas.core.frame import DataFrame

def parseCSV(filepath, header=0, skiprows=None, indexCol=0,
             na_values=None):
    """
    Parse CSV file into a DataFrame object. Try to parse dates if possible.
    """
    import csv
    f = open(filepath,'U')
    reader = csv.reader(f, dialect='excel')

    if skiprows is not None:
        skiprows = set(skiprows)
        lines = [l for i, l in enumerate(reader) if i not in skiprows]
    else:
        lines = [l for l in reader]
    f.close()
    return simpleParser(lines, header=header, indexCol=indexCol,
                        na_values=na_values)

def read_table(path, header=0, index_col=0, delimiter=','):
    data = np.genfromtxt(path, delimiter=delimiter,
                         names=header is not None,
                         dtype=object)

    columns = data.dtype.names

def parseText(filepath, sep='\t', header=0, indexCol=0, colNames=None):
    """
    Parse whitespace separated file into a DataFrame object.
    Try to parse dates if possible.
    """
    lines = [re.split(sep, l.rstrip()) for l in open(filepath,'rb').readlines()]
    return simpleParser(lines, header=header, indexCol=indexCol,
                        colNames = colNames)



def simpleParser(lines, colNames=None, header=0, indexCol=0,
                 na_values=None):
    """
    Workhorse function for processing nested list into DataFrame

    Should be replaced by np.genfromtxt eventually?
    """
    data = {}
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
            columns = list(string.ascii_uppercase[:len(lines[0])])
            # columns = ['X.%d' % (i + 1) for i in range(len(lines[0]))]
        else:
            columns = colNames
        content = lines

    data = dict(izip(columns, izip(*content)))
    if indexCol is not None:
        index_name = columns[indexCol]
        # try to parse dates
        index = _try_parse_dates(data.pop(index_name))
    else:
        index = np.arange(len(data.values()[0]))

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
        na_values = set(list(na_values))

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

def _try_parse_dates(values):
    try:
        from dateutil import parser
        parse_date = parser.parse
    except ImportError:
        def parse_date(s):
            try:
                return datetime.strptime(s, '%m/%d/%Y')
            except Exception:
                return s

    try:
        # easier to ask forgiveness than permission
        return [parse_date(val) for val in values]
    except Exception:
        # failed
        return values

#===============================================================================
# Excel tools
#===============================================================================


class ExcelFile(object):
    """
    Class for parsing tabular .xls sheets into DataFrame objects, uses xlrd

    Parameters
    ----------
    path : string
        Path to xls file
    """

    def __init__(self, path):
        import xlrd
        self.path = path
        self.book = xlrd.open_workbook(path)

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
        return simpleParser(data, header=header, indexCol=index_col)

    def parse(self, sheetname, header=None, skiprows=None, index_col=0,
              na_values=None):
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
        return simpleParser(data, header=header, indexCol=index_col,
                            na_values=na_values)

def parseExcel(filepath, header=None, indexCol=0, sheetname=None, **kwds):
    """

    """
    excel_file = ExcelFile(filepath)
    return excel_file.parse(sheetname, header=header, index_col=indexCol)

