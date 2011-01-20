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

def parseCSV(filepath, header=0, indexCol=0):
    """
    Parse CSV file into a DataFrame object. Try to parse dates if possible.
    """
    import csv
    f = open(filepath,'rb')
    reader = csv.reader(f, dialect='excel')
    lines = [l for l in reader]
    f.close()
    return simpleParser(lines, header=header, indexCol=indexCol)

def read_table(path, header=0, index_col=0, delimiter=','):
    data = np.genfromtext(path, delimiter=delimiter,
                          names=header is not None,
                          dtype=object)

    columns = data.dtype.names

def parseText(filepath, sep='\t', header=0, indexCol=0, colNames = None):
    """
    Parse whitespace separated file into a DataFrame object.
    Try to parse dates if possible.
    """
    lines = [re.split(sep, l.rstrip()) for l in open(filepath,'rb').readlines()]
    return simpleParser(lines, header=header, indexCol=indexCol,
                        colNames = colNames)

def simpleParser(lines, colNames=None, header=0, indexCol=0):
    """
    Workhorse function for processing nested list into DataFrame

    Should be replaced by np.genfromtxt eventually?
    """
    data = {}
    if header is not None:
        columns = []
        for i, c in enumerate(lines[header]):
            if c == '':
                columns.append('Unnamed: ' + string.ascii_uppercase[i])
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
            columns = string.ascii_uppercase[:len(lines[0])]
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

    data = _floatify(data)
    data = _convert_to_ndarrays(data)
    return DataFrame(data=data, index=Index(index))

NA_VALUES = set(['-1.#IND', '1.#QNAN', '1.#IND',
                 '-1.#QNAN','1.#INF','-1.#INF', '1.#INF000000',
                 'NA', '#NA', 'NULL', 'NaN', 'nan', ''])

def _floatify(data_dict):
    result = {}
    for col, values in data_dict.iteritems():
        result[col] = [_convert_float(val) for val in values]

    return result

def _convert_float(val):
    if val in NA_VALUES:
        return np.nan
    else:
        try:
            parsed = np.float64(val)
            if np.isinf(parsed):
                return val
            return parsed
        except Exception:
            return val

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



def parseExcel(filepath, header = None, indexCol = 0, dateCol = 0,
               sheetname = None):
    from pandas.core.datetools import ole2datetime
    try:
        import xlrd
    except:
        raise Exception('Sorry, you do not have xlrd.')
    book = xlrd.open_workbook(filepath)
    sheet = book.sheet_by_name(sheetname)
    data = [sheet.row_values(i) for i in range(sheet.nrows)]
    if dateCol is not None:
        for row in data:
            try:
                row[dateCol] = ole2datetime(row[dateCol])
            except Exception:
                pass
    return simpleParser(data, header = header, indexCol = indexCol)
