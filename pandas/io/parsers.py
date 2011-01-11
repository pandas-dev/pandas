"""
Module contains tools for processing files into DataFrames or other objects
"""

from datetime import datetime, timedelta
from itertools import izip
import string

import numpy as np

from pandas.core.index import Index
from pandas.core.frame import DataFrame

NA_VALUES = set(['-1.#IND', '1.#QNAN', '1.#IND',
                 '-1.#QNAN','1.#INF','-1.#INF', '1.#INF000000',
                 'NA', 'NULL', 'NaN', 'nan', ''])

def simpleParser(nestedList, colNames=None, header=0, indexCol=0):
    """
    Workhorse function for processing nested list into DataFrame

    Should be replaced by np.genfromtxt
    """
    try:
        from dateutil import parser
        parse_date = parser.parse
    except ImportError:
        def parse_date(s):
            try:
                return datetime.strptime(s, '%m/%d/%Y')
            except Exception:
                return s
    lines = nestedList
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

    index_name = columns[indexCol]

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

    for i, (name, values) in enumerate(izip(columns, izip(*content))):
        if i == indexCol:
            data[name] = values
            continue
        data[name] = [_convert_float(val) for val in values]

    # try to parse dates
    try:
        # easier to ask forgiveness than permission
        result = [parse_date(idx) for idx in data[index_name]]
        data[index_name] = result
    except Exception:
        pass

    for c, values in data.iteritems():
        try:
            data[c] = np.array(values, dtype=float)
        except Exception:
            data[c] = np.array(values, dtype=object)

    if indexCol is not None:
        index = Index(data.pop(index_name))
        return DataFrame(data=data, index=index)
    else:
        index = np.arange(len(data.values()[0]))
        return DataFrame(data=data, index=index)

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

def parseText(filepath, sep='\t', header=0, indexCol=0, colNames = None):
    """
    Parse whitespace separated file into a DataFrame object.
    Try to parse dates if possible.
    """
    lines = [l.rstrip().split(sep) for l in open(filepath,'rb').readlines()]
    return simpleParser(lines, header=header, indexCol=indexCol,
                        colNames = colNames)

#===============================================================================
# Excel tools
#===============================================================================

OLE_TIME_ZERO = datetime(1899, 12, 30, 0, 0, 0)
def ole2datetime(oledt):
    """function for converting excel date to normal date format"""
    return OLE_TIME_ZERO + timedelta(days=float(oledt))


def parseExcel(filepath, header = None, indexCol = 0, dateCol = 0,
               sheetname = None):
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
