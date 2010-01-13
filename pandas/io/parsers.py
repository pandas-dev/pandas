"""
Module contains tools for processing files into DataFrames or other objects
"""

from datetime import datetime, timedelta
from itertools import izip
import string

from dateutil import parser
import numpy as np

from pandas.core.index import Index
from pandas.core.frame import DataFrame

def simpleParser(nestedList, colNames=None, header=0, indexCol=0):
    """
    Workhorse function for processing nested list into DataFrame

    Should be replaced by np.genfromtxt
    """
    naValues = set(['-1.#IND', '1.#QNAN', '1.#IND',
                    '-1.#QNAN','1.#INF','-1.#INF', '1.#INF000000',
                    'NA', 'NULL', 'NaN', 'nan', ''])
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

    for i, (c, col) in enumerate(izip(columns, izip(*content))):
        if i == indexCol:
            data[c] = col
            continue
        data[c] = []
        for val in col:
            if val in naValues:
                val = np.nan
            else:
                try:
                    tmp = val
                    val = np.float64(val)
                    if np.isinf(val):
                        val = tmp
                except Exception:
                    pass
            data[c].append(val)

    if header is not None:
        if 'date' in columns[0].lower() or 'Unnamed' in columns[0]:
            dates = []
            for s in data[columns[0]]:
                try:
                    dates.append(parser.parse(s))
                except Exception:
                    dates.append(s)
            data[columns[0]] = dates
    for c, values in data.iteritems():
        try:
            data[c] = np.array(values, dtype = np.float64)
        except Exception:
            data[c] = np.array(values, dtype = np.object_)
    if indexCol is not None:
        index = Index(data[columns[indexCol]])
        frameData = dict([(col, data[col]) for col in columns \
                        if col != columns[indexCol]])
        return DataFrame(data=frameData, index=index)
    else:
        index = np.arange(len(data.values()[0]))
        frameData = dict([(col, data[col]) for col in columns])
        return DataFrame(data=frameData, index=index)

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
