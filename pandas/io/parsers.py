"""
Module contains tools for processing files into DataFrames or other objects
"""

from pandas.core.index import Index
from pandas.core.frame import DataFrame
from pandas.core.matrix import DataMatrix
from pandas.core.series import Series

from datetime import datetime, timedelta

try:
    from dateutil import parser
except ImportError:
    # just a little hack for now
    class parser(object):
        @classmethod
        def parse(cls, val):
            try:
                return datetime.strptime(val, '%m/%d/%Y')
            except:
                return val

from itertools import izip
import numpy as np
import string

def simpleParser(nestedList, forceFloat=True, colNames=None,
                 header=0, indexCol=0):
    """
    Workhorse function for processing nested list into DataFrame
    """
    naValues = set(['-1.#IND', '1.#QNAN', '1.#IND', 
                    '-1.#QNAN','1.#INF','-1.#INF', '1.#INF000000',
                    'NA', 'NULL', 'NaN', 'nan', ''])
    lines = nestedList
    data = {}
    if header is not None:
        columns = lines[header]
        columns = [c if c != '' else 'Unnamed: ' + string.ascii_uppercase[i] 
                   for i, c in enumerate(columns)]
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
                    if isinf(val):
                        val = tmp
                except:
                    pass
            data[c].append(val)

    if header is not None:
        if 'date' in columns[0].lower() or 'Unnamed' in columns[0]:
            dates = []
            for s in data[columns[0]]:
                try:
                    dates.append(parser.parse(s))
                except:
                    dates.append(s)
            data[columns[0]] = dates
    for c, values in data.iteritems():
        try:
            data[c] = np.array(values, dtype = np.float64)
        except:
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
    except ImportError:
        raise ImportError('Sorry, you do not have xlrd.')
    book = xlrd.open_workbook(filepath)
    sheet = book.sheet_by_name(sheetname)
    data = [sheet.row_values(i) for i in range(sheet.nrows)]
    for row in data:
        try:
            row[0] = ole2datetime(row[0])
        except:
            pass
    return simpleParser(data, header = header, indexCol = indexCol)
