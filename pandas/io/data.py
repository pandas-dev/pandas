"""
Module contains tools for collecting data from various remote sources


"""

import numpy as np
import datetime as dt
import urllib
import urllib2

from zipfile import ZipFile
from StringIO import StringIO

from pandas import DataFrame, read_csv

def DataReader(name, data_source=None, start=None, end=None):
    """
    Imports data from a number of online sources.

    Currently supports Yahoo! finance, St. Louis FED (FRED), and Kenneth
    French's data library.

    Parameters
    ----------
    name : str
        the name of the dataset
    data_source: str
        the data source ("yahoo", "fred", or "ff")
    start : {datetime, None}
        left boundary for range (defaults to 1/1/2010)
    end : {datetime, None}
        right boundary for range (defaults to today)

    Examples
    ----------

    # Data from Yahoo!
    gs = DataReader("GS", "yahoo")

    # Data from FRED
    vix = DataReader("VIXCLS", "fred")

    # Data from Fama/French
    ff = DataReader("F-F_Research_Data_Factors", "famafrench")
    ff = DataReader("F-F_Research_Data_Factors_weekly", "famafrench")
    ff = DataReader("6_Portfolios_2x3", "famafrench")
    ff = DataReader("F-F_ST_Reversal_Factor", "famafrench")
    """
    start, end = _sanitize_dates(start, end)

    if(data_source == "yahoo"):
        return get_data_yahoo(name=name, start=start, end=end)
    elif(data_source == "fred"):
        return get_data_fred(name=name, start=start, end=end)
    elif(data_source == "famafrench"):
        return get_data_famafrench(name=name)

def _sanitize_dates(start, end):
    from pandas.core.datetools import to_datetime
    start = to_datetime(start)
    end = to_datetime(end)
    if start is None:
        start = dt.datetime.today() - dt.timedelta(365)
    if end is None:
        end = dt.datetime.today()
    return start, end

def get_quote_yahoo(symbols):
    """ 
    Get current yahoo quote
    
    Returns a DataFrame
    """
    if not isinstance(symbols,list):
        raise TypeError, "symbols must be a list"
    # for codes see: http://www.gummy-stuff.org/Yahoo-data.htm
    codes = {'symbol':'s','last':'l1','change_pct':'p2','PE':'r','time':'t1','short_ratio':'s7'}
    request = str.join('',codes.values()) # code request string
    header = codes.keys()   
    
    data = dict(zip(codes.keys(),[[] for i in range(len(codes))]))
    
    urlStr = 'http://finance.yahoo.com/d/quotes.csv?s=%s&f=%s' % (str.join('+',symbols), request)
    
    try:
        lines = urllib2.urlopen(urlStr).readlines()
    except Exception, e:
        s = "Failed to download:\n{0}".format(e);
        print s
        return None
        
    for line in lines:
        fields = line.strip().split(',')
        #print fields
        for i,field in enumerate(fields):
            if field[0] == '"':
                data[header[i]].append( field.strip('"'))
            else:
                try:
                    data[header[i]].append(float(field))
                except ValueError:
                    data[header[i]].append(np.nan)

    idx = data.pop('symbol')
    
    return DataFrame(data,index=idx)

def get_data_yahoo(name=None, start=None, end=None):
    """
    Get historical data for the given name from yahoo.
    Date format is datetime

    Returns a DataFrame.
    """
    start, end = _sanitize_dates(start, end)

    if(name is None):
        print "Need to provide a name"
        return None

    yahoo_URL = 'http://ichart.yahoo.com/table.csv?'

    url = yahoo_URL + 's=%s' % name + \
      '&a=%s' % (start.month - 1) + \
      '&b=%s' % start.day + \
      '&c=%s' % start.year + \
      '&d=%s' % (end.month - 1) + \
      '&e=%s' % end.day + \
      '&f=%s' % end.year + \
      '&g=d' + \
      '&ignore=.csv'

    lines = urllib.urlopen(url).read()
    return read_csv(StringIO(lines), index_col=0, parse_dates=True)[::-1]

def get_data_fred(name=None, start=dt.datetime(2010, 1, 1),
                  end=dt.datetime.today()):
    """
    Get data for the given name from the St. Louis FED (FRED).
    Date format is datetime

    Returns a DataFrame.
    """
    start, end = _sanitize_dates(start, end)

    if(name is None):
        print "Need to provide a name"
        return None

    fred_URL = "http://research.stlouisfed.org/fred2/series/"

    url = fred_URL + '%s' % name + \
      '/downloaddata/%s' % name + '.csv'
    data = read_csv(urllib.urlopen(url), index_col=0, parse_dates=True)
    return data.truncate(start, end)

def get_data_famafrench(name, start=None, end=None):
    start, end = _sanitize_dates(start, end)

    # path of zip files
    zipFileURL = "http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/"

    url = urllib.urlopen(zipFileURL + name + ".zip")
    zipfile = ZipFile(StringIO(url.read()))
    data = zipfile.open(name + ".txt").readlines()

    file_edges = np.where(np.array([len(d) for d in data]) == 2)[0]

    datasets = {}
    for i in range(len(file_edges)-1):
        dataset = [d.split() for d in data[(file_edges[i] + 1):file_edges[i+1]]]
        if(len(dataset) > 10):
            ncol = np.median(np.array([len(d) for d in dataset]))
            header_index = np.where(np.array([len(d) for d in dataset]) == (ncol-1))[0][-1]
            header = dataset[header_index]
            # to ensure the header is unique
            header = [str(j + 1) + " " + header[j] for j in range(len(header))]
            index = np.array([d[0] for d in dataset[(header_index + 1):]], dtype=int)
            dataset = np.array([d[1:] for d in dataset[(header_index + 1):]], dtype=float)
            datasets[i] = DataFrame(dataset, index, columns=header)

    return datasets

