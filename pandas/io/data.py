"""
Module contains tools for collecting data from various http sources

Examples:

vix = get_data_fred("VIXCLS")
gs = get_data_yahoo("GS")

"""

import numpy as np
import datetime as dt
import urllib

from pandas import DataFrame, Index

def get_data_yahoo(symbol=None, start_date=dt.datetime(2010, 1, 1), end_date=dt.datetime.today()):
    """
    Get historical data for the given symbol from yahoo.
    Date format is datetime
    
    Returns a DataFrame.
    """
    
    if(symbol is None):
    	print "Need to provide a symbol"
    	return None
    
    yahoo_URL = 'http://ichart.yahoo.com/table.csv?'
    
    url = yahoo_URL + 's=%s' % symbol + \
          '&a=%s' % start_date.month + \
          '&b=%s' % start_date.day + \
          '&c=%s' % start_date.year + \
          '&d=%s' % end_date.month + \
          '&e=%s' % end_date.day + \
          '&f=%s' % end_date.year + \
          '&g=d' + \
          '&ignore=.csv'

    days = urllib.urlopen(url).readlines()
    
    data = np.array([day[:-2].split(',') for day in days])
    header = [str.lower(name) for name in data[0]]
    index = Index([dt.datetime.strptime(row[0], "%Y-%m-%d") for row in data[1:]])
    data = np.array([[row[1], row[2], row[3], row[4], int(row[5]), row[6]] for row in data[1:]], dtype=float)
    
    data = DataFrame(data, index, columns=header[1:]).sort_index()
    
    return data
    

def get_data_fred(symbol=None, start_date=dt.datetime(2010, 1, 1), end_date=dt.datetime.today()):
    """
    Get data for the given symbol from the St. Louis FED (FRED).
    Date format is datetime
    
    Returns a DataFrame.
    """
    
    if(symbol is None):
    	print "Need to provide a symbol"
    	return None
    
    fred_URL = "http://research.stlouisfed.org/fred2/series/"

    url = fred_URL + '%s' % symbol + \
          '/downloaddata/%s' % symbol + '.csv'

    days = urllib.urlopen(url).readlines()
    
    data = np.array([day[:-2].split(',') for day in days])
    header = [str.lower(name) for name in data[0]]
    index = Index([dt.datetime.strptime(row[0], "%Y-%m-%d") for row in data[1:]])
    data = np.array([[row[1]] for row in data[1:]], dtype=float)
    
    data = DataFrame(data, index, columns=header[1:]).sort_index().truncate(start_date, end_date)
    
    return data
    

