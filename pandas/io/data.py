"""
Module contains tools for collecting data from various remote sources

Examples:

# Data from FRED
vix = get_data_fred("VIXCLS")

# Data from Yahoo!
gs = get_data_yahoo("GS")

# Data from Fama/French
ff = get_data_famafrench("F-F_Research_Data_Factors")
ff = get_data_famafrench("F-F_Research_Data_Factors_weekly")
ff = get_data_famafrench("6_Portfolios_2x3")

"""

import numpy as np
import datetime as dt
import urllib
import zipfile

from pandas import DataFrame, Index

def get_data_yahoo(name=None, start_date=dt.datetime(2010, 1, 1), end_date=dt.datetime.today()):
    """
    Get historical data for the given name from yahoo.
    Date format is datetime
    
    Returns a DataFrame.
    """
    
    if(name is None):
    	print "Need to provide a name"
    	return None
    
    yahoo_URL = 'http://ichart.yahoo.com/table.csv?'
    
    url = yahoo_URL + 's=%s' % name + \
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
    

def get_data_fred(name=None, start_date=dt.datetime(2010, 1, 1), end_date=dt.datetime.today()):
    """
    Get data for the given name from the St. Louis FED (FRED).
    Date format is datetime
    
    Returns a DataFrame.
    """
    
    if(name is None):
    	print "Need to provide a name"
    	return None
    
    fred_URL = "http://research.stlouisfed.org/fred2/series/"

    url = fred_URL + '%s' % name + \
          '/downloaddata/%s' % name + '.csv'

    days = urllib.urlopen(url).readlines()
    
    data = np.array([day[:-2].split(',') for day in days])
    header = [str.lower(name) for name in data[0]]
    index = np.array([dt.datetime.strptime(row[0], "%Y-%m-%d") for row in data[1:]])
    data = np.array([[row[1]] for row in data[1:]], dtype=float)
    
    data = DataFrame(data, index, columns=header[1:]).sort_index().truncate(start_date, end_date)
    
    return data
    

def get_data_famafrench(name):

	# path of zip files
	zipFileURL = "http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/"

	from zipfile import ZipFile

	url = urllib.urlopen(zipFileURL + name + ".zip")
	zipfile = ZipFile(StringIO(url.read()))
	data = zipfile.open(name + ".txt").readlines()
	
	file_edges = np.where(np.array([len(d) for d in data]) == 2)[0]
	
	datasets = {}
	for i in range(len(file_edges)-1):
		dataset = np.array([d.split() for d in data[(file_edges[i] + 1):file_edges[i+1]]])
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

