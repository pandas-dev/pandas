"""
Module contains tools for collecting data from various remote sources


"""

import numpy as np
import datetime as dt
import urllib

from zipfile import ZipFile
from StringIO import StringIO

from pandas import DataFrame, read_csv

class DataReader(list):

    """
    Imports data from a number of online sources.

    Currently supports Yahoo! finance, St. Louis FED (FRED), and Kenneth French's data library.

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
    def __new__(cls, name, data_source=None, start=None, end=None, **kwds):

        if(start is None):
            start = dt.datetime(2010, 1, 1)
        if(end is None):
            end = dt.datetime.today()

        self = super(DataReader, cls).__new__(cls)

        if(data_source == "yahoo"):
            return self.get_data_yahoo(name=name, start=start, end=end)
        elif(data_source == "fred"):
            return self.get_data_fred(name=name, start=start, end=end)
        elif(data_source == "famafrench"):
            return self.get_data_famafrench(name=name)

    def get_data_yahoo(self, name=None, start=None, end=None):
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
          '&a=%s' % start.month + \
          '&b=%s' % start.day + \
          '&c=%s' % start.year + \
          '&d=%s' % end.month + \
          '&e=%s' % end.day + \
          '&f=%s' % end.year + \
          '&g=d' + \
          '&ignore=.csv'

        lines = urllib.urlopen(url).read()
        return read_csv(StringIO(lines), index_col=0, parse_dates=True)

    def get_data_fred(self, name=None, start=dt.datetime(2010, 1, 1),
                      end=dt.datetime.today()):
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

        data = DataFrame(data, index, columns=header[1:]).sort().truncate(start, end)

        return data


    def get_data_famafrench(self, name):

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

