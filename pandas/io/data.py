"""
Module contains tools for collecting data from various remote sources


"""

import numpy as np
import datetime as dt
import urllib
import urllib2
import time

from zipfile import ZipFile
from pandas.util.py3compat import StringIO, BytesIO, bytes_to_str

from pandas import DataFrame, read_csv, concat
from pandas.io.parsers import TextParser


def DataReader(name, data_source=None, start=None, end=None,
        retry_count=3, pause=0):
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
        return get_data_yahoo(name=name, start=start, end=end,
                   retry_count=retry_count, pause=pause)
    elif(data_source == "fred"):
        return get_data_fred(name=name, start=start, end=end)
    elif(data_source == "famafrench"):
        return get_data_famafrench(name=name)


def _sanitize_dates(start, end):
    from pandas.core.datetools import to_datetime
    start = to_datetime(start)
    end = to_datetime(end)
    if start is None:
        start = dt.datetime(2010, 1, 1)
    if end is None:
        end = dt.datetime.today()
    return start, end


def get_quote_yahoo(symbols):
    """
    Get current yahoo quote

    Returns a DataFrame
    """
    if not isinstance(symbols, list):
        raise TypeError, "symbols must be a list"
    # for codes see: http://www.gummy-stuff.org/Yahoo-data.htm
    codes = {'symbol':'s','last':'l1','change_pct':'p2','PE':'r','time':'t1','short_ratio':'s7'}
    request = str.join('',codes.values()) # code request string
    header = codes.keys()

    data = dict(zip(codes.keys(), [[] for i in range(len(codes))]))

    urlStr = 'http://finance.yahoo.com/d/quotes.csv?s=%s&f=%s' % (str.join('+',symbols), request)

    try:
        lines = urllib2.urlopen(urlStr).readlines()
    except Exception, e:
        s = "Failed to download:\n{0}".format(e)
        print s
        return None

    for line in lines:
        fields = line.strip().split(',')
        for i, field in enumerate(fields):
            if field[-2:] == '%"':
                data[header[i]].append(float(field.strip('"%')))
            elif field[0] == '"':
                data[header[i]].append(field.strip('"'))
            else:
                try:
                    data[header[i]].append(float(field))
                except ValueError:
                    data[header[i]].append(np.nan)

    idx = data.pop('symbol')

    return DataFrame(data, index=idx)


def get_data_yahoo(name=None, start=None, end=None, retry_count=3, pause=0):
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

    for _ in range(retry_count):
        resp = urllib2.urlopen(url)
        if resp.code == 200:
            lines = resp.read()
            rs = read_csv(StringIO(bytes_to_str(lines)), index_col=0,
                          parse_dates=True)[::-1]

            # Yahoo! Finance sometimes does this awesome thing where they
            # return 2 rows for the most recent business day
            if len(rs) > 2 and rs.index[-1] == rs.index[-2]:  # pragma: no cover
                rs = rs[:-1]

            return rs

        time.sleep(pause)

    raise Exception("after %d tries, Yahoo did not "
                    "return a 200 for url %s" % (pause, url))


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
    data = read_csv(urllib.urlopen(url), index_col=0, parse_dates=True, header=None,
                    skiprows=1, names=["DATE", name])
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
    for i in range(len(file_edges) - 1):
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

cur_month = dt.datetime.now().month
cur_year = dt.datetime.now().year


def _unpack(row, kind='td'):
    return [val.text for val in row.findAll(kind)]


def _parse_options_data(table):
    rows = table.findAll('tr')
    header = _unpack(rows[0], kind='th')
    data = [_unpack(r) for r in rows[1:]]
    return TextParser(data, names=header).get_chunk()


class Options(object):
    """
    This class fetches call/put data for a given stock/exipry month.

    It is instantiated with a string representing the ticker symbol.

    The class has the following methods:
        get_options_data:(month, year)
        get_call_data:(month, year)
        get_put_data: (month, year)
        get_near_stock_price(opt_frame, above_below)
        get_forward_data(months, call, put)

    Examples
    --------
    # Instantiate object with ticker
    >>> aapl = Options('aapl')

    # Fetch September 2012 call data
    >>> calls = aapl.get_call_data(9, 2012)

    # Fetch September 2012 put data
    >>> puts = aapl.get_put_data(9, 2012)

    # cut down the call data to be 3 below and 3 above the stock price.
    >>> cut_calls = aapl.get_near_stock_price(calls, above_below=3)

    # Fetch call and put data with expiry from now to 8 months out
    >>> forward_calls, forward_puts = aapl.get_forward_data(8,
        ...                                        call=True, put=True)
    """

    def __init__(self, symbol):
        """ Instantiates options_data with a ticker saved as symbol """
        self.symbol = str(symbol).upper()

    def get_options_data(self, month=cur_month, year=cur_year):
        """
        Gets call/put data for the stock with the expiration data in the
        given month and year

        Parameters
        ----------
        month: number, int
            The month the options expire.

        year: number, int
            The year the options expire.


        Returns
        -------
        call_data: pandas.DataFrame
            A DataFrame with call options data.

        put_data: pandas.DataFrame
            A DataFrame with call options data.
        """
        from BeautifulSoup import BeautifulSoup

        mon_in = month if len(str(month)) == 2 else str('0' + str(month))

        url = str('http://finance.yahoo.com/q/op?s=' + self.symbol + '&m=' +
                  str(year) + '-' + str(mon_in))

        buf = urllib2.urlopen(url)
        soup = BeautifulSoup(buf)
        body = soup.body

        tables = body.findAll('table')
        calls = tables[9]
        puts = tables[13]

        call_data = _parse_options_data(calls)
        put_data = _parse_options_data(puts)

        return [call_data, put_data]

    def get_call_data(self, month=cur_month, year=cur_year):
        """
        Gets call/put data for the stock with the expiration data in the
        given month and year

        Parameters
        ----------
        month: number, int
            The month the options expire.

        year: number, int
            The year the options expire.

        Returns
        -------
        call_data: pandas.DataFrame
            A DataFrame with call options data.
        """
        from BeautifulSoup import BeautifulSoup

        mon_in = month if len(str(month)) == 2 else str('0' + str(month))

        url = str('http://finance.yahoo.com/q/op?s=' + self.symbol + '&m=' +
                  str(year) + '-' + str(mon_in))

        buf = urllib2.urlopen(url)
        soup = BeautifulSoup(buf)
        body = soup.body

        tables = body.findAll('table')
        calls = tables[9]

        call_data = _parse_options_data(calls)

        return call_data

    def get_put_data(self, month=cur_month, year=cur_year):
        """
        Gets put data for the stock with the expiration data in the
        given month and year

        Parameters
        ----------
        month: number, int
            The month the options expire.

        year: number, int
            The year the options expire.

        Returns
        -------
        put_data: pandas.DataFrame
            A DataFrame with call options data.
        """
        from BeautifulSoup import BeautifulSoup

        mon_in = month if len(str(month)) == 2 else str('0' + str(month))

        url = str('http://finance.yahoo.com/q/op?s=' + self.symbol + '&m=' +
                  str(year) + '-' + str(mon_in))

        buf = urllib2.urlopen(url)
        soup = BeautifulSoup(buf)
        body = soup.body

        tables = body.findAll('table')
        puts = tables[13]

        put_data = _parse_options_data(puts)

        return put_data

    def get_near_stock_price(self, opt_df, above_below=2):
        """
        Cuts the data frame opt_df that is passed in to only take
        options that are near the current stock price.

        Parameters
        ----------
        opt_df: DataFrame
            The DataFrame that will be passed in to be cut down.

        above_below: number, int, optional (default=2)
            The number of strike prices above and below the stock price that
            should be taken

        Returns
        -------
        chopped: DataFrame
            The resultant DataFrame chopped down to be 2 * above_below + 1 rows
            desired. If there isn't data as far out as the user has asked for
            then
        """
        price = get_quote_yahoo([self.symbol])['last']
        start_index = np.where(opt_df['Strike'] > price)[0][0]

        get_range = range(start_index - above_below,
                          start_index + above_below + 1)

        chopped = opt_df.ix[get_range, :]

        chopped = chopped.dropna()
        chopped = chopped.reset_index()

        return chopped

    def get_forward_data(self, months, call=True, put=False):
        """
        Gets either call, put, or both data for months starting in the current
        month and going out in the future a spcified amount of time.

        Parameters
        ----------
        months: number, int
            How many months to go out in the collection of the data. This is
            inclusive.

        call: bool, optional (default=True)
            Whether or not to collect data for call options

        put: bool, optional (default=False)
            Whether or not to collect data for put options.

        Returns
        -------
        all_calls: DataFrame
            If asked for, a DataFrame containing call data from the current
            month to the current month plus months.

        all_puts: DataFrame
            If asked for, a DataFrame containing put data from the current
            month to the current month plus months.
        """
        in_months = range(cur_month, cur_month + months + 1)
        in_years = [cur_year] * months

        # Figure out how many items in in_months go past 12
        to_change = 0
        for i in range(months):
            if in_months[i] > 12:
                in_months[i] -= 12
                to_change += 1

        # Change the corresponding items in the in_years list.
        for i in range(1, to_change + 1):
            in_years[-i] += 1

        if call:
            all_calls = DataFrame()
            for mon in range(months):
                try:  # This catches cases when there isn't data for a month
                    call_frame = self.get_call_data(in_months[mon],
                                                    in_years[mon])
                    tick = str(call_frame.ix[0, 1])
                    start = len(self.symbol)
                    year = tick[start: start + 2]
                    month = tick[start + 2: start + 4]
                    day = tick[start + 4: start + 6]
                    expiry = str(month + '-' + day + '-' + year)
                    call_frame['Expiry'] = expiry
                    if mon == 0:
                        all_calls = all_calls.join(call_frame, how='right')
                    else:
                        all_calls = concat([all_calls, call_frame])
                except:
                    pass

        if put:
            all_puts = DataFrame()
            for mon in range(months):
                try:  # This catches cases when there isn't data for a month
                    put_frame = self.get_put_data(in_months[mon],
                                                  in_years[mon])

                    # Add column with expiry data to this frame.
                    tick = str(put_frame.ix[0, 1])
                    start = len(self.symbol)
                    year = tick[start: start + 2]
                    month = tick[start + 2: start + 4]
                    day = tick[start + 4: start + 6]
                    expiry = str(month + '-' + day + '-' + year)
                    put_frame['Expiry'] = expiry

                    if mon == 0:
                        all_puts = all_puts.join(put_frame, how='right')
                    else:
                        all_puts = concat([all_puts, put_frame])
                except:
                    pass

        if call and put:
            return [all_calls, all_puts]
        else:
            if call:
                return all_calls
            else:
                return all_puts
