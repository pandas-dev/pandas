"""
Module contains tools for collecting data from various remote sources


"""
import warnings

import numpy as np
import datetime as dt
import urllib
import urllib2
import time

from zipfile import ZipFile
from pandas.util.py3compat import StringIO, BytesIO, bytes_to_str

from pandas import Panel, DataFrame, Series, read_csv, concat
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
        return get_data_yahoo(symbols=name, start=start, end=end,
                              adjust_price=False, chunk=25,
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


def _in_chunks(seq, size):
    """
    Return sequence in 'chunks' of size defined by size
    """
    return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))


def get_quote_yahoo(symbols):
    """
    Get current yahoo quote

    Returns a DataFrame
    """
    if isinstance(symbols, str):
        sym_list = symbols
    elif not isinstance(symbols, Series):
        symbols  = Series(symbols)
        sym_list = str.join('+', symbols)
    else:
        sym_list = str.join('+', symbols)

    # for codes see: http://www.gummy-stuff.org/Yahoo-data.htm
    codes = {'symbol': 's', 'last': 'l1', 'change_pct': 'p2', 'PE': 'r',
             'time': 't1', 'short_ratio': 's7'}
    request = str.join('', codes.values())  # code request string
    header = codes.keys()

    data = dict(zip(codes.keys(), [[] for i in range(len(codes))]))

    urlStr = 'http://finance.yahoo.com/d/quotes.csv?s=%s&f=%s' % (
        sym_list, request)

    try:
        lines = urllib2.urlopen(urlStr).readlines()
    except Exception, e:
        s = "Failed to download:\n{0}".format(e)
        print s
        return None

    for line in lines:
        fields = line.decode('utf-8').strip().split(',')
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


def _get_hist_yahoo(sym=None, start=None, end=None, retry_count=3,
                    pause=0, **kwargs):
    """
    Get historical data for the given name from yahoo.
    Date format is datetime

    Returns a DataFrame.
    """
    if(sym is None):
        warnings.warn("Need to provide a name.")
        return None

    start, end = _sanitize_dates(start, end)

    yahoo_URL = 'http://ichart.yahoo.com/table.csv?'

    url = yahoo_URL + 's=%s' % sym + \
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


def _adjust_prices(hist_data, price_list=['Open', 'High', 'Low', 'Close']):
    """
    Return modifed DataFrame or Panel with adjusted prices based on
    'Adj Close' price. Adds 'Adj_Ratio' column.
    """
    adj_ratio = hist_data['Adj Close'] / hist_data['Close']

    data = hist_data.copy()
    for item in price_list:
        data[item] = hist_data[item] * adj_ratio
    data['Adj_Ratio'] = adj_ratio
    del data['Adj Close']
    return data


def _calc_return_index(price_df):
    """
    Return a returns index from a input price df or series. Initial value
    (typically NaN) is set to 1.
    """
    df = price_df.pct_change().add(1).cumprod()
    mask = ~df.ix[1].isnull() & df.ix[0].isnull()
    df.ix[0][mask] = 1

    #Check for first stock listings after starting date of index in ret_index
    #If True, find first_valid_index and set previous entry to 1.
    if(~mask).any():
        for sym in mask.index[~mask]:
            tstamp = df[sym].first_valid_index()
            t_idx = df.index.get_loc(tstamp) - 1
            df[sym].ix[t_idx] = 1

    ret_index = df
    return ret_index


def get_components_yahoo(idx_sym):
    """
    Returns DataFrame containing list of component information for
    index represented in idx_sym from yahoo. Includes component symbol
    (ticker), exchange, and name.

    Parameters
    ----------
    idx_sym : str
        Stock index symbol
        Examples:
        '^DJI' (Dow Jones Industrial Average)
        '^NYA' (NYSE Composite)
        '^IXIC' (NASDAQ Composite)

        See: http://finance.yahoo.com/indices for other index symbols

    Returns
    -------
    idx_df : DataFrame
    """
    stats = 'snx'
    #URL of form:
    #http://download.finance.yahoo.com/d/quotes.csv?s=@%5EIXIC&f=snxl1d1t1c1ohgv
    url = 'http://download.finance.yahoo.com/d/quotes.csv?s={0}&f={1}' \
          '&e=.csv&h={2}'

    idx_mod = idx_sym.replace('^', '@%5E')
    urlStr = url.format(idx_mod, stats, 1)

    idx_df = DataFrame()
    mask = [True]
    comp_idx = 1

    #LOOP across component index structure,
    #break when no new components are found
    while (True in mask):
        urlStr = url.format(idx_mod, stats,  comp_idx)
        lines = (urllib.urlopen(urlStr).read().decode('utf-8').strip().
                 strip('"').split('"\r\n"'))

        lines = [line.strip().split('","') for line in lines]

        temp_df = DataFrame(lines, columns=['ticker', 'name', 'exchange'])
        temp_df = temp_df.drop_duplicates()
        temp_df = temp_df.set_index('ticker')
        mask = ~temp_df.index.isin(idx_df.index)

        comp_idx = comp_idx + 50
        idx_df = idx_df.append(temp_df[mask])

    return idx_df


def get_data_yahoo(symbols=None, start=None, end=None, retry_count=3, pause=0,
                   adjust_price=False, ret_index=False, chunksize=25,
                   **kwargs):
    """
    Returns DataFrame/Panel of historical stock prices from symbols, over date
    range, start to end. To avoid being penalized by Yahoo! Finance servers,
    pauses between downloading 'chunks' of symbols can be specified.

    Parameters
    ----------
    symbols : string, array-like object (list, tuple, Series), or DataFrame
        Single stock symbol (ticker), array-like object of symbols or
        DataFrame with index containing stock symbols.
    start : string, (defaults to '1/1/2010')
        Starting date, timestamp. Parses many different kind of date
        representations (e.g., 'JAN-01-2010', '1/1/10', 'Jan, 1, 1980')
    end : string, (defaults to today)
        Ending date, timestamp. Same format as starting date.
    retry_count : int, default 3
        Number of times to retry query request.
    pause : int, default 0
        Time, in seconds, to pause between consecutive queries of chunks. If
        single value given for symbol, represents the pause between retries.
    adjust_price : bool, default False
        If True, adjusts all prices in hist_data ('Open', 'High', 'Low', 'Close')
        based on 'Adj Close' price. Adds 'Adj_Ratio' column and drops
        'Adj Close'.
    ret_index : bool, default False
        If True, includes a simple return index 'Ret_Index' in hist_data.
    chunksize : int, default 25
        Number of symbols to download consecutively before intiating pause.

    Returns
    -------
    hist_data : DataFrame (str) or Panel (array-like object, DataFrame)
    """

    def dl_mult_symbols(symbols):
        stocks = {}
        for sym_group in _in_chunks(symbols, chunksize):
            for sym in sym_group:
                try:
                    stocks[sym] = _get_hist_yahoo(sym, start=start,
                                                  end=end, **kwargs)
                except:
                    warnings.warn('Error with sym: ' + sym + '... skipping.')

            time.sleep(pause)

        return Panel(stocks).swapaxes('items', 'minor')

    if 'name' in kwargs:
        warnings.warn("Arg 'name' is deprecated, please use 'symbols' instead.",
                      FutureWarning)
        symbols = kwargs['name']

    #If a single symbol, (e.g., 'GOOG')
    if isinstance(symbols, (str, int)):
        sym = symbols
        hist_data = _get_hist_yahoo(sym, start=start, end=end)
    #Or multiple symbols, (e.g., ['GOOG', 'AAPL', 'MSFT'])
    elif isinstance(symbols, DataFrame):
        try:
            hist_data = dl_mult_symbols(Series(symbols.index))
        except ValueError:
            raise
    else: #Guess a Series
        try:
            hist_data = dl_mult_symbols(symbols)
        except TypeError:
            hist_data = dl_mult_symbols(Series(symbols))

    if(ret_index):
        hist_data['Ret_Index'] = _calc_return_index(hist_data['Adj Close'])
    if(adjust_price):
        hist_data = _adjust_prices(hist_data)

    return hist_data


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
    data = read_csv(urllib.urlopen(url), index_col=0, parse_dates=True,
                    header=None, skiprows=1, names=["DATE", name])
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
        dataset = [d.split() for d in data[(file_edges[i] + 1):
                                           file_edges[i + 1]]]
        if(len(dataset) > 10):
            ncol = np.median(np.array([len(d) for d in dataset]))
            header_index = np.where(
                np.array([len(d) for d in dataset]) == (ncol - 1))[0][-1]
            header = dataset[header_index]
            # to ensure the header is unique
            header = [str(j + 1) + " " + header[j] for j in range(len(header))]
            index = np.array(
                [d[0] for d in dataset[(header_index + 1):]], dtype=int)
            dataset = np.array(
                [d[1:] for d in dataset[(header_index + 1):]], dtype=float)
            datasets[i] = DataFrame(dataset, index, columns=header)

    return datasets

# Items needed for options class
cur_month = dt.datetime.now().month
cur_year = dt.datetime.now().year


def _unpack(row, kind='td'):
    els = row.findall('.//%s' % kind)
    return[val.text_content() for val in els]


def _parse_options_data(table):
    rows = table.findall('.//tr')
    header = _unpack(rows[0], kind='th')
    data = [_unpack(r) for r in rows[1:]]
    # Use ',' as a thousands separator as we're pulling from the US site.
    return TextParser(data, names=header, na_values=['N/A'],
                      thousands=',').get_chunk()


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

    # Can now access aapl.calls instance variable
    >>> aapl.calls

    # Fetch September 2012 put data
    >>> puts = aapl.get_put_data(9, 2012)

    # Can now access aapl.puts instance variable
    >>> aapl.puts

    # cut down the call data to be 3 below and 3 above the stock price.
    >>> cut_calls = aapl.get_near_stock_price(calls, above_below=3)

    # Fetch call and put data with expiry from now to 8 months out
    >>> forward_calls, forward_puts = aapl.get_forward_data(8,
    ...                                              call=True, put=True)

    """

    def __init__(self, symbol):
        """ Instantiates options_data with a ticker saved as symbol """
        self.symbol = str(symbol).upper()

    def get_options_data(self, month=None, year=None):
        """
        Gets call/put data for the stock with the expiration data in the
        given month and year

        Parameters
        ----------
        month: number, int, optional(default=None)
            The month the options expire. This should be either 1 or 2
            digits.

        year: number, int, optional(default=None)
            The year the options expire. This sould be a 4 digit int.


        Returns
        -------
        call_data: pandas.DataFrame
            A DataFrame with call options data.

        put_data: pandas.DataFrame
            A DataFrame with call options data.


        Notes
        -----
        When called, this function will add instance variables named
        calls and puts. See the following example:

            >>> aapl = Options('aapl')  # Create object
            >>> aapl.calls  # will give an AttributeError
            >>> aapl.get_options_data()  # Get data and set ivars
            >>> aapl.calls  # Doesn't throw AttributeError

        Also note that aapl.calls and appl.puts will always be the calls
        and puts for the next expiry. If the user calls this method with
        a different month or year, the ivar will be named callsMMYY or
        putsMMYY where MM and YY are, repsectively, two digit
        representations of the month and year for the expiry of the
        options.
        """
        from lxml.html import parse

        if month and year:  # try to get specified month from yahoo finance
            m1 = month if len(str(month)) == 2 else '0' + str(month)
            m2 = month

            if m1 != cur_month and m2 != cur_month:  # if this month use other url
                url = str('http://finance.yahoo.com/q/op?s=' + self.symbol +
                          '&m=' + str(year) + '-' + str(m1))

            else:
                url = str('http://finance.yahoo.com/q/op?s=' + self.symbol +
                                                            '+Options')

        else:  # Default to current month
            url = str('http://finance.yahoo.com/q/op?s=' + self.symbol +
                                                            '+Options')

        parsed = parse(urllib2.urlopen(url))
        doc = parsed.getroot()
        tables = doc.findall('.//table')
        calls = tables[9]
        puts = tables[13]

        call_data = _parse_options_data(calls)
        put_data = _parse_options_data(puts)

        if month:
            c_name = 'calls' + str(m1) + str(year)[2:]
            p_name = 'puts' + str(m1) + str(year)[2:]
            self.__setattr__(c_name, call_data)
            self.__setattr__(p_name, put_data)
        else:
            self.calls = call_data
            self.calls = put_data

        return [call_data, put_data]

    def get_call_data(self, month=None, year=None):
        """
        Gets call/put data for the stock with the expiration data in the
        given month and year

        Parameters
        ----------
        month: number, int, optional(default=None)
            The month the options expire. This should be either 1 or 2
            digits.

        year: number, int, optional(default=None)
            The year the options expire. This sould be a 4 digit int.

        Returns
        -------
        call_data: pandas.DataFrame
            A DataFrame with call options data.

        Notes
        -----
        When called, this function will add instance variables named
        calls and puts. See the following example:

            >>> aapl = Options('aapl')  # Create object
            >>> aapl.calls  # will give an AttributeError
            >>> aapl.get_call_data()  # Get data and set ivars
            >>> aapl.calls  # Doesn't throw AttributeError

        Also note that aapl.calls will always be the calls for the next
        expiry. If the user calls this method with a different month
        or year, the ivar will be named callsMMYY where MM and YY are,
        repsectively, two digit representations of the month and year
        for the expiry of the options.
        """
        from lxml.html import parse

        if month and year:  # try to get specified month from yahoo finance
            m1 = month if len(str(month)) == 2 else '0' + str(month)
            m2 = month

            if m1 != cur_month and m2 != cur_month:  # if this month use other url
                url = str('http://finance.yahoo.com/q/op?s=' + self.symbol +
                          '&m=' + str(year) + '-' + str(m1))

            else:
                url = str('http://finance.yahoo.com/q/op?s=' + self.symbol +
                                                            '+Options')

        else:  # Default to current month
            url = str('http://finance.yahoo.com/q/op?s=' + self.symbol +
                                                            '+Options')

        parsed = parse(urllib2.urlopen(url))
        doc = parsed.getroot()
        tables = doc.findall('.//table')
        calls = tables[9]

        call_data = _parse_options_data(calls)

        if month:
            name = 'calls' + str(m1) + str(year)[2:]
            self.__setattr__(name, call_data)
        else:
            self.calls = call_data

        return call_data

    def get_put_data(self, month=None, year=None):
        """
        Gets put data for the stock with the expiration data in the
        given month and year

        Parameters
        ----------
        month: number, int, optional(default=None)
            The month the options expire. This should be either 1 or 2
            digits.

        year: number, int, optional(default=None)
            The year the options expire. This sould be a 4 digit int.

        Returns
        -------
        put_data: pandas.DataFrame
            A DataFrame with call options data.

        Notes
        -----
        When called, this function will add instance variables named
        puts. See the following example:

            >>> aapl = Options('aapl')  # Create object
            >>> aapl.puts  # will give an AttributeError
            >>> aapl.get_put_data()  # Get data and set ivars
            >>> aapl.puts  # Doesn't throw AttributeError

                    return self.__setattr__(self, str(str(x) + str(y)))

        Also note that aapl.puts will always be the puts for the next
        expiry. If the user calls this method with a different month
        or year, the ivar will be named putsMMYY where MM and YY are,
        repsectively, two digit representations of the month and year
        for the expiry of the options.
        """
        from lxml.html import parse

        if month and year:  # try to get specified month from yahoo finance
            m1 = month if len(str(month)) == 2 else '0' + str(month)
            m2 = month

            if m1 != cur_month and m2 != cur_month:  # if this month use other url
                url = str('http://finance.yahoo.com/q/op?s=' + self.symbol +
                          '&m=' + str(year) + '-' + str(m1))

            else:
                url = str('http://finance.yahoo.com/q/op?s=' + self.symbol +
                                                            '+Options')

        else:  # Default to current month
            url = str('http://finance.yahoo.com/q/op?s=' + self.symbol +
                                                            '+Options')

        parsed = parse(urllib2.urlopen(url))
        doc = parsed.getroot()
        tables = doc.findall('.//table')
        puts = tables[13]

        put_data = _parse_options_data(puts)

        if month:
            name = 'puts' + str(m1) + str(year)[2:]
            self.__setattr__(name, put_data)
        else:
            self.puts = put_data

        return put_data

    def get_near_stock_price(self, above_below=2, call=True, put=False,
                             month=None, year=None):
        """
        Cuts the data frame opt_df that is passed in to only take
        options that are near the current stock price.

        Parameters
        ----------
        above_below: number, int, optional (default=2)
            The number of strike prices above and below the stock price that
            should be taken

        call: bool
            Tells the function weather or not it should be using
            self.calls

        put: bool
            Tells the function weather or not it should be using
            self.puts

        month: number, int, optional(default=None)
            The month the options expire. This should be either 1 or 2
            digits.

        year: number, int, optional(default=None)
            The year the options expire. This sould be a 4 digit int.

        Returns
        -------
        chopped: DataFrame
            The resultant DataFrame chopped down to be 2 * above_below + 1 rows
            desired. If there isn't data as far out as the user has asked for
            then
        """
        price = float(get_quote_yahoo([self.symbol])['last'])

        if call:
            try:
                if month:
                    m1 = month if len(str(month)) == 2 else '0' + str(month)
                    name = 'calls' + str(m1) + str(year)[2:]
                    df_c = self.__getattribute__(name)
                else:
                    df_c = self.calls
            except AttributeError:
                df_c = self.get_call_data(month, year)

            # NOTE: For some reason the put commas in all values >1000. We remove
            #       them here
            df_c.Strike = df_c.Strike.astype(str).apply(lambda x: \
                                                        x.replace(',', ''))
            # Now make sure Strike column has dtype float
            df_c.Strike = df_c.Strike.astype(float)

            start_index = np.where(df_c['Strike'] > price)[0][0]

            get_range = range(start_index - above_below,
                              start_index + above_below + 1)

            chop_call = df_c.ix[get_range, :]

            chop_call = chop_call.dropna()
            chop_call = chop_call.reset_index()

        if put:
            try:
                if month:
                    m1 = month if len(str(month)) == 2 else '0' + str(month)
                    name = 'puts' + str(m1) + str(year)[2:]
                    df_p = self.__getattribute__(name)
                else:
                    df_p = self.puts
            except AttributeError:
                df_p = self.get_put_data(month, year)

            # NOTE: For some reason the put commas in all values >1000. We remove
            #       them here
            df_p.Strike = df_p.Strike.astype(str).apply(lambda x: \
                                                        x.replace(',', ''))
            # Now make sure Strike column has dtype float
            df_p.Strike = df_p.Strike.astype(float)

            start_index = np.where(df_p.Strike > price)[0][0]

            get_range = range(start_index - above_below,
                              start_index + above_below + 1)

            chop_put = df_p.ix[get_range, :]

            chop_put = chop_put.dropna()
            chop_put = chop_put.reset_index()

        if call and put:
            return [chop_call, chop_put]
        else:
            if call:
                return chop_call
            else:
                return chop_put

    def get_forward_data(self, months, call=True, put=False, near=False,
                         above_below=2):
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

        near: bool, optional (default=False)
            Whether this function should get only the data near the
            current stock price. Uses Options.get_near_stock_price

        above_below: number, int, optional (default=2)
            The number of strike prices above and below the stock price that
            should be taken if the near option is set to True

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
        in_years = [cur_year] * (months + 1)

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
                m2 = in_months[mon]
                y2 = in_years[mon]
                try:  # This catches cases when there isn't data for a month
                    if not near:
                        try:  # Try to access the ivar if already instantiated

                            m1 = m2 if len(str(m2)) == 2 else '0' + str(m2)
                            name = 'calls' + str(m1) + str(y2)[2:]
                            call_frame = self.__getattribute__(name)
                        except:
                            call_frame = self.get_call_data(in_months[mon],
                                                        in_years[mon])

                    else:
                        call_frame = self.get_near_stock_price(call=True,
                                                               put=False,
                                                    above_below=above_below,
                                                    month=m2, year=y2)

                    tick = str(call_frame.Symbol[0])
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
                m2 = in_months[mon]
                y2 = in_years[mon]
                try:  # This catches cases when there isn't data for a month
                    if not near:
                        try:  # Try to access the ivar if already instantiated
                            m1 = m2 if len(str(m2)) == 2 else '0' + str(m2)
                            name = 'puts' + str(m1) + str(y2)[2:]
                            put_frame = self.__getattribute__(name)
                        except:
                            put_frame = self.get_call_data(in_months[mon],
                                                        in_years[mon])

                    else:
                        put_frame = self.get_near_stock_price(call=False,
                                                              put=True,
                                                    above_below=above_below,
                                                    month=m2, year=y2)

                    # Add column with expiry data to this frame.
                    tick = str(put_frame.Symbol[0])
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
