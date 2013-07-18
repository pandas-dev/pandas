"""
Module contains tools for collecting data from various remote sources


"""
import warnings
import tempfile
import itertools
import datetime as dt
import urllib
import time

from collections import defaultdict

import numpy as np

from pandas.util.py3compat import StringIO, bytes_to_str
from pandas import Panel, DataFrame, Series, read_csv, concat
from pandas.core.common import PandasError
from pandas.io.parsers import TextParser
from pandas.io.common import urlopen, ZipFile
from pandas.util.testing import _network_error_classes


class SymbolWarning(UserWarning):
    pass


class RemoteDataError(PandasError, IOError):
    pass


def DataReader(name, data_source=None, start=None, end=None,
               retry_count=3, pause=0.001):
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

    if data_source == "yahoo":
        return get_data_yahoo(symbols=name, start=start, end=end,
                              adjust_price=False, chunksize=25,
                              retry_count=retry_count, pause=pause)
    elif data_source == "google":
        return get_data_google(symbols=name, start=start, end=end,
                               adjust_price=False, chunksize=25,
                               retry_count=retry_count, pause=pause)
    elif data_source == "fred":
        return get_data_fred(name, start, end)
    elif data_source == "famafrench":
        return get_data_famafrench(name)


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


_yahoo_codes = {'symbol': 's', 'last': 'l1', 'change_pct': 'p2', 'PE': 'r',
                'time': 't1', 'short_ratio': 's7'}

def get_quote_yahoo(symbols):
    """
    Get current yahoo quote

    Returns a DataFrame
    """
    if isinstance(symbols, basestring):
        sym_list = symbols
    else:
        sym_list = '+'.join(symbols)

    # for codes see: http://www.gummy-stuff.org/Yahoo-data.htm
    request = ''.join(_yahoo_codes.itervalues())  # code request string
    header = _yahoo_codes.keys()

    data = defaultdict(list)

    url_str = 'http://finance.yahoo.com/d/quotes.csv?s=%s&f=%s' % (sym_list,
                                                                   request)

    with urlopen(url_str) as url:
        lines = url.readlines()

    for line in lines:
        fields = line.decode('utf-8').strip().split(',')
        for i, field in enumerate(fields):
            if field[-2:] == '%"':
                v = float(field.strip('"%'))
            elif field[0] == '"':
                v = field.strip('"')
            else:
                try:
                    v = float(field)
                except ValueError:
                    v = np.nan
            data[header[i]].append(v)

    idx = data.pop('symbol')
    return DataFrame(data, index=idx)


def get_quote_google(symbols):
    raise NotImplementedError("Google Finance doesn't have this functionality")


def _retry_read_url(url, retry_count, pause, name):
    for _ in xrange(retry_count):
        time.sleep(pause)

        # kludge to close the socket ASAP
        try:
            with urlopen(url) as resp:
                lines = resp.read()
        except _network_error_classes:
            pass
        else:
            rs = read_csv(StringIO(bytes_to_str(lines)), index_col=0,
                          parse_dates=True)[::-1]
            # Yahoo! Finance sometimes does this awesome thing where they
            # return 2 rows for the most recent business day
            if len(rs) > 2 and rs.index[-1] == rs.index[-2]:  # pragma: no cover
                rs = rs[:-1]
            return rs

    raise IOError("after %d tries, %s did not "
                  "return a 200 for url %r" % (retry_count, name, url))


def _get_hist_yahoo(sym, start, end, retry_count, pause):
    """
    Get historical data for the given name from yahoo.
    Date format is datetime

    Returns a DataFrame.
    """
    start, end = _sanitize_dates(start, end)
    yahoo_url = 'http://ichart.yahoo.com/table.csv?'
    url = (yahoo_url + 's=%s' % sym +
           '&a=%s' % (start.month - 1) +
           '&b=%s' % start.day +
           '&c=%s' % start.year +
           '&d=%s' % (end.month - 1) +
           '&e=%s' % end.day +
           '&f=%s' % end.year +
           '&g=d' +
           '&ignore=.csv')
    return _retry_read_url(url, retry_count, pause, 'Yahoo!')


def _get_hist_google(sym, start, end, retry_count, pause):
    """
    Get historical data for the given name from google.
    Date format is datetime

    Returns a DataFrame.
    """
    start, end = _sanitize_dates(start, end)
    google_URL = 'http://www.google.com/finance/historical?'

    # www.google.com/finance/historical?q=GOOG&startdate=Jun+9%2C+2011&enddate=Jun+8%2C+2013&output=csv
    url = google_URL + urllib.urlencode({"q": sym,
                                         "startdate": start.strftime('%b %d, '
                                                                     '%Y'),
                                         "enddate": end.strftime('%b %d, %Y'),
                                         "output": "csv"})
    return _retry_read_url(url, retry_count, pause, 'Google')


def _adjust_prices(hist_data, price_list=None):
    """
    Return modifed DataFrame or Panel with adjusted prices based on
    'Adj Close' price. Adds 'Adj_Ratio' column.
    """
    if price_list is None:
        price_list = 'Open', 'High', 'Low', 'Close'
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
    mask = df.ix[1].notnull() & df.ix[0].isnull()
    df.ix[0][mask] = 1

    # Check for first stock listings after starting date of index in ret_index
    # If True, find first_valid_index and set previous entry to 1.
    if (~mask).any():
        for sym in mask.index[~mask]:
            tstamp = df[sym].first_valid_index()
            t_idx = df.index.get_loc(tstamp) - 1
            df[sym].ix[t_idx] = 1

    return df


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
    # URL of form:
    # http://download.finance.yahoo.com/d/quotes.csv?s=@%5EIXIC&f=snxl1d1t1c1ohgv
    url = ('http://download.finance.yahoo.com/d/quotes.csv?s={0}&f={1}'
           '&e=.csv&h={2}')

    idx_mod = idx_sym.replace('^', '@%5E')
    url_str = url.format(idx_mod, stats, 1)

    idx_df = DataFrame()
    mask = [True]
    comp_idx = 1

    # LOOP across component index structure,
    # break when no new components are found
    while True in mask:
        url_str = url.format(idx_mod, stats,  comp_idx)
        with urlopen(url_str) as resp:
            raw = resp.read()
        lines = raw.decode('utf-8').strip().strip('"').split('"\r\n"')
        lines = [line.strip().split('","') for line in lines]

        temp_df = DataFrame(lines, columns=['ticker', 'name', 'exchange'])
        temp_df = temp_df.drop_duplicates()
        temp_df = temp_df.set_index('ticker')
        mask = ~temp_df.index.isin(idx_df.index)

        comp_idx = comp_idx + 50
        idx_df = idx_df.append(temp_df[mask])

    return idx_df


def _dl_mult_symbols(symbols, start, end, chunksize, retry_count, pause,
                     method):
    stocks = {}
    for sym_group in _in_chunks(symbols, chunksize):
        for sym in sym_group:
            try:
                stocks[sym] = method(sym, start, end, retry_count, pause)
            except IOError:
                warnings.warn('Failed to read symbol: {0!r}, replacing with '
                              'NaN.'.format(sym), SymbolWarning)
                stocks[sym] = np.nan

    try:
        return Panel(stocks).swapaxes('items', 'minor')
    except AttributeError:
        # cannot construct a panel with just 1D nans indicating no data
        raise RemoteDataError("No data fetched using "
                              "{0!r}".format(method.__name__))


_source_functions = {'google': _get_hist_google, 'yahoo': _get_hist_yahoo}

def _get_data_from(symbols, start, end, retry_count, pause, adjust_price,
                   ret_index, chunksize, source, name):
    if name is not None:
        warnings.warn("Arg 'name' is deprecated, please use 'symbols' "
                      "instead.", FutureWarning)
        symbols = name

    src_fn = _source_functions[source]

    # If a single symbol, (e.g., 'GOOG')
    if isinstance(symbols, (basestring, int)):
        hist_data = src_fn(symbols, start, end, retry_count, pause)
    # Or multiple symbols, (e.g., ['GOOG', 'AAPL', 'MSFT'])
    elif isinstance(symbols, DataFrame):
        hist_data = _dl_mult_symbols(symbols.index, start, end, chunksize,
                                     retry_count, pause, src_fn)
    else:
        hist_data = _dl_mult_symbols(symbols, start, end, chunksize,
                                     retry_count, pause, src_fn)
    if source.lower() == 'yahoo':
        if ret_index:
            hist_data['Ret_Index'] = _calc_return_index(hist_data['Adj Close'])
        if adjust_price:
            hist_data = _adjust_prices(hist_data)

    return hist_data


def get_data_yahoo(symbols=None, start=None, end=None, retry_count=3,
                   pause=0.001, adjust_price=False, ret_index=False,
                   chunksize=25, name=None):
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
        If True, adjusts all prices in hist_data ('Open', 'High', 'Low',
        'Close') based on 'Adj Close' price. Adds 'Adj_Ratio' column and drops
        'Adj Close'.
    ret_index : bool, default False
        If True, includes a simple return index 'Ret_Index' in hist_data.
    chunksize : int, default 25
        Number of symbols to download consecutively before intiating pause.

    Returns
    -------
    hist_data : DataFrame (str) or Panel (array-like object, DataFrame)
    """
    return _get_data_from(symbols, start, end, retry_count, pause,
                          adjust_price, ret_index, chunksize, 'yahoo', name)


def get_data_google(symbols=None, start=None, end=None, retry_count=3,
                    pause=0.001, adjust_price=False, ret_index=False,
                    chunksize=25, name=None):
    """
    Returns DataFrame/Panel of historical stock prices from symbols, over date
    range, start to end. To avoid being penalized by Google Finance servers,
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
    chunksize : int, default 25
        Number of symbols to download consecutively before intiating pause.

    Returns
    -------
    hist_data : DataFrame (str) or Panel (array-like object, DataFrame)
    """
    return _get_data_from(symbols, start, end, retry_count, pause,
                          adjust_price, ret_index, chunksize, 'google', name)


def get_data_fred(name, start=dt.datetime(2010, 1, 1),
                  end=dt.datetime.today()):
    """
    Get data for the given name from the St. Louis FED (FRED).
    Date format is datetime

    Returns a DataFrame.
    """
    start, end = _sanitize_dates(start, end)

    fred_URL = "http://research.stlouisfed.org/fred2/series/"

    url = fred_URL + '%s' % name + '/downloaddata/%s' % name + '.csv'
    with urlopen(url) as resp:
        data = read_csv(resp, index_col=0, parse_dates=True,
                        header=None, skiprows=1, names=["DATE", name],
                        na_values='.')
    try:
        return data.truncate(start, end)
    except KeyError:
        if data.ix[3].name[7:12] == 'Error':
            raise IOError("Failed to get the data. Check that {0!r} is "
                          "a valid FRED series.".format(name))
        raise


def get_data_famafrench(name):
    # path of zip files
    zip_file_url = ('http://mba.tuck.dartmouth.edu/pages/faculty/'
                    'ken.french/ftp')
    zip_file_path = '{0}/{1}.zip'.format(zip_file_url, name)

    with urlopen(zip_file_path) as url:
        raw = url.read()

    with tempfile.TemporaryFile() as tmpf:
        tmpf.write(raw)

        with ZipFile(tmpf, 'r') as zf:
            data = zf.open(name + '.txt').readlines()

    line_lengths = np.array(map(len, data))
    file_edges = np.where(line_lengths == 2)[0]

    datasets = {}
    edges = itertools.izip(file_edges + 1, file_edges[1:])
    for i, (left_edge, right_edge) in enumerate(edges):
        dataset = [d.split() for d in data[left_edge:right_edge]]
        if len(dataset) > 10:
            ncol_raw = np.array(map(len, dataset))
            ncol = np.median(ncol_raw)
            header_index = np.where(ncol_raw == ncol - 1)[0][-1]
            header = dataset[header_index]
            ds_header = dataset[header_index + 1:]
            # to ensure the header is unique
            header = ['{0} {1}'.format(j, hj) for j, hj in enumerate(header,
                                                                     start=1)]
            index = np.array([d[0] for d in ds_header], dtype=int)
            dataset = np.array([d[1:] for d in ds_header], dtype=float)
            datasets[i] = DataFrame(dataset, index, columns=header)

    return datasets


# Items needed for options class
CUR_MONTH = dt.datetime.now().month
CUR_YEAR = dt.datetime.now().year


def _unpack(row, kind):
    els = row.xpath('.//%s' % kind)
    return [val.text_content() for val in els]


def _parse_options_data(table):
    rows = table.xpath('.//tr')
    header = _unpack(rows[0], kind='th')
    data = [_unpack(row, kind='td') for row in rows[1:]]
    # Use ',' as a thousands separator as we're pulling from the US site.
    return TextParser(data, names=header, na_values=['N/A'],
                      thousands=',').get_chunk()


def _two_char_month(s):
    return '{0:0>2}'.format(s)


class Options(object):
    """
    This class fetches call/put data for a given stock/expiry month.

    It is instantiated with a string representing the ticker symbol.

    The class has the following methods:
        get_options:(month, year)
        get_calls:(month, year)
        get_puts: (month, year)
        get_near_stock_price(opt_frame, above_below)
        get_forward(months, call, put)

    Examples
    --------
    # Instantiate object with ticker
    >>> aapl = Options('aapl', 'yahoo')

    # Fetch September 2012 call data
    >>> calls = aapl.get_calls(9, 2012)

    # Can now access aapl.calls instance variable
    >>> aapl.calls

    # Fetch September 2012 put data
    >>> puts = aapl.get_puts(9, 2012)

    # Can now access aapl.puts instance variable
    >>> aapl.puts

    # cut down the call data to be 3 below and 3 above the stock price.
    >>> cut_calls = aapl.get_near_stock_price(calls, above_below=3)

    # Fetch call and put data with expiry from now to 8 months out
    >>> forward_calls, forward_puts = aapl.get_forward_data(8,
    ...                                              call=True, put=True)

    """
    def __init__(self, symbol, data_source=None):
        """ Instantiates options_data with a ticker saved as symbol """
        self.symbol = symbol.upper()
        if data_source is None:
            warnings.warn("Options(symbol) is deprecated, use Options(symbol,"
                          " data_source) instead", FutureWarning)
            data_source = "yahoo"
        if data_source != "yahoo":
            raise NotImplementedError("currently only yahoo supported")

    def get_options_data(self, month=None, year=None, expiry=None):
        """
        Gets call/put data for the stock with the expiration data in the
        given month and year

        Parameters
        ----------
        expiry: datetime.date, optional(default=None)
            The date when options expire (defaults to current month)

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

            >>> aapl = Options('aapl', 'yahoo')  # Create object
            >>> aapl.calls  # will give an AttributeError
            >>> aapl.get_options()  # Get data and set ivars
            >>> aapl.calls  # Doesn't throw AttributeError

        Also note that aapl.calls and appl.puts will always be the calls
        and puts for the next expiry. If the user calls this method with
        a different month or year, the ivar will be named callsMMYY or
        putsMMYY where MM and YY are, repsectively, two digit
        representations of the month and year for the expiry of the
        options.
        """
        return [f(month, year, expiry) for f in (self.get_put_data,
                                                 self.get_call_data)]

    def _get_option_data(self, month, year, expiry, table_loc, name):
        year, month = self._try_parse_dates(year, month, expiry)

        url = 'http://finance.yahoo.com/q/op?s={sym}'.format(sym=self.symbol)

        if month and year:  # try to get specified month from yahoo finance
            m1, m2 = _two_char_month(month), month

            # if this month use other url
            if m1 != CUR_MONTH and m2 != CUR_MONTH:
                url += '&m={year}-{m1}'.format(year=year, m1=m1)
            else:
                url += '+Options'
        else:  # Default to current month
            url += '+Options'

        try:
            from lxml.html import parse
        except ImportError:
            raise ImportError("Please install lxml if you want to use the "
                              "{0!r} class".format(self.__class__.__name__))
        try:
            doc = parse(url)
        except _network_error_classes:
            raise RemoteDataError("Unable to parse tables from URL "
                                  "{0!r}".format(url))
        else:
            root = doc.getroot()
            if root is None:
                raise RemoteDataError("Parsed URL {0!r} has no root"
                                      "element".format(url))
            tables = root.xpath('.//table')
            ntables = len(tables)
            if table_loc - 1 > ntables:
                raise IndexError("Table location {0} invalid, {1} tables"
                                 " found".format(table_loc, ntables))

        option_data = _parse_options_data(tables[table_loc])

        if month:
            name += m1 + str(year)[-2:]
        setattr(self, name, option_data)
        return option_data

    def get_call_data(self, month=None, year=None, expiry=None):
        """
        Gets call/put data for the stock with the expiration data in the
        given month and year

        Parameters
        ----------
        expiry: datetime.date, optional(default=None)
            The date when options expire (defaults to current month)

        Returns
        -------
        call_data: pandas.DataFrame
            A DataFrame with call options data.

        Notes
        -----
        When called, this function will add instance variables named
        calls and puts. See the following example:

            >>> aapl = Options('aapl', 'yahoo')  # Create object
            >>> aapl.calls  # will give an AttributeError
            >>> aapl.get_call_data()  # Get data and set ivars
            >>> aapl.calls  # Doesn't throw AttributeError

        Also note that aapl.calls will always be the calls for the next
        expiry. If the user calls this method with a different month
        or year, the ivar will be named callsMMYY where MM and YY are,
        repsectively, two digit representations of the month and year
        for the expiry of the options.
        """
        return self._get_option_data(month, year, expiry, 9, 'calls')

    def get_put_data(self, month=None, year=None, expiry=None):
        """
        Gets put data for the stock with the expiration data in the
        given month and year

        Parameters
        ----------
        expiry: datetime.date, optional(default=None)
            The date when options expire (defaults to current month)

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
        return self._get_option_data(month, year, expiry, 13, 'puts')

    def get_near_stock_price(self, above_below=2, call=True, put=False,
                             month=None, year=None, expiry=None):
        """
        Cuts the data frame opt_df that is passed in to only take
        options that are near the current stock price.

        Parameters
        ----------
        above_below: number, int, optional (default=2)
            The number of strike prices above and below the stock price that
            should be taken

        call: bool
            Tells the function whether or not it should be using
            self.calls

        put: bool
            Tells the function weather or not it should be using
            self.puts

        expiry: datetime.date, optional(default=None)
            The date when options expire (defaults to current month)

        Returns
        -------
        chopped: DataFrame
            The resultant DataFrame chopped down to be 2 * above_below + 1 rows
            desired. If there isn't data as far out as the user has asked for
            then
        """
        year, month = self._try_parse_dates(year, month, expiry)
        price = float(get_quote_yahoo([self.symbol])['last'])

        to_ret = Series({'calls': call, 'puts': put})
        to_ret = to_ret[to_ret].index

        data = {}

        for nam in to_ret:
            if month:
                m1 = _two_char_month(month)
                name = nam + m1 + str(year)[2:]

            try:
                df = getattr(self, name)
            except AttributeError:
                meth_name = 'get_{0}_data'.format(nam[:-1])
                df = getattr(self, meth_name)(month, year)

            start_index = np.where(df['Strike'] > price)[0][0]

            get_range = slice(start_index - above_below,
                              start_index + above_below + 1)
            chop = df[get_range].dropna(how='all')
            chop.reset_index(inplace=True)
            data[nam] = chop
        return [data[nam] for nam in to_ret]

    def _try_parse_dates(self, year, month, expiry):
        if year is not None or month is not None:
            warnings.warn("month, year arguments are deprecated, use expiry"
                          " instead", FutureWarning)

        if expiry is not None:
            year = expiry.year
            month = expiry.month
        return year, month

    def get_forward_data(self, months, call=True, put=False, near=False,
                         above_below=2):
        """
        Gets either call, put, or both data for months starting in the current
        month and going out in the future a specified amount of time.

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
        data : dict of str, DataFrame
        """
        warnings.warn("get_forward_data() is deprecated", FutureWarning)
        in_months = xrange(CUR_MONTH, CUR_MONTH + months + 1)
        in_years = [CUR_YEAR] * (months + 1)

        # Figure out how many items in in_months go past 12
        to_change = 0
        for i in xrange(months):
            if in_months[i] > 12:
                in_months[i] -= 12
                to_change += 1

        # Change the corresponding items in the in_years list.
        for i in xrange(1, to_change + 1):
            in_years[-i] += 1

        to_ret = Series({'calls': call, 'puts': put})
        to_ret = to_ret[to_ret].index
        data = {}

        for name in to_ret:
            all_data = DataFrame()

            for mon in xrange(months):
                m2 = in_months[mon]
                y2 = in_years[mon]

                if not near:
                    m1 = _two_char_month(m2)
                    nam = name + str(m1) + str(y2)[2:]

                    try:  # Try to access on the instance
                        frame = getattr(self, nam)
                    except AttributeError:
                        meth_name = 'get_{0}_data'.format(name[:-1])
                        frame = getattr(self, meth_name)(m2, y2)
                else:
                    frame = self.get_near_stock_price(call=call, put=put,
                                                      above_below=above_below,
                                                      month=m2, year=y2)
                tick = str(frame.Symbol[0])
                start = len(self.symbol)
                year = tick[start:start + 2]
                month = tick[start + 2:start + 4]
                day = tick[start + 4:start + 6]
                expiry = month + '-' + day + '-' + year
                frame['Expiry'] = expiry

                if not mon:
                    all_data = all_data.join(frame, how='right')
                else:
                    all_data = concat([all_data, frame])
            data[name] = all_data
        ret = [data[k] for k in to_ret]
        if len(ret) == 1:
            return ret.pop()
        if len(ret) != 2:
            raise AssertionError("should be len 2")
        return ret
