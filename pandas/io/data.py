"""
Module contains tools for collecting data from various remote sources


"""
import warnings
import tempfile
import datetime as dt
import time

from collections import defaultdict

import numpy as np

from pandas.compat import(
    StringIO, bytes_to_str, range, lmap, zip
)
import pandas.compat as compat
from pandas import Panel, DataFrame, Series, read_csv, concat, to_datetime, DatetimeIndex, DateOffset
from pandas.core.common import is_list_like, PandasError
from pandas.io.common import urlopen, ZipFile, urlencode
from pandas.tseries.offsets import MonthEnd
from pandas.util.testing import _network_error_classes
from pandas.io.html import read_html

warnings.warn("\n"
              "The pandas.io.data module is moved to a separate package "
              "(pandas-datareader) and will be removed from pandas in a "
              "future version.\nAfter installing the pandas-datareader package "
              "(https://github.com/pydata/pandas-datareader), you can change "
              "the import ``from pandas.io import data, wb`` to "
              "``from pandas_datareader import data, wb``.",
              FutureWarning)

class SymbolWarning(UserWarning):
    pass


class RemoteDataError(PandasError, IOError):
    pass


def DataReader(name, data_source=None, start=None, end=None,
               retry_count=3, pause=0.001):
    """
    Imports data from a number of online sources.

    Currently supports Yahoo! Finance, Google Finance, St. Louis FED (FRED)
    and Kenneth French's data library.

    Parameters
    ----------
    name : str or list of strs
        the name of the dataset. Some data sources (yahoo, google, fred) will
        accept a list of names.
    data_source: str, default: None
        the data source ("yahoo", "google", "fred", or "ff")
    start : datetime, default: None
        left boundary for range (defaults to 1/1/2010)
    end : datetime, default: None
        right boundary for range (defaults to today)
    retry_count : int, default 3
        Number of times to retry query request.
    pause : numeric, default 0.001
        Time, in seconds, to pause between consecutive queries of chunks. If
        single value given for symbol, represents the pause between retries.

    Examples
    ----------

    # Data from Yahoo! Finance
    gs = DataReader("GS", "yahoo")

    # Data from Google Finance
    aapl = DataReader("AAPL", "google")

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
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))


_yahoo_codes = {'symbol': 's', 'last': 'l1', 'change_pct': 'p2', 'PE': 'r',
                'time': 't1', 'short_ratio': 's7'}


_YAHOO_QUOTE_URL = 'http://finance.yahoo.com/d/quotes.csv?'


def get_quote_yahoo(symbols):
    """
    Get current yahoo quote

    Returns a DataFrame
    """
    if isinstance(symbols, compat.string_types):
        sym_list = symbols
    else:
        sym_list = '+'.join(symbols)

    # for codes see: http://www.gummy-stuff.org/Yahoo-data.htm
    request = ''.join(compat.itervalues(_yahoo_codes))  # code request string
    header = list(_yahoo_codes.keys())

    data = defaultdict(list)

    url_str = _YAHOO_QUOTE_URL + 's=%s&f=%s' % (sym_list, request)

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
                    v = field
            data[header[i]].append(v)

    idx = data.pop('symbol')
    return DataFrame(data, index=idx)


def get_quote_google(symbols):
    raise NotImplementedError("Google Finance doesn't have this functionality")


def _retry_read_url(url, retry_count, pause, name):
    for _ in range(retry_count):
        time.sleep(pause)

        # kludge to close the socket ASAP
        try:
            with urlopen(url) as resp:
                lines = resp.read()
        except _network_error_classes:
            pass
        else:
            rs = read_csv(StringIO(bytes_to_str(lines)), index_col=0,
                          parse_dates=True, na_values='-')[::-1]
            # Yahoo! Finance sometimes does this awesome thing where they
            # return 2 rows for the most recent business day
            if len(rs) > 2 and rs.index[-1] == rs.index[-2]:  # pragma: no cover
                rs = rs[:-1]

            #Get rid of unicode characters in index name.
            try:
                rs.index.name = rs.index.name.decode('unicode_escape').encode('ascii', 'ignore')
            except AttributeError:
                #Python 3 string has no decode method.
                rs.index.name = rs.index.name.encode('ascii', 'ignore').decode()

            return rs

    raise IOError("after %d tries, %s did not "
                  "return a 200 for url %r" % (retry_count, name, url))


_HISTORICAL_YAHOO_URL = 'http://ichart.finance.yahoo.com/table.csv?'


def _get_hist_yahoo(sym, start, end, interval, retry_count, pause):
    """
    Get historical data for the given name from yahoo.
    Date format is datetime

    Returns a DataFrame.
    """
    start, end = _sanitize_dates(start, end)
    url = (_HISTORICAL_YAHOO_URL + 's=%s' % sym +
           '&a=%s' % (start.month - 1) +
           '&b=%s' % start.day +
           '&c=%s' % start.year +
           '&d=%s' % (end.month - 1) +
           '&e=%s' % end.day +
           '&f=%s' % end.year +
           '&g=%s' % interval +
           '&ignore=.csv')
    return _retry_read_url(url, retry_count, pause, 'Yahoo!')


_HISTORICAL_GOOGLE_URL = 'http://www.google.com/finance/historical?'


def _get_hist_google(sym, start, end, interval, retry_count, pause):
    """
    Get historical data for the given name from google.
    Date format is datetime

    Returns a DataFrame.
    """
    start, end = _sanitize_dates(start, end)

    # www.google.com/finance/historical?q=GOOG&startdate=Jun+9%2C+2011&enddate=Jun+8%2C+2013&output=csv
    url = "%s%s" % (_HISTORICAL_GOOGLE_URL,
                    urlencode({"q": sym,
                               "startdate": start.strftime('%b %d, ' '%Y'),
                               "enddate": end.strftime('%b %d, %Y'),
                               "output": "csv"}))
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


_YAHOO_COMPONENTS_URL = 'http://download.finance.yahoo.com/d/quotes.csv?'


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
    url = _YAHOO_COMPONENTS_URL + 's={0}&f={1}&e=.csv&h={2}'

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


def _dl_mult_symbols(symbols, start, end, interval, chunksize, retry_count, pause,
                     method):
    stocks = {}
    failed = []
    passed = []
    for sym_group in _in_chunks(symbols, chunksize):
        for sym in sym_group:
            try:
                stocks[sym] = method(sym, start, end, interval, retry_count, pause)
                passed.append(sym)
            except IOError:
                warnings.warn('Failed to read symbol: {0!r}, replacing with '
                              'NaN.'.format(sym), SymbolWarning)
                failed.append(sym)

    if len(passed) == 0:
        raise RemoteDataError("No data fetched using "
                              "{0!r}".format(method.__name__))
    try:
        if len(stocks) > 0 and len(failed) > 0 and len(passed) > 0:
            df_na = stocks[passed[0]].copy()
            df_na[:] = np.nan
            for sym in failed:
                stocks[sym] = df_na
        return Panel(stocks).swapaxes('items', 'minor')
    except AttributeError:
        # cannot construct a panel with just 1D nans indicating no data
        raise RemoteDataError("No data fetched using "
                              "{0!r}".format(method.__name__))

_source_functions = {'google': _get_hist_google, 'yahoo': _get_hist_yahoo}


def _get_data_from(symbols, start, end, interval, retry_count, pause, adjust_price,
                   ret_index, chunksize, source):

    src_fn = _source_functions[source]

    # If a single symbol, (e.g., 'GOOG')
    if isinstance(symbols, (compat.string_types, int)):
        hist_data = src_fn(symbols, start, end, interval, retry_count, pause)
    # Or multiple symbols, (e.g., ['GOOG', 'AAPL', 'MSFT'])
    elif isinstance(symbols, DataFrame):
        hist_data = _dl_mult_symbols(symbols.index, start, end, interval, chunksize,
                                     retry_count, pause, src_fn)
    else:
        hist_data = _dl_mult_symbols(symbols, start, end, interval, chunksize,
                                     retry_count, pause, src_fn)
    if source.lower() == 'yahoo':
        if ret_index:
            hist_data['Ret_Index'] = _calc_return_index(hist_data['Adj Close'])
        if adjust_price:
            hist_data = _adjust_prices(hist_data)

    return hist_data


def get_data_yahoo(symbols=None, start=None, end=None, retry_count=3,
                   pause=0.001, adjust_price=False, ret_index=False,
                   chunksize=25, interval='d'):
    """
    Returns DataFrame/Panel of historical stock prices from symbols, over date
    range, start to end. To avoid being penalized by Yahoo! Finance servers,
    pauses between downloading 'chunks' of symbols can be specified.

    Parameters
    ----------
    symbols : string, array-like object (list, tuple, Series), or DataFrame, default: None
        Single stock symbol (ticker), array-like object of symbols or
        DataFrame with index containing stock symbols
    start : string, (defaults to '1/1/2010')
        Starting date, timestamp. Parses many different kind of date
        representations (e.g., 'JAN-01-2010', '1/1/10', 'Jan, 1, 1980')
    end : string, (defaults to today)
        Ending date, timestamp. Same format as starting date.
    retry_count : int, default: 3
        Number of times to retry query request.
    pause : numeric, default: 0.001
        Time, in seconds, to pause between consecutive queries of chunks. If
        single value given for symbol, represents the pause between retries.
    adjust_price : bool, default: False
        If True, adjusts all prices in hist_data ('Open', 'High', 'Low',
        'Close') based on 'Adj Close' price. Adds 'Adj_Ratio' column and drops
        'Adj Close'.
    ret_index : bool, default: False
        If True, includes a simple return index 'Ret_Index' in hist_data.
    chunksize : int, default: 25
        Number of symbols to download consecutively before intiating pause.
    interval : string, default: 'd'
        Time interval code, valid values are 'd' for daily, 'w' for weekly,
        'm' for monthly and 'v' for dividend.

    Returns
    -------
    hist_data : DataFrame (str) or Panel (array-like object, DataFrame)
    """
    if interval not in ['d', 'w', 'm', 'v']:
        raise ValueError("Invalid interval: valid values are 'd', 'w', 'm' and 'v'")
    return _get_data_from(symbols, start, end, interval, retry_count, pause,
                          adjust_price, ret_index, chunksize, 'yahoo')


def get_data_google(symbols=None, start=None, end=None, retry_count=3,
                    pause=0.001, adjust_price=False, ret_index=False,
                    chunksize=25):
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
    retry_count : int, default: 3
        Number of times to retry query request.
    pause : numeric, default: 0.001
        Time, in seconds, to pause between consecutive queries of chunks. If
        single value given for symbol, represents the pause between retries.
    chunksize : int, default: 25
        Number of symbols to download consecutively before intiating pause.
    ret_index : bool, default: False
        If True, includes a simple return index 'Ret_Index' in hist_data.

    Returns
    -------
    hist_data : DataFrame (str) or Panel (array-like object, DataFrame)
    """
    return _get_data_from(symbols, start, end, None, retry_count, pause,
                          adjust_price, ret_index, chunksize, 'google')


_FRED_URL = "http://research.stlouisfed.org/fred2/series/"


def get_data_fred(name, start=dt.datetime(2010, 1, 1),
                  end=dt.datetime.today()):
    """
    Get data for the given name from the St. Louis FED (FRED).
    Date format is datetime

    Returns a DataFrame.

    If multiple names are passed for "series" then the index of the
    DataFrame is the outer join of the indicies of each series.
    """
    start, end = _sanitize_dates(start, end)

    if not is_list_like(name):
        names = [name]
    else:
        names = name

    urls = [_FRED_URL + '%s' % n + '/downloaddata/%s' % n + '.csv' for
            n in names]

    def fetch_data(url, name):
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
    df = concat([fetch_data(url, n) for url, n in zip(urls, names)],
                axis=1, join='outer')
    return df


_FAMAFRENCH_URL = 'http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp'


def get_data_famafrench(name):
    # path of zip files
    zip_file_path = '{0}/{1}_TXT.zip'.format(_FAMAFRENCH_URL, name)

    with urlopen(zip_file_path) as url:
        raw = url.read()

    with tempfile.TemporaryFile() as tmpf:
        tmpf.write(raw)

        with ZipFile(tmpf, 'r') as zf:
            data = zf.open(zf.namelist()[0]).readlines()

    line_lengths = np.array(lmap(len, data))
    file_edges = np.where(line_lengths == 2)[0]

    datasets = {}
    edges = zip(file_edges + 1, file_edges[1:])
    for i, (left_edge, right_edge) in enumerate(edges):
        dataset = [d.split() for d in data[left_edge:right_edge]]
        if len(dataset) > 10:
            ncol_raw = np.array(lmap(len, dataset))
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
CUR_DAY = dt.datetime.now().day


def _two_char(s):
    return '{0:0>2}'.format(s)


class Options(object):
    """
    ***Experimental***
    This class fetches call/put data for a given stock/expiry month.

    It is instantiated with a string representing the ticker symbol.

    The class has the following methods:
        get_options_data:(month, year, expiry)
        get_call_data:(month, year, expiry)
        get_put_data: (month, year, expiry)
        get_near_stock_price(opt_frame, above_below)
        get_all_data(call, put)
        get_forward_data(months, call, put) (deprecated)

    Examples
    --------
    # Instantiate object with ticker
    >>> aapl = Options('aapl', 'yahoo')

    # Fetch next expiry call data
    >>> calls = aapl.get_call_data()

    # Can now access aapl.calls instance variable
    >>> aapl.calls

    # Fetch next expiry put data
    >>> puts = aapl.get_put_data()

    # Can now access aapl.puts instance variable
    >>> aapl.puts

    # cut down the call data to be 3 below and 3 above the stock price.
    >>> cut_calls = aapl.get_near_stock_price(call=True, above_below=3)

    # Fetch call and put data with expiry from now to 8 months out
    >>> forward_data = aapl.get_forward_data(8, call=True, put=True)

    # Fetch all call and put data
    >>> all_data = aapl.get_all_data()
    """

    _TABLE_LOC = {'calls': 1, 'puts': 2}
    _OPTIONS_BASE_URL = 'http://finance.yahoo.com/q/op?s={sym}'
    _FINANCE_BASE_URL = 'http://finance.yahoo.com'

    def __init__(self, symbol, data_source=None):
        """ Instantiates options_data with a ticker saved as symbol """
        self.symbol = symbol.upper()
        if data_source is None:
            warnings.warn("Options(symbol) is deprecated, use Options(symbol,"
                          " data_source) instead", FutureWarning, stacklevel=2)
            data_source = "yahoo"
        if data_source != "yahoo":
            raise NotImplementedError("currently only yahoo supported")

    def get_options_data(self, month=None, year=None, expiry=None):
        """
        ***Experimental***
        Gets call/put data for the stock with the expiration data in the
        given month and year

        Parameters
        ----------
        month : number, int, optional(default=None)
            The month the options expire. This should be either 1 or 2
            digits.

        year : number, int, optional(default=None)
            The year the options expire. This should be a 4 digit int.

        expiry : date-like or convertible or list-like object, optional (default=None)
            The date (or dates) when options expire (defaults to current month)

        Returns
        -------
        pandas.DataFrame
            A DataFrame with requested options data.

            Index:
                Strike: Option strike, int
                Expiry: Option expiry, Timestamp
                Type: Call or Put, string
                Symbol: Option symbol as reported on Yahoo, string
            Columns:
                Last: Last option price, float
                Chg: Change from prior day, float
                Bid: Bid price, float
                Ask: Ask price, float
                Vol: Volume traded, int64
                Open_Int: Open interest, int64
                IsNonstandard: True if the the deliverable is not 100 shares, otherwise false
                Underlying: Ticker of the underlying security, string
                Underlying_Price: Price of the underlying security, float64
                Quote_Time: Time of the quote, Timestamp

        Notes
        -----
        Note: Format of returned data frame is dependent on Yahoo and may change.

        When called, this function will add instance variables named
        calls and puts. See the following example:

            >>> aapl = Options('aapl', 'yahoo')  # Create object
            >>> aapl.calls  # will give an AttributeError
            >>> aapl.get_options()  # Get data and set ivars
            >>> aapl.calls  # Doesn't throw AttributeError

        Also note that aapl.calls and appl.puts will always be the calls
        and puts for the next expiry. If the user calls this method with
        a different expiry, the ivar will be named callsYYMMDD or putsYYMMDD,
        where YY, MM and DD are, respectively, two digit representations of
        the year, month and day for the expiry of the options.

        """
        return concat([f(month, year, expiry)
                       for f in (self.get_put_data,
                                 self.get_call_data)]).sortlevel()

    def _get_option_frames_from_yahoo(self, expiry):
        url = self._yahoo_url_from_expiry(expiry)
        option_frames = self._option_frames_from_url(url)
        frame_name = '_frames' + self._expiry_to_string(expiry)
        setattr(self, frame_name, option_frames)
        return option_frames

    @staticmethod
    def _expiry_to_string(expiry):
        m1 = _two_char(expiry.month)
        d1 = _two_char(expiry.day)
        return str(expiry.year)[-2:] + m1 + d1

    def _yahoo_url_from_expiry(self, expiry):
        try:
            expiry_links = self._expiry_links

        except AttributeError:
            _, expiry_links = self._get_expiry_dates_and_links()

        return self._FINANCE_BASE_URL + expiry_links[expiry]

    def _option_frames_from_url(self, url):
        frames = read_html(url)
        nframes = len(frames)
        frames_req = max(self._TABLE_LOC.values())
        if nframes < frames_req:
            raise RemoteDataError("%s options tables found (%s expected)" % (nframes, frames_req))

        if not hasattr(self, 'underlying_price'):
            try:
                self.underlying_price, self.quote_time = self._underlying_price_and_time_from_url(url)
            except IndexError:
                self.underlying_price, self.quote_time = np.nan, np.nan

        calls = frames[self._TABLE_LOC['calls']]
        puts = frames[self._TABLE_LOC['puts']]

        calls = self._process_data(calls, 'call')
        puts = self._process_data(puts, 'put')

        return {'calls': calls, 'puts': puts}

    def _underlying_price_and_time_from_url(self, url):
        root = self._parse_url(url)
        underlying_price = self._underlying_price_from_root(root)
        quote_time = self._quote_time_from_root(root)
        return underlying_price, quote_time

    @staticmethod
    def _underlying_price_from_root(root):
        underlying_price = root.xpath('.//*[@class="time_rtq_ticker Fz-30 Fw-b"]')[0]\
            .getchildren()[0].text
        underlying_price = underlying_price.replace(',', '') #GH11

        try:
            underlying_price = float(underlying_price)
        except ValueError:
            underlying_price = np.nan

        return underlying_price

    @staticmethod
    def _quote_time_from_root(root):
        #Gets the time of the quote, note this is actually the time of the underlying price.
        try:
            quote_time_text = root.xpath('.//*[@class="time_rtq Fz-m"]')[0].getchildren()[1].getchildren()[0].text
            ##TODO: Enable timezone matching when strptime can match EST with %Z
            quote_time_text = quote_time_text.split(' ')[0]
            quote_time = dt.datetime.strptime(quote_time_text, "%I:%M%p")
            quote_time = quote_time.replace(year=CUR_YEAR, month=CUR_MONTH, day=CUR_DAY)
        except ValueError:
            quote_time = np.nan

        return quote_time

    def _get_option_data(self, expiry, name):
        frame_name = '_frames' + self._expiry_to_string(expiry)

        try:
            frames = getattr(self, frame_name)
        except AttributeError:
            frames = self._get_option_frames_from_yahoo(expiry)

        option_data = frames[name]
        if expiry != self.expiry_dates[0]:
            name += self._expiry_to_string(expiry)

        setattr(self, name, option_data)
        return option_data

    def get_call_data(self, month=None, year=None, expiry=None):
        """
        ***Experimental***
        Gets call/put data for the stock with the expiration data in the
        given month and year

        Parameters
        ----------
        month : number, int, optional(default=None)
            The month the options expire. This should be either 1 or 2
            digits.

        year : number, int, optional(default=None)
            The year the options expire. This should be a 4 digit int.

        expiry : date-like or convertible or list-like object, optional (default=None)
            The date (or dates) when options expire (defaults to current month)

        Returns
        -------
        call_data: pandas.DataFrame
            A DataFrame with requested options data.

            Index:
                Strike: Option strike, int
                Expiry: Option expiry, Timestamp
                Type: Call or Put, string
                Symbol: Option symbol as reported on Yahoo, string
            Columns:
                Last: Last option price, float
                Chg: Change from prior day, float
                Bid: Bid price, float
                Ask: Ask price, float
                Vol: Volume traded, int64
                Open_Int: Open interest, int64
                IsNonstandard: True if the the deliverable is not 100 shares, otherwise false
                Underlying: Ticker of the underlying security, string
                Underlying_Price: Price of the underlying security, float64
                Quote_Time: Time of the quote, Timestamp

        Notes
        -----
        Note: Format of returned data frame is dependent on Yahoo and may change.

        When called, this function will add instance variables named
        calls and puts. See the following example:

            >>> aapl = Options('aapl', 'yahoo')  # Create object
            >>> aapl.calls  # will give an AttributeError
            >>> aapl.get_call_data()  # Get data and set ivars
            >>> aapl.calls  # Doesn't throw AttributeError

        Also note that aapl.calls will always be the calls for the next
        expiry. If the user calls this method with a different month
        or year, the ivar will be named callsYYMMDD where YY, MM and DD are,
        respectively, two digit representations of the year, month and day
        for the expiry of the options.
        """
        expiry = self._try_parse_dates(year, month, expiry)
        return self._get_data_in_date_range(expiry, call=True, put=False)

    def get_put_data(self, month=None, year=None, expiry=None):
        """
        ***Experimental***
        Gets put data for the stock with the expiration data in the
        given month and year

        Parameters
        ----------
        month : number, int, optional(default=None)
            The month the options expire. This should be either 1 or 2
            digits.

        year : number, int, optional(default=None)
            The year the options expire. This should be a 4 digit int.

        expiry : date-like or convertible or list-like object, optional (default=None)
            The date (or dates) when options expire (defaults to current month)

        Returns
        -------
        put_data: pandas.DataFrame
            A DataFrame with requested options data.

            Index:
                Strike: Option strike, int
                Expiry: Option expiry, Timestamp
                Type: Call or Put, string
                Symbol: Option symbol as reported on Yahoo, string
            Columns:
                Last: Last option price, float
                Chg: Change from prior day, float
                Bid: Bid price, float
                Ask: Ask price, float
                Vol: Volume traded, int64
                Open_Int: Open interest, int64
                IsNonstandard: True if the the deliverable is not 100 shares, otherwise false
                Underlying: Ticker of the underlying security, string
                Underlying_Price: Price of the underlying security, float64
                Quote_Time: Time of the quote, Timestamp

        Notes
        -----
        Note: Format of returned data frame is dependent on Yahoo and may change.

        When called, this function will add instance variables named
        puts. See the following example:

            >>> aapl = Options('aapl')  # Create object
            >>> aapl.puts  # will give an AttributeError
            >>> aapl.get_put_data()  # Get data and set ivars
            >>> aapl.puts  # Doesn't throw AttributeError

                    return self.__setattr__(self, str(str(x) + str(y)))

        Also note that aapl.puts will always be the puts for the next
        expiry. If the user calls this method with a different month
        or year, the ivar will be named putsYYMMDD where YY, MM and DD are,
        respectively, two digit representations of the year, month and day
        for the expiry of the options.
        """
        expiry = self._try_parse_dates(year, month, expiry)
        return self._get_data_in_date_range(expiry, put=True, call=False)

    def get_near_stock_price(self, above_below=2, call=True, put=False,
                             month=None, year=None, expiry=None):
        """
        ***Experimental***
        Returns a data frame of options that are near the current stock price.

        Parameters
        ----------
        above_below : number, int, optional (default=2)
            The number of strike prices above and below the stock price that
            should be taken

        call : bool, default: True
            Tells the function whether or not it should be using calls

        put : bool, default: False
            Tells the function weather or not it should be using puts

        month : number, int, optional(default=None)
            The month the options expire. This should be either 1 or 2
            digits.

        year : number, int, optional(default=None)
            The year the options expire. This should be a 4 digit int.

        expiry : date-like or convertible or list-like object, optional (default=None)
            The date (or dates) when options expire (defaults to current month)

        Returns
        -------
        chopped: DataFrame
            The resultant DataFrame chopped down to be 2 * above_below + 1 rows
            desired. If there isn't data as far out as the user has asked for
            then

         Note: Format of returned data frame is dependent on Yahoo and may change.

        """
        expiry = self._try_parse_dates(year, month, expiry)
        data = self._get_data_in_date_range(expiry, call=call, put=put)
        return self.chop_data(data, above_below, self.underlying_price)

    def chop_data(self, df, above_below=2, underlying_price=None):
        """Returns a data frame only options that are near the current stock price."""

        if not underlying_price:
            try:
                underlying_price = self.underlying_price
            except AttributeError:
                underlying_price = np.nan

        max_strike = max(df.index.get_level_values('Strike'))
        min_strike = min(df.index.get_level_values('Strike'))

        if not np.isnan(underlying_price) and min_strike < underlying_price < max_strike:
            start_index = np.where(df.index.get_level_values('Strike')
                                   > underlying_price)[0][0]

            get_range = slice(start_index - above_below,
                              start_index + above_below + 1)
            df = df[get_range].dropna(how='all')

        return df

    def _try_parse_dates(self, year, month, expiry):
        """
        Validates dates provided by user.  Ensures the user either provided both a month and a year or an expiry.

        Parameters
        ----------
        year : int
            Calendar year

        month : int
            Calendar month

        expiry : date-like or convertible, (preferred)
            Expiry date

        Returns
        -------
        list of expiry dates (datetime.date)
        """

        #Checks if the user gave one of the month or the year but not both and did not provide an expiry:
        if (month is not None and year is None) or (month is None and year is not None) and expiry is None:
            msg = "You must specify either (`year` and `month`) or `expiry` " \
                  "or none of these options for the next expiry."
            raise ValueError(msg)

        if expiry is not None:
            if hasattr(expiry, '__iter__'):
                expiry = [self._validate_expiry(exp) for exp in expiry]
            else:
                expiry = [self._validate_expiry(expiry)]

            if len(expiry) == 0:
                raise ValueError('No expiries available for given input.')

        elif year is None and month is None:
            #No arguments passed, provide next expiry
            year = CUR_YEAR
            month = CUR_MONTH
            expiry = dt.date(year, month, 1)
            expiry = [self._validate_expiry(expiry)]

        else:
            #Year and month passed, provide all expiries in that month
            expiry = [expiry for expiry in self.expiry_dates if expiry.year == year and expiry.month == month]
            if len(expiry) == 0:
                raise ValueError('No expiries available in %s-%s' % (year, month))

        return expiry

    def _validate_expiry(self, expiry):
        """Ensures that an expiry date has data available on Yahoo
        If the expiry date does not have options that expire on that day, return next expiry"""

        expiry_dates = self.expiry_dates
        expiry = to_datetime(expiry)
        if hasattr(expiry, 'date'):
            expiry = expiry.date()

        if expiry in expiry_dates:
            return expiry
        else:
            index = DatetimeIndex(expiry_dates).order()
            return index[index.date >= expiry][0].date()

    def get_forward_data(self, months, call=True, put=False, near=False,
                         above_below=2):
        """
        ***Experimental***
        Gets either call, put, or both data for months starting in the current
        month and going out in the future a specified amount of time.

        Parameters
        ----------
        months : number, int
            How many months to go out in the collection of the data. This is
            inclusive.

        call : bool, optional (default=True)
            Whether or not to collect data for call options

        put : bool, optional (default=False)
            Whether or not to collect data for put options.

        near : bool, optional (default=False)
            Whether this function should get only the data near the
            current stock price. Uses Options.get_near_stock_price

        above_below : number, int, optional (default=2)
            The number of strike prices above and below the stock price that
            should be taken if the near option is set to True

        Returns
        -------
        pandas.DataFrame
            A DataFrame with requested options data.

            Index:
                Strike: Option strike, int
                Expiry: Option expiry, Timestamp
                Type: Call or Put, string
                Symbol: Option symbol as reported on Yahoo, string
            Columns:
                Last: Last option price, float
                Chg: Change from prior day, float
                Bid: Bid price, float
                Ask: Ask price, float
                Vol: Volume traded, int64
                Open_Int: Open interest, int64
                IsNonstandard: True if the the deliverable is not 100 shares, otherwise false
                Underlying: Ticker of the underlying security, string
                Underlying_Price: Price of the underlying security, float64
                Quote_Time: Time of the quote, Timestamp

                Note: Format of returned data frame is dependent on Yahoo and may change.

        """
        warnings.warn("get_forward_data() is deprecated", FutureWarning,
                      stacklevel=2)
        end_date = dt.date.today() + MonthEnd(months)
        dates = (date for date in self.expiry_dates if date <= end_date.date())
        data = self._get_data_in_date_range(dates, call=call, put=put)
        if near:
            data = self.chop_data(data, above_below=above_below)
        return data

    def get_all_data(self, call=True, put=True):
        """
        ***Experimental***
        Gets either call, put, or both data for all available months starting
        in the current month.

        Parameters
        ----------
        call : bool, optional (default=True)
            Whether or not to collect data for call options

        put : bool, optional (default=True)
            Whether or not to collect data for put options.

        Returns
        -------
        pandas.DataFrame
            A DataFrame with requested options data.

            Index:
                Strike: Option strike, int
                Expiry: Option expiry, Timestamp
                Type: Call or Put, string
                Symbol: Option symbol as reported on Yahoo, string
            Columns:
                Last: Last option price, float
                Chg: Change from prior day, float
                Bid: Bid price, float
                Ask: Ask price, float
                Vol: Volume traded, int64
                Open_Int: Open interest, int64
                IsNonstandard: True if the the deliverable is not 100 shares, otherwise false
                Underlying: Ticker of the underlying security, string
                Underlying_Price: Price of the underlying security, float64
                Quote_Time: Time of the quote, Timestamp

        Note: Format of returned data frame is dependent on Yahoo and may change.

        """

        try:
            expiry_dates = self.expiry_dates
        except AttributeError:
            expiry_dates, _ = self._get_expiry_dates_and_links()

        return self._get_data_in_date_range(dates=expiry_dates, call=call, put=put)

    def _get_data_in_date_range(self, dates, call=True, put=True):

        to_ret = Series({'calls': call, 'puts': put})
        to_ret = to_ret[to_ret].index
        data = []

        for name in to_ret:
            for expiry_date in dates:
                nam = name + self._expiry_to_string(expiry_date)
                try:  # Try to access on the instance
                    frame = getattr(self, nam)
                except AttributeError:
                    frame = self._get_option_data(expiry=expiry_date, name=name)
                data.append(frame)

        return concat(data).sortlevel()

    @property
    def expiry_dates(self):
        """
        Returns a list of available expiry dates
        """
        try:
            expiry_dates = self._expiry_dates
        except AttributeError:
            expiry_dates, _ = self._get_expiry_dates_and_links()
        return expiry_dates

    def _get_expiry_dates_and_links(self):
        """
        Gets available expiry dates.

        Returns
        -------
        Tuple of:
        List of datetime.date objects
        Dict of datetime.date objects as keys and corresponding links
        """

        url = self._OPTIONS_BASE_URL.format(sym=self.symbol)
        root = self._parse_url(url)

        try:
            links = root.xpath('//*[@id="options_menu"]/form/select/option')
        except IndexError:
            raise RemoteDataError('Expiry dates not available')

        expiry_dates = [dt.datetime.strptime(element.text, "%B %d, %Y").date() for element in links]
        links = [element.attrib['data-selectbox-link'] for element in links]

        if len(expiry_dates) == 0:
            raise RemoteDataError('Data not available')

        expiry_links = dict(zip(expiry_dates, links))
        self._expiry_links = expiry_links
        self._expiry_dates = expiry_dates
        return expiry_dates, expiry_links

    def _parse_url(self, url):
        """
        Downloads and parses a URL, returns xml root.

        """
        try:
            from lxml.html import parse
        except ImportError:
            raise ImportError("Please install lxml if you want to use the "
                              "{0!r} class".format(self.__class__.__name__))
        try:
            doc = parse(url)
        except _network_error_classes:
            raise RemoteDataError("Unable to parse URL "
                                  "{0!r}".format(url))
        else:
            root = doc.getroot()
            if root is None:
                raise RemoteDataError("Parsed URL {0!r} has no root"
                                      "element".format(url))
        return root

    def _process_data(self, frame, type):
        """
        Adds columns for Expiry, IsNonstandard (ie: deliverable is not 100 shares)
        and Tag (the tag indicating what is actually deliverable, None if standard).

        """
        frame.columns = ['Strike', 'Symbol', 'Last', 'Bid', 'Ask', 'Chg', 'PctChg', 'Vol', 'Open_Int', 'IV']
        frame["Rootexp"] = frame.Symbol.str[0:-9]
        frame["Root"] = frame.Rootexp.str[0:-6]
        frame["Expiry"] = to_datetime(frame.Rootexp.str[-6:])
        #Removes dashes in equity ticker to map to option ticker.
        #Ex: BRK-B to BRKB140517C00100000
        frame["IsNonstandard"] = frame['Root'] != self.symbol.replace('-', '')
        del frame["Rootexp"]
        frame["Underlying"] = self.symbol
        try:
            frame['Underlying_Price'] = self.underlying_price
            frame["Quote_Time"] = self.quote_time
        except AttributeError:
            frame['Underlying_Price'] = np.nan
            frame["Quote_Time"] = np.nan
        frame.rename(columns={'Open Int': 'Open_Int'}, inplace=True)
        frame['Type'] = type
        frame.set_index(['Strike', 'Expiry', 'Type', 'Symbol'], inplace=True)

        return frame
