from __future__ import print_function

from pandas.compat import map, reduce, range, lrange
from pandas.io.common import urlopen
from pandas.io import json
import pandas
import numpy as np


def download(country=['MX', 'CA', 'US'], indicator=['GDPPCKD', 'GDPPCKN'],
             start=2003, end=2005):
    """
    Download data series from the World Bank's World Development Indicators

    Parameters
    ----------

    indicator: string or list of strings
        taken from the ``id`` field in ``WDIsearch()``
    country: string or list of strings.
        ``all`` downloads data for all countries
        ISO-2 character codes select individual countries (e.g.``US``,``CA``)
    start: int
        First year of the data series
    end: int
        Last year of the data series (inclusive)

    Returns
    -------

    ``pandas`` DataFrame with columns: country, iso2c, year, indicator value.
    """

    # Are ISO-2 country codes valid?
    valid_countries = [
        "AG", "AL", "AM", "AO", "AR", "AT", "AU", "AZ", "BB", "BD", "BE", "BF",
        "BG", "BH", "BI", "BJ", "BO", "BR", "BS", "BW", "BY", "BZ", "CA", "CD",
        "CF", "CG", "CH", "CI", "CL", "CM", "CN", "CO", "CR", "CV", "CY", "CZ",
        "DE", "DK", "DM", "DO", "DZ", "EC", "EE", "EG", "ER", "ES", "ET", "FI",
        "FJ", "FR", "GA", "GB", "GE", "GH", "GM", "GN", "GQ", "GR", "GT", "GW",
        "GY", "HK", "HN", "HR", "HT", "HU", "ID", "IE", "IL", "IN", "IR", "IS",
        "IT", "JM", "JO", "JP", "KE", "KG", "KH", "KM", "KR", "KW", "KZ", "LA",
        "LB", "LC", "LK", "LS", "LT", "LU", "LV", "MA", "MD", "MG", "MK", "ML",
        "MN", "MR", "MU", "MW", "MX", "MY", "MZ", "NA", "NE", "NG", "NI", "NL",
        "NO", "NP", "NZ", "OM", "PA", "PE", "PG", "PH", "PK", "PL", "PT", "PY",
        "RO", "RU", "RW", "SA", "SB", "SC", "SD", "SE", "SG", "SI", "SK", "SL",
        "SN", "SR", "SV", "SY", "SZ", "TD", "TG", "TH", "TN", "TR", "TT", "TW",
        "TZ", "UA", "UG", "US", "UY", "UZ", "VC", "VE", "VN", "VU", "YE", "ZA",
        "ZM", "ZW", "all"
    ]
    if type(country) == str:
        country = [country]
    bad_countries = np.setdiff1d(country, valid_countries)
    country = np.intersect1d(country, valid_countries)
    country = ';'.join(country)
    # Work with a list of indicators
    if type(indicator) == str:
        indicator = [indicator]
    # Download
    data = []
    bad_indicators = []
    for ind in indicator:
        try:
            tmp = _get_data(ind, country, start, end)
            tmp.columns = ['country', 'iso2c', 'year', ind]
            data.append(tmp)
        except:
            bad_indicators.append(ind)
    # Warn
    if len(bad_indicators) > 0:
        print('Failed to obtain indicator(s): %s' % '; '.join(bad_indicators))
        print('The data may still be available for download at '
              'http://data.worldbank.org')
    if len(bad_countries) > 0:
        print('Invalid ISO-2 codes: %s' % ' '.join(bad_countries))
    # Merge WDI series
    if len(data) > 0:
        out = reduce(lambda x, y: x.merge(y, how='outer'), data)
        # Clean
        out = out.drop('iso2c', axis=1)
        out = out.set_index(['country', 'year'])
        out = out.convert_objects(convert_numeric=True)
        return out


def _get_data(indicator="NY.GNS.ICTR.GN.ZS", country='US',
              start=2002, end=2005):
    # Build URL for api call
    url = ("http://api.worldbank.org/countries/" + country + "/indicators/" +
           indicator + "?date=" + str(start) + ":" + str(end) +
           "&per_page=25000&format=json")
    # Download
    with urlopen(url) as response:
        data = response.read()
    # Parse JSON file
    data = json.loads(data)[1]
    country = [x['country']['value'] for x in data]
    iso2c = [x['country']['id'] for x in data]
    year = [x['date'] for x in data]
    value = [x['value'] for x in data]
    # Prepare output
    out = pandas.DataFrame([country, iso2c, year, value]).T
    return out


def get_countries():
    '''Query information about countries
    '''
    url = 'http://api.worldbank.org/countries/?per_page=1000&format=json'
    with urlopen(url) as response:
        data = response.read()
    data = json.loads(data)[1]
    data = pandas.DataFrame(data)
    data.adminregion = [x['value'] for x in data.adminregion]
    data.incomeLevel = [x['value'] for x in data.incomeLevel]
    data.lendingType = [x['value'] for x in data.lendingType]
    data.region = [x['value'] for x in data.region]
    data = data.rename(columns={'id': 'iso3c', 'iso2Code': 'iso2c'})
    return data


def get_indicators():
    '''Download information about all World Bank data series
    '''
    url = 'http://api.worldbank.org/indicators?per_page=50000&format=json'
    with urlopen(url) as response:
        data = response.read()
    data = json.loads(data)[1]
    data = pandas.DataFrame(data)
    # Clean fields
    data.source = [x['value'] for x in data.source]
    fun = lambda x: x.encode('ascii', 'ignore')
    data.sourceOrganization = data.sourceOrganization.apply(fun)
    # Clean topic field

    def get_value(x):
        try:
            return x['value']
        except:
            return ''
    fun = lambda x: [get_value(y) for y in x]
    data.topics = data.topics.apply(fun)
    data.topics = data.topics.apply(lambda x: ' ; '.join(x))
    # Clean outpu
    data = data.sort(columns='id')
    data.index = pandas.Index(lrange(data.shape[0]))
    return data


_cached_series = None


def search(string='gdp.*capi', field='name', case=False):
    """
    Search available data series from the world bank

    Parameters
    ----------

    string: string
        regular expression
    field: string
        id, name, source, sourceNote, sourceOrganization, topics
        See notes below
    case: bool
        case sensitive search?

    Notes
    -----

    The first time this function is run it will download and cache the full
    list of available series. Depending on the speed of your network
    connection, this can take time. Subsequent searches will use the cached
    copy, so they should be much faster.

    id : Data series indicator (for use with the ``indicator`` argument of
    ``WDI()``) e.g. NY.GNS.ICTR.GN.ZS"
    name: Short description of the data series
    source: Data collection project
    sourceOrganization: Data collection organization
    note:
    sourceNote:
    topics:
    """
    # Create cached list of series if it does not exist
    global _cached_series
    if type(_cached_series) is not pandas.core.frame.DataFrame:
        _cached_series = get_indicators()
    data = _cached_series[field]
    idx = data.str.contains(string, case=case)
    out = _cached_series.ix[idx].dropna()
    return out
