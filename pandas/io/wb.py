import urllib2
import warnings
import json
import pandas
import numpy as np

def download(country=['MX','CA','US'], indicator=['GDPPCKD','GDPPCKN'], 
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
    valid_countries = ["AG", "AL", "AM", "AO", "AR", "AT", "AU", "AZ", "BB",
            "BD", "BE", "BF", "BG", "BH", "BI", "BJ", "BO", "BR", "BS", "BW",
            "BY", "BZ", "CA", "CD", "CF", "CG", "CH", "CI", "CL", "CM", "CN",
            "CO", "CR", "CV", "CY", "CZ", "DE", "DK", "DM", "DO", "DZ", "EC",
            "EE", "EG", "ER", "ES", "ET", "FI", "FJ", "FR", "GA", "GB", "GE",
            "GH", "GM", "GN", "GQ", "GR", "GT", "GW", "GY", "HK", "HN", "HR",
            "HT", "HU", "ID", "IE", "IL", "IN", "IR", "IS", "IT", "JM", "JO",
            "JP", "KE", "KG", "KH", "KM", "KR", "KW", "KZ", "LA", "LB", "LC",
            "LK", "LS", "LT", "LU", "LV", "MA", "MD", "MG", "MK", "ML", "MN",
            "MR", "MU", "MW", "MX", "MY", "MZ", "NA", "NE", "NG", "NI", "NL",
            "NO", "NP", "NZ", "OM", "PA", "PE", "PG", "PH", "PK", "PL", "PT",
            "PY", "RO", "RU", "RW", "SA", "SB", "SC", "SD", "SE", "SG", "SI",
            "SK", "SL", "SN", "SR", "SV", "SY", "SZ", "TD", "TG", "TH", "TN",
            "TR", "TT", "TW", "TZ", "UA", "UG", "US", "UY", "UZ", "VC", "VE",
            "VN", "VU", "YE", "ZA", "ZM", "ZW", "all"]    
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
            tmp = _WDIdownload(ind, country, start, end)
            tmp.columns = ['country', 'iso2c', 'year', ind]
            data.append(tmp)
        except:
            bad_indicators.append(ind) 
    # Warn
    if len(bad_indicators) > 0:
        print 'Failed to obtain indicator(s): ' + '; '.join(bad_indicators)
        print 'The data may still be available for download at http://data.worldbank.org'
    if len(bad_countries) > 0:
        print 'Invalid ISO-2 codes: ' + ' '.join(bad_countries)
    # Merge WDI series
    if len(data) > 0:
        out = reduce(lambda x,y: x.merge(y, how='outer'), data)
        # Clean
        out = out.drop('iso2c', axis=1)
        out = out.set_index(['country', 'year'])
        return out

def _WDIdownload(indicator = "NY.GNS.ICTR.GN.ZS", country = 'US', 
                start = 2002, end = 2005):
    # Build URL for api call
    url = "http://api.worldbank.org/countries/" + country + "/indicators/" + \
        indicator + "?date=" + str(start) + ":" + str(end) + "&per_page=25000" + \
        "&format=json"
    # Download
    response = urllib2.urlopen(url)
    data = response.read()
    # Parse JSON file
    data = json.loads(data)[1]
    country = map(lambda x: x['country']['value'], data)
    iso2c = map(lambda x: x['country']['id'], data)
    year = map(lambda x: x['date'], data)
    value = map(lambda x: x['value'], data)
    # Prepare output
    out = pandas.DataFrame([country, iso2c, year, value]).T
    return out
     

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
        url = 'http://api.worldbank.org/indicators?per_page=25000&format=json'
        response = urllib2.urlopen(url)
        data = response.read()
        data = json.loads(data)[1]
        _cached_series = pandas.DataFrame(data)
    data = _cached_series[field]
    idx = data.str.contains(string, case=case)
    out = _cached_series.ix[idx].dropna()
    return out
