# -*- coding: utf-8 -*-

# flake8: noqa

from __future__ import print_function

from pandas.compat import map, reduce, range, lrange
from pandas.io.common import urlopen
from pandas.io import json
import pandas
import numpy as np
import warnings

warnings.warn("\n"
              "The pandas.io.wb module is moved to a separate package "
              "(pandas-datareader) and will be removed from pandas in a "
              "future version.\nAfter installing the pandas-datareader package "
              "(https://github.com/pydata/pandas-datareader), you can change "
              "the import ``from pandas.io import data, wb`` to "
              "``from pandas_datareader import data, wb``.",
              FutureWarning)


# This list of country codes was pulled from wikipedia during October 2014.
# While some exceptions do exist, it is the best proxy for countries supported
# by World Bank.  It is an aggregation of the 2-digit ISO 3166-1 alpha-2, and
# 3-digit ISO 3166-1 alpha-3, codes, with 'all', 'ALL', and 'All' appended ot
# the end.

country_codes = ['AD', 'AE', 'AF', 'AG', 'AI', 'AL', 'AM', 'AO', 'AQ', 'AR', \
                 'AS', 'AT', 'AU', 'AW', 'AX', 'AZ', 'BA', 'BB', 'BD', 'BE', \
                 'BF', 'BG', 'BH', 'BI', 'BJ', 'BL', 'BM', 'BN', 'BO', 'BQ', \
                 'BR', 'BS', 'BT', 'BV', 'BW', 'BY', 'BZ', 'CA', 'CC', 'CD', \
                 'CF', 'CG', 'CH', 'CI', 'CK', 'CL', 'CM', 'CN', 'CO', 'CR', \
                 'CU', 'CV', 'CW', 'CX', 'CY', 'CZ', 'DE', 'DJ', 'DK', 'DM', \
                 'DO', 'DZ', 'EC', 'EE', 'EG', 'EH', 'ER', 'ES', 'ET', 'FI', \
                 'FJ', 'FK', 'FM', 'FO', 'FR', 'GA', 'GB', 'GD', 'GE', 'GF', \
                 'GG', 'GH', 'GI', 'GL', 'GM', 'GN', 'GP', 'GQ', 'GR', 'GS', \
                 'GT', 'GU', 'GW', 'GY', 'HK', 'HM', 'HN', 'HR', 'HT', 'HU', \
                 'ID', 'IE', 'IL', 'IM', 'IN', 'IO', 'IQ', 'IR', 'IS', 'IT', \
                 'JE', 'JM', 'JO', 'JP', 'KE', 'KG', 'KH', 'KI', 'KM', 'KN', \
                 'KP', 'KR', 'KW', 'KY', 'KZ', 'LA', 'LB', 'LC', 'LI', 'LK', \
                 'LR', 'LS', 'LT', 'LU', 'LV', 'LY', 'MA', 'MC', 'MD', 'ME', \
                 'MF', 'MG', 'MH', 'MK', 'ML', 'MM', 'MN', 'MO', 'MP', 'MQ', \
                 'MR', 'MS', 'MT', 'MU', 'MV', 'MW', 'MX', 'MY', 'MZ', 'NA', \
                 'NC', 'NE', 'NF', 'NG', 'NI', 'NL', 'NO', 'NP', 'NR', 'NU', \
                 'NZ', 'OM', 'PA', 'PE', 'PF', 'PG', 'PH', 'PK', 'PL', 'PM', \
                 'PN', 'PR', 'PS', 'PT', 'PW', 'PY', 'QA', 'RE', 'RO', 'RS', \
                 'RU', 'RW', 'SA', 'SB', 'SC', 'SD', 'SE', 'SG', 'SH', 'SI', \
                 'SJ', 'SK', 'SL', 'SM', 'SN', 'SO', 'SR', 'SS', 'ST', 'SV', \
                 'SX', 'SY', 'SZ', 'TC', 'TD', 'TF', 'TG', 'TH', 'TJ', 'TK', \
                 'TL', 'TM', 'TN', 'TO', 'TR', 'TT', 'TV', 'TW', 'TZ', 'UA', \
                 'UG', 'UM', 'US', 'UY', 'UZ', 'VA', 'VC', 'VE', 'VG', 'VI', \
                 'VN', 'VU', 'WF', 'WS', 'YE', 'YT', 'ZA', 'ZM', 'ZW', \
                 'ABW', 'AFG', 'AGO', 'AIA', 'ALA', 'ALB', 'AND', 'ARE', \
                 'ARG', 'ARM', 'ASM', 'ATA', 'ATF', 'ATG', 'AUS', 'AUT', \
                 'AZE', 'BDI', 'BEL', 'BEN', 'BES', 'BFA', 'BGD', 'BGR', \
                 'BHR', 'BHS', 'BIH', 'BLM', 'BLR', 'BLZ', 'BMU', 'BOL', \
                 'BRA', 'BRB', 'BRN', 'BTN', 'BVT', 'BWA', 'CAF', 'CAN', \
                 'CCK', 'CHE', 'CHL', 'CHN', 'CIV', 'CMR', 'COD', 'COG', \
                 'COK', 'COL', 'COM', 'CPV', 'CRI', 'CUB', 'CUW', 'CXR', \
                 'CYM', 'CYP', 'CZE', 'DEU', 'DJI', 'DMA', 'DNK', 'DOM', \
                 'DZA', 'ECU', 'EGY', 'ERI', 'ESH', 'ESP', 'EST', 'ETH', \
                 'FIN', 'FJI', 'FLK', 'FRA', 'FRO', 'FSM', 'GAB', 'GBR', \
                 'GEO', 'GGY', 'GHA', 'GIB', 'GIN', 'GLP', 'GMB', 'GNB', \
                 'GNQ', 'GRC', 'GRD', 'GRL', 'GTM', 'GUF', 'GUM', 'GUY', \
                 'HKG', 'HMD', 'HND', 'HRV', 'HTI', 'HUN', 'IDN', 'IMN', \
                 'IND', 'IOT', 'IRL', 'IRN', 'IRQ', 'ISL', 'ISR', 'ITA', \
                 'JAM', 'JEY', 'JOR', 'JPN', 'KAZ', 'KEN', 'KGZ', 'KHM', \
                 'KIR', 'KNA', 'KOR', 'KWT', 'LAO', 'LBN', 'LBR', 'LBY', \
                 'LCA', 'LIE', 'LKA', 'LSO', 'LTU', 'LUX', 'LVA', 'MAC', \
                 'MAF', 'MAR', 'MCO', 'MDA', 'MDG', 'MDV', 'MEX', 'MHL', \
                 'MKD', 'MLI', 'MLT', 'MMR', 'MNE', 'MNG', 'MNP', 'MOZ', \
                 'MRT', 'MSR', 'MTQ', 'MUS', 'MWI', 'MYS', 'MYT', 'NAM', \
                 'NCL', 'NER', 'NFK', 'NGA', 'NIC', 'NIU', 'NLD', 'NOR', \
                 'NPL', 'NRU', 'NZL', 'OMN', 'PAK', 'PAN', 'PCN', 'PER', \
                 'PHL', 'PLW', 'PNG', 'POL', 'PRI', 'PRK', 'PRT', 'PRY', \
                 'PSE', 'PYF', 'QAT', 'REU', 'ROU', 'RUS', 'RWA', 'SAU', \
                 'SDN', 'SEN', 'SGP', 'SGS', 'SHN', 'SJM', 'SLB', 'SLE', \
                 'SLV', 'SMR', 'SOM', 'SPM', 'SRB', 'SSD', 'STP', 'SUR', \
                 'SVK', 'SVN', 'SWE', 'SWZ', 'SXM', 'SYC', 'SYR', 'TCA', \
                 'TCD', 'TGO', 'THA', 'TJK', 'TKL', 'TKM', 'TLS', 'TON', \
                 'TTO', 'TUN', 'TUR', 'TUV', 'TWN', 'TZA', 'UGA', 'UKR', \
                 'UMI', 'URY', 'USA', 'UZB', 'VAT', 'VCT', 'VEN', 'VGB', \
                 'VIR', 'VNM', 'VUT', 'WLF', 'WSM', 'YEM', 'ZAF', 'ZMB', \
                 'ZWE', 'all', 'ALL', 'All']

def download(country=None, indicator=None,
             start=2003, end=2005,errors='warn'):
    """
    Download data series from the World Bank's World Development Indicators

    Parameters
    ----------

    indicator: string or list of strings
        taken from the ``id`` field in ``WDIsearch()``

    country: string or list of strings.
        ``all`` downloads data for all countries
        2 or 3 character ISO country codes select individual
        countries (e.g.``US``,``CA``) or (e.g.``USA``,``CAN``).  The codes
        can be mixed.

        The two ISO lists of countries, provided by wikipedia, are hardcoded
        into pandas as of 11/10/2014.

    start: int
        First year of the data series

    end: int
        Last year of the data series (inclusive)

    errors: str {'ignore', 'warn', 'raise'}, default 'warn'
        Country codes are validated against a hardcoded list.  This controls
        the outcome of that validation, and attempts to also apply
        to the results from world bank.

        errors='raise', will raise a ValueError on a bad country code.

    Returns
    -------

    ``pandas`` DataFrame with columns: country, iso_code, year,
    indicator value.

    """
    if country is None:
        country = ['MX', 'CA', 'US']
    if indicator is None:
        indicator = ['NY.GDP.MKTP.CD', 'NY.GNS.ICTR.ZS']

    if type(country) == str:
        country = [country]

    bad_countries = np.setdiff1d(country, country_codes)

    # Validate the input
    if len(bad_countries) > 0:
        tmp = ", ".join(bad_countries)
        if errors == 'raise':
            raise ValueError("Invalid Country Code(s): %s" % tmp)
        if errors == 'warn':
            warnings.warn('Non-standard ISO country codes: %s' % tmp)

    # Work with a list of indicators
    if type(indicator) == str:
        indicator = [indicator]

    # Download
    data = []
    bad_indicators = {}
    for ind in indicator:
        one_indicator_data,msg = _get_data(ind, country, start, end)
        if msg == "Success":
            data.append(one_indicator_data)
        else:
            bad_indicators[ind] = msg

    if len(bad_indicators.keys()) > 0:
        bad_ind_msgs = [i + " : " + m for i,m in bad_indicators.items()]
        bad_ind_msgs = "\n\n".join(bad_ind_msgs)
        bad_ind_msgs = "\n\nInvalid Indicators:\n\n%s" % bad_ind_msgs
        if errors == 'raise':
            raise ValueError(bad_ind_msgs)
        if errors == 'warn':
            warnings.warn(bad_ind_msgs)

    # Confirm we actually got some data, and build Dataframe
    if len(data) > 0:
        out = reduce(lambda x, y: x.merge(y, how='outer'), data)
        out = out.drop('iso_code', axis=1)
        out = out.set_index(['country', 'year'])
        out = out._convert(datetime=True, numeric=True)
        return out
    else:
        msg = "No indicators returned data."
        if errors == 'ignore':
            msg += "  Set errors='warn' for more information."
        raise ValueError(msg)

def _get_data(indicator="NY.GNS.ICTR.GN.ZS", country='US',
              start=2002, end=2005):

    if type(country) == str:
        country = [country]

    countries = ';'.join(country)

    # Build URL for api call
    url = ("http://api.worldbank.org/countries/" + countries + "/indicators/" +
           indicator + "?date=" + str(start) + ":" + str(end) +
           "&per_page=25000&format=json")

    # Download
    with urlopen(url) as response:
        data = response.read()

    # Check to see if there is a possible problem
    possible_message = json.loads(data)[0]
    if 'message' in possible_message.keys():
        msg = possible_message['message'][0]
        try:
            msg = msg['key'].split() + ["\n "] + msg['value'].split()
            wb_err = ' '.join(msg)
        except:
            wb_err = ""
            if 'key' in msg.keys():
                wb_err = msg['key'] + "\n "
            if 'value' in msg.keys():
                wb_err += msg['value']
        error_msg = "Problem with a World Bank Query \n %s"
        return None, error_msg % wb_err

    if 'total' in possible_message.keys():
        if possible_message['total'] == 0:
            return None, "No results from world bank."

    # Parse JSON file
    data = json.loads(data)[1]
    country = [x['country']['value'] for x in data]
    iso_code = [x['country']['id'] for x in data]
    year = [x['date'] for x in data]
    value = [x['value'] for x in data]
    # Prepare output
    out = pandas.DataFrame([country, iso_code, year, value]).T
    out.columns = ['country', 'iso_code', 'year', indicator]
    return out,"Success"

def get_countries():
    """Query information about countries
    """
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
    """Download information about all World Bank data series
    """
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
