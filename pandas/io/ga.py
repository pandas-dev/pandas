"""
1. Goto https://code.google.com/apis/console
2. Create new project
3. Goto APIs and register for OAuth2.0 for installed applications
4. Download JSON secret file and move into same directory as this file
"""
from datetime import datetime
import re
from pandas import compat
import numpy as np
from pandas import DataFrame
import pandas as pd
import pandas.io.parsers as psr
import pandas.lib as lib
from pandas.io.date_converters import generic_parser
import pandas.io.auth as auth
from pandas.util.decorators import Appender, Substitution

from apiclient.errors import HttpError
from oauth2client.client import AccessTokenRefreshError
from pandas.compat import zip, u

# GH11038
import warnings
warnings.warn("The pandas.io.ga module is deprecated and will be "
              "removed in a future version.",
              FutureWarning, stacklevel=2)

TYPE_MAP = {u('INTEGER'): int, u('FLOAT'): float, u('TIME'): int}

NO_CALLBACK = auth.OOB_CALLBACK_URN
DOC_URL = auth.DOC_URL

_QUERY_PARAMS = """metrics : list of str
    Un-prefixed metric names (e.g., 'visitors' and not 'ga:visitors')
dimensions : list of str
    Un-prefixed dimension variable names
start_date : str/date/datetime
end_date : str/date/datetime, optional, default is None but internally set as today
segment : list of str, optional, default: None
filters : list of str, optional, default: None
start_index : int, default 1
max_results : int, default 10000
    If >10000, must specify chunksize or ValueError will be raised"""

_QUERY_DOC = """
Construct a google analytics query using given parameters
Metrics and dimensions do not need the 'ga:' prefix

Parameters
----------
profile_id : str
%s
""" % _QUERY_PARAMS

_GA_READER_DOC = """Given query parameters, return a DataFrame with all the
data or an iterator that returns DataFrames containing chunks of the data

Parameters
----------
%s
sort : bool/list, default True
    Sort output by index or list of columns
chunksize : int, optional
    If max_results >10000, specifies the number of rows per iteration
index_col : str/list of str/dict, optional, default: None
    If unspecified then dimension variables are set as index
parse_dates : bool/list/dict, default: True
keep_date_col : boolean, default: False
date_parser : optional, default: None
na_values : optional, default: None
converters : optional, default: None
dayfirst : bool, default False
    Informs date parsing
account_name : str, optional, default: None
account_id : str, optional, default: None
property_name : str, optional, default: None
property_id : str, optional, default: None
profile_name : str, optional, default: None
profile_id : str, optional, default: None
%%(extras)s
Returns
-------
data : DataFrame or DataFrame yielding iterator
""" % _QUERY_PARAMS

_AUTH_PARAMS = """secrets : str, optional
    File path to the secrets file
scope : str, optional
    Authentication scope
token_file_name : str, optional
    Path to token storage
redirect : str, optional
    Local host redirect if unspecified
"""


def reset_token_store():
    """
    Deletes the default token store
    """
    auth.reset_default_token_store()


@Substitution(extras=_AUTH_PARAMS)
@Appender(_GA_READER_DOC)
def read_ga(metrics, dimensions, start_date, **kwargs):
    lst = ['secrets', 'scope', 'token_file_name', 'redirect']
    reader_kwds = dict((p, kwargs.pop(p)) for p in lst if p in kwargs)
    reader = GAnalytics(**reader_kwds)
    return reader.get_data(metrics=metrics, start_date=start_date,
                           dimensions=dimensions, **kwargs)


class OAuthDataReader(object):
    """
    Abstract class for handling OAuth2 authentication using the Google
    oauth2client library
    """
    def __init__(self, scope, token_file_name, redirect):
        """
        Parameters
        ----------
        scope : str
            Designates the authentication scope
        token_file_name : str
            Location of cache for authenticated tokens
        redirect : str
            Redirect URL
        """
        self.scope = scope
        self.token_store = auth.make_token_store(token_file_name)
        self.redirect_url = redirect

    def authenticate(self, secrets):
        """
        Run the authentication process and return an authorized
        http object

        Parameters
        ----------
        secrets : str
            File name for client secrets

        Notes
        -----
        See google documention for format of secrets file
        %s
        """ % DOC_URL
        flow = self._create_flow(secrets)
        return auth.authenticate(flow, self.token_store)

    def _create_flow(self, secrets):
        """
        Create an authentication flow based on the secrets file

        Parameters
        ----------
        secrets : str
            File name for client secrets

        Notes
        -----
        See google documentation for format of secrets file
        %s
        """ % DOC_URL
        return auth.get_flow(secrets, self.scope, self.redirect_url)


class GDataReader(OAuthDataReader):
    """
    Abstract class for reading data from google APIs using OAuth2
    Subclasses must implement create_query method
    """
    def __init__(self, scope=auth.DEFAULT_SCOPE,
                 token_file_name=auth.DEFAULT_TOKEN_FILE,
                 redirect=NO_CALLBACK, secrets=auth.DEFAULT_SECRETS):
        super(GDataReader, self).__init__(scope, token_file_name, redirect)
        self._service = self._init_service(secrets)

    @property
    def service(self):
        """The authenticated request service object"""
        return self._service

    def _init_service(self, secrets):
        """
        Build an authenticated google api request service using the given
        secrets file
        """
        http = self.authenticate(secrets)
        return auth.init_service(http)

    def get_account(self, name=None, id=None, **kwargs):
        """ Retrieve an account that matches the name, id, or some account
        attribute specified in **kwargs

        Parameters
        ----------
        name : str, optional, default: None
        id : str, optional, default: None
        """
        accounts = self.service.management().accounts().list().execute()
        return _get_match(accounts, name, id, **kwargs)

    def get_web_property(self, account_id=None, name=None, id=None, **kwargs):
        """
        Retrieve a web property given and account and property name, id, or
        custom attribute

        Parameters
        ----------
        account_id : str, optional, default: None
        name : str, optional, default: None
        id : str, optional, default: None
        """
        prop_store = self.service.management().webproperties()
        kwds = {}
        if account_id is not None:
            kwds['accountId'] = account_id
        prop_for_acct = prop_store.list(**kwds).execute()
        return _get_match(prop_for_acct, name, id, **kwargs)

    def get_profile(self, account_id=None, web_property_id=None, name=None,
                    id=None, **kwargs):

        """
        Retrieve the right profile for the given account, web property, and
        profile attribute (name, id, or arbitrary parameter in kwargs)

        Parameters
        ----------
        account_id : str, optional, default: None
        web_property_id : str, optional, default: None
        name : str, optional, default: None
        id : str, optional, default: None
        """
        profile_store = self.service.management().profiles()
        kwds = {}
        if account_id is not None:
            kwds['accountId'] = account_id
        if web_property_id is not None:
            kwds['webPropertyId'] = web_property_id
        profiles = profile_store.list(**kwds).execute()
        return _get_match(profiles, name, id, **kwargs)

    def create_query(self, *args, **kwargs):
        raise NotImplementedError()

    @Substitution(extras='')
    @Appender(_GA_READER_DOC)
    def get_data(self, metrics, start_date, end_date=None,
                 dimensions=None, segment=None, filters=None, start_index=1,
                 max_results=10000, index_col=None, parse_dates=True,
                 keep_date_col=False, date_parser=None, na_values=None,
                 converters=None, sort=True, dayfirst=False,
                 account_name=None, account_id=None, property_name=None,
                 property_id=None, profile_name=None, profile_id=None,
                 chunksize=None):
        if chunksize is None and max_results > 10000:
            raise ValueError('Google API returns maximum of 10,000 rows, '
                             'please set chunksize')

        account = self.get_account(account_name, account_id)
        web_property = self.get_web_property(account.get('id'), property_name,
                                             property_id)
        profile = self.get_profile(account.get('id'), web_property.get('id'),
                                   profile_name, profile_id)

        profile_id = profile.get('id')

        if index_col is None and dimensions is not None:
            if isinstance(dimensions, compat.string_types):
                dimensions = [dimensions]
            index_col = _clean_index(list(dimensions), parse_dates)

        def _read(start, result_size):
            query = self.create_query(profile_id, metrics, start_date,
                                      end_date=end_date, dimensions=dimensions,
                                      segment=segment, filters=filters,
                                      start_index=start,
                                      max_results=result_size)

            try:
                rs = query.execute()
                rows = rs.get('rows', [])
                col_info = rs.get('columnHeaders', [])
                return self._parse_data(rows, col_info, index_col,
                                        parse_dates=parse_dates,
                                        keep_date_col=keep_date_col,
                                        date_parser=date_parser,
                                        dayfirst=dayfirst,
                                        na_values=na_values,
                                        converters=converters, sort=sort)
            except HttpError as inst:
                raise ValueError('Google API error %s: %s' % (inst.resp.status,
                                 inst._get_reason()))

        if chunksize is None:
            return _read(start_index, max_results)

        def iterator():
            curr_start = start_index

            while curr_start < max_results:
                yield _read(curr_start, chunksize)
                curr_start += chunksize
        return iterator()

    def _parse_data(self, rows, col_info, index_col, parse_dates=True,
                    keep_date_col=False, date_parser=None, dayfirst=False,
                    na_values=None, converters=None, sort=True):
        # TODO use returned column types
        col_names = _get_col_names(col_info)
        df = psr._read(rows, dict(index_col=index_col, parse_dates=parse_dates,
                                  date_parser=date_parser, dayfirst=dayfirst,
                                  na_values=na_values,
                                  keep_date_col=keep_date_col,
                                  converters=converters,
                                  header=None, names=col_names))

        if isinstance(sort, bool) and sort:
            return df.sort_index()
        elif isinstance(sort, (compat.string_types, list, tuple, np.ndarray)):
            return df.sort_index(by=sort)

        return df


class GAnalytics(GDataReader):

    @Appender(_QUERY_DOC)
    def create_query(self, profile_id, metrics, start_date, end_date=None,
                     dimensions=None, segment=None, filters=None,
                     start_index=None, max_results=10000, **kwargs):
        qry = format_query(profile_id, metrics, start_date, end_date=end_date,
                           dimensions=dimensions, segment=segment,
                           filters=filters, start_index=start_index,
                           max_results=max_results, **kwargs)
        try:
            return self.service.data().ga().get(**qry)
        except TypeError as error:
            raise ValueError('Error making query: %s' % error)


def format_query(ids, metrics, start_date, end_date=None, dimensions=None,
                 segment=None, filters=None, sort=None, start_index=None,
                 max_results=10000, **kwargs):
    if isinstance(metrics, compat.string_types):
        metrics = [metrics]
    met = ','.join(['ga:%s' % x for x in metrics])

    start_date = pd.to_datetime(start_date).strftime('%Y-%m-%d')
    if end_date is None:
        end_date = datetime.today()
    end_date = pd.to_datetime(end_date).strftime('%Y-%m-%d')

    qry = dict(ids='ga:%s' % str(ids),
               metrics=met,
               start_date=start_date,
               end_date=end_date)
    qry.update(kwargs)

    names = ['dimensions', 'filters', 'sort']
    lst = [dimensions, filters, sort]
    [_maybe_add_arg(qry, n, d) for n, d in zip(names, lst)]

    if isinstance(segment, compat.string_types):
        if re.match("^[a-zA-Z0-9\-\_]+$", segment):
            _maybe_add_arg(qry, 'segment', segment, 'gaid:')
        else:
            _maybe_add_arg(qry, 'segment', segment, 'dynamic::ga')
    elif isinstance(segment, int):
        _maybe_add_arg(qry, 'segment', segment, 'gaid:')
    elif segment:
        raise ValueError("segment must be string for dynamic and int ID")

    if start_index is not None:
        qry['start_index'] = str(start_index)

    if max_results is not None:
        qry['max_results'] = str(max_results)

    return qry


def _maybe_add_arg(query, field, data, prefix='ga'):
    if data is not None:
        if isinstance(data, (compat.string_types, int)):
            data = [data]
        data = ','.join(['%s:%s' % (prefix, x) for x in data])
        query[field] = data


def _get_match(obj_store, name, id, **kwargs):
    key, val = None, None
    if len(kwargs) > 0:
        key = list(kwargs.keys())[0]
        val = list(kwargs.values())[0]

    if name is None and id is None and key is None:
        return obj_store.get('items')[0]

    name_ok = lambda item: name is not None and item.get('name') == name
    id_ok = lambda item: id is not None and item.get('id') == id
    key_ok = lambda item: key is not None and item.get(key) == val

    match = None
    if obj_store.get('items'):
        # TODO look up gapi for faster lookup
        for item in obj_store.get('items'):
            if name_ok(item) or id_ok(item) or key_ok(item):
                return item


def _clean_index(index_dims, parse_dates):
    _should_add = lambda lst: pd.Index(lst).isin(index_dims).all()
    to_remove = []
    to_add = []

    if isinstance(parse_dates, (list, tuple, np.ndarray)):
        for lst in parse_dates:
            if isinstance(lst, (list, tuple, np.ndarray)):
                if _should_add(lst):
                    to_add.append('_'.join(lst))
                to_remove.extend(lst)
    elif isinstance(parse_dates, dict):
        for name, lst in compat.iteritems(parse_dates):
            if isinstance(lst, (list, tuple, np.ndarray)):
                if _should_add(lst):
                    to_add.append(name)
                to_remove.extend(lst)

    index_dims = pd.Index(index_dims)
    to_remove = pd.Index(set(to_remove))
    to_add = pd.Index(set(to_add))

    return index_dims.difference(to_remove).union(to_add)


def _get_col_names(header_info):
    return [x['name'][3:] for x in header_info]


def _get_column_types(header_info):
    return [(x['name'][3:], x['columnType']) for x in header_info]


def _get_dim_names(header_info):
    return [x['name'][3:] for x in header_info
            if x['columnType'] == u('DIMENSION')]


def _get_met_names(header_info):
    return [x['name'][3:] for x in header_info
            if x['columnType'] == u('METRIC')]


def _get_data_types(header_info):
    return [(x['name'][3:], TYPE_MAP.get(x['dataType'], object))
            for x in header_info]
