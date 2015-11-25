import warnings
from datetime import datetime
import json
import logging
from time import sleep
import uuid

import numpy as np

from distutils.version import StrictVersion
from pandas import compat
from pandas.core.api import DataFrame
from pandas.tools.merge import concat
from pandas.core.common import PandasError
from pandas.util.decorators import deprecate
from pandas.compat import lzip, bytes_to_str

def _check_google_client_version():

    try:
        import pkg_resources

    except ImportError:
        raise ImportError('Could not import pkg_resources (setuptools).')

    if compat.PY3:
        google_api_minimum_version = '1.4.1'
    else:
        google_api_minimum_version = '1.2.0'

    _GOOGLE_API_CLIENT_VERSION = pkg_resources.get_distribution('google-api-python-client').version

    if StrictVersion(_GOOGLE_API_CLIENT_VERSION) < StrictVersion(google_api_minimum_version):
        raise ImportError("pandas requires google-api-python-client >= {0} for Google BigQuery support, "
                          "current version {1}".format(google_api_minimum_version, _GOOGLE_API_CLIENT_VERSION))

logger = logging.getLogger('pandas.io.gbq')
logger.setLevel(logging.ERROR)


class AccessDenied(PandasError, ValueError):
    """
    Raised when invalid credentials are provided, or tokens have expired.
    """
    pass


class DatasetCreationError(PandasError, ValueError):
    """
    Raised when the create dataset method fails
    """
    pass


class GenericGBQException(PandasError, ValueError):
    """
    Raised when an unrecognized Google API Error occurs.
    """
    pass


class InvalidColumnOrder(PandasError, ValueError):
    """
    Raised when the provided column order for output
    results DataFrame does not match the schema
    returned by BigQuery.
    """
    pass


class InvalidPageToken(PandasError, ValueError):
    """
    Raised when Google BigQuery fails to return,
    or returns a duplicate page token.
    """
    pass


class InvalidSchema(PandasError, ValueError):
    """
    Raised when the provided DataFrame does
    not match the schema of the destination
    table in BigQuery.
    """
    pass


class NotFoundException(PandasError, ValueError):
    """
    Raised when the project_id, table or dataset provided in the query could not be found.
    """
    pass


class StreamingInsertError(PandasError, ValueError):
    """
    Raised when BigQuery reports a streaming insert error.
    For more information see `Streaming Data Into BigQuery
    <https://cloud.google.com/bigquery/streaming-data-into-bigquery>`__
    """


class TableCreationError(PandasError, ValueError):
    """
    Raised when the create table method fails
    """
    pass


class GbqConnector(object):

    def __init__(self, project_id, reauth=False):
        self.test_google_api_imports()
        self.project_id = project_id
        self.reauth = reauth
        self.credentials = self.get_credentials()
        self.service = self.get_service(self.credentials)

    def test_google_api_imports(self):
        try:
            import httplib2
            from apiclient.discovery import build
            from apiclient.errors import HttpError
            from oauth2client.client import AccessTokenRefreshError
            from oauth2client.client import OAuth2WebServerFlow
            from oauth2client.file import Storage
            from oauth2client.tools import run_flow, argparser
        except ImportError as e:
            raise ImportError("Missing module required for Google BigQuery support: {0}".format(str(e)))

    def get_credentials(self):
        from oauth2client.client import OAuth2WebServerFlow
        from oauth2client.file import Storage
        from oauth2client.tools import run_flow, argparser

        _check_google_client_version()

        flow = OAuth2WebServerFlow(client_id='495642085510-k0tmvj2m941jhre2nbqka17vqpjfddtd.apps.googleusercontent.com',
                                   client_secret='kOc9wMptUtxkcIFbtZCcrEAc',
                                   scope='https://www.googleapis.com/auth/bigquery',
                                   redirect_uri='urn:ietf:wg:oauth:2.0:oob')

        storage = Storage('bigquery_credentials.dat')
        credentials = storage.get()

        if credentials is None or credentials.invalid or self.reauth:
            credentials = run_flow(flow, storage, argparser.parse_args([]))

        return credentials

    @staticmethod
    def get_service(credentials):
        import httplib2
        from apiclient.discovery import build

        _check_google_client_version()

        http = httplib2.Http()
        http = credentials.authorize(http)
        bigquery_service = build('bigquery', 'v2', http=http)

        return bigquery_service

    @staticmethod
    def process_http_error(ex):
        # See `BigQuery Troubleshooting Errors <https://cloud.google.com/bigquery/troubleshooting-errors>`__

        status = json.loads(bytes_to_str(ex.content))['error']
        errors = status.get('errors', None)

        if errors:
            for error in errors:
                reason = error['reason']
                message = error['message']

                raise GenericGBQException("Reason: {0}, Message: {1}".format(reason, message))

        raise GenericGBQException(errors)

    @staticmethod
    def process_insert_errors(insert_errors, verbose):
        for insert_error in insert_errors:
            row = insert_error['index']
            errors = insert_error.get('errors', None)
            for error in errors:
                reason = error['reason']
                message = error['message']
                location = error['location']
                error_message = 'Error at Row: {0}, Reason: {1}, Location: {2}, Message: {3}'.format(row, reason, location, message)

                # Report all error messages if verbose is set
                if verbose:
                    print(error_message)
                else:
                    raise StreamingInsertError(error_message + '\nEnable verbose logging to see all errors')

        raise StreamingInsertError

    def run_query(self, query, verbose=True):
        from apiclient.errors import HttpError
        from oauth2client.client import AccessTokenRefreshError

        _check_google_client_version()

        job_collection = self.service.jobs()
        job_data = {
            'configuration': {
                'query': {
                    'query': query
                    # 'allowLargeResults', 'createDisposition', 'preserveNulls', destinationTable, useQueryCache
                }
            }
        }

        try:
            query_reply = job_collection.insert(projectId=self.project_id, body=job_data).execute()
        except AccessTokenRefreshError:
            raise AccessDenied("The credentials have been revoked or expired, please re-run the application "
                               "to re-authorize")
        except HttpError as ex:
            self.process_http_error(ex)

        job_reference = query_reply['jobReference']

        while not query_reply.get('jobComplete', False):
            if verbose:
                print('Waiting for job to complete...')
            try:
                query_reply = job_collection.getQueryResults(projectId=job_reference['projectId'],
                                                             jobId=job_reference['jobId']).execute()
            except HttpError as ex:
                self.process_http_error(ex)

        total_rows = int(query_reply['totalRows'])
        result_pages = list()
        seen_page_tokens = list()
        current_row = 0
        # Only read schema on first page
        schema = query_reply['schema']

        # Loop through each page of data
        while 'rows' in query_reply and current_row < total_rows:
            page = query_reply['rows']
            result_pages.append(page)
            current_row += len(page)
            page_token = query_reply.get('pageToken', None)

            if not page_token and current_row < total_rows:
                raise InvalidPageToken(
                    "Required pageToken was missing. Received {0} of {1} rows".format(current_row,
                                                                                      total_rows))

            elif page_token in seen_page_tokens:
                raise InvalidPageToken("A duplicate pageToken was returned")

            seen_page_tokens.append(page_token)

            try:
                query_reply = job_collection.getQueryResults(
                                projectId=job_reference['projectId'],
                                jobId=job_reference['jobId'],
                                pageToken=page_token).execute()
            except HttpError as ex:
                self.process_http_error(ex)

        if current_row < total_rows:
            raise InvalidPageToken()

        return schema, result_pages

    def load_data(self, dataframe, dataset_id, table_id, chunksize, verbose):
        from apiclient.errors import HttpError

        job_id = uuid.uuid4().hex
        rows = []
        remaining_rows = len(dataframe)

        if verbose:
            total_rows = remaining_rows
            print("\n\n")

        for index, row in dataframe.reset_index(drop=True).iterrows():
            row_dict = dict()
            row_dict['json'] = json.loads(row.to_json(force_ascii=False,
                                                      date_unit='s',
                                                      date_format='iso'))
            row_dict['insertId'] = job_id + str(index)
            rows.append(row_dict)
            remaining_rows -= 1

            if (len(rows) % chunksize == 0) or (remaining_rows == 0):
                if verbose:
                    print("\rStreaming Insert is {0}% Complete".format(((total_rows - remaining_rows) * 100) / total_rows))

                body = {'rows': rows}

                try:
                    response = self.service.tabledata().insertAll(
                                                            projectId = self.project_id,
                                                            datasetId = dataset_id,
                                                            tableId = table_id,
                                                            body = body).execute()
                except HttpError as ex:
                    self.process_http_error(ex)

                # For streaming inserts, even if you receive a success HTTP response code, you'll need to check the
                # insertErrors property of the response to determine if the row insertions were successful, because
                # it's possible that BigQuery was only partially successful at inserting the rows.
                # See the `Success HTTP Response Codes <https://cloud.google.com/bigquery/streaming-data-into-bigquery#troubleshooting>`__
                # section

                insert_errors = response.get('insertErrors', None)
                if insert_errors:
                    self.process_insert_errors(insert_errors, verbose)

                sleep(1)  # Maintains the inserts "per second" rate per API
                rows = []

        if verbose:
            print("\n")

    def verify_schema(self, dataset_id, table_id, schema):
        from apiclient.errors import HttpError

        try:
            return (self.service.tables().get(
                projectId=self.project_id,
                datasetId=dataset_id,
                tableId=table_id
            ).execute()['schema']) == schema

        except HttpError as ex:
            self.process_http_error(ex)

    def delete_and_recreate_table(self, dataset_id, table_id, table_schema, verbose):
        delay = 0

        # Changes to table schema may take up to 2 minutes as of May 2015
        # See `Issue 191 <https://code.google.com/p/google-bigquery/issues/detail?id=191>`__
        # Compare previous schema with new schema to determine if there should be a 120 second delay

        if not self.verify_schema(dataset_id, table_id, table_schema):
            if verbose:
                print('The existing table has a different schema. Please wait 2 minutes. See Google BigQuery issue #191')
            delay = 120

        table = _Table(self.project_id, dataset_id)
        table.delete(table_id)
        table.create(table_id, table_schema)
        sleep(delay)


def _parse_data(schema, rows):
    # see: http://pandas.pydata.org/pandas-docs/dev/missing_data.html#missing-data-casting-rules-and-indexing
    dtype_map = {'INTEGER': np.dtype(float),
                 'FLOAT': np.dtype(float),
                 'TIMESTAMP': 'M8[ns]'}     # This seems to be buggy without nanosecond indicator

    fields = schema['fields']
    col_types = [field['type'] for field in fields]
    col_names = [str(field['name']) for field in fields]
    col_dtypes = [dtype_map.get(field['type'], object) for field in fields]
    page_array = np.zeros((len(rows),),
                          dtype=lzip(col_names, col_dtypes))

    for row_num, raw_row in enumerate(rows):
        entries = raw_row.get('f', [])
        for col_num, field_type in enumerate(col_types):
            field_value = _parse_entry(entries[col_num].get('v', ''),
                                       field_type)
            page_array[row_num][col_num] = field_value

    return DataFrame(page_array, columns=col_names)


def _parse_entry(field_value, field_type):
    if field_value is None or field_value == 'null':
        return None
    if field_type == 'INTEGER' or field_type == 'FLOAT':
        return float(field_value)
    elif field_type == 'TIMESTAMP':
        timestamp = datetime.utcfromtimestamp(float(field_value))
        return np.datetime64(timestamp)
    elif field_type == 'BOOLEAN':
        return field_value == 'true'
    return field_value


def read_gbq(query, project_id=None, index_col=None, col_order=None, reauth=False, verbose=True):
    """Load data from Google BigQuery.

    THIS IS AN EXPERIMENTAL LIBRARY

    The main method a user calls to execute a Query in Google BigQuery and read results
    into a pandas DataFrame using the v2 Google API client for Python.  Documentation for
    the API is available at https://developers.google.com/api-client-library/python/.
    Authentication to the Google BigQuery service is via OAuth 2.0 using the product name
    'pandas GBQ'.

    Parameters
    ----------
    query : str
        SQL-Like Query to return data values
    project_id : str
        Google BigQuery Account project ID.
    index_col : str (optional)
        Name of result column to use for index in results DataFrame
    col_order : list(str) (optional)
        List of BigQuery column names in the desired order for results
        DataFrame
    reauth : boolean (default False)
        Force Google BigQuery to reauthenticate the user. This is useful
        if multiple accounts are used.
    verbose : boolean (default True)
        Verbose output

    Returns
    -------
    df: DataFrame
        DataFrame representing results of query

    """

    if not project_id:
        raise TypeError("Missing required parameter: project_id")

    connector = GbqConnector(project_id, reauth=reauth)
    schema, pages = connector.run_query(query, verbose=verbose)
    dataframe_list = []
    while len(pages) > 0:
        page = pages.pop()
        dataframe_list.append(_parse_data(schema, page))

    if len(dataframe_list) > 0:
        final_df = concat(dataframe_list, ignore_index=True)
    else:
        final_df = _parse_data(schema, [])

    # Reindex the DataFrame on the provided column
    if index_col is not None:
        if index_col in final_df.columns:
            final_df.set_index(index_col, inplace=True)
        else:
            raise InvalidColumnOrder(
                'Index column "{0}" does not exist in DataFrame.'
                .format(index_col)
            )

    # Change the order of columns in the DataFrame based on provided list
    if col_order is not None:
        if sorted(col_order) == sorted(final_df.columns):
            final_df = final_df[col_order]
        else:
            raise InvalidColumnOrder(
                'Column order does not match this DataFrame.'
            )

    # Downcast floats to integers and objects to booleans
    # if there are no NaN's. This is presently due to a
    # limitation of numpy in handling missing data.
    final_df._data = final_df._data.downcast(dtypes='infer')
    return final_df


def to_gbq(dataframe, destination_table, project_id, chunksize=10000,
           verbose=True, reauth=False, if_exists='fail'):
    """Write a DataFrame to a Google BigQuery table.

    THIS IS AN EXPERIMENTAL LIBRARY

    Parameters
    ----------
    dataframe : DataFrame
        DataFrame to be written
    destination_table : string
        Name of table to be written, in the form 'dataset.tablename'
    project_id : str
        Google BigQuery Account project ID.
    chunksize : int (default 10000)
        Number of rows to be inserted in each chunk from the dataframe.
    verbose : boolean (default True)
        Show percentage complete
    reauth : boolean (default False)
        Force Google BigQuery to reauthenticate the user. This is useful
        if multiple accounts are used.
    if_exists : {'fail', 'replace', 'append'}, default 'fail'
        'fail': If table exists, do nothing.
        'replace': If table exists, drop it, recreate it, and insert data.
        'append': If table exists, insert data. Create if does not exist.
    """

    if if_exists not in ('fail', 'replace', 'append'):
        raise ValueError("'{0}' is not valid for if_exists".format(if_exists))

    if '.' not in destination_table:
        raise NotFoundException("Invalid Table Name. Should be of the form 'datasetId.tableId' ")

    connector = GbqConnector(project_id, reauth=reauth)
    dataset_id, table_id = destination_table.rsplit('.', 1)

    table = _Table(project_id, dataset_id, reauth=reauth)

    table_schema = _generate_bq_schema(dataframe)

    # If table exists, check if_exists parameter
    if table.exists(table_id):
        if if_exists == 'fail':
            raise TableCreationError("Could not create the table because it already exists. "
                                     "Change the if_exists parameter to append or replace data.")
        elif if_exists == 'replace':
            connector.delete_and_recreate_table(dataset_id, table_id, table_schema, verbose)
        elif if_exists == 'append':
            if not connector.verify_schema(dataset_id, table_id, table_schema):
                raise InvalidSchema("Please verify that the column order, structure and data types in the DataFrame "
                                    "match the schema of the destination table.")
    else:
        table.create(table_id, table_schema)

    connector.load_data(dataframe, dataset_id, table_id, chunksize, verbose)


def generate_bq_schema(df, default_type='STRING'):

    # deprecation TimeSeries, #11121
    warnings.warn("generate_bq_schema is deprecated and will be removed in a future version",
                  FutureWarning, stacklevel=2)

    return _generate_bq_schema(df, default_type=default_type)

def _generate_bq_schema(df, default_type='STRING'):
    """ Given a passed df, generate the associated Google BigQuery schema.

    Parameters
    ----------
    df : DataFrame
    default_type : string
        The default big query type in case the type of the column
        does not exist in the schema.
    """

    type_mapping = {
        'i': 'INTEGER',
        'b': 'BOOLEAN',
        'f': 'FLOAT',
        'O': 'STRING',
        'S': 'STRING',
        'U': 'STRING',
        'M': 'TIMESTAMP'
    }

    fields = []
    for column_name, dtype in df.dtypes.iteritems():
        fields.append({'name': column_name,
                       'type': type_mapping.get(dtype.kind, default_type)})

    return {'fields': fields}

class _Table(GbqConnector):

    def __init__(self, project_id, dataset_id, reauth=False):
        from apiclient.errors import HttpError
        self.test_google_api_imports()
        self.project_id = project_id
        self.reauth = reauth
        self.credentials = self.get_credentials()
        self.service = self.get_service(self.credentials)
        self.http_error = HttpError
        self.dataset_id = dataset_id

    def exists(self, table_id):
        """ Check if a table exists in Google BigQuery

        .. versionadded:: 0.17.0

        Parameters
        ----------
        table : str
            Name of table to be verified

        Returns
        -------
        boolean
            true if table exists, otherwise false
        """

        try:
            self.service.tables().get(
                 projectId=self.project_id,
                 datasetId=self.dataset_id,
                 tableId=table_id).execute()
            return True
        except self.http_error as ex:
            if ex.resp.status == 404:
                return False
            else:
                self.process_http_error(ex)

    def create(self, table_id, schema):
        """ Create a table in Google BigQuery given a table and schema

        .. versionadded:: 0.17.0

        Parameters
        ----------
        table : str
            Name of table to be written
        schema : str
            Use the generate_bq_schema to generate your table schema from a dataframe.
        """

        if self.exists(table_id):
            raise TableCreationError("The table could not be created because it already exists")

        if not _Dataset(self.project_id).exists(self.dataset_id):
            _Dataset(self.project_id).create(self.dataset_id)

        body = {
            'schema': schema,
            'tableReference': {
                'tableId': table_id,
                'projectId': self.project_id,
                'datasetId': self.dataset_id
            }
        }

        try:
            self.service.tables().insert(
                 projectId=self.project_id,
                 datasetId=self.dataset_id,
                 body=body).execute()
        except self.http_error as ex:
            self.process_http_error(ex)

    def delete(self, table_id):
        """ Delete a table in Google BigQuery

        .. versionadded:: 0.17.0

        Parameters
        ----------
        table : str
            Name of table to be deleted
        """

        if not self.exists(table_id):
            raise NotFoundException("Table does not exist")

        try:
            self.service.tables().delete(
                 datasetId=self.dataset_id,
                 projectId=self.project_id,
                 tableId=table_id).execute()
        except self.http_error as ex:
            self.process_http_error(ex)


class _Dataset(GbqConnector):

    def __init__(self, project_id, reauth=False):
        from apiclient.errors import HttpError
        self.test_google_api_imports()
        self.project_id = project_id
        self.reauth = reauth
        self.credentials = self.get_credentials()
        self.service = self.get_service(self.credentials)
        self.http_error = HttpError

    def exists(self, dataset_id):
        """ Check if a dataset exists in Google BigQuery

        .. versionadded:: 0.17.0

        Parameters
        ----------
        dataset_id : str
            Name of dataset to be verified

        Returns
        -------
        boolean
            true if dataset exists, otherwise false
        """

        try:
            self.service.datasets().get(
                 projectId=self.project_id,
                 datasetId=dataset_id).execute()
            return True
        except self.http_error as ex:
            if ex.resp.status == 404:
                return False
            else:
                self.process_http_error(ex)

    def datasets(self):
        """ Return a list of datasets in Google BigQuery

        .. versionadded:: 0.17.0

        Parameters
        ----------
        None

        Returns
        -------
        list
            List of datasets under the specific project
        """

        try:
            list_dataset_response = self.service.datasets().list(
                                         projectId=self.project_id).execute().get('datasets', None)

            if not list_dataset_response:
                return []

            dataset_list = list()

            for row_num, raw_row in enumerate(list_dataset_response):
                dataset_list.append(raw_row['datasetReference']['datasetId'])

            return dataset_list
        except self.http_error as ex:
            self.process_http_error(ex)

    def create(self, dataset_id):
        """ Create a dataset in Google BigQuery

        .. versionadded:: 0.17.0

        Parameters
        ----------
        dataset : str
            Name of dataset to be written
        """

        if self.exists(dataset_id):
            raise DatasetCreationError("The dataset could not be created because it already exists")

        body = {
            'datasetReference': {
                'projectId': self.project_id,
                'datasetId': dataset_id
            }
        }

        try:
            self.service.datasets().insert(
                 projectId=self.project_id,
                 body=body).execute()
        except self.http_error as ex:
            self.process_http_error(ex)

    def delete(self, dataset_id):
        """ Delete a dataset in Google BigQuery

        .. versionadded:: 0.17.0

        Parameters
        ----------
        dataset : str
            Name of dataset to be deleted
        """

        if not self.exists(dataset_id):
            raise NotFoundException("Dataset {0} does not exist".format(dataset_id))

        try:
            self.service.datasets().delete(
                 datasetId=dataset_id,
                 projectId=self.project_id).execute()

        except self.http_error as ex:
            self.process_http_error(ex)

    def tables(self, dataset_id):
        """ List tables in the specific dataset in Google BigQuery

        .. versionadded:: 0.17.0

        Parameters
        ----------
        dataset : str
            Name of dataset to list tables for

        Returns
        -------
        list
            List of tables under the specific dataset
        """

        try:
            list_table_response = self.service.tables().list(
                                       projectId=self.project_id,
                                       datasetId=dataset_id).execute().get('tables', None)

            if not list_table_response:
                return []

            table_list = list()

            for row_num, raw_row in enumerate(list_table_response):
                table_list.append(raw_row['tableReference']['tableId'])

            return table_list
        except self.http_error as ex:
            self.process_http_error(ex)
