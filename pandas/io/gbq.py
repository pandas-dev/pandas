"""
Pandas module to interface with Google BigQuery.
"""
import os
import sys
import tempfile
import csv
import logging
from datetime import datetime

import pandas as pd
import numpy as np

from pandas import DataFrame, concat
from pandas.core.common import PandasError

try:
    import bq
    import bigquery_client
    import gflags as flags
    _BQ_INSTALLED = True
except ImportError:
    _BQ_INSTALLED = False


# Setup the logger
logger = logging.getLogger('pandas.io.gbq')

# These are some custom exceptions that the
# to_gbq() method can throw

class SchemaMissing(PandasError,IOError):
    """
    Raised when attempting to write a DataFrame to
    a new table in Google BigQuery without specifying
    a schema describing the DataFrame.
    """
    pass

class InvalidSchema(PandasError,IOError):
    """
    Raised when attempting to write a DataFrame to
    Google BigQuery with an invalid table schema.
    """
    pass

class TableExistsFail(PandasError,IOError):
    """
    Raised when attempting to write a DataFrame to
    an existing Google BigQuery table without specifying
    that a replace/update action be taken.
    """
    pass

class InvalidColumnOrder(PandasError,IOError):
    """
    Raised when the provided column order for output
    results DataFrame does not match the schema
    returned by BigQuery.
    """
    pass


def _authenticate():
    """
    For testing, we abstract the authentication to BigQuery API.
    Presently this is implemented using the bq.py Client.Get()
    method. Any exceptions raised are considered fatal, so we
    do not process them.

    Returns
    -------
    BigqueryClient : Configured connection to Google BigQuery
    """
    return bq.Client.Get()

def _parse_entry(field_value, field_type):
    """
    Given a value and the corresponding BigQuery data type,
    perform any operations needed and return in a format
    appropriate for a numpy record dictionary

    Parameters
    ----------
    field_value : Source object to be transformed
    field_type : String representation of Google BigQuery
                 data type (per schema)

    Returns
    -------
    field_value : object or primitive of type corresponding
                  to field_type
    """

    # Avoid any casting problems
    if field_value is None or field_value == 'null':
        return None
    if field_type == 'INTEGER' or field_type == 'FLOAT':
        field_value = float(field_value)
    elif field_type == 'TIMESTAMP':
        timestamp = datetime.utcfromtimestamp(float(field_value))
        field_value = np.datetime64(timestamp)
    elif field_type == 'BOOLEAN':
        field_value = field_value == 'true'
    else:
        field_value = str(field_value)
    return field_value


def _parse_page(raw_page, col_names, col_types, col_dtypes):
    """
    Given a list of rows produced by the client.apiclient.tabledata().list(),
    build a numpy array with proper dtypes and column names as specified
    by the arguments.

    Parameters
    ----------
    raw_page : Resulting list of rows from a page retrieved via
        bigquery API
        client.apiclient.tabledata().list().execute()['rows']
    col_names: An ordered list of names for the columns
    col_types: String representation of the BigQuery DataType for that
        column
    col_dtypes: Target numpy.dtype for the column

    Returns
    -------
    page_array : numpy record array corresponding
        to the page data
    """

    # Should be at most 100,000 per the API, but this could
    # be increased in the future. Should only be less than
    # this for the last page to reduce API calls
    page_row_count = len(raw_page)

    # Place to hold the results for a page of data
    page_array = np.zeros(
                    (page_row_count,),
                    dtype=zip(col_names,col_dtypes)
                 )
    for row_num, raw_row in enumerate(raw_page):
        entries = raw_row.get('f', [])
         # Iterate over each entry - setting proper field types
        for col_num, field_type in enumerate(col_types):
            # Process the field's types using schema
            field_value = _parse_entry(entries[col_num].get('v', ''),
                                       field_type)
            # Fill the value into the final array
            page_array[row_num][col_num] = field_value

    return page_array

def _parse_data(client, job, index_col=None, col_order=None):
    """
    Iterate through the query results and piece together the
    final DataFrame. Builds a DataFrame for each page of
    results, then concatenates them together when finished.
    To save memory, we use numpy record arrays to build these
    DataFrames.

    Parameters
    ----------
    client: An instance of bq.Client
    job: An array containing the job info for a completed query
    index_col: str (optional)
        Name of result column to use for index in results DataFrame
    col_order: list() (optional)
        List of BigQuery column names in the desired order for results
        DataFrame

    Returns
    -------
    df: pandas DataFrame
        DataFrame representing results of query

    Raises:
    ------
    InvalidColumnOrder:
        Raised if 'col_order' parameter doesn't match returned DataFrame
    BigqueryError:
        Raised by bigquery_client if a Google API error is encountered


    Notes:
    -----
    This script relies on Google being consistent with their 
    pagination API. We are using the most flexible iteration method
    that we could find in the bq.py/bigquery_client.py API's, but 
    these have undergone large amounts of change recently.

    We have encountered bugs with this functionality, see:
    http://stackoverflow.com/questions/19145587/bq-py-not-paging-results
    """

    # dtype Map -
    # see: http://pandas.pydata.org/pandas-docs/dev/missing_data.html#missing-data-casting-rules-and-indexing
    dtype_map = {'INTEGER': np.dtype(float),
                 'FLOAT': np.dtype(float),
                 'TIMESTAMP': 'M8[ns]'}     # This seems to be buggy without nanosecond indicator

    # We first need the schema to get information about the columns of
    # our dataframe. 

    table_dict = job['configuration']['query']['destinationTable']
    fields = client.GetTableSchema(table_dict)['fields']

    # Get the schema into a format useable to create our
    # dataframe
    col_dtypes = []
    col_types = []
    col_names = []

    # TODO: Do this in one clean step
    for field in fields:
        col_types.append(field['type'])
        # Note the encoding... numpy doesn't like titles that are UTF8, which is the return
        # type from the API
        col_names.append(field['name'].encode('ascii', 'ignore'))
        # Note, it would be nice to use 'str' types, but BigQuery doesn't have a fixed length
        # in mind - just maxes out at 64k
        col_dtypes.append(dtype_map.get(field['type'],object))

  
    # How many columns are there
    num_columns = len(col_names)
    
    # Iterate over the result rows.
    # Since Google's API now requires pagination of results,
    # we do that here. The following is repurposed from 
    # bigquery_client.py  :: Client.ReadTableRows()

    # Initially, no page token is set
    page_token = None

    # Most of Google's client API's allow one to set total_rows in case
    # the user only wants the first 'n' results from a query. Typically
    # they set this to sys.maxint by default, but this caused problems
    # during testing - specifically on OS X. It appears that at some
    # point in bigquery_client.py, there is an attempt to cast this value
    # to an unsigned integer. Depending on the python install,
    # sys.maxint may exceed the limitations of unsigned integers.
    #
    # See:
    # https://code.google.com/p/google-bigquery-tools/issues/detail?id=14

    # This is hardcoded value for 32bit sys.maxint per
    # the above note. Theoretically, we could simply use
    # 100,000 (or whatever the current max page size is),
    # but this is more flexible in the event of an API change
    total_rows = 2147483647

    # Keep track of rows read
    row_count = 0

    # Keep our page DataFrames until the end when we 
    # concatentate them
    dataframe_list = list()

    # Iterate over all rows
    while row_count < total_rows:
        data = client.apiclient.tabledata().list(maxResults=total_rows - row_count,
                                                 pageToken=page_token,
                                                 **table_dict).execute()

        # If there are more results than will fit on a page, 
        # you will recieve a token for the next page.
        page_token = data.get('pageToken', None)

        # How many rows are there across all pages?
        total_rows = min(total_rows, int(data['totalRows'])) # Changed to use get(data[rows],0)
        raw_page = data.get('rows', [])
        page_array = _parse_page(raw_page, col_names, col_types, col_dtypes)

        row_count += len(page_array)
        if total_rows > 0:
            completed = (100 * row_count) / total_rows
            logger.info('Remaining Rows: ' + str(total_rows - row_count) + '(' + str(completed) + '% Complete)') 
        else:
            logger.info('No Rows')

        dataframe_list.append(DataFrame(page_array))

        # Handle any exceptions that might have occured
        if not page_token and row_count != total_rows:
            raise bigquery_client.BigqueryInterfaceError(
                'PageToken missing for %r' % (
                    bigquery_client.ApiClientHelper.TableReference.Create(**table_dict),))
        if not raw_page and row_count != total_rows:
            raise bigquery_client.BigqueryInterfaceError(
                'Not enough rows returned by server for %r' % (
                    bigquery_client.ApiClientHelper.TableReference.Create(**table_dict),))

    # Build final dataframe
    final_df = concat(dataframe_list, ignore_index=True)

    # Reindex the DataFrame on the provided column
    if index_col is not None:
        if index_col in col_names:
            final_df.set_index(index_col, inplace=True)
            col_names.remove(index_col)
        else:
            raise InvalidColumnOrder('Index column "{0}" does not exist in DataFrame.'.format(index_col))

    # Change the order of columns in the DataFrame based on provided list
    if col_order is not None:
        if sorted(col_order) == sorted(col_names):
            final_df = final_df[col_order]
        else:
            raise InvalidColumnOrder('Column order does not match this DataFrame.')

    # Downcast floats to integers and objects to booleans
    # if there are no NaN's. This is presently due to a
    # limitation of numpy in handling missing data.
    final_df._data = final_df._data.downcast(dtypes='infer')
    return final_df

def to_gbq(dataframe, destination_table, schema=None, col_order=None, if_exists='fail', **kwargs):
    """Write a DataFrame to a Google BigQuery table. 
    
    If the table exists, the DataFrame will be appended. If not, a new table 
    will be created, in which case the schema will have to be specified. By default,
    rows will be written in the order they appear in the DataFrame, though
    the user may specify an alternative order.

    Parameters
    ----------
    dataframe : DataFrame
        DataFrame to be written
    destination_table : string
         name of table to be written, in the form 'dataset.tablename'
    schema : sequence (optional)
         list of column types in order for data to be inserted, e.g. ['INTEGER', 'TIMESTAMP', 'BOOLEAN']
    col_order : sequence (optional)
         order which columns are to be inserted, e.g. ['primary_key', 'birthday', 'username']
    if_exists : {'fail', 'replace', 'append'} (optional)
        - fail: If table exists, do nothing.
        - replace: If table exists, drop it, recreate it, and insert data.
        - append: If table exists, insert data. Create if does not exist.
    kwargs are passed to the Client constructor

    Raises
    ------
    SchemaMissing :
        Raised if the 'if_exists' parameter is set to 'replace', but no schema is specified
    TableExists :
        Raised if the specified 'destination_table' exists but the 'if_exists' parameter is set to 'fail' (the default)
    InvalidSchema :
        Raised if the 'schema' parameter does not match the provided DataFrame
    """

    if not _BQ_INSTALLED:
        if sys.version_info >= (3, 0):
            raise NotImplementedError('gbq module does not support Python 3 yet')
        else:
            raise ImportError('Could not import Google BigQuery Client.')

    ALLOWED_TYPES = ['STRING', 'INTEGER', 'FLOAT', 'BOOLEAN', 'TIMESTAMP', 'RECORD']

    if if_exists == 'replace' and schema is None:
        raise SchemaMissing('Cannot replace a table without specifying the data schema')
    else:
        client = _authenticate()
        table_reference = client.GetTableReference(destination_table)
        if client.TableExists(table_reference):
            if if_exists == 'fail':
                raise TableExistsFail('Cannot overwrite existing tables if \'if_exists="fail"\'')
            else:
                # Build up a string representation of the 
                # table's schema. Since the table already
                # exists, we ask ask the API for it, which
                # is returned in a list of dictionaries
                # describing column data. Iterate over these
                # and build up a string of form:
                # "col_name1 : col_type1, col_name2 : col_type2..."
                schema_full = client.GetTableSchema(dict(table_reference))['fields']
                schema = ''
                for count, row in enumerate(schema_full):
                    if count > 0:
                        schema += ', '
                    schema += row['name'] + ':' + row['type']
        else:
            logger.info('Creating New Table')
            if schema is None:
                raise SchemaMissing('Cannot create a new table without specifying the data schema')
            else:
                columns = dataframe.columns
                if len(schema) != len(columns):
                    raise InvalidSchema('Incorrect number of columns in schema')
                else:
                    schema_string = ''
                    for count, name in enumerate(columns):
                        if count > 0:
                            schema_string += ', '
                        column_type = schema[count].upper()
                        if column_type in ALLOWED_TYPES:
                            schema_string += name + ':' + schema[count].lower()
                        else:
                            raise InvalidSchema('Invalid Type: ' + column_type + ". Must be one of: " + str(ALLOWED_TYPES))
                    schema = schema_string

    opts = kwargs
    opts['sync'] = True
    opts['skip_leading_rows'] = 1
    opts['encoding'] = 'UTF-8'
    opts['max_bad_records'] = 0

    # See: https://developers.google.com/bigquery/docs/reference/v2/jobs
    if if_exists == 'replace':
        opts['write_disposition'] = 'WRITE_TRUNCATE'
    elif if_exists == 'append':
        opts['write_disposition'] = 'WRITE_APPEND'

    with tempfile.NamedTemporaryFile() as csv_file:
        dataframe.to_csv(csv_file.name, index=False, encoding='utf-8')
        job = client.Load(table_reference, csv_file.name, schema=schema, **opts)

def read_gbq(query, project_id = None, destination_table = None, index_col=None, col_order=None, **kwargs):
    """Load data from Google BigQuery.
    
    The main method a user calls to load data from Google BigQuery into a pandas DataFrame. 
    This is a simple wrapper for Google's bq.py and bigquery_client.py, which we use
    to get the source data. Because of this, this script respects the user's bq settings 
    file, '~/.bigqueryrc', if it exists. Such a file can be generated using 'bq init'. Further,
    additional parameters for the query can be specified as either ``**kwds`` in the command,
    or using FLAGS provided in the 'gflags' module. Particular options can be found in
    bigquery_client.py.

    Parameters
    ----------
    query : str
        SQL-Like Query to return data values
    project_id : str (optional)
        Google BigQuery Account project ID. Optional, since it may be
        located in ~/.bigqueryrc
    index_col : str (optional)
        Name of result column to use for index in results DataFrame
    col_order : list(str) (optional)
        List of BigQuery column names in the desired order for results
        DataFrame
    destination_table : string (optional)
        If provided, send the results to the given table.
    **kwargs : 
        To be passed to bq.Client.Create(). Particularly: 'trace', 
        'sync', 'api', 'api_version'

    Returns
    -------
    df: DataFrame
        DataFrame representing results of query

    """
    if not _BQ_INSTALLED:
        if sys.version_info >= (3, 0):
            raise NotImplementedError('gbq module does not support Python 3 yet')
        else:
            raise ImportError('Could not import Google BigQuery Client.')

    query_args = kwargs
    query_args['project_id'] = project_id
    query_args['query'] = query
    query_args['destination_table'] = destination_table
    query_args['sync'] = True

    client = _authenticate()

    job = client.Query(**query_args)
    
    return _parse_data(client, job, index_col=index_col, col_order=col_order)
