""" Google BigQuery support """
import datetime


def _try_import():
    # since pandas is a dependency of pandas-gbq
    # we need to import on first use
    try:
        import pandas_gbq
    except ImportError:

        # give a nice error message
        raise ImportError("Load data from Google BigQuery\n"
                          "\n"
                          "the pandas-gbq package is not installed\n"
                          "see the docs: https://pandas-gbq.readthedocs.io\n"
                          "\n"
                          "you can install via pip or conda:\n"
                          "pip install pandas-gbq\n"
                          "conda install pandas-gbq -c conda-forge\n")

    return pandas_gbq


def read_gbq(query, project_id=None, index_col=None, col_order=None,
             reauth=False, verbose=True, private_key=None, dialect='legacy',
             allow_large_results=False, query_dataset='query_dataset',
             query_tableid=None,
             **kwargs):
    r"""Load data from Google BigQuery.

    The main method a user calls to execute a Query in Google BigQuery
    and read results into a pandas DataFrame.

    This function requires the `pandas-gbq package
    <https://pandas-gbq.readthedocs.io>`__.

    Authentication to the Google BigQuery service is via OAuth 2.0.

    - If "private_key" is not provided:

      By default "application default credentials" are used.

      If default application credentials are not found or are restrictive,
      user account credentials are used. In this case, you will be asked to
      grant permissions for product name 'pandas GBQ'.

    - If "private_key" is provided:

      Service account credentials will be used to authenticate.

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
    private_key : str (optional)
        Service account private key in JSON format. Can be file path
        or string contents. This is useful for remote server
        authentication (eg. jupyter iPython notebook on remote host)

    dialect : {'legacy', 'standard'}, default 'legacy'
        'legacy' : Use BigQuery's legacy SQL dialect.
        'standard' : Use BigQuery's standard SQL, which is
        compliant with the SQL 2011 standard. For more information
        see `BigQuery SQL Reference
        <https://cloud.google.com/bigquery/sql-reference/>`__

    allow_large_results : boolean (default False)
        Allows large queries greater than quota limit - Use when a GenericGBQException
        error is thrown with "Reason: responseTooLarge" 
        See: https://cloud.google.com/bigquery/docs/writing-results#large-results
        
    query_dataset : str (optional)
        Google BigQuery dataset to which the results of a large query will
        be saved
    query_tableid : str (optional)
        Google BigQuery tableid to which the results will be saved

    `**kwargs` : Arbitrary keyword arguments
        configuration (dict): query config parameters for job processing.
        For example:

            configuration = {'query': {'useQueryCache': False}}

        For more information see `BigQuery SQL Reference
        <https://cloud.google.com/bigquery/docs/reference/rest/v2/jobs#configuration.query>`__

    Returns
    -------
    df: DataFrame
        DataFrame representing results of query

    """
    pandas_gbq = _try_import()
    kwargs = update_read_gbq_kwargs(kwargs, project_id, allow_large_results, query_dataset, query_tableid)
    return pandas_gbq.read_gbq(
        query, project_id=project_id,
        index_col=index_col, col_order=col_order,
        reauth=reauth, verbose=verbose,
        private_key=private_key,
        dialect=dialect,
        **kwargs)


def to_gbq(dataframe, destination_table, project_id, chunksize=10000,
           verbose=True, reauth=False, if_exists='fail', private_key=None):
    pandas_gbq = _try_import()
    pandas_gbq.to_gbq(dataframe, destination_table, project_id,
                      chunksize=chunksize,
                      verbose=verbose, reauth=reauth,
                      if_exists=if_exists, private_key=private_key)


def update_read_gbq_kwargs(read_gbq_kwargs, project_id, allow_large_results, query_dataset,
                            query_tableid, write_disposition="WRITE_TRUNCATE"):
    """
    
    Parameters
    ----------
    read_gbq_kwargs: dict
        query config parameters for job processing passed to the read_gbq function 
        For example:

            configuration = {'query': {'useQueryCache': False}}

        For more information see `BigQuery SQL Reference
        <https://cloud.google.com/bigquery/docs/reference/rest/v2/jobs#configuration.query>`__
    project_id : str
        Google BigQuery Account project ID.    
    allow_large_results : boolean (default False)
        Allows large queries greater than quota limit - Use when a GenericGBQException
        error is thrown with "Reason: responseTooLarge" 
        See: https://cloud.google.com/bigquery/docs/writing-results#large-results
    intermediate_dataset : str
        Google BigQuery dataset to which the results of a large query will
        be saved
    query_tableid : str
        Google BigQuery tableid to which the results will be saved
    write_disposition : str
        Use "WRITE_TRUNCATE" to overwrite old intermediate BigQuery query data



    Returns
    -------
    updated kwargs

    """

    # Create tableId if left as None
    if query_tableid is None:
        # Generic name with timestamp unique down to the microsecond:
        query_tableid = 'intermediate_query_results_' + datetime.datetime.utcnow().strftime("%Y%M%d-%H%m-%f")

    # New configuration for allowing large queries:
    updated_query_config = {'query': {
        'allowLargeResults': allow_large_results,
        'destinationTable': {
            'projectId': project_id,
            'datasetId': query_dataset,
            'tableId': query_tableid
        },
        'writeDisposition': write_disposition
    }}

    # Append to predefined configuration:
    if 'configuration' not in read_gbq_kwargs:
        read_gbq_kwargs['configuration'] = updated_query_config
    else:
        # Append new configuration to user prescribed query configuration:
        read_gbq_kwargs['configuration']
        if 'query' not in read_gbq_kwargs['configuration']:
            read_gbq_kwargs['configuration']['query'] = updated_query_config['query']
        else:
            for updated_key in updated_query_config:
                if updated_key not in read_gbq_kwargs['configuration']['query']:
                    read_gbq_kwargs['configuration']['query'][updated_key] = updated_query_config['query'][updated_key]

    return read_gbq_kwargs
