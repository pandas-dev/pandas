""" Google BigQuery support """


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


def read_gbq(
        query, project_id=None, index_col=None, col_order=None, reauth=False,
        verbose=True, private_key=None, dialect='legacy', **kwargs):
    """
    Load data from Google BigQuery.

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
        SQL-Like Query to return data values.
    project_id : str
        Google BigQuery Account project ID.
    index_col : str (optional)
        Name of result column to use for index in results DataFrame.
    col_order : list(str) (optional)
        List of BigQuery column names in the desired order for results
        DataFrame.
    reauth : boolean (default False)
        Force Google BigQuery to reauthenticate the user. This is useful
        if multiple accounts are used.
    verbose : boolean (default True)
        Verbose output.
    private_key : str (optional)
        Service account private key in JSON format. Can be file path
        or string contents. This is useful for remote server
        authentication (eg. Jupyter/IPython notebook on remote host).
    dialect : {'legacy', 'standard'}, default 'legacy'
        SQL syntax dialect to use.
        'legacy' : Use BigQuery's legacy SQL dialect.
        'standard' : Use BigQuery's standard SQL, which is
        compliant with the SQL 2011 standard. For more information
        see `BigQuery SQL Reference
        <https://cloud.google.com/bigquery/sql-reference/>`__.
    kwargs : dict
        Arbitrary keyword arguments.
        configuration (dict): query config parameters for job processing.
        For example:

            configuration = {'query': {'useQueryCache': False}}

        For more information see `BigQuery SQL Reference
        <https://cloud.google.com/bigquery/docs/reference/rest/v2/jobs#configuration.query>`__

    Returns
    -------
    df: DataFrame
        DataFrame representing results of query.

    See Also
    --------
    pandas_gbq.read_gbq
    """
    pandas_gbq = _try_import()
    return pandas_gbq.read_gbq(
        query, project_id=project_id,
        index_col=index_col, col_order=col_order,
        reauth=reauth, verbose=verbose,
        private_key=private_key,
        dialect=dialect,
        **kwargs)


def to_gbq(
        dataframe, destination_table, project_id, chunksize=10000,
        verbose=True, reauth=False, if_exists='fail', private_key=None,
        **kwargs):
    """
    Write a DataFrame to a Google BigQuery table.

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
    dataframe : DataFrame
        DataFrame to be written.
    destination_table : string
        Name of table to be written, in the form 'dataset.tablename'.
    project_id : str
        Google BigQuery Account project ID.
    chunksize : int (default 10000)
        Number of rows to be inserted in each chunk from the dataframe.
        Set to ``None`` to load the whole dataframe at once.
    verbose : boolean (default True)
        Show percentage complete.
    reauth : boolean (default False)
        Force Google BigQuery to reauthenticate the user. This is useful
        if multiple accounts are used.
    if_exists : {'fail', 'replace', 'append'}, default 'fail'
        Behavior when the destination table exists.
        'fail': If table exists, do nothing.
        'replace': If table exists, drop it, recreate it, and insert data.
        'append': If table exists, insert data. Create if does not exist.
    private_key : str (optional)
        Service account private key in JSON format. Can be file path
        or string contents. This is useful for remote server
        authentication (eg. Jupyter/IPython notebook on remote host).
    kwargs : dict
        Arbitrary keyword arguments.

        auth_local_webserver (boolean): default False
            Use the [local webserver flow] instead of the [console flow] when
            getting user credentials.

            .. [local webserver flow]
                http://google-auth-oauthlib.readthedocs.io/en/latest/reference/google_auth_oauthlib.flow.html#google_auth_oauthlib.flow.InstalledAppFlow.run_local_server
            .. [console flow]
                http://google-auth-oauthlib.readthedocs.io/en/latest/reference/google_auth_oauthlib.flow.html#google_auth_oauthlib.flow.InstalledAppFlow.run_console
            .. versionadded:: pandas-gbq 0.2.0
        table_schema (list of dicts):
            List of BigQuery table fields to which according DataFrame columns
            conform to, e.g. `[{'name': 'col1', 'type': 'STRING'},...]`. If
            schema is not provided, it will be generated according to dtypes
            of DataFrame columns. See BigQuery API documentation on available
            names of a field.
            .. versionadded:: pandas-gbq 0.3.1

    See Also
    --------
    pandas_gbq.to_gbq
    pandas.DataFrame.to_gbq
    """
    pandas_gbq = _try_import()
    pandas_gbq.to_gbq(
        dataframe, destination_table, project_id, chunksize=chunksize,
        verbose=verbose, reauth=reauth, if_exists=if_exists,
        private_key=private_key, **kwargs)
