""" Google BigQuery support """

import warnings


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
             reauth=False, private_key=None, auth_local_webserver=False,
             dialect=None, location=None, configuration=None,
             verbose=None):
    """
    Load data from Google BigQuery.

    This function requires the `pandas-gbq package
    <https://pandas-gbq.readthedocs.io>`__.

    See the `How to authenticate with Google BigQuery
    <https://pandas-gbq.readthedocs.io/en/latest/howto/authentication.html>`__
    guide for authentication instructions.

    Parameters
    ----------
    query : str
        SQL-Like Query to return data values.
    project_id : str, optional
        Google BigQuery Account project ID. Optional when available from
        the environment.
    index_col : str, optional
        Name of result column to use for index in results DataFrame.
    col_order : list(str), optional
        List of BigQuery column names in the desired order for results
        DataFrame.
    reauth : boolean, default False
        Force Google BigQuery to re-authenticate the user. This is useful
        if multiple accounts are used.
    private_key : str, optional
        Service account private key in JSON format. Can be file path
        or string contents. This is useful for remote server
        authentication (eg. Jupyter/IPython notebook on remote host).
    auth_local_webserver : boolean, default False
        Use the `local webserver flow`_ instead of the `console flow`_
        when getting user credentials.

        .. _local webserver flow:
            http://google-auth-oauthlib.readthedocs.io/en/latest/reference/google_auth_oauthlib.flow.html#google_auth_oauthlib.flow.InstalledAppFlow.run_local_server
        .. _console flow:
            http://google-auth-oauthlib.readthedocs.io/en/latest/reference/google_auth_oauthlib.flow.html#google_auth_oauthlib.flow.InstalledAppFlow.run_console

        *New in version 0.2.0 of pandas-gbq*.
    dialect : str, default 'legacy'
        Note: The default value is changing to 'standard' in a future verion.

        SQL syntax dialect to use. Value can be one of:

        ``'legacy'``
            Use BigQuery's legacy SQL dialect. For more information see
            `BigQuery Legacy SQL Reference
            <https://cloud.google.com/bigquery/docs/reference/legacy-sql>`__.
        ``'standard'``
            Use BigQuery's standard SQL, which is
            compliant with the SQL 2011 standard. For more information
            see `BigQuery Standard SQL Reference
            <https://cloud.google.com/bigquery/docs/reference/standard-sql/>`__.

        .. versionchanged:: 0.24.0
    location : str, optional
        Location where the query job should run. See the `BigQuery locations
        documentation
        <https://cloud.google.com/bigquery/docs/dataset-locations>`__ for a
        list of available locations. The location must match that of any
        datasets used in the query.

        *New in version 0.5.0 of pandas-gbq*.
    configuration : dict, optional
        Query config parameters for job processing.
        For example:

            configuration = {'query': {'useQueryCache': False}}

        For more information see `BigQuery REST API Reference
        <https://cloud.google.com/bigquery/docs/reference/rest/v2/jobs#configuration.query>`__.
    verbose : None, deprecated
        Deprecated in Pandas-GBQ 0.4.0. Use the `logging module
        to adjust verbosity instead
        <https://pandas-gbq.readthedocs.io/en/latest/intro.html#logging>`__.

    Returns
    -------
    df: DataFrame
        DataFrame representing results of query.

    See Also
    --------
    pandas_gbq.read_gbq : This function in the pandas-gbq library.
    pandas.DataFrame.to_gbq : Write a DataFrame to Google BigQuery.
    """
    pandas_gbq = _try_import()

    if dialect is None:
        dialect = "legacy"
        warnings.warn(
            'The default value for dialect is changing to "standard" in a '
            'future version of pandas-gbq. Pass in dialect="legacy" to '
            "disable this warning.",
            FutureWarning,
            stacklevel=2,
        )

    return pandas_gbq.read_gbq(
        query, project_id=project_id, index_col=index_col,
        col_order=col_order, reauth=reauth, verbose=verbose,
        private_key=private_key, auth_local_webserver=auth_local_webserver,
        dialect=dialect, location=location, configuration=configuration)


def to_gbq(dataframe, destination_table, project_id=None, chunksize=None,
           verbose=None, reauth=False, if_exists='fail', private_key=None,
           auth_local_webserver=False, table_schema=None, location=None,
           progress_bar=True):
    pandas_gbq = _try_import()
    return pandas_gbq.to_gbq(
        dataframe, destination_table, project_id=project_id,
        chunksize=chunksize, verbose=verbose, reauth=reauth,
        if_exists=if_exists, private_key=private_key,
        auth_local_webserver=auth_local_webserver,
        table_schema=table_schema, location=location,
        progress_bar=progress_bar)
