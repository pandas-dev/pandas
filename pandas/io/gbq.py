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


def read_gbq(*args, **kwargs):
    """
    Load data from Google BigQuery.

    This function requires the `pandas-gbq package
    <https://pandas-gbq.readthedocs.io>`__.

    See :meth:`pandas_gbq.read_gbq`.

    Returns
    -------
    df: DataFrame
        DataFrame representing results of query
    """
    pandas_gbq = _try_import()
    return pandas_gbq.read_gbq(*args, **kwargs)


def to_gbq(*args, **kwargs):
    """
    Write a DataFrame to a Google BigQuery table.

    This function requires the `pandas-gbq package
    <https://pandas-gbq.readthedocs.io>`__.

    See :meth:`pandas_gbq.to_gbq`.
    """
    pandas_gbq = _try_import()
    pandas_gbq.to_gbq(*args, **kwargs)
