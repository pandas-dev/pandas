""" Google BigQuery support """

from pandas.util.decorators import docstring_wrapper


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
                          "you can install via:\n"
                          "pip install pandas-gbq\n")

    return pandas_gbq


def read_gbq(query, project_id=None, index_col=None, col_order=None,
             reauth=False, verbose=True, private_key=None, dialect='legacy',
             **kwargs):
    pandas_gbq = _try_import()
    return pandas_gbq.read_gbq(
        query, project_id=project_id,
        index_col=index_col, col_order=col_order,
        reauth=reauth, verbose=verbose,
        private_key=private_key,
        dialect=dialect,
        **kwargs)


read_gbq = docstring_wrapper(read_gbq,
                             lambda: _try_import().read_gbq.__doc__)


def to_gbq(dataframe, destination_table, project_id, chunksize=10000,
           verbose=True, reauth=False, if_exists='fail', private_key=None):
    pandas_gbq = _try_import()
    pandas_gbq.to_gbq(dataframe, destination_table, project_id,
                      chunksize=chunksize,
                      verbose=verbose, reauth=reauth,
                      if_exists=if_exists, private_key=private_key)


to_gbq = docstring_wrapper(to_gbq,
                           lambda: _try_import().to_gbq.__doc__)
