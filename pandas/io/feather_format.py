""" feather-format compat """

from distutils.version import LooseVersion

from pandas.util._decorators import deprecate_kwarg

from pandas import DataFrame, Int64Index, RangeIndex

from pandas.io.common import _stringify_path


def _try_import():
    # since pandas is a dependency of pyarrow
    # we need to import on first use
    try:
        import pyarrow
        from pyarrow import feather
    except ImportError:
        # give a nice error message
        raise ImportError("pyarrow is not installed\n\n"
                          "you can install via conda\n"
                          "conda install pyarrow -c conda-forge\n"
                          "or via pip\n"
                          "pip install -U pyarrow\n")

    if LooseVersion(pyarrow.__version__) < LooseVersion('0.9.0'):
        raise ImportError("pyarrow >= 0.9.0 required for feather support\n\n"
                          "you can install via conda\n"
                          "conda install pyarrow -c conda-forge"
                          "or via pip\n"
                          "pip install -U pyarrow\n")

    return feather, pyarrow


def to_feather(df, path):
    """
    Write a DataFrame to the feather-format

    Parameters
    ----------
    df : DataFrame
    path : string file path, or file-like object

    """
    path = _stringify_path(path)
    if not isinstance(df, DataFrame):
        raise ValueError("feather only support IO with DataFrames")

    feather = _try_import()[0]
    valid_types = {'string', 'unicode'}

    # validate index
    # --------------

    # validate that we have only a default index
    # raise on anything else as we don't serialize the index

    if not isinstance(df.index, Int64Index):
        raise ValueError("feather does not support serializing {} "
                         "for the index; you can .reset_index()"
                         "to make the index into column(s)".format(
                             type(df.index)))

    if not df.index.equals(RangeIndex.from_range(range(len(df)))):
        raise ValueError("feather does not support serializing a "
                         "non-default index for the index; you "
                         "can .reset_index() to make the index "
                         "into column(s)")

    if df.index.name is not None:
        raise ValueError("feather does not serialize index meta-data on a "
                         "default index")

    # validate columns
    # ----------------

    # must have value column names (strings only)
    if df.columns.inferred_type not in valid_types:
        raise ValueError("feather must have string column names")

    feather.write_feather(df, path)


@deprecate_kwarg(old_arg_name='nthreads', new_arg_name='use_threads')
def read_feather(path, columns=None, use_threads=True):
    """
    Load a feather-format object from the file path

    .. versionadded 0.20.0

    Parameters
    ----------
    path : string file path, or file-like object
    columns : sequence, default None
        If not provided, all columns are read

        .. versionadded 0.24.0
    nthreads : int, default 1
        Number of CPU threads to use when reading to pandas.DataFrame

       .. versionadded 0.21.0
       .. deprecated 0.24.0
    use_threads : bool, default True
        Whether to parallelize reading using multiple threads

       .. versionadded 0.24.0

    Returns
    -------
    type of object stored in file
    """

    feather, pyarrow = _try_import()
    path = _stringify_path(path)

    if LooseVersion(pyarrow.__version__) < LooseVersion('0.11.0'):
        int_use_threads = int(use_threads)
        if int_use_threads < 1:
            int_use_threads = 1
        return feather.read_feather(path, columns=columns,
                                    nthreads=int_use_threads)

    return feather.read_feather(path, columns=columns,
                                use_threads=bool(use_threads))
