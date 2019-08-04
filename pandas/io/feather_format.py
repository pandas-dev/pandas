""" feather-format compat """

from distutils.version import LooseVersion

from pandas.compat._optional import import_optional_dependency
from pandas.util._decorators import deprecate_kwarg

from pandas import DataFrame, Int64Index, RangeIndex

from pandas.io.common import _stringify_path


def to_feather(df, path):
    """
    Write a DataFrame to the feather-format

    Parameters
    ----------
    df : DataFrame
    path : string file path, or file-like object

    """
    import_optional_dependency("pyarrow")
    from pyarrow import feather

    path = _stringify_path(path)

    if not isinstance(df, DataFrame):
        raise ValueError("feather only support IO with DataFrames")

    valid_types = {"string", "unicode"}

    # validate index
    # --------------

    # validate that we have only a default index
    # raise on anything else as we don't serialize the index

    if not isinstance(df.index, Int64Index):
        raise ValueError(
            "feather does not support serializing {} "
            "for the index; you can .reset_index()"
            "to make the index into column(s)".format(type(df.index))
        )

    if not df.index.equals(RangeIndex.from_range(range(len(df)))):
        raise ValueError(
            "feather does not support serializing a "
            "non-default index for the index; you "
            "can .reset_index() to make the index "
            "into column(s)"
        )

    if df.index.name is not None:
        raise ValueError(
            "feather does not serialize index meta-data on a " "default index"
        )

    # validate columns
    # ----------------

    # must have value column names (strings only)
    if df.columns.inferred_type not in valid_types:
        raise ValueError("feather must have string column names")

    feather.write_feather(df, path)


@deprecate_kwarg(old_arg_name="nthreads", new_arg_name="use_threads")
def read_feather(path, columns=None, use_threads=True):
    """
    Load a feather-format object from the file path.

    .. versionadded 0.20.0

    Parameters
    ----------
    path : str, path object or file-like object
        Any valid string path is acceptable. The string could be a URL. Valid
        URL schemes include http, ftp, s3, and file. For file URLs, a host is
        expected. A local file could be:
        ``file://localhost/path/to/table.feather``.

        If you want to pass in a path object, pandas accepts any
        ``os.PathLike``.

        By file-like object, we refer to objects with a ``read()`` method,
        such as a file handler (e.g. via builtin ``open`` function)
        or ``StringIO``.
    columns : sequence, default None
        If not provided, all columns are read.

        .. versionadded 0.24.0
    nthreads : int, default 1
        Number of CPU threads to use when reading to pandas.DataFrame.

       .. versionadded 0.21.0
       .. deprecated 0.24.0
    use_threads : bool, default True
        Whether to parallelize reading using multiple threads.

       .. versionadded 0.24.0

    Returns
    -------
    type of object stored in file
    """
    pyarrow = import_optional_dependency("pyarrow")
    from pyarrow import feather

    path = _stringify_path(path)

    if LooseVersion(pyarrow.__version__) < LooseVersion("0.11.0"):
        int_use_threads = int(use_threads)
        if int_use_threads < 1:
            int_use_threads = 1
        return feather.read_feather(path, columns=columns, nthreads=int_use_threads)

    return feather.read_feather(path, columns=columns, use_threads=bool(use_threads))
