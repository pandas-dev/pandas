""" pickle compat """
import pickle
from typing import Any, Optional
import warnings

from pandas._typing import FilePathOrBuffer
from pandas.compat import pickle_compat as pc

from pandas.io.common import get_filepath_or_buffer, get_handle


def to_pickle(
    obj: Any,
    filepath_or_buffer: FilePathOrBuffer,
    compression: Optional[str] = "infer",
    protocol: int = pickle.HIGHEST_PROTOCOL,
):
    """
    Pickle (serialize) object to file.

    Parameters
    ----------
    obj : any object
        Any python object.
    filepath_or_buffer : str, path object or file-like object
        File path, URL, or buffer where the pickled object will be stored.

        .. versionchanged:: 1.0.0
           Accept URL. URL has to be of S3 or GCS.

    compression : {'infer', 'gzip', 'bz2', 'zip', 'xz', None}, default 'infer'
        If 'infer' and 'path_or_url' is path-like, then detect compression from
        the following extensions: '.gz', '.bz2', '.zip', or '.xz' (otherwise no
        compression) If 'infer' and 'path_or_url' is not path-like, then use
        None (= no decompression).
    protocol : int
        Int which indicates which protocol should be used by the pickler,
        default HIGHEST_PROTOCOL (see [1], paragraph 12.1.2). The possible
        values for this parameter depend on the version of Python. For Python
        2.x, possible values are 0, 1, 2. For Python>=3.0, 3 is a valid value.
        For Python >= 3.4, 4 is a valid value. A negative value for the
        protocol parameter is equivalent to setting its value to
        HIGHEST_PROTOCOL.

        .. [1] https://docs.python.org/3/library/pickle.html
        .. versionadded:: 0.21.0

    See Also
    --------
    read_pickle : Load pickled pandas object (or any object) from file.
    DataFrame.to_hdf : Write DataFrame to an HDF5 file.
    DataFrame.to_sql : Write DataFrame to a SQL database.
    DataFrame.to_parquet : Write a DataFrame to the binary parquet format.

    Examples
    --------
    >>> original_df = pd.DataFrame({"foo": range(5), "bar": range(5, 10)})
    >>> original_df
       foo  bar
    0    0    5
    1    1    6
    2    2    7
    3    3    8
    4    4    9
    >>> pd.to_pickle(original_df, "./dummy.pkl")

    >>> unpickled_df = pd.read_pickle("./dummy.pkl")
    >>> unpickled_df
       foo  bar
    0    0    5
    1    1    6
    2    2    7
    3    3    8
    4    4    9

    >>> import os
    >>> os.remove("./dummy.pkl")
    """
    fp_or_buf, _, compression, should_close = get_filepath_or_buffer(
        filepath_or_buffer, compression=compression, mode="wb"
    )
    if not isinstance(fp_or_buf, str) and compression == "infer":
        compression = None
    f, fh = get_handle(fp_or_buf, "wb", compression=compression, is_text=False)
    if protocol < 0:
        protocol = pickle.HIGHEST_PROTOCOL
    try:
        f.write(pickle.dumps(obj, protocol=protocol))
    finally:
        f.close()
        for _f in fh:
            _f.close()
        if should_close:
            try:
                fp_or_buf.close()
            except ValueError:
                pass


def read_pickle(
    filepath_or_buffer: FilePathOrBuffer, compression: Optional[str] = "infer"
):
    """
    Load pickled pandas object (or any object) from file.

    .. warning::

       Loading pickled data received from untrusted sources can be
       unsafe. See `here <https://docs.python.org/3/library/pickle.html>`__.

    Parameters
    ----------
    filepath_or_buffer : str, path object or file-like object
        File path, URL, or buffer where the pickled object will be loaded from.

        .. versionchanged:: 1.0.0
           Accept URL. URL is not limited to S3 and GCS.

    compression : {'infer', 'gzip', 'bz2', 'zip', 'xz', None}, default 'infer'
        If 'infer' and 'path_or_url' is path-like, then detect compression from
        the following extensions: '.gz', '.bz2', '.zip', or '.xz' (otherwise no
        compression) If 'infer' and 'path_or_url' is not path-like, then use
        None (= no decompression).

    Returns
    -------
    unpickled : same type as object stored in file

    See Also
    --------
    DataFrame.to_pickle : Pickle (serialize) DataFrame object to file.
    Series.to_pickle : Pickle (serialize) Series object to file.
    read_hdf : Read HDF5 file into a DataFrame.
    read_sql : Read SQL query or database table into a DataFrame.
    read_parquet : Load a parquet object, returning a DataFrame.

    Notes
    -----
    read_pickle is only guaranteed to be backwards compatible to pandas 0.20.3.

    Examples
    --------
    >>> original_df = pd.DataFrame({"foo": range(5), "bar": range(5, 10)})
    >>> original_df
       foo  bar
    0    0    5
    1    1    6
    2    2    7
    3    3    8
    4    4    9
    >>> pd.to_pickle(original_df, "./dummy.pkl")

    >>> unpickled_df = pd.read_pickle("./dummy.pkl")
    >>> unpickled_df
       foo  bar
    0    0    5
    1    1    6
    2    2    7
    3    3    8
    4    4    9

    >>> import os
    >>> os.remove("./dummy.pkl")
    """
    fp_or_buf, _, compression, should_close = get_filepath_or_buffer(
        filepath_or_buffer, compression=compression
    )
    if not isinstance(fp_or_buf, str) and compression == "infer":
        compression = None
    f, fh = get_handle(fp_or_buf, "rb", compression=compression, is_text=False)

    # 1) try standard library Pickle
    # 2) try pickle_compat (older pandas version) to handle subclass changes
    # 3) try pickle_compat with latin-1 encoding upon a UnicodeDecodeError

    try:
        excs_to_catch = (AttributeError, ImportError, ModuleNotFoundError)
        try:
            with warnings.catch_warnings(record=True):
                # We want to silence any warnings about, e.g. moved modules.
                warnings.simplefilter("ignore", Warning)
                return pickle.load(f)
        except excs_to_catch:
            # e.g.
            #  "No module named 'pandas.core.sparse.series'"
            #  "Can't get attribute '__nat_unpickle' on <module 'pandas._libs.tslib"
            return pc.load(f, encoding=None)
    except UnicodeDecodeError:
        # e.g. can occur for files written in py27; see GH#28645 and GH#31988
        return pc.load(f, encoding="latin-1")
    finally:
        f.close()
        for _f in fh:
            _f.close()
        if should_close:
            try:
                fp_or_buf.close()
            except ValueError:
                pass
