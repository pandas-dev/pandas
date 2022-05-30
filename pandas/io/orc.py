""" orc compat """
from __future__ import annotations

import io
from types import ModuleType
from typing import (
    TYPE_CHECKING,
    Literal,
)

from pandas._typing import (
    FilePath,
    ReadBuffer,
    WriteBuffer,
)
from pandas.compat._optional import import_optional_dependency

from pandas.io.common import get_handle

if TYPE_CHECKING:
    from pandas import DataFrame


def read_orc(
    path: FilePath | ReadBuffer[bytes], columns: list[str] | None = None, **kwargs
) -> DataFrame:
    """
    Load an ORC object from the file path, returning a DataFrame.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    path : str, path object, or file-like object
        String, path object (implementing ``os.PathLike[str]``), or file-like
        object implementing a binary ``read()`` function. The string could be a URL.
        Valid URL schemes include http, ftp, s3, and file. For file URLs, a host is
        expected. A local file could be:
        ``file://localhost/path/to/table.orc``.
    columns : list, default None
        If not None, only these columns will be read from the file.
    **kwargs
        Any additional kwargs are passed to pyarrow.

    Returns
    -------
    DataFrame

    Notes
    -----
    Before using this function you should read the :ref:`user guide about ORC <io.orc>`
    and :ref:`install optional dependencies <install.warn_orc>`.
    """
    # we require a newer version of pyarrow than we support for parquet

    orc = import_optional_dependency("pyarrow.orc")

    with get_handle(path, "rb", is_text=False) as handles:
        orc_file = orc.ORCFile(handles.handle)
        return orc_file.read(columns=columns, **kwargs).to_pandas()


def to_orc(
    df: DataFrame,
    path: FilePath | WriteBuffer[bytes] | None = None,
    *,
    engine: Literal["pyarrow"] = "pyarrow",
    index: bool | None = None,
    **kwargs,
) -> bytes | None:
    """
    Write a DataFrame to the ORC format.

    .. versionadded:: 1.5.0

    Parameters
    ----------
    df : DataFrame
        The dataframe to be written to ORC. Raises NotImplementedError
        if dtype of one or more columns is category, unsigned integers,
        intervals, periods or sparse.
    path : str, file-like object or None, default None
        If a string, it will be used as Root Directory path
        when writing a partitioned dataset. By file-like object,
        we refer to objects with a write() method, such as a file handle
        (e.g. via builtin open function). If path is None,
        a bytes object is returned.
    engine : {{'pyarrow'}}, default 'pyarrow'
        ORC library to use, or library it self, checked with 'pyarrow' name
        and version >= 7.0.0. Raises ValueError if it is anything but
        'pyarrow'.
    index : bool, optional
        If ``True``, include the dataframe's index(es) in the file output. If
        ``False``, they will not be written to the file.
        If ``None``, similar to ``infer`` the dataframe's index(es)
        will be saved. However, instead of being saved as values,
        the RangeIndex will be stored as a range in the metadata so it
        doesn't require much space and is faster. Other indexes will
        be included as columns in the file output.
    **kwargs
        Additional keyword arguments passed to the engine

    Returns
    -------
    bytes if no path argument is provided else None

    Raises
    ------
    NotImplementedError
        * Dtype of one or more columns is unsigned integers, intervals,
            periods, sparse or categorical.
    ValueError
        * engine is not pyarrow.

    Notes
    -----
    * Before using this function you should read the
    :ref:`user guide about ORC <io.orc>` and
    :ref:`install optional dependencies <install.warn_orc>`.
    * This function requires `pyarrow <https://arrow.apache.org/docs/python/>`
    _ library.
    * Unsigned integers, intervals, periods, sparse and categorical Dtypes
    are not supported yet.
    * Currently timezones in datetime columns are not preserved when a
    dataframe is converted into ORC files.
    """
    if index is None:
        index = df.index.names[0] is not None

    # If unsupported dtypes are found raise NotImplementedError
    for dtype in df.dtypes:
        dtype_str = dtype.__str__().lower()
        if (
            "category" in dtype_str
            or "interval" in dtype_str
            or "sparse" in dtype_str
            or "period" in dtype_str
            or "uint" in dtype_str
        ):
            raise NotImplementedError(
                """The dtype of one or more columns is unsigned integers,
intervals, periods, sparse or categorical which is not supported yet."""
            )

    if engine != "pyarrow":
        raise ValueError("engine must be 'pyarrow'")
    engine = import_optional_dependency(engine, min_version="7.0.0")
    orc = import_optional_dependency("pyarrow.orc")

    was_none = path is None
    if was_none:
        path = io.BytesIO()
    assert path is not None  # For mypy
    with get_handle(path, "wb", is_text=False) as handles:
        assert isinstance(engine, ModuleType)  # For mypy
        orc.write_table(
            engine.Table.from_pandas(df, preserve_index=index), handles.handle, **kwargs
        )

    if was_none:
        assert isinstance(path, io.BytesIO)  # For mypy
        return path.getvalue()
    return None
