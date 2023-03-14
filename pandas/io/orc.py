""" orc compat """
from __future__ import annotations

import io
from types import ModuleType
from typing import (
    Any,
    Literal,
)

from pandas._libs import lib
from pandas._typing import (
    DtypeBackend,
    FilePath,
    ReadBuffer,
    WriteBuffer,
)
from pandas.compat._optional import import_optional_dependency
from pandas.util._validators import check_dtype_backend

from pandas.core.dtypes.common import (
    is_categorical_dtype,
    is_interval_dtype,
    is_period_dtype,
    is_unsigned_integer_dtype,
)

import pandas as pd
from pandas.core.frame import DataFrame

from pandas.io.common import get_handle


def read_orc(
    path: FilePath | ReadBuffer[bytes],
    columns: list[str] | None = None,
    dtype_backend: DtypeBackend | lib.NoDefault = lib.no_default,
    **kwargs,
) -> DataFrame:
    """
    Load an ORC object from the file path, returning a DataFrame.

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
        Output always follows the ordering of the file and not the columns list.
        This mirrors the original behaviour of
        :external+pyarrow:py:meth:`pyarrow.orc.ORCFile.read`.
    dtype_backend : {"numpy_nullable", "pyarrow"}, defaults to NumPy backed DataFrames
        Which dtype_backend to use, e.g. whether a DataFrame should have NumPy
        arrays, nullable dtypes are used for all dtypes that have a nullable
        implementation when "numpy_nullable" is set, pyarrow is used for all
        dtypes if "pyarrow" is set.

        The dtype_backends are still experimential.

        .. versionadded:: 2.0

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

    check_dtype_backend(dtype_backend)

    with get_handle(path, "rb", is_text=False) as handles:
        orc_file = orc.ORCFile(handles.handle)
        pa_table = orc_file.read(columns=columns, **kwargs)
    if dtype_backend is not lib.no_default:
        if dtype_backend == "pyarrow":
            df = pa_table.to_pandas(types_mapper=pd.ArrowDtype)
        else:
            from pandas.io._util import _arrow_dtype_mapping

            mapping = _arrow_dtype_mapping()
            df = pa_table.to_pandas(types_mapper=mapping.get)
        return df
    else:
        return pa_table.to_pandas()


def to_orc(
    df: DataFrame,
    path: FilePath | WriteBuffer[bytes] | None = None,
    *,
    engine: Literal["pyarrow"] = "pyarrow",
    index: bool | None = None,
    engine_kwargs: dict[str, Any] | None = None,
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
    engine : str, default 'pyarrow'
        ORC library to use. Pyarrow must be >= 7.0.0.
    index : bool, optional
        If ``True``, include the dataframe's index(es) in the file output. If
        ``False``, they will not be written to the file.
        If ``None``, similar to ``infer`` the dataframe's index(es)
        will be saved. However, instead of being saved as values,
        the RangeIndex will be stored as a range in the metadata so it
        doesn't require much space and is faster. Other indexes will
        be included as columns in the file output.
    engine_kwargs : dict[str, Any] or None, default None
        Additional keyword arguments passed to :func:`pyarrow.orc.write_table`.

    Returns
    -------
    bytes if no path argument is provided else None

    Raises
    ------
    NotImplementedError
        Dtype of one or more columns is category, unsigned integers, interval,
        period or sparse.
    ValueError
        engine is not pyarrow.

    Notes
    -----
    * Before using this function you should read the
      :ref:`user guide about ORC <io.orc>` and
      :ref:`install optional dependencies <install.warn_orc>`.
    * This function requires `pyarrow <https://arrow.apache.org/docs/python/>`_
      library.
    * For supported dtypes please refer to `supported ORC features in Arrow
      <https://arrow.apache.org/docs/cpp/orc.html#data-types>`__.
    * Currently timezones in datetime columns are not preserved when a
      dataframe is converted into ORC files.
    """
    if index is None:
        index = df.index.names[0] is not None
    if engine_kwargs is None:
        engine_kwargs = {}

    # If unsupported dtypes are found raise NotImplementedError
    # In Pyarrow 9.0.0 this check will no longer be needed
    for dtype in df.dtypes:
        if (
            is_categorical_dtype(dtype)
            or is_interval_dtype(dtype)
            or is_period_dtype(dtype)
            or is_unsigned_integer_dtype(dtype)
        ):
            raise NotImplementedError(
                "The dtype of one or more columns is not supported yet."
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
        try:
            orc.write_table(
                engine.Table.from_pandas(df, preserve_index=index),
                handles.handle,
                **engine_kwargs,
            )
        except TypeError as e:
            raise NotImplementedError(
                "The dtype of one or more columns is not supported yet."
            ) from e

    if was_none:
        assert isinstance(path, io.BytesIO)  # For mypy
        return path.getvalue()
    return None
