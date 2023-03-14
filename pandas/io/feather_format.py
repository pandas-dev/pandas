""" feather-format compat """
from __future__ import annotations

from typing import (
    Hashable,
    Sequence,
)

from pandas._libs import lib
from pandas._typing import (
    DtypeBackend,
    FilePath,
    ReadBuffer,
    StorageOptions,
    WriteBuffer,
)
from pandas.compat._optional import import_optional_dependency
from pandas.util._decorators import doc
from pandas.util._validators import check_dtype_backend

import pandas as pd
from pandas.core.api import (
    DataFrame,
    RangeIndex,
)
from pandas.core.shared_docs import _shared_docs

from pandas.io.common import get_handle


@doc(storage_options=_shared_docs["storage_options"])
def to_feather(
    df: DataFrame,
    path: FilePath | WriteBuffer[bytes],
    storage_options: StorageOptions = None,
    **kwargs,
) -> None:
    """
    Write a DataFrame to the binary Feather format.

    Parameters
    ----------
    df : DataFrame
    path : str, path object, or file-like object
    {storage_options}

        .. versionadded:: 1.2.0

    **kwargs :
        Additional keywords passed to `pyarrow.feather.write_feather`.

        .. versionadded:: 1.1.0
    """
    import_optional_dependency("pyarrow")
    from pyarrow import feather

    if not isinstance(df, DataFrame):
        raise ValueError("feather only support IO with DataFrames")

    valid_types = {"string", "unicode"}

    # validate index
    # --------------

    # validate that we have only a default index
    # raise on anything else as we don't serialize the index

    if not df.index.dtype == "int64":
        typ = type(df.index)
        raise ValueError(
            f"feather does not support serializing {typ} "
            "for the index; you can .reset_index() to make the index into column(s)"
        )

    if not df.index.equals(RangeIndex.from_range(range(len(df)))):
        raise ValueError(
            "feather does not support serializing a non-default index for the index; "
            "you can .reset_index() to make the index into column(s)"
        )

    if df.index.name is not None:
        raise ValueError(
            "feather does not serialize index meta-data on a default index"
        )

    # validate columns
    # ----------------

    # must have value column names (strings only)
    if df.columns.inferred_type not in valid_types:
        raise ValueError("feather must have string column names")

    with get_handle(
        path, "wb", storage_options=storage_options, is_text=False
    ) as handles:
        feather.write_feather(df, handles.handle, **kwargs)


@doc(storage_options=_shared_docs["storage_options"])
def read_feather(
    path: FilePath | ReadBuffer[bytes],
    columns: Sequence[Hashable] | None = None,
    use_threads: bool = True,
    storage_options: StorageOptions = None,
    dtype_backend: DtypeBackend | lib.NoDefault = lib.no_default,
):
    """
    Load a feather-format object from the file path.

    Parameters
    ----------
    path : str, path object, or file-like object
        String, path object (implementing ``os.PathLike[str]``), or file-like
        object implementing a binary ``read()`` function. The string could be a URL.
        Valid URL schemes include http, ftp, s3, and file. For file URLs, a host is
        expected. A local file could be: ``file://localhost/path/to/table.feather``.
    columns : sequence, default None
        If not provided, all columns are read.
    use_threads : bool, default True
        Whether to parallelize reading using multiple threads.
    {storage_options}

        .. versionadded:: 1.2.0

    dtype_backend : {{"numpy_nullable", "pyarrow"}}, defaults to NumPy backed DataFrames
        Which dtype_backend to use, e.g. whether a DataFrame should have NumPy
        arrays, nullable dtypes are used for all dtypes that have a nullable
        implementation when "numpy_nullable" is set, pyarrow is used for all
        dtypes if "pyarrow" is set.

        The dtype_backends are still experimential.

        .. versionadded:: 2.0

    Returns
    -------
    type of object stored in file
    """
    import_optional_dependency("pyarrow")
    from pyarrow import feather

    check_dtype_backend(dtype_backend)

    with get_handle(
        path, "rb", storage_options=storage_options, is_text=False
    ) as handles:
        if dtype_backend is lib.no_default:
            return feather.read_feather(
                handles.handle, columns=columns, use_threads=bool(use_threads)
            )

        pa_table = feather.read_table(
            handles.handle, columns=columns, use_threads=bool(use_threads)
        )

        if dtype_backend == "numpy_nullable":
            from pandas.io._util import _arrow_dtype_mapping

            return pa_table.to_pandas(types_mapper=_arrow_dtype_mapping().get)

        elif dtype_backend == "pyarrow":
            return pa_table.to_pandas(types_mapper=pd.ArrowDtype)
