""" feather-format compat """
from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Hashable,
    Sequence,
)

from pandas._config import using_nullable_dtypes

from pandas._libs import lib
from pandas.compat._optional import import_optional_dependency
from pandas.util._decorators import doc

from pandas import (
    arrays,
    get_option,
)
from pandas.core.api import DataFrame
from pandas.core.shared_docs import _shared_docs

from pandas.io.common import get_handle

if TYPE_CHECKING:
    from pandas._typing import (
        FilePath,
        ReadBuffer,
        StorageOptions,
        WriteBuffer,
    )


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
    use_nullable_dtypes: bool | lib.NoDefault = lib.no_default,
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

    use_nullable_dtypes : bool = False
        Whether or not to use nullable dtypes as default when reading data. If
        set to True, nullable dtypes are used for all dtypes that have a nullable
        implementation, even if no nulls are present.

        .. note::

            The nullable dtype implementation can be configured by calling
            ``pd.set_option("mode.dtype_backend", "pandas")`` to use
            numpy-backed nullable dtypes or
            ``pd.set_option("mode.dtype_backend", "pyarrow")`` to use
            pyarrow-backed nullable dtypes (using ``pd.ArrowDtype``).

        .. versionadded:: 2.0

    Returns
    -------
    type of object stored in file
    """
    import_optional_dependency("pyarrow")
    from pyarrow import feather

    use_nullable_dtypes = (
        use_nullable_dtypes
        if use_nullable_dtypes is not lib.no_default
        else using_nullable_dtypes()
    )

    with get_handle(
        path, "rb", storage_options=storage_options, is_text=False
    ) as handles:
        if not use_nullable_dtypes:
            return feather.read_feather(
                handles.handle, columns=columns, use_threads=bool(use_threads)
            )

        dtype_backend = get_option("mode.dtype_backend")

        pa_table = feather.read_table(
            handles.handle, columns=columns, use_threads=bool(use_threads)
        )

        if dtype_backend == "pandas":
            from pandas.io._util import _arrow_dtype_mapping

            return pa_table.to_pandas(types_mapper=_arrow_dtype_mapping().get)

        elif dtype_backend == "pyarrow":
            return DataFrame(
                {
                    col_name: arrays.ArrowExtensionArray(pa_col)
                    for col_name, pa_col in zip(
                        pa_table.column_names, pa_table.itercolumns()
                    )
                }
            )
