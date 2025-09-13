"""feather-format compat"""

from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
)
import warnings

from pandas._config import using_string_dtype

from pandas._libs import lib
from pandas.compat._optional import import_optional_dependency
from pandas.errors import Pandas4Warning
from pandas.util._decorators import doc
from pandas.util._validators import check_dtype_backend

from pandas.core.api import DataFrame
from pandas.core.shared_docs import _shared_docs

from pandas.io._util import arrow_table_to_pandas
from pandas.io.common import get_handle

if TYPE_CHECKING:
    from collections.abc import (
        Hashable,
        Sequence,
    )

    from pandas._typing import (
        DtypeBackend,
        FilePath,
        ReadBuffer,
        StorageOptions,
        WriteBuffer,
    )


@doc(storage_options=_shared_docs["storage_options"])
def to_feather(
    df: DataFrame,
    path: FilePath | WriteBuffer[bytes],
    storage_options: StorageOptions | None = None,
    **kwargs: Any,
) -> None:
    """
    Write a DataFrame to the binary Feather format.

    Parameters
    ----------
    df : DataFrame
    path : str, path object, or file-like object
    {storage_options}
    **kwargs :
        Additional keywords passed to `pyarrow.feather.write_feather`.

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
    storage_options: StorageOptions | None = None,
    dtype_backend: DtypeBackend | lib.NoDefault = lib.no_default,
) -> DataFrame:
    """
    Load a feather-format object from the file path.

    Feather is particularly useful for scenarios that require efficient
    serialization and deserialization of tabular data. It supports
    schema preservation, making it a reliable choice for use cases
    such as sharing data between Python and R, or persisting intermediate
    results during data processing pipelines. This method provides additional
    flexibility with options for selective column reading, thread parallelism,
    and choosing the backend for data types.

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

    dtype_backend : {{'numpy_nullable', 'pyarrow'}}
        Back-end data type applied to the resultant :class:`DataFrame`
        (still experimental). If not specified, the default behavior
        is to not use nullable data types. If specified, the behavior
        is as follows:

        * ``"numpy_nullable"``: returns nullable-dtype-backed :class:`DataFrame`.
        * ``"pyarrow"``: returns pyarrow-backed nullable
          :class:`ArrowDtype` :class:`DataFrame`

        .. versionadded:: 2.0

    Returns
    -------
    type of object stored in file
        DataFrame object stored in the file.

    See Also
    --------
    read_csv : Read a comma-separated values (csv) file into a pandas DataFrame.
    read_excel : Read an Excel file into a pandas DataFrame.
    read_spss : Read an SPSS file into a pandas DataFrame.
    read_orc : Load an ORC object into a pandas DataFrame.
    read_sas : Read SAS file into a pandas DataFrame.

    Examples
    --------
    >>> df = pd.read_feather("path/to/file.feather")  # doctest: +SKIP
    """
    import_optional_dependency("pyarrow")
    from pyarrow import feather

    # import utils to register the pyarrow extension types
    import pandas.core.arrays.arrow.extension_types  # pyright: ignore[reportUnusedImport] # noqa: F401

    check_dtype_backend(dtype_backend)

    with get_handle(
        path, "rb", storage_options=storage_options, is_text=False
    ) as handles:
        if dtype_backend is lib.no_default and not using_string_dtype():
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    "make_block is deprecated",
                    Pandas4Warning,
                )

                return feather.read_feather(
                    handles.handle, columns=columns, use_threads=bool(use_threads)
                )

        pa_table = feather.read_table(
            handles.handle, columns=columns, use_threads=bool(use_threads)
        )
        return arrow_table_to_pandas(pa_table, dtype_backend=dtype_backend)
