"""
Vortex format support for pandas.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any

    from pandas import DataFrame
    from pandas._typing import (
        FilePath,
        ReadBuffer,
        StorageOptions,
        WriteBuffer,
    )


def read_vortex(
    path: FilePath | ReadBuffer[bytes],
    columns: Sequence[str] | None = None,
    storage_options: StorageOptions | None = None,
    **kwargs: Any,
) -> DataFrame:
    """
    Load a Vortex file from the file path, returning a DataFrame.

    .. versionadded:: 3.0.0

    Parameters
    ----------
    path : str, path object, or file-like object
        String or path object (implementing ``os.PathLike[str]``).
    columns : list, optional
        If not None, only these columns will be read from the file.
    storage_options : dict, optional
        Extra options that make sense for a particular storage connection.
        Currently not supported for Vortex.
    **kwargs
        Any additional kwargs are passed to ``vortex.open`` and ``scan``.

    Returns
    -------
    DataFrame

    See Also
    --------
    DataFrame.to_vortex : Write a DataFrame to the Vortex format.

    Examples
    --------
    >>> df = pd.read_vortex("path/to/file.vortex")  # doctest: +SKIP
    """
    from pandas.compat._optional import import_optional_dependency

    vortex = import_optional_dependency("vortex", extra="vortex is required.")
    from pathlib import Path

    if isinstance(path, Path):
        path = str(path)

    v_file = vortex.open(path)
    scan = v_file.scan(projection=columns, **kwargs)
    result_array = scan.read_all()

    if hasattr(result_array, "to_arrow_table"):
        arrow_table = result_array.to_arrow_table()
    elif hasattr(result_array, "to_arrow"):
        import pyarrow as pa

        arrow_obj = result_array.to_arrow()
        arrow_table = arrow_obj if isinstance(arrow_obj, pa.Table) else pa.Table.from_batches([arrow_obj])
    else:
        raise ValueError("Unable to convert Vortex result to Arrow format.")

    return arrow_table.to_pandas()


def to_vortex(
    df: DataFrame,
    path: FilePath | WriteBuffer[bytes],
    *,
    storage_options: StorageOptions | None = None,
    **kwargs: Any,
) -> None:
    """
    Write a DataFrame to the Vortex binary format.

    .. versionadded:: 3.0.0

    Parameters
    ----------
    df : DataFrame
    path : str or path object
    storage_options : dict, optional
        Currently not supported for Vortex.
    **kwargs
        Additional arguments passed to ``vortex.io.write``.

    See Also
    --------
    read_vortex : Read a Vortex file.

    Examples
    --------
    >>> df.to_vortex("output.vortex")  # doctest: +SKIP
    """
    from pandas.compat._optional import import_optional_dependency

    vortex = import_optional_dependency("vortex", extra="vortex is required.")
    pa = import_optional_dependency("pyarrow", extra="pyarrow is required.")
    from pathlib import Path

    if isinstance(path, Path):
        path = str(path)

    table = pa.Table.from_pandas(df)
    v_array = vortex.array(table)
    vortex.io.write(v_array, path, **kwargs)
