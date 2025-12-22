"""
Vortex format support for pandas.
"""
from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
)

from pandas.compat._optional import import_optional_dependency

if TYPE_CHECKING:
    from pandas import DataFrame
    from pandas._typing import (
        FilePath,
        ReadBuffer,
        StorageOptions,
        WriteBuffer,
    )


def read_vortex(
    path: FilePath | ReadBuffer[bytes],
    columns: list[str] | None = None,
    storage_options: StorageOptions | None = None,
    **kwargs: Any,
) -> DataFrame:
    """
    Load a Vortex file from the file path, returning a DataFrame.

    Parameters
    ----------
    path : str, path object, or file-like object
        String or path object (implementing ``os.PathLike[str]``).
        The string could be a URL. Valid URL schemes include http, ftp, s3, 
        gs, and file.
    columns : list, optional
        If not None, only these columns will be read from the file.
    storage_options : dict, optional
        Extra options that make sense for a particular storage connection, e.g.
        host, port, username, password, etc. Currently not supported for Vortex.

        .. versionadded:: 3.0.0
    **kwargs
        Any additional kwargs are passed to ``vortex.open`` and ``scan``.

    Returns
    -------
    DataFrame
        DataFrame containing the data from the Vortex file.

    See Also
    --------
    DataFrame.to_vortex : Write a DataFrame to the Vortex format.
    read_parquet : Read a Parquet file.
    read_feather : Read a Feather file.
    read_orc : Read an ORC file.

    Examples
    --------
    >>> df = pd.read_vortex("path/to/file.vortex")  # doctest: +SKIP

    Read only certain columns:

    >>> df = pd.read_vortex(
    ...     "path/to/file.vortex",
    ...     columns=["col1", "col2"]
    ... )  # doctest: +SKIP
    """
    vortex = import_optional_dependency(
        "vortex", extra="vortex is required for Vortex support."
    )

    # Convert Path object to string if necessary
    from pathlib import Path
    if isinstance(path, Path):
        path = str(path)

    # Open the Vortex file
    v_file = vortex.open(path)

    # Perform scan with optional column projection
    scan = v_file.scan(projection=columns, **kwargs)

    # Read all data and convert to Arrow format
    result_array = scan.read_all()

    # Convert Vortex result to Arrow Table
    if hasattr(result_array, "to_arrow_table"):
        arrow_table = result_array.to_arrow_table()
    elif hasattr(result_array, "to_arrow"):
        arrow_obj = result_array.to_arrow()
        import pyarrow as pa

        if isinstance(arrow_obj, pa.Table):
            arrow_table = arrow_obj
        else:
            # Construct Table from batches if needed
            arrow_table = pa.Table.from_batches([arrow_obj])
    else:
        raise ValueError(
            " vortex is properly installed."
        )

    # Convert Arrow Table to pandas DataFrame
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

    Parameters
    ----------
    df : DataFrame
        The DataFrame to write.
    path : str or path object
        String or path object (implementing ``os.PathLike[str]``).
    storage_options : dict, optional
        Extra options that make sense for a particular storage connection, e.g.
        host, port, username, password, etc. Currently not supported for Vortex.

        .. versionadded:: 3.0.0
    **kwargs
        Additional arguments passed to ``vortex.io.write``.

    See Also
    --------
    read_vortex : Read a Vortex file.
    DataFrame.to_parquet : Write a DataFrame to the Parquet format.
    DataFrame.to_feather : Write a DataFrame to the Feather format.
    DataFrame.to_orc : Write a DataFrame to the ORC format.

    Notes
    -----
    This function writes the DataFrame to the Vortex columnar storage format,
    which is optimized for analytical workloads.

    Examples
    --------
    >>> df = pd.DataFrame({"A": [1, 2, 3], "B": ["x", "y", "z"]})
    >>> df.to_vortex("output.vortex")  # doctest: +SKIP
    """
    vortex = import_optional_dependency(
        "vortex", extra="vortex is required for Vortex support."
    )
    pa = import_optional_dependency(
        "pyarrow", extra="pyarrow is required for Vortex support."
    )

    # Convert Path object to string if necessary
    from pathlib import Path
    if isinstance(path, Path):
        path = str(path)

    # Convert DataFrame to Arrow Table
    table = pa.Table.from_pandas(df)

    # Convert Arrow Table to Vortex Array
    v_array = vortex.array(table)

    # Write to file
    vortex.io.write(v_array, path, **kwargs)
