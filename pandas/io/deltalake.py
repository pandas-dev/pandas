""" Delta Lake support """
from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
)

from pandas.compat._optional import import_optional_dependency

if TYPE_CHECKING:
    import pyarrow.fs as pa_fs

    from pandas import DataFrame


def _try_import():
    # since pandas is a dependency of deltalake
    # we need to import on first use
    msg = (
        "deltalake is required to load data from Delta Lake. "
        "See the docs: https://delta-io.github.io/delta-rs/python."
    )
    deltalake = import_optional_dependency("deltalake", extra=msg)
    return deltalake


def read_deltalake(
    table_uri: str,
    version: int | None = None,
    storage_options: dict[str, str] | None = None,
    without_files: bool = False,
    partitions: list[tuple[str, str, Any]] | None = None,
    columns: list[str] | None = None,
    filesystem: str | pa_fs.FileSystem | None = None,
) -> DataFrame:
    """
    Load data from Deltalake.

    This function requires the `deltalake package
    <https://delta-io.github.io/delta-rs/python>`__.

    See the `How to load a Delta table
    <https://delta-io.github.io/delta-rs/python/usage.html#loading-a-delta-table>`__
    guide for loading instructions.

    Parameters
    ----------
    table_uri: str
        The path of the DeltaTable.
    version: int, optional
        The version of the DeltaTable.
    storage_options: Dict[str, str], optional
        A dictionary of the options to use for the storage backend.
    without_files: bool, default False
        If True, will load table without tracking files.
        Some append-only applications might have no need of tracking any files.
        So, the DeltaTable will be loaded with a significant memory reduction.
    partitions: List[Tuple[str, str, Any], optional
        A list of partition filters, see help(DeltaTable.files_by_partitions)
        for filter syntax.
    columns: List[str], optional
        The columns to project. This can be a list of column names to include
        (order and duplicates will be preserved).
    filesystem: Union[str, pa_fs.FileSystem], optional
        A concrete implementation of the Pyarrow FileSystem or
        a fsspec-compatible interface. If None, the first file path will be used
        to determine the right FileSystem.

    Returns
    -------
    df: DataFrame
        DataFrame including the results.

    See Also
    --------
    deltalake.DeltaTable : Create a DeltaTable instance with the deltalake library.
    """
    deltalake = _try_import()

    table = deltalake.DeltaTable(
        table_uri=table_uri,
        version=version,
        storage_options=storage_options,
        without_files=without_files,
    )
    return table.to_pandas(
        partitions=partitions, columns=columns, filesystem=filesystem
    )
