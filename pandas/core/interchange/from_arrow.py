"""
Reconstruct a pandas :class:`DataFrame` from a ``pyarrow.Table``.

This is the pandas-side home for the "pandas metadata" handling that has
historically lived in pyarrow's ``pandas_compat.py`` (GH#59780). The low-level
conversion of each column's memory is still delegated to pyarrow; what lives
here is the pandas-specific reconstruction of the index, the column labels and
(eventually) the dtypes from the ``b"pandas"`` schema metadata, so that pandas
controls the arrow -> pandas conversion.

This is an initial, deliberately small implementation: it covers the common
cases (a default/serialized :class:`Index`, a :class:`RangeIndex` stored as
metadata, and string column labels) and leaves the more involved cases
(``MultiIndex`` columns, full ``numpy_type`` dtype restoration, ``attrs``) as
follow-ups, see the ``TODO`` markers below.
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pyarrow as pa

    from pandas import (
        DataFrame,
        Index,
    )

_INDEX_LEVEL_PATTERN = re.compile(r"^__index_level_\d+__$")


def _index_name(field_name: str) -> str | None:
    """
    User-facing name for a serialized index level.

    Auto-generated ``__index_level_N__`` names correspond to an unnamed index,
    so they map back to ``None``.
    """
    if _INDEX_LEVEL_PATTERN.match(field_name):
        return None
    return field_name


def _reconstruct_index(
    table: pa.Table,
    index_descriptors: list,
    columns: dict[str, pa.ChunkedArray],
) -> tuple[Index, set[str]]:
    """
    Build the row :class:`Index` from the pandas metadata index descriptors.

    Returns the reconstructed index and the set of field names that were
    consumed as index levels (so the caller can drop them from the data).
    """
    from pandas import (
        Index,
        MultiIndex,
        RangeIndex,
    )

    arrays: list = []
    names: list[str | None] = []
    consumed: set[str] = set()

    for descriptor in index_descriptors:
        if isinstance(descriptor, str):
            # a serialized index level stored as an actual column in the table
            level = columns[descriptor].to_pandas()
            # clear the arrow field name so it does not leak into the index name
            # (the user-facing name is derived from the descriptor below)
            level.name = None
            arrays.append(level)
            names.append(_index_name(descriptor))
            consumed.add(descriptor)
        elif isinstance(descriptor, dict) and descriptor.get("kind") == "range":
            # a RangeIndex stored purely as metadata (no backing column)
            return (
                RangeIndex(
                    start=descriptor["start"],
                    stop=descriptor["stop"],
                    step=descriptor["step"],
                    name=descriptor["name"],
                ),
                consumed,
            )

    if not arrays:
        return RangeIndex(table.num_rows), consumed
    if len(arrays) == 1:
        return Index(arrays[0], name=names[0]), consumed
    return MultiIndex.from_arrays(arrays, names=names), consumed


def table_to_dataframe(table: pa.Table) -> DataFrame:
    """
    Convert a ``pyarrow.Table`` to a pandas :class:`DataFrame`.

    Parameters
    ----------
    table : pyarrow.Table

    Returns
    -------
    DataFrame
    """
    from pandas import DataFrame

    columns = {name: table.column(i) for i, name in enumerate(table.column_names)}
    metadata = table.schema.pandas_metadata

    if metadata is None:
        # No pandas metadata: straight column-wise conversion with a default index.
        data = {name: col.to_pandas() for name, col in columns.items()}
        return DataFrame(data)

    index_descriptors = metadata.get("index_columns", [])
    index, index_fields = _reconstruct_index(table, index_descriptors, columns)

    # data columns are everything that was not consumed as an index level
    data = {
        name: col.to_pandas()
        for name, col in columns.items()
        if name not in index_fields
    }
    # to_pandas() yields Series on a default RangeIndex; assemble then attach the
    # reconstructed index (all columns share that same default index, so this is
    # an alignment-free assignment).
    result = DataFrame(data)
    result.index = index

    # TODO(GH#59780): restore original (possibly non-string / MultiIndex) column
    # labels and dtypes from ``metadata["columns"]`` / ``metadata["column_indexes"]``,
    # and restore ``DataFrame.attrs`` from ``metadata["attributes"]``.
    return result
