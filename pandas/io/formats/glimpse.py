"""
Formatting for DataFrame.glimpse.
"""

from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
)

import numpy as np

from pandas._libs import lib

from pandas.io.formats.printing import (
    get_adjustment,
    pprint_thing,
)

if TYPE_CHECKING:
    from collections.abc import Hashable

    from pandas import DataFrame


def _format_value(value: Any, missing: bool) -> str:
    if (lib.is_float(value) or lib.is_complex(value)) and np.isnan(value):
        # match the DataFrame repr, which shows "NaN" rather than "nan"
        return "NaN"
    if missing:
        return str(value)
    if isinstance(value, str):
        # quote so that whitespace and numeric-looking strings stand out
        return repr(str(value))
    # escape so that e.g. a nested DataFrame stays on its column's line
    return pprint_thing(str(value), escape_chars=("\t", "\r", "\n"))


def glimpse_columns(
    frame: DataFrame, max_items_per_column: int
) -> list[tuple[Hashable, str, list[str], list[str | None]]]:
    """
    Collect each column's label, dtype and formatted first values.

    The values are returned twice: as displayed, and with None for missing.
    """
    head = frame.iloc[:max_items_per_column]
    result = []
    for label, ser in head.items():
        mask = ser.isna().tolist()
        # iterate the array, not the Series, to keep numpy scalars (float32
        # formats as 0.1, not 0.10000000149011612)
        texts = [
            _format_value(val, na) for val, na in zip(ser.array, mask, strict=True)
        ]
        values = [None if na else text for text, na in zip(texts, mask, strict=True)]
        result.append((label, str(ser.dtype), texts, values))
    return result


def format_glimpse(
    shape: tuple[int, int],
    columns: list[tuple[Hashable, str, list[str], list[str | None]]],
    max_colname_length: int,
) -> str:
    """
    Render the output of ``glimpse_columns`` with one line per column.
    """
    lines = [f"Rows: {shape[0]}", f"Columns: {shape[1]}"]
    if columns:
        adj = get_adjustment()
        names = []
        for label, *_ in columns:
            name = pprint_thing(label, escape_chars=("\t", "\r", "\n"))
            if len(name) > max_colname_length:
                name = name[: max_colname_length - 1] + "…"
            names.append(name)
        dtypes = [f"<{dtype}>" for _, dtype, _, _ in columns]
        names = adj.justify(names, max(map(adj.len, names)), mode="left")
        dtypes = adj.justify(dtypes, max(map(adj.len, dtypes)), mode="right")
        for name, dtype, (_, _, texts, _) in zip(names, dtypes, columns, strict=True):
            lines.append(f"$ {name} {dtype} {', '.join(texts)}".rstrip())
    return "\n".join(lines) + "\n"
