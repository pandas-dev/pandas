"""Rich render helpers for b2view."""

from __future__ import annotations

from pprint import pformat
from textwrap import wrap
from typing import Any

import numpy as np


def make_metadata_renderable(info):
    """Return a Rich renderable for ObjectInfo metadata."""
    from rich.table import Table

    table = Table(show_header=False, box=None, expand=True)
    table.add_column("key", style="bold cyan", no_wrap=True)
    table.add_column("value")
    table.add_row("path", info.path)
    table.add_row("kind", info.kind)
    for key, value in info.metadata.items():
        table.add_row(str(key), _format_metadata_value(value))
    return table


def make_preview_renderable(preview: Any):
    """Return a single Rich renderable for a preview object."""
    _, body = make_preview_renderables(preview)
    return body


def make_preview_renderables(preview: Any):
    """Return ``(header, body)`` Rich renderables for a preview object.

    CTable previews get a separate header renderable so the UI can keep column
    titles fixed while only the row body scrolls. Other preview kinds return
    ``None`` for the header.
    """
    from rich.pretty import Pretty
    from rich.table import Table
    from rich.text import Text

    if isinstance(preview, np.ndarray):
        return None, Text(np.array2string(preview, threshold=200, edgeitems=5), no_wrap=False)

    if isinstance(preview, dict) and "data" in preview and "columns" in preview:
        widths = _preview_column_widths(preview)
        header = _make_ctable_header(preview, widths)
        body = Table(expand=True, show_header=False, show_lines=False)
        for name in preview["columns"]:
            body.add_column(name, width=widths[name], overflow="fold")
        nrows = preview["stop"] - preview["start"]
        for i in range(nrows):
            body.add_row(*[_format_cell(preview["data"][name][i]) for name in preview["columns"]])
        if preview.get("hidden_columns", 0):
            body.caption = f"{preview['hidden_columns']} columns hidden"
        return header, body

    if isinstance(preview, dict) and "message" in preview:
        return None, Text(str(preview["message"]))

    return None, Pretty(preview)


def _make_ctable_header(preview: dict[str, Any], widths: dict[str, int]):
    from rich.align import Align
    from rich.console import Group
    from rich.text import Text

    title = Align.center(Text(f"rows {preview['start']}:{preview['stop']} of {preview['nrows']}"))
    wrapped_columns = []
    for name in preview["columns"]:
        width = widths[name]
        parts = wrap(name, width=width, break_long_words=True, break_on_hyphens=False) or [""]
        wrapped_columns.append(parts)
    height = max(len(parts) for parts in wrapped_columns) if wrapped_columns else 0
    lines = []
    for row in range(height):
        cells = []
        for name, parts in zip(preview["columns"], wrapped_columns, strict=True):
            width = widths[name]
            text = parts[row] if row < len(parts) else ""
            cells.append(f" {text:<{width}} ")
        lines.append("│".join(cells))
    return Group(title, Text("\n".join(lines)))


def _preview_column_widths(preview: dict[str, Any], *, max_width: int = 40) -> dict[str, int]:
    widths = {}
    nrows = preview["stop"] - preview["start"]
    for name in preview["columns"]:
        values = preview["data"][name]
        width = len(name)
        for i in range(nrows):
            width = max(width, min(max_width, len(_format_cell(values[i]))))
        widths[name] = min(max_width, max(4, width))
    return widths


def _format_metadata_value(value: Any) -> str:
    if isinstance(value, dict):
        return "\n".join(f"{key}: {val}" for key, val in value.items()) or "{}"
    if isinstance(value, (list, tuple)):
        return repr(value)
    return str(value)


def column_float_decimals(values: Any) -> int | None:
    """Return a uniform decimal count for a float column, or None.

    The count derives from the column's maximum magnitude so every cell fits
    the same ~9 character budget that _fmt_float uses per value: digits move
    from the fraction to the integer part as magnitudes grow, but uniformly
    for the whole column, keeping the decimal points aligned.

    Returns None when *values* is not a float column or when its magnitude
    calls for scientific notation (handled per value by _fmt_float).
    """
    arr = np.asarray(values)
    if arr.dtype.kind != "f" or arr.size == 0:
        return None
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return None
    largest = float(np.max(np.abs(finite)))
    if largest == 0:
        return 1  # an all-zero column reads best as plain 0.0
    if largest >= 1e9 or largest < 1e-6:
        return None
    int_digits = max(1, int(np.floor(np.log10(largest))) + 1)
    # 9-char budget: sign/pad + int digits + decimal point + decimals
    return max(0, 7 - int_digits)


def format_cell(value: Any, *, float_decimals: int | None = None) -> str:
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, np.ndarray):
        text = np.array2string(value, threshold=20, formatter={"float_kind": lambda x: _fmt_float(x)})
    elif isinstance(value, (list, tuple, dict)):
        text = pformat(value, compact=True, width=80)
    elif isinstance(value, float):
        text = _fmt_float(value) if float_decimals is None else f"{value:9.{float_decimals}f}"
    else:
        text = str(value)
    text = " ".join(text.splitlines())
    return text if len(text) <= 200 else text[:197] + "..."


def _fmt_float(x: float) -> str:
    """Show floats with a fixed width of 9 characters and up to 6 decimal digits, right-aligned."""
    if abs(x) >= 1e9 or (abs(x) < 1e-6 and abs(x) > 0):
        return f"{x: .6e}"
    if abs(x) == 0:
        return " 0.0"
    abs_x = abs(x)
    # Choose format to keep total width ~9 chars including leading space for sign
    if abs_x < 10:
        return f"{x:9.6f}"[:9]
    if abs_x < 1000:
        return f"{x:9.3f}"[:9]
    if abs_x < 1e6:
        return f"{x:9.0f}"[:9]
    return f"{x:9.0f}"[:9]


_format_cell = format_cell
