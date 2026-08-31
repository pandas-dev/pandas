"""
Thin adapter around the JSON serializer used by ``to_json`` / ``read_json``.

pandas owns all of the semantics (orients, date formatting, missing values,
label stringification); the engine only turns JSON-native Python objects into
text and back. ``orjson`` is used when installed, the standard library
``json`` module otherwise. Both produce the same values, and the same text
apart from the exponent notation of very small or very large floats.
"""

from __future__ import annotations

import json
import math
from typing import Any

import numpy as np

from pandas._libs import writers
from pandas.compat._optional import import_optional_dependency

orjson = import_optional_dependency("orjson", errors="ignore")


def native_nan() -> bool:
    """Whether the engine writes NaN / infinity as ``null`` itself."""
    return orjson is not None


def native_numpy() -> bool:
    """Whether the engine serializes numpy arrays of int, uint, float64, bool."""
    return orjson is not None


def native_fragments() -> bool:
    """Whether the engine splices pre-encoded JSON text (see ``fragment``)."""
    return orjson is not None


class Fragment:
    """
    Already-encoded JSON text, copied verbatim into the output by ``dumps``.
    """

    __slots__ = ("data",)

    def __init__(self, data: bytes) -> None:
        self.data = data


def fragment(text: bytes) -> Fragment:
    return Fragment(text)


def _orjson_default(obj: Any) -> Any:
    if isinstance(obj, Fragment):
        return orjson.Fragment(obj.data)
    raise TypeError(f"Type is not JSON serializable: {type(obj).__name__}")


def dumps_bytes(obj: Any) -> bytes:
    """
    Compact serialization of ``obj`` as UTF-8 bytes, used to build fragments.
    """
    if isinstance(obj, Fragment):
        return obj.data
    if orjson is not None:
        try:
            return orjson.dumps(
                obj, default=_orjson_default, option=orjson.OPT_SERIALIZE_NUMPY
            )
        except TypeError:
            obj = _sanitize_for_stdlib(obj)
    return json.dumps(obj, ensure_ascii=False, separators=(",", ":")).encode("utf-8")


def dumps(obj: Any, *, indent: int = 0) -> str:
    """
    Serialize ``obj``, which holds only JSON-native Python objects.

    ``obj`` may additionally hold NaN / infinity floats and numpy arrays when
    ``native_nan()`` / ``native_numpy()`` are True.
    """
    if isinstance(obj, Fragment) and not indent:
        return obj.data.decode("utf-8")
    if orjson is not None:
        option = orjson.OPT_SERIALIZE_NUMPY
        if indent:
            option |= orjson.OPT_INDENT_2
        try:
            return orjson.dumps(obj, default=_orjson_default, option=option).decode(
                "utf-8"
            )
        except TypeError:
            # orjson rejects a few things the stdlib encoder accepts: integers
            # beyond 64 bits, nesting deeper than 254 levels, and lone
            # surrogates. Those are rare, so take the slower path for them.
            obj = _sanitize_for_stdlib(obj)
    return json.dumps(
        obj,
        ensure_ascii=False,
        indent=2 if indent else None,
        separators=(",", ": ") if indent else (",", ":"),
    )


def _sanitize_for_stdlib(obj: Any) -> Any:
    """
    Convert numpy arrays to lists and NaN / infinity to None ahead of the
    stdlib encoder, which cannot handle either.
    """
    if isinstance(obj, np.ndarray):
        if obj.dtype.kind == "f":
            mask = ~np.isfinite(obj)
            obj = obj.astype(object)
            obj[mask] = None
        return obj.tolist()
    if isinstance(obj, Fragment):
        return json.loads(obj.data)
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, dict):
        return {key: _sanitize_for_stdlib(value) for key, value in obj.items()}
    if isinstance(obj, list):
        return [_sanitize_for_stdlib(value) for value in obj]
    return obj


def _parse_constant(name: str) -> float | None:
    # ``NaN`` decodes to None so that it becomes a missing value in any dtype
    if name == "NaN":
        return None
    return float(name)


def loads(text: str | bytes) -> Any:
    """
    Deserialize JSON text.

    Beyond strict JSON this accepts the ``NaN``, ``Infinity`` and
    ``-Infinity`` literals. ``NaN`` is read as ``None``; it is rewritten to
    ``null`` ahead of orjson, which rejects it. ``Infinity`` and floats that
    overflow are read by the stdlib decoder, since orjson cannot produce
    infinity. Integers beyond 64 bits raise ``ValueError``.
    """
    data = text.encode("utf-8") if isinstance(text, str) else text
    pos, sign, nan_count, needs_stdlib = writers.json_scan(data)
    if sign > 0:
        raise ValueError(f"Value is too big! at position {pos}")
    if sign < 0:
        raise ValueError(f"Value is too small at position {pos}")
    if orjson is not None and not needs_stdlib:
        if nan_count:
            data = writers.json_nan_to_null(data, nan_count)
        try:
            return orjson.loads(data)
        except orjson.JSONDecodeError:
            # e.g. lone surrogate escapes, which the stdlib decoder accepts
            pass
    return json.loads(data, parse_constant=_parse_constant)
