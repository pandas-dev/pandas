#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# This source code is licensed under a BSD-style license (found in the
# LICENSE file in the root directory of this source tree)
#######################################################################

"""Schema compiler: turns a dataclass row definition into a CompiledSchema."""

from __future__ import annotations

import base64
import copy
import dataclasses
import datetime
import typing
from dataclasses import MISSING
from typing import Any

import numpy as np

from blosc2.schema import (
    BLOSC2_FIELD_METADATA_KEY,
    DictionarySpec,
    ListSpec,
    NDArraySpec,
    ObjectSpec,
    SchemaSpec,
    StructSpec,
    Utf8Spec,
    VLBytesSpec,
    VLStringSpec,
    complex64,
    complex128,
    float32,
    float64,
    int8,
    int16,
    int32,
    int64,
    string,
    timestamp,
    uint8,
    uint16,
    uint32,
    uint64,
)
from blosc2.schema import (
    bool as b2_bool,
)
from blosc2.schema import (
    bytes as b2_bytes,
)

# Maps the "kind" string used in serialized dicts back to spec constructors.
_KIND_TO_SPEC: dict[str, type[SchemaSpec]] = {
    # signed integers
    "int8": int8,
    "int16": int16,
    "int32": int32,
    "int64": int64,
    # unsigned integers
    "uint8": uint8,
    "uint16": uint16,
    "uint32": uint32,
    "uint64": uint64,
    # floats
    "float32": float32,
    "float64": float64,
    # complex
    "complex64": complex64,
    "complex128": complex128,
    # bool / string / bytes / varlen
    "bool": b2_bool,
    "string": string,
    "bytes": b2_bytes,
    "vlstring": VLStringSpec,
    "vlbytes": VLBytesSpec,
    "utf8": Utf8Spec,
    "object": ObjectSpec,
    "timestamp": timestamp,
    # dictionary
    "dictionary": DictionarySpec,
    # fixed-shape N-D arrays
    "ndarray": NDArraySpec,
}

# ---------------------------------------------------------------------------
# Display-width helper (used by CTable.__str__ / info())
# ---------------------------------------------------------------------------

_DTYPE_DISPLAY_WIDTH: dict[str, int] = {
    "int8": 6,
    "int16": 8,
    "int32": 10,
    "int64": 12,
    "uint8": 6,
    "uint16": 8,
    "uint32": 10,
    "uint64": 12,
    "float32": 12,
    "float64": 15,
    "bool": 6,
    "complex64": 20,
    "complex128": 25,
    "timestamp": 26,
}


def compute_display_width(spec: SchemaSpec) -> int:
    """Return a reasonable terminal display width for *spec*'s column."""
    if isinstance(spec, DictionarySpec):
        return 32
    if isinstance(spec, (VLStringSpec, VLBytesSpec, ObjectSpec, Utf8Spec)):
        return 40
    if isinstance(spec, NDArraySpec):
        return max(20, len(spec.display_label()) + 4)
    if isinstance(spec, (ListSpec, StructSpec)):
        return max(40, len(spec.display_label()) + 4)
    if isinstance(spec, timestamp):
        return _DTYPE_DISPLAY_WIDTH["timestamp"]
    dtype = spec.dtype
    if dtype is None:
        return 20
    if dtype.kind == "U":  # fixed-width unicode (string spec)
        return max(10, min(dtype.itemsize // 4, 50))
    if dtype.kind == "S":  # fixed-width bytes
        return max(10, min(dtype.itemsize, 50))
    return _DTYPE_DISPLAY_WIDTH.get(dtype.name, 20)


# ---------------------------------------------------------------------------
# Mapping from Python primitive annotations to default spec constructors.
# Keys are the actual builtin types (bool before int because bool <: int).
# ---------------------------------------------------------------------------
_ANNOTATION_TO_SPEC: dict[type, type[SchemaSpec]] = {
    bool: b2_bool,  # must come before int (bool is a subclass of int)
    int: int64,
    float: float64,
    complex: complex128,
    str: string,
    bytes: b2_bytes,
}


# ---------------------------------------------------------------------------
# Compiled representations
# ---------------------------------------------------------------------------


@dataclasses.dataclass(slots=True)
class ColumnConfig:
    """Per-column NDArray storage options."""

    cparams: dict[str, Any] | None
    dparams: dict[str, Any] | None
    chunks: tuple[int, ...] | None
    blocks: tuple[int, ...] | None


@dataclasses.dataclass(slots=True)
class CompiledColumn:
    """All compile-time information about a single CTable column."""

    name: str
    py_type: Any
    spec: SchemaSpec
    dtype: np.dtype | None
    default: Any  # MISSING means no default declared
    config: ColumnConfig
    display_width: int = 20  # terminal column width for __str__ / info()


@dataclasses.dataclass(slots=True)
class CompiledSchema:
    """Compiled representation of a CTable row schema.

    Built once per row class by :func:`compile_schema` and cached on the
    ``CTable`` instance.  Drives NDArray creation, row validation, and
    future schema serialization.
    """

    row_cls: type[Any]
    columns: list[CompiledColumn]
    columns_by_name: dict[str, CompiledColumn]
    validator_model: type[Any] | None = None  # filled in by schema_validation
    metadata: dict[str, Any] = dataclasses.field(default_factory=dict)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def get_blosc2_field_metadata(dc_field: dataclasses.Field) -> dict[str, Any] | None:
    """Return the ``blosc2`` metadata dict stored on a dataclass field, or ``None``."""
    return dc_field.metadata.get(BLOSC2_FIELD_METADATA_KEY)


def infer_spec_from_annotation(annotation: Any) -> SchemaSpec:
    """Build a default :class:`SchemaSpec` from a plain Python type annotation.

    Supports ``bool``, ``int``, ``float``, ``str``, and ``bytes``.

    Raises
    ------
    TypeError
        If the annotation is not one of the supported primitive types.
    """
    spec_cls = _ANNOTATION_TO_SPEC.get(annotation)
    if spec_cls is None:
        raise TypeError(
            f"Cannot infer a Blosc2 schema spec from annotation {annotation!r}. "
            f"Use b2.field(b2.<type>(...)) to declare this column explicitly."
        )
    return spec_cls()


def validate_annotation_matches_spec(name: str, annotation: Any, spec: SchemaSpec) -> None:
    """Raise :exc:`TypeError` if *annotation* is incompatible with *spec*."""
    if isinstance(spec, NDArraySpec):
        if annotation not in (np.ndarray, object):
            raise TypeError(
                f"Column {name!r}: annotation {annotation!r} is incompatible with "
                "NDArraySpec (expected np.ndarray or object)."
            )
        return

    if isinstance(spec, ListSpec):
        origin = typing.get_origin(annotation)
        if origin not in (list, list):
            raise TypeError(
                f"Column {name!r}: annotation {annotation!r} is incompatible with list spec; "
                "expected list[T]."
            )
        args = typing.get_args(annotation)
        if len(args) != 1:
            raise TypeError(f"Column {name!r}: list annotations must specify exactly one item type.")
        item_annotation = args[0]
        expected = spec.item_spec.python_type
        if item_annotation is not expected:
            raise TypeError(
                f"Column {name!r}: list item annotation {item_annotation!r} is incompatible with "
                f"item spec {type(spec.item_spec).__name__!r} (expected {expected.__name__!r})."
            )
        return

    if isinstance(spec, DictionarySpec):
        if annotation is not str:
            raise TypeError(
                f"Column {name!r}: annotation {annotation!r} is incompatible with "
                f"DictionarySpec (expected str)."
            )
        return

    if isinstance(spec, timestamp):
        if annotation in (object, np.datetime64, datetime.datetime, str, int):
            return
        raise TypeError(
            f"Column {name!r}: annotation {annotation!r} is incompatible with "
            "timestamp spec (expected object, np.datetime64, datetime.datetime, str, or int)."
        )

    # VLStringSpec and VLBytesSpec are only reachable via blosc2.field(blosc2.vlstring()/vlbytes()),
    # so annotation must match python_type (str or bytes respectively).
    expected = spec.python_type
    if annotation is not expected:
        raise TypeError(
            f"Column {name!r}: annotation {annotation!r} is incompatible with "
            f"spec {type(spec).__name__!r} (expected Python type {expected.__name__!r})."
        )


# ---------------------------------------------------------------------------
# Public compiler entry point
# ---------------------------------------------------------------------------


_RESERVED_COLUMN_NAMES: frozenset[str] = frozenset({"_meta", "_valid_rows", "_cols", "_indexes"})


def _validate_column_name(name: str) -> None:
    """Raise :exc:`ValueError` if *name* is not a legal CTable column name.

    Rules (enforced for both in-memory and persistent tables so that an
    in-memory schema can always be persisted without surprises):

    * must be a non-empty string
    * must not start with ``_``  (reserved for internal table layout)
    * must not be one of the reserved internal names

    Literal ``/`` characters are allowed in logical names; persistent CTable
    storage percent-encodes path segments before writing under ``_cols``.
    """
    if not name:
        raise ValueError("Column name cannot be empty.")
    if name.startswith("_"):
        raise ValueError(f"Column name cannot start with '_' (reserved for internal use): {name!r}")
    if name in _RESERVED_COLUMN_NAMES:
        raise ValueError(f"Column name {name!r} is reserved for internal CTable use.")


def compile_schema(row_cls: type[Any]) -> CompiledSchema:
    """Compile *row_cls* (a dataclass) into a :class:`CompiledSchema`.

    Parameters
    ----------
    row_cls:
        A class decorated with ``@dataclass``.  Each field must either carry a
        ``b2.field(...)`` default or use a supported plain annotation
        (``int``, ``float``, ``bool``, ``str``, ``bytes``).

    Returns
    -------
    CompiledSchema

    Raises
    ------
    TypeError
        If *row_cls* is not a dataclass, if a field spec is incompatible with
        its annotation, or if an unsupported annotation is encountered.
    ValueError
        If any column name violates the naming rules.
    """
    if not dataclasses.is_dataclass(row_cls) or not isinstance(row_cls, type):
        raise TypeError(
            f"{row_cls!r} is not a dataclass type. CTable row schemas must be defined with @dataclass."
        )

    # Resolve string annotations (handles `from __future__ import annotations`)
    try:
        hints = typing.get_type_hints(row_cls)
    except Exception as exc:
        raise TypeError(f"Could not resolve type hints for {row_cls!r}: {exc}") from exc

    columns: list[CompiledColumn] = []

    for dc_field in dataclasses.fields(row_cls):
        name = dc_field.name
        _validate_column_name(name)
        annotation = hints.get(name, dc_field.type)
        meta = get_blosc2_field_metadata(dc_field)

        if meta is not None:
            # Explicit b2.field(...) path
            spec = copy.deepcopy(meta["spec"])
            if not isinstance(spec, SchemaSpec):
                raise TypeError(
                    f"Column {name!r}: b2.field() requires a SchemaSpec as its first "
                    f"argument, got {type(spec)!r}."
                )
            validate_annotation_matches_spec(name, annotation, spec)
            config = ColumnConfig(
                cparams=meta.get("cparams"),
                dparams=meta.get("dparams"),
                chunks=meta.get("chunks"),
                blocks=meta.get("blocks"),
            )
        else:
            # Inferred shorthand: plain annotation without b2.field()
            spec = infer_spec_from_annotation(annotation)
            config = ColumnConfig(cparams=None, dparams=None, chunks=None, blocks=None)

        # Resolve default value
        if dc_field.default is not MISSING:
            default = dc_field.default
        elif dc_field.default_factory is not MISSING:  # type: ignore[misc]
            default = dc_field.default_factory
        else:
            default = MISSING

        columns.append(
            CompiledColumn(
                name=name,
                py_type=object if isinstance(spec, (timestamp, NDArraySpec)) else annotation,
                spec=spec,
                dtype=getattr(spec, "dtype", None),
                default=default,
                config=config,
                display_width=compute_display_width(spec),
            )
        )

    return CompiledSchema(
        row_cls=row_cls,
        columns=columns,
        columns_by_name={col.name: col for col in columns},
    )


# ---------------------------------------------------------------------------
# Schema serialization helpers  (Step 12 — persistence groundwork)
# ---------------------------------------------------------------------------


def _bytes_to_json(value: bytes) -> dict[str, Any]:
    return {"__bytes__": True, "base64": base64.b64encode(value).decode("ascii")}


def _json_to_bytes(value: dict[str, Any]) -> bytes:
    return base64.b64decode(value["base64"])


def _default_to_json(value: Any) -> Any:
    """Convert a declared field default to a JSON-compatible value."""
    if value is MISSING:
        raise ValueError("Cannot serialize MISSING as a declared default.")
    if isinstance(value, complex):
        return {"__complex__": True, "real": value.real, "imag": value.imag}
    if isinstance(value, bytes):
        return _bytes_to_json(value)
    return value


def _default_from_json(value: Any) -> Any:
    """Reverse of :func:`_default_to_json`."""
    if isinstance(value, dict) and value.get("__complex__"):
        return complex(value["real"], value["imag"])
    if isinstance(value, dict) and value.get("__bytes__"):
        return _json_to_bytes(value)
    return value


def spec_from_metadata_dict(data: dict[str, Any]) -> SchemaSpec:
    """Reconstruct one SchemaSpec from serialized metadata."""
    data = dict(data)
    kind = data.pop("kind")
    if isinstance(data.get("null_value"), dict) and data["null_value"].get("__bytes__"):
        data["null_value"] = _json_to_bytes(data["null_value"])
    if kind == "list":
        item_spec = spec_from_metadata_dict(data.pop("item"))
        return ListSpec(item_spec, **data)
    if kind == "struct":
        return StructSpec.from_metadata_dict({"fields": data.pop("fields"), **data})
    if kind == "dictionary":
        index_type = spec_from_metadata_dict(data.pop("index_type"))
        value_type = spec_from_metadata_dict(data.pop("value_type"))
        return DictionarySpec(index_type=index_type, value_type=value_type, **data)
    if kind == "ndarray":
        return NDArraySpec(
            item_shape=tuple(data.pop("item_shape")), dtype=np.dtype(data.pop("dtype_str")), **data
        )
    spec_cls = _KIND_TO_SPEC.get(kind)
    if spec_cls is None:
        raise ValueError(f"Unknown column kind {kind!r}")
    return spec_cls(**data)


def schema_to_dict(schema: CompiledSchema) -> dict[str, Any]:
    """Serialize *schema* to a JSON-compatible dict.

    The result is self-contained: it can be stored as table metadata and
    later passed to :func:`schema_from_dict` to reconstruct the schema
    without the original Python dataclass.

    Example output::

        {
            "version": 1,
            "row_cls": "Row",
            "columns": [
                {"name": "id", "kind": "int64", "ge": 0},
                {"name": "score", "kind": "float64", "ge": 0, "le": 100, "default": 0.0},
                {"name": "active", "kind": "bool", "default": true},
            ]
        }
    """
    cols = []
    for col in schema.columns:
        entry: dict[str, Any] = {"name": col.name}
        entry.update(col.spec.to_metadata_dict())  # adds "kind" + constraints
        if isinstance(entry.get("null_value"), bytes):
            entry["null_value"] = _bytes_to_json(entry["null_value"])
        if col.default is not MISSING:
            entry["default"] = _default_to_json(col.default)
        if col.config.cparams is not None:
            entry["cparams"] = col.config.cparams
        if col.config.dparams is not None:
            entry["dparams"] = col.config.dparams
        if col.config.chunks is not None:
            entry["chunks"] = list(col.config.chunks)
        if col.config.blocks is not None:
            entry["blocks"] = list(col.config.blocks)
        cols.append(entry)

    schema_version = 2 if schema.metadata.get("nested") is not None else 1
    result = {
        "version": schema_version,
        "columns": cols,
    }
    if schema.metadata:
        result["metadata"] = schema.metadata
    return result


def schema_from_dict(data: dict[str, Any]) -> CompiledSchema:
    """Reconstruct a :class:`CompiledSchema` from a dict produced by
    :func:`schema_to_dict`.

    The original Python dataclass is *not* required.  ``row_cls`` on the
    returned schema will be ``None``.

    Raises
    ------
    ValueError
        If *data* uses an unknown schema version or an unknown column kind.
    """
    version = data.get("version", 1)
    if version not in (1, 2):
        raise ValueError(f"Unsupported schema version {version!r}")

    columns: list[CompiledColumn] = []
    for entry in data["columns"]:
        entry = dict(entry)  # don't mutate caller's data
        name = entry.pop("name")
        kind = entry.pop("kind")
        default = _default_from_json(entry.pop("default")) if "default" in entry else MISSING
        cparams = entry.pop("cparams", None)
        dparams = entry.pop("dparams", None)
        chunks = tuple(entry.pop("chunks")) if "chunks" in entry else None
        blocks = tuple(entry.pop("blocks")) if "blocks" in entry else None

        spec = spec_from_metadata_dict({"kind": kind, **entry})

        columns.append(
            CompiledColumn(
                name=name,
                py_type=spec.python_type,
                spec=spec,
                dtype=getattr(spec, "dtype", None),
                default=default,
                config=ColumnConfig(cparams=cparams, dparams=dparams, chunks=chunks, blocks=blocks),
                display_width=compute_display_width(spec),
            )
        )

    return CompiledSchema(
        row_cls=None,
        columns=columns,
        columns_by_name={col.name: col for col in columns},
        metadata=dict(data.get("metadata", {})),
    )
