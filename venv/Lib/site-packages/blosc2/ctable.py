#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# This source code is licensed under a BSD-style license (found in the
# LICENSE file in the root directory of this source tree)
#######################################################################

"""CTable: a columnar compressed table built on top of blosc2.NDArray."""

from __future__ import annotations

import ast
import contextlib
import contextvars
import copy
import dataclasses
import json
import operator
import os
import pprint
import re
import shutil
from collections import deque, namedtuple
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import MISSING, dataclass
from dataclasses import field as dataclass_field
from textwrap import TextWrapper
from typing import TYPE_CHECKING, Any, Generic, Literal, TypeVar

import numpy as np

import blosc2
from blosc2 import compute_chunks_blocks
from blosc2.ctable_indexing import _CTableIndexingMixin
from blosc2.ctable_storage import (
    FileTableStorage,
    InMemoryTableStorage,
    TableStorage,
    TreeStoreTableStorage,
    join_field_path,
    split_field_path,
)
from blosc2.info import InfoReporter, format_nbytes_human, format_nbytes_info
from blosc2.list_array import ListArray, coerce_list_cell
from blosc2.scalar_array import _ScalarVarLenArray

if TYPE_CHECKING:
    from blosc2.dictionary_column import DictionaryColumn
from blosc2.schema import (
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
from blosc2.schema_compiler import (
    ColumnConfig,
    CompiledColumn,
    CompiledSchema,
    _validate_column_name,
    compile_schema,
    compute_display_width,
    get_blosc2_field_metadata,
    schema_from_dict,
    schema_to_dict,
)


def _is_arrow_string_type(pa, pa_type) -> bool:
    """True for any Arrow string-like type, including the view-based layout.

    ``string_view`` (Arrow's variable-length view layout, e.g. Polars'
    default export type for string columns via the PyCapsule interface) is
    not one of ``string``/``large_string``/``utf8``/``large_utf8`` but is
    handled identically everywhere those are.
    """
    if pa_type in (pa.string(), pa.large_string(), pa.utf8(), pa.large_utf8()):
        return True
    is_string_view = getattr(pa.types, "is_string_view", None)
    return bool(is_string_view is not None and is_string_view(pa_type))


def _is_arrow_binary_type(pa, pa_type) -> bool:
    """True for any Arrow binary-like type, including the view-based layout."""
    if pa.types.is_binary(pa_type) or pa.types.is_large_binary(pa_type):
        return True
    is_binary_view = getattr(pa.types, "is_binary_view", None)
    return bool(is_binary_view is not None and is_binary_view(pa_type))


@dataclass(frozen=True)
class NullPolicy:
    """Default sentinels for inferred CTable scalar nulls.

    CTable nullable scalar columns are represented with per-column sentinel
    values. This policy is used when CTable has to infer those sentinels, such
    as when importing nullable scalar Arrow or Parquet columns without an
    explicit column-level null sentinel. The selected sentinel is stored in the
    resulting CTable schema, so existing tables remain self-describing.

    Examples
    --------
    Use :func:`blosc2.null_policy` to apply a policy while creating a CTable
    from data with nullable scalar columns::

        policy = blosc2.NullPolicy(
            signed_int_strategy="max",
            string_value="<NULL>",
            column_null_values={"user_id": -1, "country": "NA"},
        )

        with blosc2.null_policy(policy):
            table = blosc2.CTable.from_parquet("data.parquet")

    The same policy is used for explicit nullable schema specs::

        @dataclass
        class Row:
            user_id: int = blosc2.field(blosc2.int64(nullable=True))
            country: str = blosc2.field(blosc2.string(nullable=True))

        with blosc2.null_policy(policy):
            table = blosc2.CTable(Row)

    ``column_null_values`` takes precedence over the type-wide defaults in the
    policy.  This is useful when a particular column needs a sentinel that is
    known not to collide with its real values.
    """

    string_value: str = "__BLOSC2_NULL__"
    bytes_value: bytes = b"__BLOSC2_NULL__"
    float_value: float = float("nan")
    bool_value: int = 255
    signed_int_strategy: Literal["min", "max"] = "min"
    unsigned_int_strategy: Literal["min", "max"] = "max"
    timestamp_value: int = int(np.iinfo(np.int64).min)
    column_null_values: Mapping[str, Any] = dataclass_field(default_factory=dict)

    def sentinel_for_arrow_type(self, pa, pa_type):
        """Return the default sentinel for *pa_type*, or ``None`` if unsupported."""
        signed_ints = [
            (pa.int8(), np.int8),
            (pa.int16(), np.int16),
            (pa.int32(), np.int32),
            (pa.int64(), np.int64),
        ]
        unsigned_ints = [
            (pa.uint8(), np.uint8),
            (pa.uint16(), np.uint16),
            (pa.uint32(), np.uint32),
            (pa.uint64(), np.uint64),
        ]
        for arrow_type, dtype in signed_ints:
            if pa_type == arrow_type:
                info = np.iinfo(dtype)
                return info.min if self.signed_int_strategy == "min" else info.max
        for arrow_type, dtype in unsigned_ints:
            if pa_type == arrow_type:
                info = np.iinfo(dtype)
                return info.min if self.unsigned_int_strategy == "min" else info.max
        if pa_type in (pa.float32(), pa.float64()):
            return self.float_value
        if pa_type == pa.bool_():
            return self.bool_value
        if _is_arrow_string_type(pa, pa_type):
            return self.string_value
        if _is_arrow_binary_type(pa, pa_type):
            return self.bytes_value
        if pa.types.is_timestamp(pa_type):
            return self.timestamp_value
        return None


DEFAULT_NULL_POLICY = NullPolicy()
_NULL_POLICY = contextvars.ContextVar("blosc2_null_policy", default=DEFAULT_NULL_POLICY)
# Sentinel for set_printoptions params whose valid value includes ``None``
# (so ``None`` can be set explicitly rather than meaning "leave unchanged").
_UNSET = object()

_CTABLE_PRINT_OPTIONS: dict[str, Any] = {
    "display_index": True,
    "display_rows": 60,
    "display_width": None,
    "display_precision": 6,
    "fancy": False,
}
_SMALL_NROWS_LIMIT = 10_000_000
_SMALL_SORT_MATERIALIZE_LIMIT = _SMALL_NROWS_LIMIT
_MAX_GROWTH_ROWS = 1_048_576


def get_null_policy() -> NullPolicy:
    """Return the current default null policy."""
    return _NULL_POLICY.get()


def set_printoptions(
    *,
    display_index: bool | None = None,
    display_rows: int | None = None,
    display_width: int | None = _UNSET,
    display_precision: int | None = None,
    fancy: bool | None = None,
) -> None:
    """Set global display options for :class:`CTable` string representations.

    These options affect ``str(ctable)``/``repr(ctable)``/``print(ctable)`` (the
    interactive, truncated view).  They do *not* affect :meth:`CTable.to_string`,
    which renders everything by default.

    Parameters
    ----------
    display_index:
        Whether the display should include a pandas-like logical row index
        column.  ``None`` leaves the current setting unchanged.
    display_rows:
        Maximum number of rows shown before truncating to a compact head/tail
        view (five first and five last rows, when possible).  ``-1`` shows all
        rows, ``0`` shows none.  ``None`` leaves the current setting unchanged.
    display_width:
        Character budget used to decide how many columns fit before truncating
        the middle ones with ``...``.  ``None`` (the default) auto-detects the
        terminal width, ``-1`` shows all columns, a positive int sets a fixed
        budget.  Omit the argument to leave the current setting unchanged.
    display_precision:
        Number of digits after the decimal point for floating-point values in
        table displays.  Trailing zeros are trimmed.  ``None`` leaves the
        current setting unchanged.
    fancy:
        Whether to use the more decorated table display, including separator
        rules and a detailed footer.  ``False`` (default) uses a simpler
        pandas-like footer such as ``[726017 rows x 5 columns]`` and omits
        separator rules.  ``None`` leaves the current setting unchanged.
    """
    if display_index is not None:
        if not isinstance(display_index, bool):
            raise TypeError("display_index must be a bool or None")
        _CTABLE_PRINT_OPTIONS["display_index"] = display_index
    if display_rows is not None:
        if not isinstance(display_rows, int) or isinstance(display_rows, bool) or display_rows < -1:
            raise TypeError("display_rows must be -1 (all), a non-negative int, or None")
        _CTABLE_PRINT_OPTIONS["display_rows"] = display_rows
    if display_width is not _UNSET:
        if not (
            display_width is None
            or (
                isinstance(display_width, int)
                and not isinstance(display_width, bool)
                and display_width >= -1
            )
        ):
            raise TypeError("display_width must be None (auto), -1 (all), or a non-negative int")
        _CTABLE_PRINT_OPTIONS["display_width"] = display_width
    if display_precision is not None:
        if (
            not isinstance(display_precision, int)
            or isinstance(display_precision, bool)
            or display_precision < 0
        ):
            raise TypeError("display_precision must be a non-negative int or None")
        _CTABLE_PRINT_OPTIONS["display_precision"] = display_precision
    if fancy is not None:
        if not isinstance(fancy, bool):
            raise TypeError("fancy must be a bool or None")
        _CTABLE_PRINT_OPTIONS["fancy"] = fancy


def get_printoptions() -> dict[str, Any]:
    """Return a copy of the global :class:`CTable` display options."""
    return dict(_CTABLE_PRINT_OPTIONS)


@contextlib.contextmanager
def printoptions(**kwargs: Any):
    """Temporarily set :class:`CTable` display options, restored on exit.

    Accepts the same keyword arguments as :func:`set_printoptions`.  Handy for a
    one-off full dump, e.g.::

        with blosc2.printoptions(display_rows=-1, display_width=-1):
            print(ctable)
    """
    saved = dict(_CTABLE_PRINT_OPTIONS)
    try:
        set_printoptions(**kwargs)
        yield
    finally:
        _CTABLE_PRINT_OPTIONS.clear()
        _CTABLE_PRINT_OPTIONS.update(saved)


@contextlib.contextmanager
def null_policy(policy: NullPolicy):
    """Temporarily set the default policy for CTable null sentinel inference."""
    token = _NULL_POLICY.set(policy)
    try:
        yield
    finally:
        _NULL_POLICY.reset(token)


# ---------------------------------------------------------------------------
# Index proxy
# ---------------------------------------------------------------------------


_DTYPE_SPEC_FACTORIES = {
    np.dtype(np.int8): int8,
    np.dtype(np.int16): int16,
    np.dtype(np.int32): int32,
    np.dtype(np.int64): int64,
    np.dtype(np.uint8): uint8,
    np.dtype(np.uint16): uint16,
    np.dtype(np.uint32): uint32,
    np.dtype(np.uint64): uint64,
    np.dtype(np.float32): float32,
    np.dtype(np.float64): float64,
    np.dtype(np.complex64): complex64,
    np.dtype(np.complex128): complex128,
    np.dtype(np.bool_): b2_bool,
}


class _CTableInfoReporter(InfoReporter):
    """Info reporter that also preserves the historic ``t.info()`` call style."""

    def __len__(self) -> int:
        return len(self.obj.info_items)

    def __repr__(self) -> str:
        items = self.obj.info_items
        max_key_len = max(len(k) for k, _ in items)
        parts = []
        for key, value in items:
            if isinstance(value, dict):
                parts.append(f"{key.ljust(max_key_len)} :")
                pretty = pprint.pformat(value, sort_dicts=False)
                parts.extend(f" {line}" for line in pretty.splitlines())
                continue

            wrapper = TextWrapper(
                width=96,
                initial_indent=key.ljust(max_key_len) + " : ",
                subsequent_indent=" " * max_key_len + " : ",
            )
            parts.append(wrapper.fill(str(value)))
        return "\n".join(parts) + "\n"

    def __call__(self) -> None:
        print(repr(self), end="")


class _InfoLiteral:
    """Pretty-printer helper for unquoted literal values inside info dicts."""

    def __init__(self, text: str) -> None:
        self.text = text

    def __repr__(self) -> str:
        return self.text


# RowT is intentionally left unbound so CTable works with both dataclasses
# and legacy Pydantic models during the transition period.
RowT = TypeVar("RowT")

# Arrays larger than this threshold use blosc2.arange instead of np.arange to
# avoid large transient allocations when mapping logical to physical row positions.
_BLOSC2_ARANGE_THRESHOLD = 1_000_000


def _arange(start, stop=None, step=1) -> blosc2.NDArray | np.ndarray:
    """Return a range array, using blosc2 for large n to save memory."""
    if stop is None:
        start, stop = 0, start
    n = len(range(start, stop, step))
    return (
        blosc2.arange(start, stop, step) if n >= _BLOSC2_ARANGE_THRESHOLD else np.arange(start, stop, step)
    )


# ---------------------------------------------------------------------------
# Legacy Pydantic-compat helpers
# Keep these so existing code that uses Annotated[type, NumpyDtype(...)] or
# Annotated[str, MaxLen(...)] on a pydantic.BaseModel continues to work.
# ---------------------------------------------------------------------------


class NumpyDtype:
    """Metadata tag for Pydantic-based schemas (legacy)."""

    def __init__(self, dtype):
        self.dtype = dtype


class MaxLen:
    """Metadata tag for fixed-width string/bytes columns in Pydantic-based schemas (legacy)."""

    def __init__(self, length: int):
        self.length = int(length)


def _default_display_width(origin) -> int:
    """Return a sensible display column width for a given Python type (legacy)."""
    return {int: 12, float: 15, bool: 6, complex: 25}.get(origin, 20)


def _resolve_field_dtype(field) -> tuple[np.dtype, int]:
    """Return (numpy dtype, display_width) for a Pydantic model field (legacy).

    Extracts dtype from NumpyDtype metadata when present (same class), otherwise
    falls back to a sensible default for each Python primitive type.
    """
    annotation = field.annotation
    origin = getattr(annotation, "__origin__", annotation)

    # str / bytes → look for MaxLen metadata, build fixed-width dtype
    if origin in (str, bytes) or annotation in (str, bytes):
        is_bytes = origin is bytes or annotation is bytes
        max_len = 32
        if hasattr(annotation, "__metadata__"):
            for meta in annotation.__metadata__:
                if isinstance(meta, MaxLen):
                    max_len = meta.length
                    break
        kind = "S" if is_bytes else "U"
        dt = np.dtype(f"{kind}{max_len}")
        display_width = max(10, min(max_len, 50))
        return dt, display_width

    # Check for explicit NumpyDtype metadata (same class as defined here)
    if hasattr(annotation, "__metadata__"):
        for meta in annotation.__metadata__:
            if isinstance(meta, NumpyDtype):
                dt = np.dtype(meta.dtype)
                display_width = _default_display_width(origin)
                return dt, display_width

    # Primitive defaults
    _PRIMITIVE_MAP = {
        int: (np.int64, 12),
        float: (np.float64, 15),
        bool: (np.bool_, 6),
        complex: (np.complex128, 25),
    }
    if origin in _PRIMITIVE_MAP:
        dt_raw, display_width = _PRIMITIVE_MAP[origin]
        return np.dtype(dt_raw), display_width

    return np.dtype(np.object_), 20


class _LegacySpec(SchemaSpec):
    """Internal compatibility spec wrapping a dtype extracted from a Pydantic schema."""

    def __init__(self, dtype: np.dtype):
        self.dtype = np.dtype(dtype)
        self.python_type = object

    def to_pydantic_kwargs(self) -> dict[str, Any]:
        return {}

    def to_metadata_dict(self) -> dict[str, Any]:
        return {"kind": "legacy", "dtype": str(self.dtype)}


def _compile_pydantic_schema(row_cls: type) -> CompiledSchema:
    """Compatibility adapter: build a CompiledSchema from a Pydantic BaseModel subclass."""
    columns: list[CompiledColumn] = []
    for name, pyd_field in row_cls.model_fields.items():
        dtype, display_width = _resolve_field_dtype(pyd_field)
        spec = _LegacySpec(dtype)
        col = CompiledColumn(
            name=name,
            py_type=object,
            spec=spec,
            dtype=dtype,
            default=MISSING,
            config=ColumnConfig(cparams=None, dparams=None, chunks=None, blocks=None),
            display_width=display_width,
        )
        columns.append(col)
    return CompiledSchema(
        row_cls=row_cls,
        columns=columns,
        columns_by_name={col.name: col for col in columns},
    )


# ---------------------------------------------------------------------------
# ColumnViewIndexer
# ---------------------------------------------------------------------------


class ColumnViewIndexer:
    """Returned by :attr:`Column.view`; indexing returns a Column sub-view.

    Use ``t.price.view[2:10]`` to obtain a writable logical sub-view for
    chained operations (``sum()``, ``[:] = values``, …).
    Use ``t.price[2:10]`` to materialise values as a NumPy array.
    """

    def __init__(self, column: Column) -> None:
        self._column = column

    def __getitem__(self, key) -> Column:
        return self._column._view_from_key(key)

    def __repr__(self) -> str:
        return f"<ColumnViewIndexer col={self._column._col_name!r}>"


# ---------------------------------------------------------------------------
# Internal row/indexing helpers (unchanged)
# ---------------------------------------------------------------------------


def _find_physical_index(arr: blosc2.NDArray, logical_key: int) -> int:
    """Translate a logical (valid-row) index into a physical array index.

    Iterates chunk metadata of the boolean *arr* (valid_rows) to locate the
    *logical_key*-th True value without fully decompressing the array.

    Returns
    -------
    int
        Physical position in the underlying storage array.

    Raises
    ------
    IndexError
        If the logical index is out of range or the array is inconsistent.
    """
    count = 0
    chunk_size = arr.chunks[0]

    for info in arr.iterchunks_info():
        actual_size = min(chunk_size, arr.shape[0] - info.nchunk * chunk_size)
        chunk_start = info.nchunk * chunk_size

        if info.special == blosc2.SpecialValue.ZERO:
            continue

        if info.special == blosc2.SpecialValue.VALUE:
            val = np.frombuffer(info.repeated_value, dtype=arr.dtype)[0]
            if not val:
                continue
            if count + actual_size <= logical_key:
                count += actual_size
                continue
            return chunk_start + (logical_key - count)

        chunk_data = arr[chunk_start : chunk_start + actual_size]
        n_true = int(np.count_nonzero(chunk_data))
        if count + n_true <= logical_key:
            count += n_true
            continue

        return chunk_start + int(np.flatnonzero(chunk_data)[logical_key - count])

    raise IndexError("Unexpected error finding physical index.")


def _make_namedtuple_row_type(col_names: tuple[str, ...]):
    base = namedtuple("CTableRow", col_names, rename=True)
    field_name_map = dict(zip(col_names, base._fields, strict=True))

    class CTableRow(base):
        __slots__ = ()
        _field_name_map = field_name_map
        _original_fields = col_names

        def __getitem__(self, key):
            if isinstance(key, str):
                try:
                    return getattr(self, self._field_name_map[key])
                except KeyError:
                    pass
                # Not a top-level field under its literal spelling: `key` may be
                # an escaped/dotted logical path (e.g. "trip.sec" for a "trip"
                # struct column with a "sec" leaf, or "trip\.info" escaping a
                # literal dot in a top-level name), as reported by
                # schema_dict()/col_names. Walk the path through nested dicts.
                parts = split_field_path(key)
                if parts and parts[0] in self._field_name_map:
                    value = getattr(self, self._field_name_map[parts[0]])
                    try:
                        for part in parts[1:]:
                            value = value[part]
                        return value
                    except (KeyError, TypeError, IndexError):
                        pass
                raise KeyError(f"No field named {key!r}. Available: {list(self._original_fields)}")
            return tuple.__getitem__(self, key)

        def as_dict(self) -> dict[str, Any]:
            return {name: self[name] for name in self._original_fields}

    return CTableRow


# ---------------------------------------------------------------------------
# Column
# ---------------------------------------------------------------------------


class RowTransformer:
    """Row-wise transformer for fixed-shape ndarray columns.

    A row transformer sees one table row at a time.  For a source column with
    physical shape ``(nrows, *item_shape)``, axes passed to reductions are axes
    within ``item_shape`` (so they are shifted by one for batch evaluation).
    """

    def __init__(
        self,
        source: str,
        *,
        selection=(),
        op: str | None = None,
        axis=None,
        ord=None,
    ) -> None:
        self.source = source
        self.selection = tuple(selection)
        self.op = op
        self.axis = axis
        self.ord = ord
        self.kind = "row_transformer"
        self.source_columns = [source]

    def __getitem__(self, key):
        if not isinstance(key, tuple):
            key = (key,)
        return RowTransformer(
            self.source,
            selection=(*self.selection, *key),
            op=self.op,
            axis=self.axis,
            ord=self.ord,
        )

    def _with_op(self, op: str, *, axis=None, ord=None):
        return RowTransformer(self.source, selection=self.selection, op=op, axis=axis, ord=ord)

    def sum(self, *, axis=None):
        return self._with_op("sum", axis=axis)

    def mean(self, *, axis=None):
        return self._with_op("mean", axis=axis)

    def min(self, *, axis=None):
        return self._with_op("min", axis=axis)

    def max(self, *, axis=None):
        return self._with_op("max", axis=axis)

    def argmin(self, *, axis=None):
        return self._with_op("argmin", axis=axis)

    def argmax(self, *, axis=None):
        return self._with_op("argmax", axis=axis)

    def norm(self, *, axis=None, ord=None):
        return self._with_op("norm", axis=axis, ord=ord)

    @staticmethod
    def _serialize_selector(selector):
        if isinstance(selector, slice):
            return {"kind": "slice", "start": selector.start, "stop": selector.stop, "step": selector.step}
        if selector is Ellipsis:
            return {"kind": "ellipsis"}
        if selector is None:
            return {"kind": "newaxis"}
        if isinstance(selector, (int, np.integer)):
            return {"kind": "int", "value": int(selector)}
        raise TypeError(f"Unsupported row-transformer selector {selector!r}")

    @staticmethod
    def _deserialize_selector(data):
        kind = data["kind"]
        if kind == "slice":
            return slice(data.get("start"), data.get("stop"), data.get("step"))
        if kind == "ellipsis":
            return Ellipsis
        if kind == "newaxis":
            return None
        if kind == "int":
            return int(data["value"])
        raise ValueError(f"Unsupported row-transformer selector kind {kind!r}")

    @staticmethod
    def _serialize_axis(axis):
        if isinstance(axis, tuple):
            return list(axis)
        return axis

    @staticmethod
    def _deserialize_axis(axis):
        if isinstance(axis, list):
            return tuple(axis)
        return axis

    def to_metadata(self) -> dict:
        meta = {
            "kind": "row_transformer",
            "source": self.source,
            "selection": [self._serialize_selector(s) for s in self.selection],
        }
        if self.op is not None:
            meta["op"] = self.op
            meta["axis"] = self._serialize_axis(self.axis)
            if self.ord is not None:
                meta["ord"] = self.ord
        return meta

    @classmethod
    def from_metadata(cls, meta: dict):
        return cls(
            meta["source"],
            selection=tuple(cls._deserialize_selector(s) for s in meta.get("selection", ())),
            op=meta.get("op"),
            axis=cls._deserialize_axis(meta.get("axis")),
            ord=meta.get("ord"),
        )

    def _row_axis_to_batch_axis(self, ndim: int, *, none_means_all_item: bool = False):
        axis = self.axis
        item_ndim = max(0, ndim - 1)
        if axis is None:
            return tuple(range(1, ndim)) if none_means_all_item and item_ndim else None

        def one(ax):
            ax = int(ax)
            if ax < 0:
                ax += item_ndim
            if not 0 <= ax < item_ndim:
                raise ValueError(f"axis {ax} is out of bounds for row item with {item_ndim} dimensions")
            return ax + 1

        if isinstance(axis, tuple):
            return tuple(one(ax) for ax in axis)
        return one(axis)

    def _apply_selection(self, arr: np.ndarray) -> np.ndarray:
        if not self.selection:
            return arr
        return arr[(slice(None), *self.selection)]

    def evaluate_batch(self, raw_columns: Mapping[str, Any]) -> np.ndarray:
        arr = np.asarray(raw_columns[self.source])
        if arr.ndim == 0:
            arr = arr.reshape((1,))
        arr = self._apply_selection(arr)
        if self.op is None:
            return np.asarray(arr)
        axis = self._row_axis_to_batch_axis(arr.ndim, none_means_all_item=True)
        if self.op == "sum":
            return np.asarray(np.sum(arr, axis=axis))
        if self.op == "mean":
            return np.asarray(np.mean(arr, axis=axis))
        if self.op == "min":
            return np.asarray(np.min(arr, axis=axis))
        if self.op == "max":
            return np.asarray(np.max(arr, axis=axis))
        if self.op == "argmin":
            if self.axis is None:
                return np.asarray(np.argmin(arr.reshape((arr.shape[0], -1)), axis=1), dtype=np.int64)
            return np.asarray(np.argmin(arr, axis=axis), dtype=np.int64)
        if self.op == "argmax":
            if self.axis is None:
                return np.asarray(np.argmax(arr.reshape((arr.shape[0], -1)), axis=1), dtype=np.int64)
            return np.asarray(np.argmax(arr, axis=axis), dtype=np.int64)
        if self.op == "norm":
            if self.axis is None:
                return np.asarray(np.linalg.norm(arr.reshape((arr.shape[0], -1)), ord=self.ord, axis=1))
            return np.asarray(np.linalg.norm(arr, ord=self.ord, axis=axis))
        raise ValueError(f"Unsupported row-transformer op {self.op!r}")

    def evaluate_row(self, row: Mapping[str, Any]):
        arr = np.asarray(row[self.source])
        if self.selection:
            arr = arr[self.selection]
        if self.op is None:
            return arr
        if self.op == "sum":
            return np.sum(arr, axis=self.axis)
        if self.op == "mean":
            return np.mean(arr, axis=self.axis)
        if self.op == "min":
            return np.min(arr, axis=self.axis)
        if self.op == "max":
            return np.max(arr, axis=self.axis)
        if self.op == "argmin":
            return np.asarray(np.argmin(arr, axis=self.axis), dtype=np.int64)
        if self.op == "argmax":
            return np.asarray(np.argmax(arr, axis=self.axis), dtype=np.int64)
        if self.op == "norm":
            return np.linalg.norm(arr, ord=self.ord, axis=self.axis)
        raise ValueError(f"Unsupported row-transformer op {self.op!r}")

    def evaluate_existing(self, table: CTable) -> np.ndarray:
        return self.evaluate_batch({self.source: table[self.source][:]})


class NullableExpr:
    """Lazy result of arithmetic involving nullable columns.

    Arithmetic on nullable int/timestamp columns promotes to float64 with
    NaN marking the null rows (nullable float columns already use NaN), so
    NaN is the single null representation for every derived expression.  This wrapper carries that fact plus the
    owning table, so reductions (``sum``/``mean``/``min``/``max``/``std``)
    skip nulls and dead physical rows exactly like the corresponding
    :class:`Column` reductions — instead of NaN-poisoning the way a plain
    :class:`blosc2.LazyExpr` reduction would.

    Everything else (``compute()``, slicing, use as an operand) behaves like
    the wrapped expression; further arithmetic keeps the wrapper, since NaN
    propagates through it.

    ``null_pred`` is the boolean lazy predicate (over the raw physical
    columns) marking the null rows. It is carried along instead of being
    re-derived as ``isnan(expr)`` because applying further lazy operations
    on top of a ``where()``-carrying expression is unreliable, and because
    it keeps nulls distinct from NaNs the arithmetic itself may produce
    (e.g. ``0/0``): those are values, not nulls, and are not skipped.
    """

    def __init__(self, expr, table, null_pred):
        self._expr = expr
        self._table = table
        self._null_pred = null_pred

    def __getattr__(self, name):
        return getattr(self._expr, name)

    def __getitem__(self, key):
        return self._expr[key]

    def __repr__(self):
        return f"NullableExpr({self._expr!r})"

    def compute(self, **kwargs):
        return self._expr.compute(**kwargs)

    # ---- chaining: arithmetic keeps the NaN-null channel ----

    def _wrap(self, expr, other=None):
        pred = self._null_pred
        if isinstance(other, NullableExpr):
            pred = pred | other._null_pred
        return NullableExpr(expr, self._table, pred)

    @staticmethod
    def _operand(other):
        return other._expr if isinstance(other, NullableExpr) else other

    @staticmethod
    def _defer(other) -> bool:
        # Column operands re-enter through Column.__r<op>__, which applies
        # the column's own sentinel-null rewrite (operating on its raw
        # storage here would leak sentinel values into the result).
        return isinstance(other, Column)

    def __add__(self, other):
        if self._defer(other):
            return NotImplemented
        return self._wrap(self._expr + self._operand(other), other)

    def __radd__(self, other):
        return self._wrap(other + self._expr)

    def __sub__(self, other):
        if self._defer(other):
            return NotImplemented
        return self._wrap(self._expr - self._operand(other), other)

    def __rsub__(self, other):
        return self._wrap(other - self._expr)

    def __mul__(self, other):
        if self._defer(other):
            return NotImplemented
        return self._wrap(self._expr * self._operand(other), other)

    def __rmul__(self, other):
        return self._wrap(other * self._expr)

    def __truediv__(self, other):
        if self._defer(other):
            return NotImplemented
        return self._wrap(self._expr / self._operand(other), other)

    def __rtruediv__(self, other):
        return self._wrap(other / self._expr)

    def __floordiv__(self, other):
        if self._defer(other):
            return NotImplemented
        return self._wrap(self._expr // self._operand(other), other)

    def __rfloordiv__(self, other):
        return self._wrap(other // self._expr)

    def __mod__(self, other):
        if self._defer(other):
            return NotImplemented
        return self._wrap(self._expr % self._operand(other), other)

    def __rmod__(self, other):
        return self._wrap(other % self._expr)

    def __pow__(self, other):
        if self._defer(other):
            return NotImplemented
        # nan ** 0 == 1.0 would silently resurrect a null as a real value.
        pred = self._null_pred
        if isinstance(other, NullableExpr):
            pred = pred | other._null_pred
        return self._wrap(blosc2.where(pred, np.nan, self._expr ** self._operand(other)), other)

    def __rpow__(self, other):
        # 1 ** nan == 1.0, same hazard as __pow__.
        return self._wrap(blosc2.where(self._null_pred, np.nan, other**self._expr))

    def __neg__(self):
        return self._wrap(-self._expr)

    # ---- comparisons: IEEE NaN already gives SQL False, except for != ----

    def __lt__(self, other):
        return NotImplemented if self._defer(other) else self._expr < self._operand(other)

    def __le__(self, other):
        return NotImplemented if self._defer(other) else self._expr <= self._operand(other)

    def __gt__(self, other):
        return NotImplemented if self._defer(other) else self._expr > self._operand(other)

    def __ge__(self, other):
        return NotImplemented if self._defer(other) else self._expr >= self._operand(other)

    def __eq__(self, other):
        return NotImplemented if self._defer(other) else self._expr == self._operand(other)

    def __ne__(self, other):
        # IEEE says nan != x is True; SQL null semantics say a null satisfies
        # no comparison. Guard both sides (Column and NullableExpr operands
        # carry their own null predicates).
        result = (self._expr != Column._unwrap_operand(other)) & ~self._null_pred
        other_pred = None
        if isinstance(other, Column):
            other_pred = other._raw_null_pred()
        elif isinstance(other, NullableExpr):
            other_pred = other._null_pred
        if other_pred is not None:
            result = result & ~other_pred
        return result

    # ---- reductions: skip nulls and dead physical rows ----

    def _reduction_mask(self, where=None):
        """Live, non-null rows in physical coordinates — the same recipe as
        ``Column._lazy_nonnull_mask``, with the carried null predicate in
        place of the per-column sentinel check."""
        t = self._table
        n_rows = t._known_n_rows()
        mask = None if (n_rows is not None and n_rows == len(t._valid_rows)) else t._valid_rows
        if where is not None:
            mask = where if mask is None else mask & where
        nonnull = ~self._null_pred
        return nonnull if mask is None else mask & nonnull

    def sum(self, dtype=None, *, where=None):
        """Sum of live, non-null values; zero when every value is null."""
        return float(self._expr.sum(where=self._reduction_mask(where), dtype=dtype or np.float64))

    def mean(self, *, where=None):
        """Mean of live, non-null values; NaN when every value is null."""
        try:
            return float(self._expr.mean(where=self._reduction_mask(where), dtype=np.float64))
        except ValueError:
            return float("nan")

    def std(self, ddof: int = 0, *, where=None):
        """Standard deviation of live, non-null values; NaN when every value is null."""
        try:
            return float(self._expr.std(where=self._reduction_mask(where), dtype=np.float64, ddof=ddof))
        except ValueError:
            return float("nan")

    def _minmax(self, op: str, where):
        mask = self._reduction_mask(where)
        count = int(mask.where(blosc2.ones(self._expr.shape, dtype=np.int64), 0).sum(dtype=np.int64))
        if count == 0:
            raise ValueError(f"{op}() called on an expression where all values are null.")
        return getattr(self._expr, op)(where=mask)

    def min(self, *, where=None):
        """Minimum of live, non-null values; raises ``ValueError`` if all are null."""
        return self._minmax("min", where)

    def max(self, *, where=None):
        """Maximum of live, non-null values; raises ``ValueError`` if all are null."""
        return self._minmax("max", where)


class Column:
    """Column view for a :class:`CTable`, with vectorized operations and reductions."""

    _REPR_PREVIEW_ITEMS = 8

    def __init__(self, table: CTable, col_name: str, mask=None):
        self._table = table
        self._col_name = col_name
        self._mask = mask

    @property
    def _raw_col(self):
        cc = self._table._computed_cols.get(self._col_name)
        if cc is not None:
            return self._table._build_computed_lazy(cc)
        return self._table._cols[self._col_name]

    @property
    def is_computed(self) -> bool:
        """True if this column is a virtual computed column (read-only)."""
        return self._col_name in self._table._computed_cols

    @property
    def is_generated(self) -> bool:
        """True if this column is a stored generated/materialized column."""
        return self._col_name in self._table._root_table._materialized_cols

    @property
    def is_stale(self) -> bool:
        """True if this generated column needs to be refreshed before use."""
        meta = self._table._root_table._materialized_cols.get(self._col_name)
        return bool(meta and meta.get("stale", False))

    def _ensure_not_stale(self) -> None:
        if self.is_stale:
            raise ValueError(
                f"Generated column {self._col_name!r} is stale because one or more source columns were "
                f"modified. Call refresh_generated_column({self._col_name!r}) before reading it, or use "
                f"t[{self._col_name!r}].read_stale() to explicitly read the last stored stale values."
            )

    def read_stale(self, key=slice(None)):
        """Read stored values even when this generated column is marked stale.

        This is an explicit escape hatch for inspecting the last materialized
        values.  Normal reads raise for stale generated columns so outdated
        values are not used accidentally.
        """
        return self._values_from_key(key, check_stale=False)

    @property
    def is_list(self) -> bool:
        col = self._table._schema.columns_by_name.get(self._col_name)
        return col is not None and isinstance(col.spec, ListSpec)

    @property
    def is_varlen_scalar(self) -> bool:
        """True if this column holds variable-length scalar strings or bytes."""
        col = self._table._schema.columns_by_name.get(self._col_name)
        return col is not None and isinstance(
            col.spec, (VLStringSpec, VLBytesSpec, StructSpec, ObjectSpec, Utf8Spec)
        )

    @property
    def is_utf8(self) -> bool:
        """True if this column stores variable-length UTF-8 strings (offsets + bytes)."""
        col = self._table._schema.columns_by_name.get(self._col_name)
        return col is not None and isinstance(col.spec, Utf8Spec)

    @property
    def is_dictionary(self) -> bool:
        """True if this column is a dictionary-encoded string column."""
        col = self._table._schema.columns_by_name.get(self._col_name)
        return col is not None and isinstance(col.spec, DictionarySpec)

    @property
    def is_ndarray(self) -> bool:
        """True if this column stores fixed-shape N-D array values per row."""
        col = self._table._schema.columns_by_name.get(self._col_name)
        return col is not None and isinstance(col.spec, NDArraySpec)

    @property
    def raw(self):
        """The underlying storage container for this column, without null-value processing.

        Returns the raw :class:`blosc2.NDArray`, :class:`~blosc2.ListArray`,
        :class:`~blosc2.DictionaryColumn`, or scalar varlen array directly.
        Unlike :meth:`__getitem__`, which always materializes NumPy arrays,
        this is the column as a blosc2-native compressed object: usable as a
        lazy-expression operand without decompressing, and exposing storage
        details such as ``schunk``, ``chunks``, ``cparams`` or
        ``iterchunks_info()``.

        This is a physical view of the column: fixed-width containers are
        over-allocated to chunk capacity for appends, so their first axis is
        longer than ``len(column)`` and positions of rows deleted from the
        table still hold their old values.  No validity-mask or null-sentinel
        processing is applied; use the :class:`Column` interface for logical
        reads.

        Raises :exc:`AttributeError` for computed (virtual) columns, which have
        no backing storage.
        """
        if self.is_computed:
            raise AttributeError(
                f"Column {self._col_name!r} is a computed column and has no underlying array"
            )
        return self._raw_col

    @property
    def _valid_rows(self):
        if self._mask is None:
            return self._table._valid_rows

        return (self._table._valid_rows & self._mask).compute()

    def _lazy_valid_rows(self):
        """Return this column's visible-row mask without forcing lazy evaluation."""
        if self._mask is None:
            return self._table._valid_rows
        return self._table._valid_rows & self._mask

    def _resolve_live_positions(self) -> np.ndarray:
        """Physical positions for all live rows, respecting sorted-view order."""
        slp = getattr(self._table, "_cached_live_positions", None)
        if slp is not None and self._table.base is not None:
            return slp
        return np.where(self._valid_rows[:])[0]

    def _has_identity_positions(self) -> bool:
        """True when logical row ``k`` maps to physical row ``k`` for every row.

        Holds for a base table with no column mask, no sorted/filtered view, and
        no deletions.  In that case a logical slice is a physical slice, so it
        can be read straight from the underlying NDArray instead of resolving
        and gathering explicit live positions.  All checks are O(1) (cached
        counts / lengths) — no validity scan is triggered.
        """
        t = self._table
        if self._mask is not None or t.base is not None:
            return False
        if getattr(t, "_cached_live_positions", None) is not None:
            return False
        n = t._known_n_rows()  # cached live-row count, may be None
        return n is not None and n == len(t._valid_rows)

    def __getitem__(self, key: int | slice | list | np.ndarray):
        """Return values for the given logical index.

        - ``int``              → scalar
        - ``slice``            → :class:`numpy.ndarray`
        - ``list / np.ndarray`` → :class:`numpy.ndarray`
        - ``bool np.ndarray``  → :class:`numpy.ndarray`

        For a writable logical sub-view use :attr:`view`.
        """
        return self._values_from_key(key)

    def _values_from_key(self, key, *, check_stale: bool = True):  # noqa: C901
        """Materialise values for a logical index key."""
        if check_stale:
            self._ensure_not_stale()
        if isinstance(key, tuple) and self.is_ndarray:
            if len(key) == 0:
                raise IndexError("empty tuple index is not valid for Column")
            row_key, inner_key = key[0], key[1:]
            values = self._values_from_key(row_key, check_stale=False)
            if not inner_key:
                return values
            if isinstance(row_key, (int, np.integer)) and not isinstance(row_key, (bool, np.bool_)):
                return values[inner_key]
            return values[(slice(None), *inner_key)]

        if isinstance(key, int):
            n_rows = len(self)
            if key < 0:
                key += n_rows
            if not (0 <= key < n_rows):
                raise IndexError(f"index {key} is out of bounds for column with size {n_rows}")
            # For sorted views the cached positions hold rows in sorted (not
            # physical-ascending) order — use the key-th entry directly.
            _slp = getattr(self._table, "_cached_live_positions", None)
            if _slp is not None and self._table.base is not None:
                pos_true = int(_slp[key])
            else:
                pos_true = _find_physical_index(self._valid_rows, key)
            if self.is_dictionary:
                return self._raw_col[int(pos_true)]
            return self._maybe_decode_timestamp_values(self._raw_col[int(pos_true)])

        elif isinstance(key, slice):
            # Identity fast path: when logical positions equal physical ones, a
            # logical slice is a physical slice.  Read it straight from the
            # underlying NDArray, skipping the O(nrows) live-position scan and
            # letting NDArray's strided-gather fast path handle coarse steps.
            # Plain stored columns only; everything else falls through to the
            # position-gather path below.  utf8 is a varlen-scalar kind but
            # Utf8Array slices itself efficiently (offsets+bytes span read),
            # so it takes the fast path too instead of the index-gather one.
            if (
                not (
                    self.is_computed
                    or self.is_list
                    or self.is_dictionary
                    or (self.is_varlen_scalar and not self.is_utf8)
                )
                and self._has_identity_positions()
            ):
                return self._maybe_decode_timestamp_values(np.asarray(self._raw_col[key]))
            real_pos = self._resolve_live_positions()
            # Apply the slice straight to the physical positions so that all
            # slice semantics (including negative steps) follow NumPy.
            selected_pos = real_pos[key]
            if selected_pos.size == 0:
                if self.is_utf8:
                    return self._raw_col[selected_pos]
                if self.is_list or self.is_varlen_scalar or self.is_dictionary:
                    return []
                if self.is_ndarray:
                    spec = self._table._schema.columns_by_name[self._col_name].spec
                    return np.empty((0, *spec.item_shape), dtype=self.dtype)
                return np.array([], dtype=self.dtype)
            if self.is_computed:
                lo, hi = int(selected_pos.min()), int(selected_pos.max())
                chunk = np.asarray(self._raw_col[lo : hi + 1])
                return chunk[selected_pos - lo]
            if self.is_list or self.is_varlen_scalar or self.is_dictionary:
                return self._raw_col[selected_pos]
            return self._maybe_decode_timestamp_values(np.asarray(self._raw_col[selected_pos]))

        elif isinstance(key, np.ndarray) and key.dtype == np.bool_:
            n_live = len(self)
            if len(key) != n_live:
                raise IndexError(
                    f"Boolean mask length {len(key)} does not match number of live rows {n_live}."
                )
            all_pos = self._resolve_live_positions()
            phys_indices = all_pos[key]
            if self.is_computed:
                raw_np = np.asarray(self._raw_col[:])
                return raw_np[phys_indices]
            if self.is_list or self.is_varlen_scalar or self.is_dictionary:
                return self._raw_col[phys_indices]
            return self._maybe_decode_timestamp_values(self._raw_col[phys_indices])

        elif isinstance(key, (list, tuple, np.ndarray)):
            real_pos = self._resolve_live_positions()
            phys_indices = np.array([real_pos[i] for i in key], dtype=np.int64)
            if self.is_computed:
                raw_np = np.asarray(self._raw_col[:])
                return raw_np[phys_indices]
            if self.is_list or self.is_varlen_scalar or self.is_dictionary:
                return self._raw_col[phys_indices]
            return self._maybe_decode_timestamp_values(self._raw_col[phys_indices])

        raise TypeError(f"Invalid index type: {type(key)}")

    def _view_from_key(self, key) -> Column:
        """Build a Column sub-view for the given logical index key.

        Called by :class:`ColumnViewIndexer`.  Supports slice, boolean mask,
        and integer list / array keys.  The returned :class:`Column` shares
        the underlying physical storage and writes through to the table.
        """
        if isinstance(key, slice):
            valid = self._valid_rows
            real_pos = blosc2.where(valid, _arange(len(valid))).compute()
            start, stop, step = key.indices(len(real_pos))
            mask = blosc2.zeros(len(self._table._valid_rows), dtype=np.bool_)
            if start < stop:
                if step == 1:
                    phys_start = real_pos[start]
                    phys_stop = real_pos[stop - 1]
                    mask[phys_start : phys_stop + 1] = True
                else:
                    lindices = _arange(start, stop, step)
                    phys_indices = real_pos[lindices]
                    mask[phys_indices[:]] = True
            return Column(self._table, self._col_name, mask=mask)

        elif isinstance(key, np.ndarray) and key.dtype == np.bool_:
            n_live = len(self)
            if len(key) != n_live:
                raise IndexError(
                    f"Boolean mask length {len(key)} does not match number of live rows {n_live}."
                )
            all_pos = np.where(self._valid_rows[:])[0]
            phys_indices = all_pos[key]
            mask_np = np.zeros(len(self._table._valid_rows), dtype=np.bool_)
            mask_np[phys_indices] = True
            return Column(self._table, self._col_name, mask=blosc2.asarray(mask_np))

        elif isinstance(key, (list, tuple, np.ndarray)):
            real_pos = blosc2.where(self._valid_rows, _arange(len(self._valid_rows))).compute()
            phys_indices = np.array([real_pos[i] for i in key], dtype=np.int64)
            mask_np = np.zeros(len(self._table._valid_rows), dtype=np.bool_)
            mask_np[phys_indices] = True
            return Column(self._table, self._col_name, mask=blosc2.asarray(mask_np))

        raise TypeError(
            f"Column.view[] does not support key type {type(key).__name__!r}. "
            "Supported: slice, boolean array, list / integer array."
        )

    @property
    def view(self) -> ColumnViewIndexer:
        """Return a :class:`ColumnViewIndexer` for creating logical sub-views.

        Examples
        --------
        Read a sub-view for chained aggregates::

            sub = t.price.view[2:10]
            sub.sum()

        Bulk write through a sub-view::

            t.price.view[0:5][:] = np.zeros(5)
        """
        return ColumnViewIndexer(self)

    def take(self, indices, /) -> Column:
        """Return a column containing values at the requested logical positions.

        Indices are relative to the live values visible through this column
        (including any column view mask).  The result preserves the order of
        ``indices`` and any duplicates.
        """
        if self.is_computed:
            raise ValueError("Column.take is not supported for computed columns yet.")
        table_view = self._table.view(self._valid_rows).select([self._col_name])
        return table_view.take(indices)[self._col_name]

    def __setitem__(self, key: int | slice | list | np.ndarray, value):  # noqa: C901
        """Set one or more live column values; accepts the same index forms as :meth:`__getitem__`.

        Raises ``ValueError`` if the table is read-only or is a view; use
        :meth:`CTable.take` or :meth:`CTable.copy` to get an independent,
        writable table first.
        """
        if self._table._read_only:
            raise ValueError("Table is read-only (opened with mode='r').")
        if self._table.base is not None:
            raise ValueError(
                "Cannot assign values through a view. Use .take(indices) or .copy() "
                "to get an independent, writable table first."
            )
        if self.is_computed:
            raise ValueError(f"Column {self._col_name!r} is a computed column and cannot be written to.")
        # In-place column mutation invalidates any incremental summary accumulator
        # for this column (the close-time builder falls back to a full rescan).
        self._table._root_table._invalidate_summary_accumulator(self._col_name)
        if isinstance(key, int):
            n_rows = len(self)
            if key < 0:
                key += n_rows
            if not (0 <= key < n_rows):
                raise IndexError(f"index {key} is out of bounds for column with size {n_rows}")
            pos_true = _find_physical_index(self._valid_rows, key)
            if self.is_ndarray:
                spec = self._table._schema.columns_by_name[self._col_name].spec
                value = CTable._coerce_ndarray_value(self._col_name, spec, value)
            self._raw_col[int(pos_true)] = value

        elif isinstance(key, np.ndarray) and key.dtype == np.bool_:
            n_live = len(self)
            if len(key) != n_live:
                raise IndexError(
                    f"Boolean mask length {len(key)} does not match number of live rows {n_live}."
                )
            all_pos = np.where(self._valid_rows[:])[0]
            phys_indices = all_pos[key]
            if self.is_list or self.is_varlen_scalar:
                if len(value) != len(phys_indices):
                    raise ValueError("Length mismatch in list-column assignment")
                for pos, cell in zip(phys_indices, value, strict=True):
                    self._raw_col[int(pos)] = cell
            else:
                if self.is_ndarray:
                    spec = self._table._schema.columns_by_name[self._col_name].spec
                    value = CTable._coerce_ndarray_batch(self._col_name, spec, value, len(phys_indices))
                elif isinstance(value, (list, tuple)):
                    value = np.array(value, dtype=self._raw_col.dtype)
                self._raw_col[phys_indices] = value

        elif isinstance(key, (slice, list, tuple, np.ndarray)):
            # Fast path: slice of a blosc2.NDArray into a scalar or ndarray column when
            # there are no deleted rows (physical positions == logical positions).
            # Skips the O(n) validity-mask gather and decompresses one chunk at a time.
            # Guards: must not be a view (views have _last_pos=_n_rows=None, so
            # None==None would be True and bypass the physical-position remapping);
            # _resolve_last_pos() handles disk-opened tables where _last_pos starts None.
            _tbl = self._table
            if (
                isinstance(key, slice)
                and isinstance(value, blosc2.NDArray)
                and not self.is_list
                and not self.is_varlen_scalar
                and _tbl.base is None
                and _tbl._resolve_last_pos() == _tbl._n_rows
            ):
                n_live = _tbl._n_rows
                start, stop, step = key.indices(n_live)
                chunk_size = value.chunks[0] if value.chunks else 65536
                if self.is_ndarray:
                    spec = self._table._schema.columns_by_name[self._col_name].spec

                    def _coerce(v, n):
                        return CTable._coerce_ndarray_batch(self._col_name, spec, v, n)
                else:
                    tgt_dtype = self._raw_col.dtype

                    def _coerce(v, n):
                        return np.ascontiguousarray(v, dtype=tgt_dtype)

                if step == 1:
                    n_selected = stop - start
                    for c in range(0, n_selected, chunk_size):
                        c_end = min(c + chunk_size, n_selected)
                        chunk = _coerce(value[c:c_end], c_end - c)
                        self._raw_col[start + c : start + c_end] = chunk
                else:
                    # Non-unit step: build phys_indices but still decompress value in chunks.
                    logi_indices = list(range(start, stop, step))
                    n_selected = len(logi_indices)
                    phys_indices = np.array(logi_indices, dtype=np.int64)
                    for c in range(0, n_selected, chunk_size):
                        c_end = min(c + chunk_size, n_selected)
                        chunk = _coerce(value[c:c_end], c_end - c)
                        self._raw_col[phys_indices[c:c_end]] = chunk
            else:
                real_pos = blosc2.where(self._valid_rows, _arange(len(self._valid_rows))).compute()
                if isinstance(key, slice):
                    lindices = range(*key.indices(len(real_pos)))
                    phys_indices = np.array([real_pos[i] for i in lindices], dtype=np.int64)
                else:
                    phys_indices = np.array([real_pos[i] for i in key], dtype=np.int64)

                if self.is_list or self.is_varlen_scalar:
                    if len(value) != len(phys_indices):
                        raise ValueError("Length mismatch in list-column assignment")
                    for pos, cell in zip(phys_indices, value, strict=True):
                        self._raw_col[int(pos)] = cell
                else:
                    if self.is_ndarray:
                        spec = self._table._schema.columns_by_name[self._col_name].spec
                        value = CTable._coerce_ndarray_batch(self._col_name, spec, value, len(phys_indices))
                    elif isinstance(value, (list, tuple)):
                        value = np.array(value, dtype=self._raw_col.dtype)
                    # Decompress value in chunks to bound peak memory when it is a blosc2.NDArray.
                    if isinstance(value, blosc2.NDArray):
                        chunk_size = value.chunks[0] if value.chunks else 65536
                        n = len(phys_indices)
                        tgt_dtype = self._raw_col.dtype
                        for c in range(0, n, chunk_size):
                            c_end = min(c + chunk_size, n)
                            chunk = np.ascontiguousarray(value[c:c_end], dtype=tgt_dtype)
                            self._raw_col[phys_indices[c:c_end]] = chunk
                    else:
                        self._raw_col[phys_indices] = value

        else:
            raise TypeError(f"Invalid index type: {type(key)}")
        root = self._table._root_table
        root._mark_generated_columns_stale(self._col_name)
        root._mark_all_indexes_stale()

    def __iter__(self):
        """Iterate over live column values in row order, skipping deleted rows.

        For an ordered view (sorted view, position view, or a reordering slice
        such as ``t[::-1]``), rows are yielded in the view's order rather than
        physical-ascending order.  Physical-order iteration below stays chunked.
        """
        self._ensure_not_stale()
        slp = getattr(self._table, "_cached_live_positions", None)
        if slp is not None and self._table.base is not None:
            yield from self._values_from_key(slice(None), check_stale=False)
            return
        if self.is_computed:
            yield from self._iter_chunks_computed(size=None)
            return
        if self.is_list or self.is_varlen_scalar:
            yield from self._raw_col[np.where(self._valid_rows[:])[0]]
            return
        arr = self._valid_rows
        chunk_size = arr.chunks[0]

        for info in arr.iterchunks_info():
            actual_size = min(chunk_size, arr.shape[0] - info.nchunk * chunk_size)
            chunk_start = info.nchunk * chunk_size

            if info.special == blosc2.SpecialValue.ZERO:
                continue

            if info.special == blosc2.SpecialValue.VALUE:
                val = np.frombuffer(info.repeated_value, dtype=arr.dtype)[0]
                if not val:
                    continue
                yield from self._raw_col[chunk_start : chunk_start + actual_size]
                continue

            mask_chunk = arr[chunk_start : chunk_start + actual_size]
            data_chunk = self._raw_col[chunk_start : chunk_start + actual_size]
            yield from data_chunk[mask_chunk]

    @staticmethod
    def _format_array_value(value) -> str:
        arr = np.asarray(value)
        if arr.ndim == 1 and arr.size <= 6:
            return np.array2string(arr, separator=", ", max_line_width=10_000)
        return f"ndarray(shape={arr.shape}, dtype={arr.dtype})"

    def __repr__(self) -> str:
        preview_len = self._REPR_PREVIEW_ITEMS + 1
        if self.is_list:
            label = self._table._dtype_info_label(
                self.dtype, self._table._schema.columns_by_name[self._col_name].spec
            )
            preview_values = [f"<{label}>"] * min(len(self), preview_len)
        else:
            # Honor an ordered view (sorted/position view, reordering slice) so
            # the preview shows the view's leading rows, not physical-first ones.
            slp = getattr(self._table, "_cached_live_positions", None)
            if slp is not None and self._table.base is not None:
                preview_pos = np.asarray(slp[:preview_len])
            else:
                preview_pos = np.where(self._valid_rows[:])[0][:preview_len]
            if self.is_dictionary or self.is_varlen_scalar:
                preview_values = self._raw_col[preview_pos]
            elif len(preview_pos) == 0:
                preview_values = []
            else:
                preview_values = self._maybe_decode_timestamp_values(self._raw_col[preview_pos]).tolist()
        truncated = len(preview_values) > self._REPR_PREVIEW_ITEMS
        if truncated:
            preview_values = preview_values[: self._REPR_PREVIEW_ITEMS]

        if self.is_ndarray and preview_values:
            preview_items = [self._format_array_value(value) for value in preview_values]
            if truncated:
                preview_items.append("...")
            preview = ", ".join(preview_items)
        elif self.dtype is not None and self.dtype.kind in "biufc" and preview_values:
            arr = np.asarray(preview_values, dtype=self.dtype)
            preview = np.array2string(arr, separator=", ", max_line_width=10_000)[1:-1]
            if truncated:
                preview = f"{preview}, ..." if preview else "..."
        else:
            preview_items = []
            for value in preview_values:
                if isinstance(value, np.generic):
                    value = value.item()
                preview_items.append(repr(value))
            if truncated:
                preview_items.append("...")
            preview = ", ".join(preview_items)

        return f"Column({self._col_name!r}, dtype={self.dtype}, len={len(self)}, values=[{preview}])"

    def __len__(self):
        """Return the number of live (non-deleted) values in this column."""
        return blosc2.count_nonzero(self._valid_rows)

    @property
    def shape(self) -> tuple[int, ...]:
        """Logical shape of the live column values."""
        if self.is_ndarray:
            spec = self._table._schema.columns_by_name[self._col_name].spec
            return (len(self), *spec.item_shape)
        return (len(self),)

    def summary(self) -> str:
        """Return and print a compact summary for this column.

        For fixed-shape ndarray columns this includes logical shape, storage, and
        row-norm statistics when numeric.  Scalar columns fall back to ``info``.
        """
        if not self.is_ndarray:
            text = str(self.info)
            print(text)
            return text
        raw = self._raw_col
        rows = len(self)
        capacity = raw.shape[0] if hasattr(raw, "shape") else len(self._table._valid_rows)
        lines = [
            f"ndarray column {self._col_name!r}",
            f"  rows       : {rows:,} live / {capacity:,} capacity",
            f"  item_shape : {self.item_shape}",
            f"  dtype      : {self.dtype}",
            f"  storage    : NDArray shape={getattr(raw, 'shape', None)}, chunks={getattr(raw, 'chunks', None)}, blocks={getattr(raw, 'blocks', None)}",
        ]
        cbytes = getattr(raw, "cbytes", None)
        if cbytes is not None:
            lines.append(f"  cbytes     : {format_nbytes_info(cbytes)}")
        if rows and self.dtype is not None and self.dtype.kind in "biufc":
            flat = np.asarray(self[:]).reshape(rows, -1)
            norms = np.linalg.norm(flat, axis=1)
            lines.append(
                "  row stats  : "
                f"min(norm(axis=1))={norms.min():.6g}, "
                f"mean(norm(axis=1))={norms.mean():.6g}, "
                f"max(norm(axis=1))={norms.max():.6g}"
            )
        text = "\n".join(lines)
        print(text)
        return text

    @property
    def info(self) -> _CTableInfoReporter:
        """Get information about this column.

        The report includes both logical/live-row details and, when available,
        the physical storage details used internally by lazy predicates.

        Examples
        --------
        >>> print(t["score"].info)
        >>> t["score"].info()
        """
        return _CTableInfoReporter(self)

    @property
    def info_items(self) -> list[tuple[str, object]]:
        """Structured summary items used by :attr:`info`."""
        raw = self._raw_col
        table = self._table
        col_meta = table._schema.columns_by_name.get(self._col_name)
        spec = col_meta.spec if col_meta is not None else None
        chunks = getattr(raw, "chunks", None)
        blocks = getattr(raw, "blocks", None)

        if self.is_list:
            backend = "list"
        elif self.is_varlen_scalar:
            backend = "variable-length scalar"
        elif self.is_dictionary:
            backend = "dictionary"
        else:
            backend = "NDArray" if isinstance(raw, blosc2.NDArray) else type(raw).__name__

        # Virtual computed columns are not stored; otherwise report the table's
        # storage kind, mirroring CTable.info's persistent/in-memory wording.
        if self.is_computed:
            storage = "computed"
        elif isinstance(table._storage, FileTableStorage):
            storage = "persistent"
        else:
            storage = "in-memory"

        # Block order mirrors CTable.info: identity, shape/grid, sizes, content,
        # then compression params.
        items: list[tuple[str, object]] = [
            ("type", self.__class__.__name__),
            ("name", self._col_name),
            ("dtype", table._dtype_info_label(self.dtype, spec)),
            ("backend", backend),
            ("storage", storage),
        ]

        items.append(("nrows", len(self)))
        items.append(("shape", self.shape))
        if chunks is not None:
            items.append(("chunks", chunks))
        if blocks is not None:
            items.append(("blocks", blocks))

        nbytes = getattr(raw, "nbytes", None)
        cbytes = getattr(raw, "cbytes", None)
        cratio = getattr(raw, "cratio", None)
        if nbytes is not None:
            items.append(("nbytes", format_nbytes_info(nbytes)))
        if cbytes is not None:
            items.append(("cbytes", format_nbytes_info(cbytes)))
        if cratio is not None:
            items.append(("cratio", f"{cratio:.2f}x"))

        items.append(("nullable", self.null_value is not None or getattr(spec, "nullable", False)))
        if self.is_dictionary:
            items.append(("dictionary_size", len(raw.dictionary)))

        cparams = getattr(raw, "cparams", None)
        dparams = getattr(raw, "dparams", None)
        if cparams is not None:
            items.append(("cparams", cparams))
        if dparams is not None:
            items.append(("dparams", dparams))
        return items

    @property
    def item_shape(self) -> tuple[int, ...]:
        """Per-row item shape; ``()`` for scalar columns."""
        if self.is_ndarray:
            return tuple(self._table._schema.columns_by_name[self._col_name].spec.item_shape)
        return ()

    @property
    def item_ndim(self) -> int:
        """Number of per-row item dimensions."""
        return len(self.item_shape)

    @property
    def item_size(self) -> int:
        """Number of scalar values stored in each row item."""
        return int(np.prod(self.item_shape, dtype=np.int64)) if self.item_shape else 1

    @property
    def ndim(self) -> int:
        """Number of logical dimensions."""
        return 1 + self.item_ndim

    @property
    def size(self) -> int:
        """Number of live scalar values in the logical column array."""
        return len(self) * self.item_size

    @property
    def row_transformer(self) -> RowTransformer:
        """Build row-wise projections/reductions for generated columns."""
        if not self.is_ndarray:
            raise TypeError(f"Column {self._col_name!r} is not a fixed-shape ndarray column.")
        return RowTransformer(self._col_name)

    def _ensure_queryable(self) -> None:
        self._ensure_not_stale()
        if self.is_utf8:
            raise NotImplementedError(
                f"Column {self._col_name!r} is a variable-length utf8 column; "
                "only comparisons (==, !=, <, <=, >, >=) are supported, not arithmetic "
                "or bitwise operations."
            )
        if self.is_varlen_scalar:
            raise NotImplementedError(
                f"Column {self._col_name!r} is a vlstring/vlbytes column; "
                "lazy expressions and vectorized comparisons are not supported yet."
            )
        if self.is_dictionary:
            raise NotImplementedError(
                f"Column {self._col_name!r} is a dictionary column; "
                "use == and isin() for dictionary column comparisons."
            )

    def _raise_ndarray_compare(self) -> None:
        raise TypeError(
            f"Cannot compare ndarray column {self._col_name!r} directly; the result would not be a "
            "1-D row mask. Use an element projection like t.embedding[:, 0] > 0.5 or an "
            "axis-aware reduction like t.embedding.max(axis=1) > 0.5."
        )

    def _ensure_comparable(self) -> None:
        self._ensure_queryable()
        if self.is_ndarray:
            self._raise_ndarray_compare()

    @staticmethod
    def _unwrap_operand(other):
        if isinstance(other, Column):
            other._ensure_queryable()
            return other._raw_col
        if isinstance(other, NullableExpr):
            return other._expr
        return other

    @property
    def _is_nullable_bool(self) -> bool:
        col = self._table._schema.columns_by_name.get(self._col_name)
        return (
            col is not None
            and col.spec.to_metadata_dict().get("kind") == "bool"
            and getattr(col.spec, "null_value", None) is not None
        )

    @property
    def _timestamp_spec(self):
        col = self._table._schema.columns_by_name.get(self._col_name)
        return col.spec if col is not None and isinstance(col.spec, timestamp) else None

    def _maybe_decode_timestamp_values(self, values):
        spec = self._timestamp_spec
        if spec is None:
            return values
        if np.isscalar(values):
            return np.datetime64(int(values), spec.unit)
        return np.asarray(values).astype(f"datetime64[{spec.unit}]")

    def _coerce_timestamp_operand(self, other):
        spec = self._timestamp_spec
        if isinstance(other, Column) and other.is_ndarray:
            other._raise_ndarray_compare()
        other = self._unwrap_operand(other)
        if spec is None:
            return other
        if isinstance(other, np.datetime64):
            return other.astype(f"datetime64[{spec.unit}]").astype(np.int64)
        if isinstance(other, str):
            return np.datetime64(other).astype(f"datetime64[{spec.unit}]").astype(np.int64)
        if hasattr(other, "isoformat"):
            return np.datetime64(other).astype(f"datetime64[{spec.unit}]").astype(np.int64)
        if isinstance(other, np.ndarray) and np.issubdtype(other.dtype, np.datetime64):
            return other.astype(f"datetime64[{spec.unit}]").astype(np.int64)
        return other

    def _raw_null_pred(self):
        """Boolean lazy predicate over the raw physical array, True where the
        value is this column's null sentinel.

        Returns ``None`` when there is nothing to propagate: no ``null_value``
        configured, or a fixed-shape ndarray column (whose per-item sentinel
        mask does not align 1:1 with the row-level predicates built here;
        see ``Column.is_null()`` for those instead). Dictionary and
        variable-length scalar columns never reach here because
        ``_ensure_queryable`` already rejects them for arithmetic/comparisons.
        """
        if self.is_ndarray:
            return None
        nv = self.null_value
        if nv is None:
            return None
        if isinstance(nv, (float, np.floating)) and np.isnan(nv):
            return blosc2.isnan(self._raw_col)
        return self._raw_col == nv

    def _combined_null_pred(self, other):
        """OR of self's and other's raw null predicates; ``None`` if neither
        operand is nullable (the zero-overhead case)."""
        self_pred = self._raw_null_pred()
        if isinstance(other, Column):
            other_pred = other._raw_null_pred()
        elif isinstance(other, NullableExpr):
            other_pred = other._null_pred
        else:
            other_pred = None
        if self_pred is None:
            return other_pred
        if other_pred is None:
            return self_pred
        return self_pred | other_pred

    def _null_aware_arith(self, other, raw_result):
        """Sentinel-based null propagation for arithmetic: rows where any
        nullable operand is null
        become NaN in the result, promoting integer/timestamp results to
        float64 the same way pandas' legacy int-null arithmetic does. The
        result is a :class:`NullableExpr`, so reductions on it skip nulls
        like the corresponding Column reductions. Costs nothing when neither
        operand is nullable — *raw_result* is returned unchanged.
        """
        null_pred = self._combined_null_pred(other)
        if null_pred is None:
            return raw_result
        return NullableExpr(blosc2.where(null_pred, np.nan, raw_result), self._table, null_pred)

    def _null_aware_compare(self, other, raw_result):
        """SQL ``WHERE`` semantics for comparisons: a null
        operand never satisfies any comparison, so null rows are forced to
        False. Costs nothing when neither operand is nullable.
        """
        null_pred = self._combined_null_pred(other)
        if null_pred is None:
            return raw_result
        return raw_result & ~null_pred

    def __neg__(self):
        self._ensure_queryable()
        return -self._raw_col

    def __pos__(self):
        self._ensure_queryable()
        return +self._raw_col

    def __abs__(self):
        self._ensure_queryable()
        return abs(self._raw_col)

    def __add__(self, other):
        self._ensure_queryable()
        return self._null_aware_arith(other, self._raw_col + self._unwrap_operand(other))

    def __radd__(self, other):
        self._ensure_queryable()
        return self._null_aware_arith(other, self._unwrap_operand(other) + self._raw_col)

    def __sub__(self, other):
        self._ensure_queryable()
        return self._null_aware_arith(other, self._raw_col - self._unwrap_operand(other))

    def __rsub__(self, other):
        self._ensure_queryable()
        return self._null_aware_arith(other, self._unwrap_operand(other) - self._raw_col)

    def __mul__(self, other):
        self._ensure_queryable()
        return self._null_aware_arith(other, self._raw_col * self._unwrap_operand(other))

    def __rmul__(self, other):
        self._ensure_queryable()
        return self._null_aware_arith(other, self._unwrap_operand(other) * self._raw_col)

    def __truediv__(self, other):
        self._ensure_queryable()
        return self._null_aware_arith(other, self._raw_col / self._unwrap_operand(other))

    def __rtruediv__(self, other):
        self._ensure_queryable()
        return self._null_aware_arith(other, self._unwrap_operand(other) / self._raw_col)

    def __floordiv__(self, other):
        self._ensure_queryable()
        return self._null_aware_arith(other, self._raw_col // self._unwrap_operand(other))

    def __rfloordiv__(self, other):
        self._ensure_queryable()
        return self._null_aware_arith(other, self._unwrap_operand(other) // self._raw_col)

    def __mod__(self, other):
        self._ensure_queryable()
        return self._null_aware_arith(other, self._raw_col % self._unwrap_operand(other))

    def __rmod__(self, other):
        self._ensure_queryable()
        return self._null_aware_arith(other, self._unwrap_operand(other) % self._raw_col)

    def __pow__(self, other):
        self._ensure_queryable()
        return self._null_aware_arith(other, self._raw_col ** self._unwrap_operand(other))

    def __rpow__(self, other):
        self._ensure_queryable()
        return self._null_aware_arith(other, self._unwrap_operand(other) ** self._raw_col)

    def __and__(self, other):
        self._ensure_queryable()
        return self._raw_col & self._unwrap_operand(other)

    def __rand__(self, other):
        self._ensure_queryable()
        return self._unwrap_operand(other) & self._raw_col

    def __or__(self, other):
        self._ensure_queryable()
        return self._raw_col | self._unwrap_operand(other)

    def __ror__(self, other):
        self._ensure_queryable()
        return self._unwrap_operand(other) | self._raw_col

    def __xor__(self, other):
        self._ensure_queryable()
        return self._raw_col ^ self._unwrap_operand(other)

    def __rxor__(self, other):
        self._ensure_queryable()
        return self._unwrap_operand(other) ^ self._raw_col

    def __invert__(self):
        self._ensure_queryable()
        if self._is_nullable_bool:
            return self._raw_col == 0
        return ~self._raw_col

    def __lt__(self, other):
        if self.is_utf8:
            return self._utf8_compare(np.less, other)
        self._ensure_comparable()
        return self._null_aware_compare(other, self._raw_col < self._coerce_timestamp_operand(other))

    def __le__(self, other):
        if self.is_utf8:
            return self._utf8_compare(np.less_equal, other)
        self._ensure_comparable()
        return self._null_aware_compare(other, self._raw_col <= self._coerce_timestamp_operand(other))

    def __eq__(self, other):
        if self.is_dictionary:
            return self._dictionary_eq(other)
        if self.is_utf8:
            return self._utf8_compare(np.equal, other)
        self._ensure_comparable()
        if self._is_nullable_bool and isinstance(other, (bool, np.bool_)):
            return self._null_aware_compare(other, self._raw_col == int(other))
        return self._null_aware_compare(other, self._raw_col == self._coerce_timestamp_operand(other))

    def __ne__(self, other):
        if self.is_dictionary:
            result = self._dictionary_eq(other)
            if isinstance(result, np.ndarray):
                return ~result
            return ~np.asarray(result, dtype=bool)
        if self.is_utf8:
            return self._utf8_compare(np.not_equal, other)
        self._ensure_comparable()
        if self._is_nullable_bool and isinstance(other, (bool, np.bool_)):
            return self._null_aware_compare(other, self._raw_col == int(not other))
        return self._null_aware_compare(other, self._raw_col != self._coerce_timestamp_operand(other))

    def _utf8_chunked_bool(self, fn, *, chunk_size: int = 65536) -> np.ndarray:
        """Apply ``fn(chunk, start, stop)`` over this utf8 column's logical rows.

        *fn* returns a boolean array for each ``StringDType`` chunk read from
        the underlying :class:`~blosc2.utf8_array.Utf8Array`.  Returns a
        physical-length (``_valid_rows``-length) boolean NumPy array; rows
        beyond the column's logical length are left ``False``.
        """
        arr = self._raw_col
        n_phys = len(self._table._valid_rows)
        n_logical = len(arr)
        result = np.zeros(n_phys, dtype=np.bool_)
        for start in range(0, n_logical, chunk_size):
            stop = min(start + chunk_size, n_logical)
            result[start:stop] = fn(arr[start:stop], start, stop)
        return result

    def _utf8_chunked_bytes(self, fn, *, chunk_size: int = 65536) -> np.ndarray:
        """Apply ``fn(arr, start, stop)`` over this utf8 column's logical rows.

        Like :meth:`_utf8_chunked_bool`, but *fn* operates directly on the
        underlying :class:`~blosc2.utf8_array.Utf8Array` (raw offsets/bytes)
        instead of a materialized ``StringDType`` chunk, so no per-row decode
        happens.  Returns a physical-length boolean NumPy array; rows beyond
        the column's logical length are left ``False``.
        """
        arr = self._raw_col
        n_phys = len(self._table._valid_rows)
        n_logical = len(arr)
        result = np.zeros(n_phys, dtype=np.bool_)
        for start in range(0, n_logical, chunk_size):
            stop = min(start + chunk_size, n_logical)
            result[start:stop] = fn(arr, start, stop)
        return result

    def _utf8_compare(self, numpy_op, other):
        """Comparison predicate for a utf8 column.

        Compares against a Python ``str`` scalar or another utf8
        :class:`Column` (element-wise).  Returns a physical-length boolean
        ``blosc2.NDArray``, already intersected with this column's live-row
        mask.  A null value on either side never satisfies any comparison
        (SQL ``WHERE`` semantics), matching :meth:`_null_aware_compare` for
        every other column kind.
        """
        if isinstance(other, Column):
            if not other.is_utf8:
                raise TypeError(
                    f"Column {self._col_name!r} is a utf8 column; it can only be compared with a "
                    f"str or another utf8 Column, got Column {other._col_name!r}."
                )
            return self._utf8_compare_column(numpy_op, other)
        if isinstance(other, str):
            return self._utf8_compare_scalar(numpy_op, other)
        raise TypeError(
            f"Column {self._col_name!r} is a utf8 column; it can only be compared with a str "
            f"or another utf8 Column, got {type(other).__name__!r}."
        )

    def _utf8_compare_column(self, numpy_op, other: Column):
        """Column-vs-Column comparison, evaluated chunk by chunk on decoded
        ``StringDType`` values since numexpr/miniexpr cannot operate on them.
        """
        other_arr = other._raw_col
        nv = self.null_value
        other_nv = other.null_value

        def fn(chunk, start, stop):
            rhs = other_arr[start:stop]
            res = numpy_op(chunk, rhs)
            if nv is not None:
                res &= chunk != nv
            if other_nv is not None:
                res &= rhs != other_nv
            return res

        raw = self._utf8_chunked_bool(fn)
        return blosc2.asarray(raw) & self._lazy_valid_rows()

    def _utf8_compare_scalar(self, numpy_op, value: str):
        """Scalar comparison, evaluated chunk by chunk directly on raw UTF-8
        bytes (no decode to ``StringDType``) via
        :meth:`~blosc2.utf8_array.Utf8Array.equal_mask_span` /
        :meth:`~blosc2.utf8_array.Utf8Array.order_masks_span`.
        """
        nv = self.null_value

        if numpy_op in (np.equal, np.not_equal):

            def fn(arr, start, stop):
                res = arr.equal_mask_span(value, start, stop)
                if numpy_op is np.not_equal:
                    res = ~res
                if nv is not None:
                    res &= ~arr.equal_mask_span(nv, start, stop)
                return res
        else:

            def fn(arr, start, stop):
                lt, gt = arr.order_masks_span(value, start, stop)
                if numpy_op is np.less:
                    res = lt
                elif numpy_op is np.less_equal:
                    res = ~gt
                elif numpy_op is np.greater:
                    res = gt
                else:  # np.greater_equal
                    res = ~lt
                if nv is not None:
                    res = res & ~arr.equal_mask_span(nv, start, stop)
                return res

        raw = self._utf8_chunked_bytes(fn)
        return blosc2.asarray(raw) & self._lazy_valid_rows()

    def _dictionary_eq(self, other):
        """Return a physical-slot boolean predicate for dictionary equality.

        Regular fixed-width columns build predicates against their raw physical
        arrays, whose length is the table slot capacity.  Dictionary predicates
        need to use the same coordinate system so they can be combined with
        regular predicates before aggregate/view code intersects them with
        ``_valid_rows``.
        """
        dc = self._raw_col  # DictionaryColumn
        spec = self._table._schema.columns_by_name[self._col_name].spec
        if other is None:
            target_code = spec.null_code
        elif isinstance(other, str):
            try:
                target_code = dc.value_to_code(other)
            except KeyError:
                return blosc2.zeros(len(self._table._valid_rows), dtype=np.bool_)
        else:
            raise TypeError(
                f"Dictionary column {self._col_name!r} can only be compared with str or None, "
                f"got {type(other).__name__!r}."
            )
        pred = dc.codes == np.int32(target_code)
        valid = self._lazy_valid_rows()
        if len(dc.codes) != len(self._table._valid_rows):
            physical = blosc2.zeros(len(self._table._valid_rows), dtype=np.bool_)
            physical[: len(dc.codes)] = pred
            pred = physical
        return pred & valid

    def isin(self, values) -> np.ndarray:
        """Return a boolean array True where the live value is in *values*.

        For dictionary columns this performs efficient integer-code membership
        testing (no decoding of all values).  Values absent from the
        dictionary are treated as not-present.

        For non-dictionary columns this decodes all live values and tests
        membership in a set.
        """
        if self.is_dictionary:
            return self._dictionary_isin(values)
        live_values = self[:]
        test_set = set(values)
        if isinstance(live_values, np.ndarray):
            return np.array([v in test_set for v in live_values.tolist()], dtype=bool)
        return np.array([v in test_set for v in live_values], dtype=bool)

    def _dictionary_isin(self, values) -> np.ndarray:
        """Return a boolean array for in-membership tests against a dictionary column."""
        dc = self._raw_col  # DictionaryColumn
        spec = self._table._schema.columns_by_name[self._col_name].spec
        valid = self._valid_rows
        live_pos = np.where(valid[:])[0]
        if len(live_pos) == 0:
            return np.zeros(0, dtype=bool)
        # Map requested values to codes, ignoring absent values.
        target_codes: set[int] = set()
        for v in values:
            if v is None:
                target_codes.add(spec.null_code)
            elif isinstance(v, str):
                with contextlib.suppress(KeyError):
                    target_codes.add(dc.value_to_code(v))
        if not target_codes:
            return np.zeros(len(live_pos), dtype=bool)
        live_codes = dc.codes[live_pos]
        mask = np.zeros(len(live_codes), dtype=bool)
        for code in target_codes:
            mask |= live_codes == np.int32(code)
        return mask

    def __gt__(self, other):
        if self.is_utf8:
            return self._utf8_compare(np.greater, other)
        self._ensure_comparable()
        return self._null_aware_compare(other, self._raw_col > self._coerce_timestamp_operand(other))

    def __ge__(self, other):
        if self.is_utf8:
            return self._utf8_compare(np.greater_equal, other)
        self._ensure_comparable()
        return self._null_aware_compare(other, self._raw_col >= self._coerce_timestamp_operand(other))

    @property
    def dtype(self):
        """NumPy dtype of the underlying storage.

        ``None`` for variable-length columns with no fixed element dtype
        (:func:`~blosc2.vlstring`, :func:`~blosc2.vlbytes`,
        :func:`~blosc2.list`).  :func:`~blosc2.utf8` columns report
        ``numpy.dtypes.StringDType()``, the dtype of their materialized reads.
        """
        return getattr(self._raw_col, "dtype", None)

    def iter_chunks(self, size: int = 65536):
        """Iterate over live column values in chunks of *size* rows.

        Yields numpy arrays of at most *size* elements each, skipping deleted
        rows.  The last chunk may be smaller than *size*.

        Parameters
        ----------
        size:
            Number of live rows per yielded chunk.  Defaults to 65 536.

        Yields
        ------
        numpy.ndarray
            A 1-D array of up to *size* live values with this column's dtype.

        Examples
        --------
        >>> for chunk in t["score"].iter_chunks(size=100_000):
        ...     process(chunk)
        """
        self._ensure_not_stale()
        if self.is_computed:
            yield from self._iter_chunks_computed(size=size)
            return
        if self.is_list:
            raise TypeError("Column.iter_chunks() is not supported for list columns in V1.")
        if self.is_varlen_scalar and not self.is_utf8:
            raise TypeError("Column.iter_chunks() is not supported for varlen scalar columns.")
        valid = self._valid_rows
        raw = self._raw_col
        arr_len = len(valid)
        phys_chunk = valid.chunks[0]

        pending: list[np.ndarray] = []
        pending_count = 0

        for info in valid.iterchunks_info():
            actual = min(phys_chunk, arr_len - info.nchunk * phys_chunk)
            start = info.nchunk * phys_chunk

            if info.special == blosc2.SpecialValue.ZERO:
                continue

            if info.special == blosc2.SpecialValue.VALUE:
                val = np.frombuffer(info.repeated_value, dtype=valid.dtype)[0]
                if not val:
                    continue
                segment = raw[start : start + actual]
            else:
                mask = valid[start : start + actual]
                data_part = raw[start : start + actual]
                if len(data_part) < actual:
                    # Logically-sized storage (utf8) is shorter than the
                    # capacity-sized validity mask; rows past its end are
                    # never live, so the extra mask tail is all False.
                    mask = mask[: len(data_part)]
                segment = data_part[mask]

            if len(segment) == 0:
                continue

            pending.append(segment)
            pending_count += len(segment)

            while pending_count >= size:
                combined = np.concatenate(pending)
                yield combined[:size]
                rest = combined[size:]
                pending = [rest] if len(rest) > 0 else []
                pending_count = len(rest)

        if pending:
            yield np.concatenate(pending)

    def _iter_chunks_computed(self, size):
        """Yield live values from a computed column, chunk-by-chunk.

        Evaluates the LazyExpr slice-by-slice using the physical chunk layout
        of *valid_rows* and applies the valid-rows mask before accumulating.
        When *size* is None (used by ``__iter__``), each physical chunk is
        yielded directly.
        """
        lazy = self._raw_col  # a LazyExpr
        valid = self._valid_rows
        phys_len = len(valid)
        chunk_size = valid.chunks[0]

        pending: list[np.ndarray] = []
        pending_n = 0

        for chunk_start in range(0, phys_len, chunk_size):
            chunk_end = min(chunk_start + chunk_size, phys_len)
            mask = valid[chunk_start:chunk_end]  # numpy bool array
            n_live = int(np.count_nonzero(mask))
            if n_live == 0:
                continue

            # Evaluate the expression only for this physical slice
            data_chunk = np.asarray(lazy[chunk_start:chunk_end])
            segment = data_chunk[mask] if n_live < (chunk_end - chunk_start) else data_chunk

            if size is None:
                # __iter__ path: yield each chunk directly
                yield from segment
                continue

            pending.append(segment)
            pending_n += len(segment)

            while pending_n >= size:
                combined = np.concatenate(pending)
                yield combined[:size]
                rest = combined[size:]
                pending = [rest] if len(rest) > 0 else []
                pending_n = len(rest)

        if size is not None and pending:
            yield np.concatenate(pending)

    def assign(self, data) -> None:
        """Replace all live values in this column with *data*.

        Parameters
        ----------
        data:
            List, numpy array, or any iterable.  Must have exactly as many
            elements as there are live rows in this column.  Values are
            coerced to the column's dtype if possible.

        Raises
        ------
        ValueError
            If ``len(data)`` does not match the number of live rows, the
            table is opened read-only, or the table is a view (use
            :meth:`CTable.take` or :meth:`CTable.copy` to get an
            independent, writable table first).
        TypeError
            If values cannot be coerced to the column's dtype.
        """
        if self._table._read_only:
            raise ValueError("Table is read-only (opened with mode='r').")
        if self._table.base is not None:
            raise ValueError(
                "Cannot assign values through a view. Use .take(indices) or .copy() "
                "to get an independent, writable table first."
            )
        if self.is_computed:
            raise ValueError(f"Column {self._col_name!r} is a computed column and cannot be written to.")
        if self.is_list:
            values = list(data)
            if len(values) != len(self):
                raise ValueError(f"assign() requires {len(self)} values (live rows), got {len(values)}.")
            live_pos = np.where(self._valid_rows[:])[0]
            for pos, cell in zip(live_pos, values, strict=True):
                self._raw_col[int(pos)] = cell
            root = self._table._root_table
            root._mark_generated_columns_stale(self._col_name)
            root._mark_all_indexes_stale()
            return
        n_live = len(self)
        arr = np.asarray(data)
        if len(arr) != n_live:
            raise ValueError(f"assign() requires {n_live} values (live rows), got {len(arr)}.")
        try:
            arr = arr.astype(self.dtype)
        except (ValueError, OverflowError) as exc:
            raise TypeError(f"Cannot coerce data to column dtype {self.dtype!r}: {exc}") from exc
        live_pos = np.where(self._valid_rows[:])[0]
        self._raw_col[live_pos] = arr
        root = self._table._root_table
        root._mark_generated_columns_stale(self._col_name)
        root._mark_all_indexes_stale()

    # ------------------------------------------------------------------
    # Null sentinel support
    # ------------------------------------------------------------------

    @property
    def null_value(self):
        """The sentinel value that represents NULL for this column, or ``None``."""
        col_info = self._table._schema.columns_by_name.get(self._col_name)
        if col_info is None:
            return None
        return getattr(col_info.spec, "null_value", None)

    def _null_mask_for(self, arr: np.ndarray) -> np.ndarray:
        """Return a bool array True where *arr* contains the null sentinel.

        Always returns an array of the same length as *arr*; all False when
        no null_value is configured.
        """
        nv = self.null_value
        if nv is None:
            return np.zeros(len(arr), dtype=np.bool_)
        arr = np.asarray(arr)
        if self.is_ndarray:
            if arr.ndim <= self.item_ndim:
                arr = arr.reshape((1, *arr.shape))
            if isinstance(nv, (float, np.floating)) and np.isnan(nv):
                elem_mask = np.isnan(arr)
            else:
                elem_mask = arr == nv
            inner_axes = tuple(range(1, elem_mask.ndim))
            return elem_mask.all(axis=inner_axes) if inner_axes else elem_mask.astype(np.bool_)
        if np.issubdtype(arr.dtype, np.datetime64):
            # Timestamp columns materialize with the int64 sentinel already
            # decoded into np.datetime64('NaT') (they share the same bit
            # pattern), so the sentinel value itself never appears in arr.
            return np.isnat(arr)
        if isinstance(nv, (float, np.floating)) and np.isnan(nv):
            return np.isnan(arr)
        return arr == nv

    def is_null(self) -> np.ndarray:
        """Return a boolean array True where the live value is the null sentinel.

        For varlen scalar columns (vlstring/vlbytes) nullability is represented
        as native ``None`` values, so this returns True wherever the value is
        ``None``.  For dictionary columns, returns True where the code equals
        the null_code (``-1`` by default).
        """
        if self.is_dictionary:
            return self._dictionary_eq(None)
        if self.is_varlen_scalar and not self.is_utf8:
            return np.array([v is None for v in self], dtype=np.bool_)
        return self._null_mask_for(self[:])

    def notnull(self) -> np.ndarray:
        """Return a boolean array True where the live value is *not* the null sentinel."""
        return ~self.is_null()

    def null_count(self) -> int:
        """Return the number of live rows whose value equals the null sentinel.

        Returns ``0`` in O(1) if no ``null_value`` is configured for this column
        and the column is not a varlen scalar column.
        """
        if self.is_dictionary:
            return int(self.is_null().sum())
        if self.is_varlen_scalar and not self.is_utf8:
            return sum(1 for v in self if v is None)
        if self.null_value is None:
            return 0
        return int(self.is_null().sum())

    def fillna(self, value):
        """Return live values with null sentinels replaced by *value*.

        Dictionary and variable-length scalar columns (whose nulls are
        native ``None`` cells) return a list; other columns return a NumPy
        array.
        """
        if (self.is_dictionary or self.is_varlen_scalar) and not self.is_utf8:
            return [value if v is None else v for v in self[:]]
        arr = np.array(self[:], copy=True)
        if self.null_value is not None:
            arr[self._null_mask_for(arr)] = value
        return arr

    def _nonnull_chunks(self):
        """Yield chunks of live, non-null values.

        Each yielded array has the null sentinel values removed.  If no
        null_value is configured this behaves identically to
        :meth:`iter_chunks`.
        """
        nv = self.null_value
        if nv is None:
            yield from self.iter_chunks()
            return
        is_nan_nv = isinstance(nv, float) and np.isnan(nv)
        for chunk in self.iter_chunks():
            if is_nan_nv:
                mask = ~np.isnan(chunk)
            else:
                mask = chunk != nv
            filtered = chunk[mask]
            if len(filtered) > 0:
                yield filtered

    def unique(self) -> np.ndarray:
        """Return sorted array of unique live, non-null values.

        Null sentinel values are excluded.
        Processes data in chunks — never loads the full column at once.
        """
        seen: set = set()
        for chunk in self._nonnull_chunks():
            seen.update(chunk.tolist())
        return np.array(sorted(seen), dtype=self.dtype)

    def value_counts(self) -> dict:
        """Return a ``{value: count}`` dict sorted by count descending.

        Null sentinel values are excluded.
        Processes data in chunks — never loads the full column at once.

        Example
        -------
        >>> t["active"].value_counts()
        {True: 8432, False: 1568}
        """
        counts: dict = {}
        for chunk in self._nonnull_chunks():
            for val in chunk.tolist():
                counts[val] = counts.get(val, 0) + 1
        return dict(sorted(counts.items(), key=lambda kv: -kv[1]))

    # ------------------------------------------------------------------
    # Aggregate helpers
    # ------------------------------------------------------------------

    def _require_nonempty(self, op: str) -> None:
        if len(self) == 0:
            raise ValueError(f"Column.{op}() called on an empty column.")

    def _require_kind(self, kinds: str, op: str) -> None:
        """Raise TypeError if this column's dtype is not in *kinds*."""
        self._ensure_not_stale()
        if self.dtype.kind not in kinds:
            _kind_names = {
                "b": "bool",
                "i": "signed int",
                "u": "unsigned int",
                "f": "float",
                "c": "complex",
                "U": "string",
                "S": "bytes",
            }
            raise TypeError(
                f"Column.{op}() is not supported for dtype {self.dtype!r} "
                f"({_kind_names.get(self.dtype.kind, self.dtype.kind)})."
            )

    # ------------------------------------------------------------------
    # Aggregates
    # ------------------------------------------------------------------

    def _normalize_sum_where(self, where):
        """Normalize an optional ``sum(where=...)`` predicate to a boolean array/expression."""
        if where is None:
            return None
        if isinstance(where, str):
            self._table._guard_varlen_scalar_expression(where)
            operands = self._table._where_expression_operands(where)
            where, operands = self._table._rewrite_nested_expression(where, operands)
            where = blosc2.lazyexpr(where, operands)
        if isinstance(where, np.ndarray) and where.dtype == np.bool_:
            where = blosc2.asarray(where)
        if isinstance(where, Column):
            where = where._raw_col == 1 if where._is_nullable_bool else where._raw_col
        if not (
            isinstance(where, (blosc2.NDArray, blosc2.LazyExpr))
            and getattr(where, "dtype", None) == np.bool_
        ):
            raise TypeError(f"Expected boolean blosc2.NDArray or LazyExpr, got {type(where).__name__}")
        return where

    def _lazy_nonnull_mask(self, where=None):
        """Build a lazy visible-row mask, optionally intersected with non-null values.

        When all physical rows are visible, avoid injecting ``_valid_rows`` into
        the expression.  This keeps aggregate predicates aligned with the data
        columns, which lets the miniexpr reduction fast path run for common
        no-deletes/no-filtered-view cases.
        """
        raw = self._raw_col
        if not isinstance(raw, (blosc2.NDArray, blosc2.LazyExpr)):
            return NotImplemented

        table_n_rows = self._table._known_n_rows()
        all_rows_visible = (
            self._mask is None and table_n_rows is not None and table_n_rows == len(self._table._valid_rows)
        )
        mask = None if all_rows_visible else self._lazy_valid_rows()
        if where is not None:
            mask = where if mask is None else mask & where
        nv = self.null_value
        if nv is not None:
            if isinstance(nv, (float, np.floating)) and np.isnan(nv):
                nonnull = ~blosc2.isnan(raw)
            else:
                nonnull = raw != nv
            mask = nonnull if mask is None else mask & nonnull
        return mask

    def _sum_lazy_fastpath(self, acc_dtype, where=None, *, jit=None, jit_backend=None):
        """Try to compute ``sum`` as a pushed-down lazy masked reduction."""
        if self.is_list or self.is_varlen_scalar or self.dtype is None or self.dtype.kind not in "biufc":
            return NotImplemented

        raw = self._raw_col
        if not isinstance(raw, (blosc2.NDArray, blosc2.LazyExpr)):
            return NotImplemented

        # A lazy masked reduction scans the full physical column.  For very
        # selective filtered views, the existing iterator can skip all-zero mask
        # chunks and is usually faster.  Explicit sum(where=...) is already a
        # direct pushed-down aggregate, so do not apply the density guard there.
        total_rows = len(self._table._valid_rows)
        if (
            where is None
            and self._table.base is not None
            and total_rows
            and self._table.nrows / total_rows < 0.25
        ):
            return NotImplemented

        mask = self._lazy_nonnull_mask(where=where)
        if mask is NotImplemented:
            return NotImplemented

        try:
            if mask is None:
                return raw.sum(dtype=acc_dtype, jit=jit, jit_backend=jit_backend)
            force_miniexpr = jit is True or jit_backend is not None
            if force_miniexpr and isinstance(raw, blosc2.NDArray):
                zero = blosc2.zeros(
                    raw.shape, dtype=np.dtype(acc_dtype), chunks=raw.chunks, blocks=raw.blocks
                )
            else:
                zero = acc_dtype(0)
            return blosc2.where(mask, raw, zero).sum(dtype=acc_dtype, jit=jit, jit_backend=jit_backend)
        except Exception:
            return NotImplemented

    def _ndarray_values_for_reduction(self, where=None) -> np.ndarray:
        arr = np.asarray(self[:])
        null_mask = self._null_mask_for(arr) if self.null_value is not None else None
        if null_mask is not None and null_mask.any():
            arr = arr[~null_mask]
        if where is None:
            return arr
        where = self._normalize_sum_where(where)
        mask = where.compute() if isinstance(where, blosc2.LazyExpr) else where[:]
        mask = np.asarray(mask, dtype=bool)
        if mask.ndim != 1:
            raise ValueError("Column reduction where= must be a 1-D row mask.")
        if len(mask) != len(self._table._valid_rows):
            if len(mask) != len(self):
                raise ValueError(
                    f"Column reduction where= mask length {len(mask)} does not match live rows {len(self)}."
                )
            if null_mask is not None and len(null_mask) == len(mask):
                mask = mask[~null_mask]
            return arr[mask]
        live_pos = np.where(self._valid_rows[:])[0]
        row_mask = mask[live_pos]
        if null_mask is not None and len(null_mask) == len(row_mask):
            row_mask = row_mask[~null_mask]
        return arr[row_mask]

    def _ndarray_reduce(self, op: str, *, axis=None, dtype=None, where=None, ddof: int = 0):
        arr = self._ndarray_values_for_reduction(where=where)
        if op == "sum":
            return np.sum(arr, axis=axis, dtype=dtype)
        if op == "mean":
            return np.mean(arr, axis=axis, dtype=dtype)
        if op == "min":
            return np.min(arr, axis=axis)
        if op == "max":
            return np.max(arr, axis=axis)
        if op == "argmin":
            return np.argmin(arr, axis=axis)
        if op == "argmax":
            return np.argmax(arr, axis=axis)
        if op == "std":
            return np.std(arr, axis=axis, ddof=ddof, dtype=dtype)
        raise ValueError(f"Unsupported ndarray reduction {op!r}")

    def sum(self, dtype=None, axis=None, *, where=None, jit=None, jit_backend=None):
        """Sum of all live, non-null values.

        Returns zero for an empty column or filtered view.

        Supported dtypes: bool, int, uint, float, complex.
        Bool values are counted as 0 / 1.
        Null sentinel values are skipped.

        Parameters
        ----------
        dtype:
            Optional accumulator dtype.  When omitted, float columns use
            ``np.float64``, complex columns use ``np.complex128``, and integer
            / bool columns use ``np.int64``.
        where:
            Optional boolean predicate. Only rows where the predicate is true,
            the table row is live, and this column is non-null are included.
            This enables direct filtered aggregate pushdown, avoiding creation
            of an intermediate filtered table view.
        jit:
            Optional miniexpr JIT policy passed to the lazy reduction engine.
        jit_backend:
            Optional miniexpr JIT backend. Use ``"tcc"`` or ``"cc"``.

        Examples
        --------
        Sum values matching a predicate without materializing a filtered view::

            total = t["amount"].sum(where=t.category == 3)

        Combine several column predicates::

            total = t.col2.sum(where=(t.col1 < 300) & (t.col2 < 400))

        Nullable sentinel values are skipped automatically::

            # Equivalent to summing only live rows where predicate is true and
            # t.col2 is not its configured null sentinel.
            total = t.col2.sum(where=t.col1 < 300)
        """
        if self.is_ndarray:
            self._require_kind("biufc", "sum")
            return self._ndarray_reduce("sum", axis=axis, dtype=dtype, where=where)
        if axis not in (None, 0):
            return np.sum(self[:], axis=axis, dtype=dtype)
        self._require_kind("biufc", "sum")
        where = self._normalize_sum_where(where)
        # Use a wide accumulator to reduce overflow risk
        acc_dtype = np.dtype(dtype).type if dtype is not None else None
        if acc_dtype is None:
            acc_dtype = (
                np.float64
                if self.dtype.kind == "f"
                else (
                    np.complex128
                    if self.dtype.kind == "c"
                    else np.int64
                    if self.dtype.kind in "biu"
                    else None
                )
            )

        result = self._sum_lazy_fastpath(acc_dtype, where=where, jit=jit, jit_backend=jit_backend)
        if result is NotImplemented:
            if where is not None:
                return self._table.where(where)[self._col_name].sum(
                    dtype=dtype, jit=jit, jit_backend=jit_backend
                )
            result = acc_dtype(0)
            for chunk in self._nonnull_chunks():
                result += chunk.sum(dtype=acc_dtype)

        # Return in the column's natural dtype when it fits, else keep the requested/wide dtype
        if dtype is None and self.dtype.kind in "biu":
            return int(result)
        return result

    def _lazy_aggregate_fastpath(self, op: str, *, where=None, dtype=None, ddof: int = 0):
        if self.is_list or self.is_varlen_scalar or self.dtype is None or self.dtype.kind not in "biuf":
            return NotImplemented
        raw = self._raw_col
        if not isinstance(raw, (blosc2.NDArray, blosc2.LazyExpr)):
            return NotImplemented
        mask = self._lazy_nonnull_mask(where=where)
        if mask is NotImplemented:
            return NotImplemented
        try:
            count = None
            if op in {"min", "max"} and mask is not None:
                count = int(mask.where(blosc2.ones(raw.shape, dtype=np.int64), 0).sum(dtype=np.int64))
                if count == 0:
                    raise ValueError(f"{op}() called on a column where all values are null.")
            if op == "mean":
                return float(
                    raw.mean(dtype=dtype or np.float64)
                    if mask is None
                    else raw.mean(where=mask, dtype=dtype or np.float64)
                )
            if op == "std":
                return float(
                    raw.std(dtype=dtype or np.float64, ddof=ddof)
                    if mask is None
                    else raw.std(where=mask, dtype=dtype or np.float64, ddof=ddof)
                )
            if op == "min":
                return raw.min() if mask is None else raw.min(where=mask)
            if op == "max":
                return raw.max() if mask is None else raw.max(where=mask)
        except ValueError:
            if op in {"mean", "std"}:
                return float("nan")
            raise
        except Exception:
            return NotImplemented
        return NotImplemented

    def _summary_minmax_source(self):
        """Return ``(sidecar_path, dtype, nullable)`` for a summary-readable
        ``min``/``max``, or ``None`` when the index shortcut is not provably
        correct.

        Excluded: a view (its summary describes the base table); a column kind
        without numeric/string block extrema; a leaky null sentinel (only a
        non-nullable column, or a NaN-sentinel float — whose NaNs the summary
        drops — match the nulls-skipped contract of ``min()``); and a stale,
        absent, or in-memory-only index.  Deletions/appends are covered by the
        stale flag (every mutation marks the index stale; a rebuild re-summarises
        only the live rows, and capacity padding never enters the summaries).
        """
        table = self._table
        if table.base is not None:
            return None
        if (
            self.is_computed
            or self.is_ndarray
            or self.is_list
            or self.is_varlen_scalar
            or self.is_dictionary
        ):
            return None
        dtype = self.dtype
        if dtype is None or dtype.kind not in "biufUS":
            return None
        col = table._schema.columns_by_name.get(self._col_name)
        spec = col.spec if col is not None else None
        nullable = getattr(spec, "nullable", False)
        null_value = getattr(spec, "null_value", None)
        is_nan_float = dtype.kind == "f" and isinstance(null_value, float) and np.isnan(null_value)
        if nullable and not is_nan_float:
            return None  # non-NaN sentinel leaks into the block extrema
        desc = table._root_table._get_index_catalog().get(self._col_name)
        if not desc or desc.get("stale", False):
            return None
        levels = desc.get("levels") or {}
        level = "block" if "block" in levels else next(iter(levels), None)
        if level is None:
            return None
        path = levels[level].get("path")
        if path is None:
            return None  # in-memory sidecar: the scan is already fast
        return path, dtype, nullable

    def _index_summary_minmax(self, op: str):
        """Exact ``min``/``max`` from the column index's block summaries, or
        ``NotImplemented`` when that shortcut is not applicable (see
        :meth:`_summary_minmax_source`).

        Every index kind (SUMMARY/FULL/PARTIAL/BUCKET/OPSI) persists per-block
        ``(min, max, flags)``, so reducing those is decompression-free (~240x
        faster than scanning tens of millions of rows).
        """
        source = self._summary_minmax_source()
        if source is None:
            return NotImplemented
        path, dtype, nullable = source
        try:
            from blosc2.indexing import _INDEX_MMAP_MODE, FLAG_ALL_NAN, FLAG_HAS_NAN, _open_sidecar_file

            # Read the tiny (min, max, flags) arrays and drop the handle: unlike
            # the cached _open_level_summary_handle, this releases the file
            # descriptor immediately (min()/max() must not leak one per table).
            handle = _open_sidecar_file(path, _INDEX_MMAP_MODE)
            flags = np.asarray(handle["flags"][:])
            vals = np.asarray(handle[op][:])
            del handle
        except Exception:
            return NotImplemented
        if vals.shape[0] == 0:
            return NotImplemented
        # A non-nullable float with NaN *data* makes numpy min/max return NaN,
        # but the summary dropped those NaNs — they would disagree, so bail.
        if dtype.kind == "f" and not nullable and bool((flags & (FLAG_HAS_NAN | FLAG_ALL_NAN)).any()):
            return NotImplemented
        valid = (flags & FLAG_ALL_NAN) == 0
        if not valid.any():
            return NotImplemented  # whole column null → let the scan raise
        vals = vals[valid]
        if dtype.kind in "US":
            return min(vals) if op == "min" else max(vals)
        return vals.min() if op == "min" else vals.max()

    def min(self, axis=None, *, where=None):
        """Minimum live, non-null value.

        Supported dtypes: bool, int, uint, float, string, bytes.
        Strings are compared lexicographically.
        Null sentinel values are skipped. When *where* is provided, only rows
        matching the boolean predicate are included.
        """
        if self.is_ndarray:
            self._require_kind("biuf", "min")
            return self._ndarray_reduce("min", axis=axis, where=where)
        if axis not in (None, 0):
            return np.min(self[:], axis=axis)
        self._require_kind("biufUS", "min")
        where = self._normalize_sum_where(where)
        if where is None:
            # Try the index-summary shortcut first: it returns NotImplemented for
            # an empty/all-null column, so the emptiness check (which counts live
            # rows — expensive on a nullable column) only runs on the fallback.
            fast_idx = self._index_summary_minmax("min")
            if fast_idx is not NotImplemented:
                return fast_idx
            self._require_nonempty("min")
        fast = self._lazy_aggregate_fastpath("min", where=where)
        if fast is not NotImplemented:
            return fast
        if where is not None:
            return self._table.where(where)[self._col_name].min()
        result = None
        is_str = self.dtype.kind in "US"
        for chunk in self._nonnull_chunks():
            # numpy .min()/.max() don't support string dtypes in recent NumPy;
            # fall back to Python's built-in min/max which work on any comparable type.
            chunk_min = min(chunk) if is_str else chunk.min()
            if result is None or chunk_min < result:
                result = chunk_min
        if result is None:
            raise ValueError("min() called on a column where all values are null.")
        return result

    def max(self, axis=None, *, where=None):
        """Maximum live, non-null value.

        Supported dtypes: bool, int, uint, float, string, bytes.
        Strings are compared lexicographically.
        Null sentinel values are skipped. When *where* is provided, only rows
        matching the boolean predicate are included.
        """
        if self.is_ndarray:
            self._require_kind("biuf", "max")
            return self._ndarray_reduce("max", axis=axis, where=where)
        if axis not in (None, 0):
            return np.max(self[:], axis=axis)
        self._require_kind("biufUS", "max")
        where = self._normalize_sum_where(where)
        if where is None:
            # See min(): shortcut before the live-row count.
            fast_idx = self._index_summary_minmax("max")
            if fast_idx is not NotImplemented:
                return fast_idx
            self._require_nonempty("max")
        fast = self._lazy_aggregate_fastpath("max", where=where)
        if fast is not NotImplemented:
            return fast
        if where is not None:
            return self._table.where(where)[self._col_name].max()
        result = None
        is_str = self.dtype.kind in "US"
        for chunk in self._nonnull_chunks():
            chunk_max = max(chunk) if is_str else chunk.max()
            if result is None or chunk_max > result:
                result = chunk_max
        if result is None:
            raise ValueError("max() called on a column where all values are null.")
        return result

    def argmin(self, axis=None, *, where=None):
        """Index of the minimum live, non-null value.

        For fixed-shape ndarray columns, this follows NumPy axis semantics on
        the logical array of shape ``(nrows, *item_shape)``.  For scalar
        columns, the result is the logical row position within this column (or
        filtered view).
        """
        if self.is_ndarray:
            self._require_kind("biuf", "argmin")
            return self._ndarray_reduce("argmin", axis=axis, where=where)
        if axis not in (None, 0):
            return np.argmin(self[:], axis=axis)
        self._require_kind("biuf", "argmin")
        if where is not None:
            return self._table.where(self._normalize_sum_where(where))[self._col_name].argmin()
        arr = np.asarray(self[:])
        if arr.size == 0:
            raise ValueError("argmin() called on an empty column.")
        mask = (
            self._null_mask_for(arr) if self.null_value is not None else np.zeros(len(arr), dtype=np.bool_)
        )
        if mask.all():
            raise ValueError("argmin() called on a column where all values are null.")
        positions = np.where(~mask)[0]
        return int(positions[np.argmin(arr[positions])])

    def argmax(self, axis=None, *, where=None):
        """Index of the maximum live, non-null value.

        For fixed-shape ndarray columns, this follows NumPy axis semantics on
        the logical array of shape ``(nrows, *item_shape)``.  For scalar
        columns, the result is the logical row position within this column (or
        filtered view).
        """
        if self.is_ndarray:
            self._require_kind("biuf", "argmax")
            return self._ndarray_reduce("argmax", axis=axis, where=where)
        if axis not in (None, 0):
            return np.argmax(self[:], axis=axis)
        self._require_kind("biuf", "argmax")
        if where is not None:
            return self._table.where(self._normalize_sum_where(where))[self._col_name].argmax()
        arr = np.asarray(self[:])
        if arr.size == 0:
            raise ValueError("argmax() called on an empty column.")
        mask = (
            self._null_mask_for(arr) if self.null_value is not None else np.zeros(len(arr), dtype=np.bool_)
        )
        if mask.all():
            raise ValueError("argmax() called on a column where all values are null.")
        positions = np.where(~mask)[0]
        return int(positions[np.argmax(arr[positions])])

    def mean(self, axis=None, *, where=None):
        """Arithmetic mean of all live, non-null values.

        Supported dtypes: bool, int, uint, float.
        Null sentinel values are skipped. When *where* is provided, only rows
        matching the boolean predicate are included.
        Always returns a Python float.
        """
        if self.is_ndarray:
            self._require_kind("biuf", "mean")
            return self._ndarray_reduce("mean", axis=axis, where=where)
        if axis not in (None, 0):
            return np.mean(self[:], axis=axis)
        self._require_kind("biuf", "mean")
        where = self._normalize_sum_where(where)
        if where is None and len(self) == 0:
            if self._table.base is not None:
                return float("nan")
            self._require_nonempty("mean")
        fast = self._lazy_aggregate_fastpath("mean", where=where)
        if fast is not NotImplemented:
            return fast
        if where is not None:
            return self._table.where(where)[self._col_name].mean()
        total = np.float64(0)
        count = 0
        for chunk in self._nonnull_chunks():
            total += chunk.sum(dtype=np.float64)
            count += len(chunk)
        if count == 0:
            return float("nan")
        return float(total / count)

    def std(self, ddof: int = 0, axis=None, *, where=None):
        """Standard deviation of all live, non-null values (single-pass, Welford's algorithm).

        Parameters
        ----------
        ddof:
            Delta degrees of freedom.  ``0`` (default) gives the population
            std; ``1`` gives the sample std (divides by N-1).
        where:
            Optional boolean predicate. Only rows where the predicate is true,
            the table row is live, and this column is non-null are included.

        Supported dtypes: bool, int, uint, float.
        Null sentinel values are skipped.
        Always returns a Python float.
        """
        if self.is_ndarray:
            self._require_kind("biuf", "std")
            return self._ndarray_reduce("std", axis=axis, where=where, ddof=ddof)
        if axis not in (None, 0):
            return np.std(self[:], axis=axis, ddof=ddof)
        self._require_kind("biuf", "std")
        where = self._normalize_sum_where(where)
        if where is None and len(self) == 0:
            if self._table.base is not None:
                return float("nan")
            self._require_nonempty("std")
        fast = self._lazy_aggregate_fastpath("std", where=where, ddof=ddof)
        if fast is not NotImplemented:
            return fast
        if where is not None:
            return self._table.where(where)[self._col_name].std(ddof=ddof)

        # Chan's parallel update — combines per-chunk (n, mean, M2) tuples.
        # This is numerically stable and requires only a single pass.
        n_total = np.int64(0)
        mean_total = np.float64(0)
        M2_total = np.float64(0)

        for chunk in self._nonnull_chunks():
            chunk = chunk.astype(np.float64)
            n_b = np.int64(len(chunk))
            mean_b = chunk.mean()
            M2_b = np.float64(((chunk - mean_b) ** 2).sum())

            if n_total == 0:
                n_total, mean_total, M2_total = n_b, mean_b, M2_b
            else:
                delta = mean_b - mean_total
                n_new = n_total + n_b
                mean_total = (n_total * mean_total + n_b * mean_b) / n_new
                M2_total += M2_b + delta**2 * n_total * n_b / n_new
                n_total = n_new

        divisor = n_total - ddof
        if divisor <= 0:
            return float("nan")
        return float(np.sqrt(M2_total / divisor))

    def norm(self, ord=None, axis=None, *, where=None):
        """Vector/matrix norm of a fixed-shape ndarray column.

        The column is treated as a logical array of shape ``(nrows, *item_shape)``.
        For example, ``axis=1`` computes one norm per row for a 1-D item shape.
        """
        if not self.is_ndarray:
            raise TypeError(f"Column.norm() is only supported for ndarray columns, got {self._col_name!r}.")
        self._require_kind("biuf", "norm")
        arr = self._ndarray_values_for_reduction(where=where)
        return np.linalg.norm(arr, ord=ord, axis=axis)

    def any(self) -> bool:
        """Return True if at least one live, non-null value is True.

        Supported dtypes: bool.
        Null sentinel values are skipped.
        Short-circuits on the first True found.
        """
        self._require_kind("b", "any")
        return any(chunk.any() for chunk in self._nonnull_chunks())

    def all(self) -> bool:
        """Return True if every live, non-null value is True.

        Supported dtypes: bool.
        Null sentinel values are skipped.
        Short-circuits on the first False found.
        """
        self._require_kind("b", "all")
        return all(chunk.all() for chunk in self._nonnull_chunks())


# ---------------------------------------------------------------------------
# CTable
# ---------------------------------------------------------------------------


def _fmt_bytes(n: int) -> str:
    """Human-readable byte count (e.g. '1.23 MB')."""
    if n < 1024:
        return f"{n} B"
    if n < 1024**2:
        return f"{n / 1024:.2f} KB"
    if n < 1024**3:
        return f"{n / 1024**2:.2f} MB"
    return f"{n / 1024**3:.2f} GB"


_EXPECTED_SIZE_DEFAULT = 1_048_576
_BATCH_SIZE_DEFAULT = 2048

# ---------------------------------------------------------------------------
# Computed-column definition (virtual columns backed by a LazyExpr)
# ---------------------------------------------------------------------------

# Each entry in CTable._computed_cols maps column name → this dict shape:
#   {
#     "expression": str,        # LazyExpr.expression string (for serialization)
#     "col_deps":  list[str],   # dep column names in operand order (o0=col_deps[0], …)
#     "lazy":      LazyExpr,    # the live lazy expression (holds NDArray refs)
#     "dtype":     np.dtype,    # result dtype
#   }
# We use a plain dict so that nothing extra needs to be imported.


class _StructPathColumn:
    """Virtual read-only column representing a struct prefix path.

    Values are reconstructed per row from descendant dotted leaf columns.
    """

    def __init__(self, table: CTable, prefix: str, leaves: list[str]):
        self._table = table
        self._prefix = prefix
        self._leaves = list(leaves)

    def _leaf_is_null_at_logical(self, leaf: str, idx: int) -> bool:
        col = self._table[leaf]
        v = col[idx]
        nv = col.null_value
        if nv is None:
            return v is None
        try:
            return bool(col._null_mask_for(np.asarray([v]))[0])
        except Exception:
            return v is None

    def _row_value_at_logical(self, idx: int):
        # If every descendant leaf is null at this row, represent the struct as None.
        if self._leaves and all(self._leaf_is_null_at_logical(leaf, idx) for leaf in self._leaves):
            return None
        prefix_parts = split_field_path(self._prefix)
        result: dict[str, Any] = {}
        for leaf in self._leaves:
            parts = split_field_path(leaf)
            rel_parts = parts[len(prefix_parts) :]
            if not rel_parts:
                continue
            node = result
            for part in rel_parts[:-1]:
                child = node.get(part)
                if not isinstance(child, dict):
                    child = {}
                    node[part] = child
                node = child
            node[rel_parts[-1]] = self._table._normalize_scalar_value(self._table[leaf][idx])
        return result

    def __getitem__(self, key):
        if isinstance(key, int):
            return self._row_value_at_logical(key)
        if isinstance(key, slice):
            start, stop, step = key.indices(self._table.nrows)
            return [self._row_value_at_logical(i) for i in range(start, stop, step)]
        if isinstance(key, (list, np.ndarray)):
            if len(key) == 0:
                return []
            if isinstance(key, np.ndarray) and key.dtype == np.bool_:
                idxs = np.where(key)[0]
            elif isinstance(key[0], (bool, np.bool_)):
                idxs = [i for i, v in enumerate(key) if v]
            else:
                idxs = [int(i) for i in key]
            return [self._row_value_at_logical(i) for i in idxs]
        raise TypeError(f"Invalid index type: {type(key)}")

    def __iter__(self):
        for i in range(self._table.nrows):
            yield self._row_value_at_logical(i)


class NestedColumn:
    """A read-only accessor for a nested (dotted) group of CTable columns.

    Returned by attribute access on a :class:`CTable` (or on another
    ``NestedColumn``) when the name refers to an internal node of the dotted
    column tree rather than a leaf.  For a table flattened from a
    ``struct``/``list<struct>`` schema, ``t.trip`` is a ``NestedColumn``
    grouping every leaf under the ``trip.`` prefix, while a leaf such as
    ``t.trip.sec`` (or ``t.trip.begin.lon``) is a :class:`Column`.  Drilling
    into an intermediate node (e.g. ``t.trip.begin``) yields another
    ``NestedColumn``.

    Exposes aggregate metadata over its descendant leaf columns
    (:attr:`col_names`, :attr:`nrows`, :attr:`ncols`, :attr:`nbytes`,
    :attr:`cbytes`, :attr:`cratio`) and an :attr:`info` report.

    Examples
    --------
    >>> t.trip                      # doctest: +SKIP
    <NestedColumn 'trip'>
    >>> t.trip.col_names            # doctest: +SKIP
    ['sec', 'km', 'begin.lon', ...]
    >>> t.trip.sec                  # a leaf -> Column            # doctest: +SKIP
    """

    def __init__(self, table: CTable, prefix: str):
        self._table = table
        self._prefix = prefix

    def _descendant_col_names(self) -> list[str]:
        prefix_parts = split_field_path(self._prefix)
        return [
            name
            for name in self._table.col_names
            if (parts := split_field_path(name))[: len(prefix_parts)] == prefix_parts
            and len(parts) > len(prefix_parts)
        ]

    def _relative_col_name(self, name: str) -> str:
        prefix_parts = split_field_path(self._prefix)
        return join_field_path(split_field_path(name)[len(prefix_parts) :])

    @property
    def col_names(self) -> list[str]:
        """Descendant leaf column names relative to this nested prefix."""
        return [self._relative_col_name(name) for name in self._descendant_col_names()]

    @property
    def nrows(self) -> int:
        """Number of logical rows in this nested namespace."""
        return self._table.nrows

    @property
    def ncols(self) -> int:
        """Number of descendant leaf columns in this nested namespace."""
        return len(self._descendant_col_names())

    @property
    def nbytes(self) -> int:
        """Uncompressed size in bytes for stored descendant columns."""
        return sum(
            getattr(self._table._cols[name], "nbytes", 0)
            for name in self._descendant_col_names()
            if name in self._table._cols
        )

    @property
    def cbytes(self) -> int:
        """Compressed size in bytes for stored descendant columns."""
        return sum(
            getattr(self._table._cols[name], "cbytes", 0)
            for name in self._descendant_col_names()
            if name in self._table._cols
        )

    @property
    def cratio(self) -> float:
        """Compression ratio for stored descendant columns."""
        if self.cbytes == 0:
            return float("inf")
        return self.nbytes / self.cbytes

    @property
    def info_items(self) -> list[tuple[str, object]]:
        """Structured summary items used by :attr:`info`."""
        table = self._table
        storage_type = "persistent" if isinstance(table._storage, FileTableStorage) else "in-memory"
        column_summary = {}
        for name in self._descendant_col_names():
            rel_name = self._relative_col_name(name)
            if name in table._computed_cols:
                cc = table._computed_cols[name]
                column_summary[rel_name] = _InfoLiteral(
                    f"{cc['dtype']} (computed: {table._readable_computed_expr(cc)})"
                )
            else:
                col_meta = table._schema.columns_by_name.get(name)
                dtype_label = table._dtype_info_label(
                    getattr(table._cols[name], "dtype", None), col_meta.spec if col_meta else None
                )
                cbytes = getattr(table._cols[name], "cbytes", None)
                if cbytes is not None:
                    nbytes = getattr(table._cols[name], "nbytes", None)
                    detail = f"cbytes: {format_nbytes_human(cbytes)}"
                    if nbytes is not None and cbytes:
                        detail += f", cratio: {nbytes / cbytes:.2f}x"
                    column_summary[rel_name] = _InfoLiteral(f"{dtype_label} ({detail})")
                else:
                    column_summary[rel_name] = _InfoLiteral(dtype_label)

        descendant = set(self._descendant_col_names())
        index_summary = {}
        for idx in table.indexes:
            if idx.col_name not in descendant:
                continue
            stale = " stale" if idx.stale else ""
            label = f" name={idx.name!r}" if idx.name and idx.name != "__self__" else ""
            stats = idx.storage_stats()
            if stats is None:
                suffix = "(size=n/a, sidecars not directly addressable)"
            else:
                _, cbytes, _ = stats
                suffix = f"({format_nbytes_human(cbytes)})"
            index_summary[self._relative_col_name(idx.col_name)] = f"[{idx.kind}{stale}{label}] {suffix}"

        return [
            ("type", self.__class__.__name__),
            ("storage", storage_type),
            ("nrows", self.nrows),
            ("nbytes", format_nbytes_info(self.nbytes)),
            ("cbytes", format_nbytes_info(self.cbytes)),
            ("cratio", f"{self.cratio:.2f}x"),
            ("columns", column_summary),
            ("indexes", index_summary if index_summary else "none"),
        ]

    @property
    def info(self) -> _CTableInfoReporter:
        """Get information about this nested column namespace.

        Examples
        --------
        >>> print(t.trip.info)
        >>> t.trip.info()
        """
        return _CTableInfoReporter(self)

    def __getattr__(self, name: str):
        path = join_field_path((*split_field_path(self._prefix), name))
        if path in self._table._cols or path in self._table._computed_cols:
            return Column(self._table, path)
        path_parts = split_field_path(path)
        for col_name in self._table.col_names:
            parts = split_field_path(col_name)
            if parts[: len(path_parts)] == path_parts and len(parts) > len(path_parts):
                return NestedColumn(self._table, path)
        raise AttributeError(path)

    def __repr__(self) -> str:
        return f"<NestedColumn {self._prefix!r}>"


class _LazyColumnDict(dict):
    """Dict-like column cache that opens persistent columns on first use.

    Persistent CTables can be wide, and opening every stored column eagerly is
    expensive for workloads that touch only a small subset of columns, e.g.
    ``blosc2.open(path).trip.km.sum()`` on a nested table.  Keep the public and
    internal ``_cols`` access pattern mostly unchanged while deferring each
    ``storage.open_*_column()`` call until that column is actually requested.

    Methods that logically need all materialized columns, such as ``items()``
    and ``values()``, force-load the cache for compatibility with normal
    ``dict`` usage.  Name-oriented operations, such as ``keys()``, iteration,
    ``len()``, and ``in``, operate from the schema column list without opening
    the column payloads.
    """

    def __init__(
        self,
        table: CTable,
        storage: TableStorage,
        col_names: list[str],
        *,
        source_cols: dict | None = None,
    ):
        super().__init__()
        self._table = table
        self._storage = storage
        self._col_names = list(col_names)
        self._available = set(col_names)
        # When set, columns are projected lazily from another table's ``_cols``
        # mapping (e.g. a ``select()`` view) instead of opened from storage.
        # The source mapping itself opens on demand, so the same NDArray object
        # is shared (no copy) — identical to eager projection, just deferred.
        self._source_cols = source_cols

    def _load(self, name: str):
        if name not in self._available:
            raise KeyError(name)
        if not dict.__contains__(self, name):
            if self._source_cols is not None:
                value = self._source_cols[name]
            else:
                value = self._table._open_column_from_storage(self._storage, name)
            dict.__setitem__(self, name, value)
        return dict.__getitem__(self, name)

    def _load_all(self) -> None:
        for name in self._col_names:
            self._load(name)

    def __getitem__(self, name: str):
        return self._load(name)

    def get(self, name: str, default=None):
        return self._load(name) if name in self._available else default

    def __contains__(self, name: object) -> bool:
        return name in self._available

    def __iter__(self):
        return iter(self._col_names)

    def __len__(self) -> int:
        return len(self._col_names)

    def keys(self):
        return dict.fromkeys(self._col_names).keys()

    def items(self):
        self._load_all()
        return dict.items(self)

    def values(self):
        self._load_all()
        return dict.values(self)

    def __setitem__(self, name: str, value) -> None:
        if name not in self._available:
            self._available.add(name)
            self._col_names.append(name)
        dict.__setitem__(self, name, value)

    def __delitem__(self, name: str) -> None:
        self._available.remove(name)
        self._col_names.remove(name)
        if dict.__contains__(self, name):
            dict.__delitem__(self, name)


class _ChunkAlignedWriter:
    """Buffer writes to a fixed-size NDArray and flush them chunk-aligned.

    During Arrow/Parquet import the incoming batches have variable, non
    chunk-aligned sizes, so writing each one straight to ``arr[pos:pos+m]``
    makes most writes straddle chunk boundaries — forcing a
    decompress-merge-recompress of partially filled chunks.  This buffer
    accumulates appended arrays and writes them out in exact ``chunk_len``
    blocks aligned to chunk boundaries, so every chunk is compressed once.
    Only the final flush may write a partial (sub-chunk) tail.
    """

    __slots__ = ("arr", "chunk_len", "on_write", "pending", "pending_n", "wpos")

    def __init__(self, arr, chunk_len: int, on_write=None) -> None:
        self.arr = arr
        self.chunk_len = chunk_len
        self.pending: list[np.ndarray] = []
        self.pending_n = 0
        self.wpos = 0
        # Optional callback(start_pos, block) invoked for each chunk-aligned
        # write, used to fold per-block summaries incrementally.
        self.on_write = on_write

    def append(self, block: np.ndarray) -> None:
        if len(block) == 0:
            return
        self.pending.append(block)
        self.pending_n += len(block)
        while self.pending_n >= self.chunk_len:
            self._write(self._take(self.chunk_len))

    def flush(self) -> None:
        if self.pending_n:
            self._write(self._take(self.pending_n))

    def _write(self, block: np.ndarray) -> None:
        n = len(block)
        if self.on_write is not None:
            self.on_write(self.wpos, block)
        self.arr[self.wpos : self.wpos + n] = block
        self.wpos += n

    def _take(self, n: int) -> np.ndarray:
        """Pull exactly *n* rows from the front of the pending queue."""
        # Fast path: the head array already holds exactly the requested rows.
        if len(self.pending[0]) == n:
            self.pending_n -= n
            return self.pending.pop(0)
        parts: list[np.ndarray] = []
        need = n
        while need > 0:
            head = self.pending[0]
            if len(head) <= need:
                parts.append(head)
                need -= len(head)
                self.pending.pop(0)
            else:
                parts.append(head[:need])
                self.pending[0] = head[need:]
                need = 0
        self.pending_n -= n
        return parts[0] if len(parts) == 1 else np.concatenate(parts)


class _ColumnSummaryAccumulator:
    """Incrementally accumulate per-block min/max for a SUMMARY index.

    As contiguous numpy data is appended to a scalar column during a build
    (import, ``extend``, row ``append``), this folds it into per-block min/max
    summaries while the data is still uncompressed in memory.  At close, the
    accumulated summaries are handed to the index builder, avoiding a full
    decompression pass over the column just to recompute min/max.

    The accumulator only stays valid while writes are pure forward-contiguous
    appends covering ``[0, n)``.  Any out-of-order feed, in-place update, or
    compaction marks it invalid via :meth:`invalidate`, and the index builder
    transparently falls back to the out-of-core (decompress-and-scan) path.
    """

    __slots__ = ("_carry", "_parts", "_total", "block_len", "dtype", "summary_dtype", "valid")

    def __init__(self, dtype: np.dtype, block_len: int) -> None:
        from blosc2.indexing import _summary_dtype

        self.dtype = np.dtype(dtype)
        self.block_len = int(block_len)
        self.summary_dtype = _summary_dtype(self.dtype)
        self._parts: list[np.ndarray] = []  # completed per-block summary chunks
        self._carry: np.ndarray | None = None  # trailing < block_len values
        self._total = 0  # number of elements fed (== next expected start_pos)
        self.valid = self.block_len > 0

    def invalidate(self) -> None:
        self.valid = False
        self._parts = []
        self._carry = None

    def feed(self, start_pos: int, values: np.ndarray) -> None:
        if not self.valid:
            return
        if start_pos != self._total:
            self.invalidate()
            return
        from blosc2.indexing import _fill_summaries_from_2d

        if values.ndim != 1:
            # Not a plain scalar column write (e.g. an ndarray column); summaries
            # don't apply -- disable rather than risk a wrong result.
            self.invalidate()
            return
        values = np.ascontiguousarray(values)
        m = values.shape[0]
        if m == 0:
            return
        buf = values if self._carry is None else np.concatenate([self._carry, values])
        n = buf.shape[0]
        n_complete = n // self.block_len
        if n_complete:
            data_2d = buf[: n_complete * self.block_len].reshape(n_complete, self.block_len)
            block_summ = np.empty(n_complete, dtype=self.summary_dtype)
            _fill_summaries_from_2d(data_2d, block_summ, 0, self.dtype)
            self._parts.append(block_summ)
        rem = n - n_complete * self.block_len
        self._carry = buf[n_complete * self.block_len :].copy() if rem else None
        self._total = start_pos + m

    def finalize(self, expected_size: int) -> np.ndarray | None:
        """Return the full per-block summary array, or None if unusable.

        *expected_size* is the column's physical length the index will summarize;
        the accumulator must have covered exactly that many elements.
        """
        if not self.valid or self._total != expected_size:
            return None
        from blosc2.indexing import _fill_summaries_from_2d

        parts = list(self._parts)
        if self._carry is not None and self._carry.shape[0]:
            tail = np.empty(1, dtype=self.summary_dtype)
            _fill_summaries_from_2d(self._carry.reshape(1, -1), tail, 0, self.dtype)
            parts.append(tail)
        if not parts:
            return np.empty(0, dtype=self.summary_dtype)
        return np.concatenate(parts)


class ColExpr:
    """Unbound column expression: a recipe that, given a table, evaluates
    against that table's columns.

    ``blosc2.col("x") + 1`` builds a deferred computation; passing it to
    :meth:`CTable.assign` or using it to index/filter a table binds it,
    replaying the operators on ``table["x"]`` — so all :class:`Column`
    semantics (null propagation, SQL comparison rules, dictionary/timestamp
    handling) apply identically to the bound form ``t.x + 1``.

    Only operators are supported; method calls such as ``col("x").sum()``
    are not, since there is no table to evaluate against yet. Use the bound
    form (``t.x.sum()``) for those.
    """

    def __init__(self, bind, repr_str):
        self._bind = bind  # callable: CTable -> Column/LazyExpr/NullableExpr/scalar
        self._repr = repr_str

    def __repr__(self):
        return self._repr

    def __getattr__(self, name):
        raise AttributeError(
            f"{name!r} is not supported on an unbound column expression (only operators are). "
            f"Bind it to a table first, e.g. use t.{self._repr}.{name}(...) instead of "
            f"{self._repr}.{name}(...)."
        )


def _colexpr_binop(op, symbol, reflected=False):
    def bind(self, other, t):
        left = self._bind(t)
        right = other._bind(t) if isinstance(other, ColExpr) else other
        return op(right, left) if reflected else op(left, right)

    def method(self, other):
        other_repr = other._repr if isinstance(other, ColExpr) else repr(other)
        if reflected:
            repr_str = f"({other_repr} {symbol} {self._repr})"
        else:
            repr_str = f"({self._repr} {symbol} {other_repr})"
        return ColExpr(lambda t: bind(self, other, t), repr_str)

    return method


def _colexpr_unary(op, symbol):
    def method(self):
        return ColExpr(lambda t: op(self._bind(t)), f"({symbol}{self._repr})")

    return method


_COLEXPR_BINOPS = [
    ("__add__", operator.add, "+", False),
    ("__radd__", operator.add, "+", True),
    ("__sub__", operator.sub, "-", False),
    ("__rsub__", operator.sub, "-", True),
    ("__mul__", operator.mul, "*", False),
    ("__rmul__", operator.mul, "*", True),
    ("__truediv__", operator.truediv, "/", False),
    ("__rtruediv__", operator.truediv, "/", True),
    ("__floordiv__", operator.floordiv, "//", False),
    ("__rfloordiv__", operator.floordiv, "//", True),
    ("__mod__", operator.mod, "%", False),
    ("__rmod__", operator.mod, "%", True),
    ("__pow__", operator.pow, "**", False),
    ("__rpow__", operator.pow, "**", True),
    ("__and__", operator.and_, "&", False),
    ("__rand__", operator.and_, "&", True),
    ("__or__", operator.or_, "|", False),
    ("__ror__", operator.or_, "|", True),
    ("__lt__", operator.lt, "<", False),
    ("__le__", operator.le, "<=", False),
    ("__gt__", operator.gt, ">", False),
    ("__ge__", operator.ge, ">=", False),
    ("__eq__", operator.eq, "==", False),
    ("__ne__", operator.ne, "!=", False),
]
for _dunder, _op, _symbol, _reflected in _COLEXPR_BINOPS:
    setattr(ColExpr, _dunder, _colexpr_binop(_op, _symbol, _reflected))
ColExpr.__neg__ = _colexpr_unary(operator.neg, "-")
ColExpr.__invert__ = _colexpr_unary(operator.invert, "~")
del _dunder, _op, _symbol, _reflected


def col(name: str) -> ColExpr:
    """Build an unbound column expression referencing a column by name.

    The name is resolved against a table's columns only when the expression
    is bound — passed to :meth:`CTable.assign`, or used to index/filter a
    table (``t[col("x") > 0]``, ``t.where(col("x") > 0)``). An unknown name
    therefore fails at bind time with the table's normal unknown-column
    error, not when ``col()`` is called.

    Examples
    --------
    >>> import blosc2
    >>> from blosc2 import col
    >>> t.assign(profit=col("revenue") - col("cost"))  # doctest: +SKIP
    """
    return ColExpr(lambda t: t[name], name)


class CTable(_CTableIndexingMixin, Generic[RowT]):
    """Columnar compressed table with typed columns and row-oriented access."""

    #: Ordered list of stored column names.  Computed columns are **not**
    #: included; access those via :attr:`computed_columns`.
    col_names: list[str]

    #: Parent table when this instance is a row-filter or column-projection
    #: view (created by :meth:`where`, :meth:`select`, or :meth:`view`).
    #: ``None`` for top-level tables.  Structural mutations such as
    #: :meth:`add_column` and :meth:`drop_column` are blocked on views.
    base: CTable | None

    @property
    def _n_rows(self) -> int:
        """Number of live rows, computed lazily for reopened tables."""
        n_rows = getattr(self, "_n_rows_cached", None)
        if n_rows is None:
            n_rows = int(blosc2.count_nonzero(self._valid_rows))
            self._n_rows_cached = n_rows
        return n_rows

    @_n_rows.setter
    def _n_rows(self, value: int | None) -> None:
        self._n_rows_cached = value

    def _known_n_rows(self) -> int | None:
        """Return cached live-row count without triggering a scan."""
        return getattr(self, "_n_rows_cached", None)

    def _iter_live_positions_chunks(self):
        """Yield chunks of physical positions for live rows without materialising the full mask."""
        valid_rows = self._valid_rows
        n = len(valid_rows)
        chunks = getattr(valid_rows, "chunks", None)
        chunk_len = chunks[0] if chunks else n

        for start in range(0, n, chunk_len):
            stop = min(start + chunk_len, n)
            local_pos = np.flatnonzero(valid_rows[start:stop])
            if len(local_pos):
                yield (local_pos + start).astype(np.intp, copy=False)

    def _live_positions_from_valid_rows_chunks(self) -> np.ndarray:
        """Return live physical row positions by scanning the validity NDArray chunk-wise."""
        cached = getattr(self, "_cached_live_positions", None)
        if cached is not None:
            return cached
        positions = list(self._iter_live_positions_chunks())
        if not positions:
            result = np.empty(0, dtype=np.intp)
        elif len(positions) == 1:
            result = positions[0]
        else:
            result = np.concatenate(positions).astype(np.intp, copy=False)
        if self.base is not None:
            self._cached_live_positions = result
        return result

    def __init__(
        self,
        row_type: type[RowT],
        new_data=None,
        *,
        urlpath: str | None = None,
        mode: str = "a",
        expected_size: int | None = None,
        compact: bool = False,
        validate: bool = True,
        cparams: dict[str, Any] | None = None,
        dparams: dict[str, Any] | None = None,
        create_summary_index: bool = True,
    ) -> None:
        """Create a new CTable or open an existing one.

        Parameters
        ----------
        create_summary_index:
            If ``True`` (default), SUMMARY indexes are automatically built for
            all eligible scalar columns.  These indexes are extremely cheap to
            store (< 0.1% of column size) and accelerate ``where()`` queries
            without any user action.  Set to ``False`` to disable.

            The build is triggered by :meth:`close`, not by table creation, so
            *when* it happens depends on the table's lifecycle:

            - **Persistent** tables (``urlpath=...``) are closed as part of
              normal use, so they get these indexes and reopen with them.
            - A **purely in-memory** table is never closed automatically, so it
              is *not* indexed unless you close it explicitly or use it as a
              context manager (``with blosc2.CTable(...) as t:``).  Otherwise
              call :meth:`create_index` yourself.

            Note that :meth:`to_b2z` and :meth:`save` write live rows through a
            logical copy and do **not** trigger the build; index the source
            table (or the reopened result) explicitly if you need it.
        """
        # Auto-size: if the caller didn't specify expected_size and new_data has a
        # known length, pre-allocate just enough (×2 for headroom, min 64).
        # Fall back to 1 M when new_data has no __len__ or is absent.
        if expected_size is None:
            if new_data is not None and hasattr(new_data, "__len__"):
                expected_size = max(len(new_data) * 2, 64)
            else:
                expected_size = _EXPECTED_SIZE_DEFAULT
        self._row_type = row_type
        self._validate = validate
        self._table_cparams = cparams
        self._table_dparams = dparams
        self._cols: dict[str, blosc2.NDArray | ListArray] = {}
        self._computed_cols: dict[str, dict] = {}  # virtual/computed columns
        self._materialized_cols: dict[str, dict] = {}  # stored columns auto-filled from expressions
        self._expr_index_arrays: dict[str, blosc2.NDArray] = {}
        self._cached_index_catalog: dict | None = None
        self._cached_index_catalog_revision: int | None = None
        self._cached_live_positions: np.ndarray | None = None
        self._col_widths: dict[str, int] = {}
        self.col_names: list[str] = []
        self.auto_compact = compact
        self._create_summary_index = create_summary_index
        self._summary_indexes_built = False
        self.base = None

        # Choose storage backend
        if urlpath is not None:
            if mode == "w" and os.path.exists(urlpath):
                if os.path.isdir(urlpath):
                    shutil.rmtree(urlpath)
                else:
                    os.remove(urlpath)
            storage: TableStorage = FileTableStorage(urlpath, mode)
        else:
            storage = InMemoryTableStorage()
        self._storage = storage
        self._read_only = storage.is_read_only()

        if storage.table_exists() and mode != "w":
            # ---- Open existing persistent table ----
            if new_data is not None:
                raise ValueError(
                    "Cannot pass new_data when opening an existing table. Use mode='w' to overwrite."
                )
            storage.check_kind()
            schema_dict = storage.load_schema()
            self._schema: CompiledSchema = schema_from_dict(schema_dict)
            self._schema = CompiledSchema(
                row_cls=row_type,
                columns=self._schema.columns,
                columns_by_name=self._schema.columns_by_name,
            )
            self.col_names = [c["name"] for c in schema_dict["columns"]]
            self._valid_rows = storage.open_valid_rows()
            self._cols = _LazyColumnDict(self, storage, self.col_names)
            for name in self.col_names:
                cc = self._schema.columns_by_name[name]
                self._col_widths[name] = max(len(name), cc.display_width)
            self._n_rows = None
            # Restore cached row count from saved metadata so that
            # where() can skip the _valid_rows intersection for all-valid tables.
            if "n_rows" in schema_dict:
                self._n_rows_cached = schema_dict["n_rows"]
            self._last_pos = None  # resolve lazily on first write
            # ---- Restore computed/materialized column metadata (if any) ----
            self._computed_cols = {}
            self._materialized_cols = {}
            self._expr_index_arrays = {}
            self._load_computed_cols_from_schema(schema_dict)
            self._load_materialized_cols_from_schema(schema_dict)
            # Restore auto-index preference from the schema.
            self._create_summary_index = schema_dict.get("create_summary_index", True)
            self._summary_indexes_built = schema_dict.get("summary_indexes_built", False)
        else:
            # ---- Create new table ----
            if storage.is_read_only():
                raise FileNotFoundError(f"No CTable found at {urlpath!r}")
            if urlpath is not None and mode == "a":
                raise FileNotFoundError(
                    f"No CTable found at {urlpath!r}: mode='a' opens an existing table; "
                    "use mode='w' to create a new one."
                )

            # Build compiled schema from either a dataclass or a legacy Pydantic model
            if dataclasses.is_dataclass(row_type) and isinstance(row_type, type):
                self._schema = compile_schema(row_type)
            else:
                self._schema = _compile_pydantic_schema(row_type)
            self._resolve_nullable_specs(self._schema)

            self._n_rows = 0
            self._last_pos = 0

            default_chunks, default_blocks = compute_chunks_blocks((expected_size,))
            # Compute the table-wide shared grid once so both the _valid_rows
            # mask and the fixed-size columns use it; this keeps where() on the
            # fast_eval path (the boolean mask is combined with the condition).
            shared_chunks, shared_blocks, aligned_names = self._compute_aligned_grid(
                self._schema.columns, expected_size
            )
            valid_chunks = shared_chunks if shared_chunks is not None else default_chunks
            valid_blocks = shared_blocks if shared_blocks is not None else default_blocks
            self._valid_rows = storage.create_valid_rows(
                shape=(expected_size,),
                chunks=valid_chunks,
                blocks=valid_blocks,
            )
            self._init_columns(
                expected_size,
                default_chunks,
                default_blocks,
                storage,
                aligned_grid=(shared_chunks, shared_blocks, aligned_names),
            )
            storage.save_schema(schema_to_dict(self._schema))

            if new_data is not None:
                self._load_initial_data(new_data)
                # Persist the row count so subsequent opens can skip the
                # _valid_rows intersection in where().
                self._save_n_rows_to_meta()

    def close(self) -> None:
        """Close any persistent backing store held by this table.

        On the first close of a writable root table, this also builds the
        automatic SUMMARY indexes (unless ``create_summary_index=False``); see
        the ``create_summary_index`` parameter of :class:`CTable` for how this
        interacts with in-memory vs. persistent tables.
        """
        storage = getattr(self, "_storage", None)
        # Persist row count for root tables so subsequent opens can skip
        # the _valid_rows intersection in where() for all-valid tables.
        if not self._read_only and self.base is None:
            self._save_n_rows_to_meta()
        # Persist user vlmeta if a dedicated SChunk was created
        if storage is not None:
            uv = getattr(storage, "_vlmeta", None)
            if uv is not None and hasattr(storage, "save_vlmeta"):
                storage.save_vlmeta(uv)
        try:
            self._flush_varlen_columns()
            if not self._read_only and self.base is None:
                self.trim_capacity()
            # Build SUMMARY indexes for eligible columns on first close
            # (one-time).  These are cheap (~<0.1% of column size) and
            # accelerate queries without any user action.
            if (
                not self._read_only
                and self.base is None
                and getattr(self, "_create_summary_index", True)
                and not getattr(self, "_summary_indexes_built", False)
            ):
                self._build_summary_indexes()
                # Persist that indexes have been built so subsequent opens
                # skip the catalog check.
                self._save_n_rows_to_meta()
        except Exception:
            with contextlib.suppress(Exception):
                if storage is not None and hasattr(storage, "close"):
                    storage.close()
            raise
        if storage is not None and hasattr(storage, "close"):
            storage.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    def __del__(self):
        with contextlib.suppress(Exception):
            storage = getattr(self, "_storage", None)
            if storage is not None and hasattr(storage, "discard"):
                storage.discard()
            elif storage is not None and hasattr(storage, "close"):
                storage.close()

    @staticmethod
    def _is_list_column(col: CompiledColumn) -> bool:
        return isinstance(col.spec, ListSpec)

    @staticmethod
    def _is_varlen_scalar_column(col: CompiledColumn) -> bool:
        return isinstance(col.spec, (VLStringSpec, VLBytesSpec, StructSpec, ObjectSpec, Utf8Spec))

    @staticmethod
    def _is_utf8_column(col: CompiledColumn) -> bool:
        return isinstance(col.spec, Utf8Spec)

    @staticmethod
    def _is_dictionary_column(col: CompiledColumn) -> bool:
        return isinstance(col.spec, DictionarySpec)

    def _dict_rank_index_stale(self, name: str, dict_rank_meta: dict) -> bool:
        """True if a dict-rank FULL index no longer matches the live dictionary.

        The index encodes alphabetical ranks frozen at build time; if the
        dictionary gained/changed entries the ranks are wrong, so callers must
        fall back to lexsort until the index is rebuilt.
        """
        from blosc2.ctable_indexing import _dict_rank_hash

        col = self._root_table._cols.get(name)
        if col is None:
            return True
        dictionary = list(col.dictionary)
        if len(dictionary) != dict_rank_meta.get("dict_len"):
            return True
        return _dict_rank_hash(dictionary) != dict_rank_meta.get("dict_hash")

    @staticmethod
    def _is_ndarray_column(col: CompiledColumn) -> bool:
        return isinstance(col.spec, NDArraySpec)

    @staticmethod
    def _column_physical_shape(col: CompiledColumn, capacity: int) -> tuple[int, ...]:
        if CTable._is_ndarray_column(col):
            return (capacity, *col.spec.item_shape)
        return (capacity,)

    @staticmethod
    def _ndarray_null_item(spec: NDArraySpec) -> np.ndarray:
        null_value = getattr(spec, "null_value", None)
        if null_value is None:
            raise TypeError("NDArraySpec is not nullable")
        return np.full(spec.item_shape, null_value, dtype=spec.dtype)

    @staticmethod
    def _coerce_ndarray_value(name: str, spec: NDArraySpec, val) -> np.ndarray:
        if val is None:
            if getattr(spec, "null_value", None) is None:
                raise TypeError(f"Column {name!r} is not nullable; received None.")
            return CTable._ndarray_null_item(spec)
        arr = np.asarray(val, dtype=spec.dtype)
        if arr.shape != spec.item_shape:
            raise ValueError(f"Column {name!r}: expected item shape {spec.item_shape}, got {arr.shape}")
        return np.ascontiguousarray(arr)

    @staticmethod
    def _coerce_ndarray_batch(name: str, spec: NDArraySpec, values, nrows: int) -> np.ndarray:
        if values is None:
            null_item = CTable._coerce_ndarray_value(name, spec, None)
            return np.broadcast_to(null_item, (nrows, *spec.item_shape)).copy()
        if isinstance(values, np.ndarray) and values.dtype != object:
            arr = np.ascontiguousarray(values, dtype=spec.dtype)
            if arr.ndim == len(spec.item_shape):
                arr = arr.reshape((1, *arr.shape))
            if arr.shape != (nrows, *spec.item_shape):
                raise ValueError(
                    f"Column {name!r}: expected batch shape {(nrows, *spec.item_shape)}, got {arr.shape}"
                )
            return arr
        rows = [CTable._coerce_ndarray_value(name, spec, value) for value in values]
        arr = np.ascontiguousarray(rows, dtype=spec.dtype)
        if arr.shape != (nrows, *spec.item_shape):
            raise ValueError(
                f"Column {name!r}: expected batch shape {(nrows, *spec.item_shape)}, got {arr.shape}"
            )
        return arr

    @staticmethod
    def _column_chunks_blocks(col: CompiledColumn, shape: tuple[int, ...]):
        return compute_chunks_blocks(shape, dtype=col.dtype)

    @staticmethod
    def _is_list_spec(spec: SchemaSpec) -> bool:
        return isinstance(spec, ListSpec)

    @staticmethod
    def _policy_null_value_for_spec(spec: SchemaSpec, policy: NullPolicy):
        if isinstance(spec, NDArraySpec):
            dtype = spec.dtype
            if dtype == np.dtype(np.bool_):
                return policy.bool_value
            if dtype.kind == "i":
                info = np.iinfo(dtype)
                return info.min if policy.signed_int_strategy == "min" else info.max
            if dtype.kind == "u":
                info = np.iinfo(dtype)
                return info.min if policy.unsigned_int_strategy == "min" else info.max
            if dtype.kind == "f":
                return policy.float_value
            return None
        if isinstance(spec, (int8, int16, int32, int64)):
            info = np.iinfo(spec.dtype)
            return info.min if policy.signed_int_strategy == "min" else info.max
        if isinstance(spec, (uint8, uint16, uint32, uint64)):
            info = np.iinfo(spec.dtype)
            return info.min if policy.unsigned_int_strategy == "min" else info.max
        if isinstance(spec, (float32, float64)):
            return policy.float_value
        if isinstance(spec, b2_bool):
            return policy.bool_value
        if isinstance(spec, (string, Utf8Spec)):
            return policy.string_value
        if isinstance(spec, b2_bytes):
            return policy.bytes_value
        if isinstance(spec, timestamp):
            return policy.timestamp_value
        return None

    @staticmethod
    def _validate_null_value_for_spec(name: str, spec: SchemaSpec, null_value) -> None:  # noqa: C901
        if isinstance(spec, NDArraySpec):
            dtype = spec.dtype
            if dtype == np.dtype(np.bool_):
                if null_value != 255:
                    raise ValueError(f"Null sentinel for nullable bool ndarray column {name!r} must be 255")
                return
            if dtype.kind in "iu":
                if isinstance(null_value, (bool, np.bool_)) or not isinstance(null_value, (int, np.integer)):
                    raise TypeError(f"Null sentinel for ndarray column {name!r} must be an integer")
                info = np.iinfo(dtype)
                if not info.min <= int(null_value) <= info.max:
                    raise ValueError(
                        f"Null sentinel for ndarray column {name!r}={null_value!r} is outside {dtype} range"
                    )
                return
            if dtype.kind == "f":
                if not isinstance(null_value, (int, float, np.integer, np.floating)):
                    raise TypeError(f"Null sentinel for ndarray column {name!r} must be numeric")
                return
            raise TypeError(
                f"Nullable ndarray column {name!r} has unsupported dtype {dtype!r}; "
                "use bool, integer, unsigned integer, or floating dtype."
            )
        if isinstance(spec, (int8, int16, int32, int64, uint8, uint16, uint32, uint64, timestamp)):
            if isinstance(null_value, (bool, np.bool_)) or not isinstance(null_value, (int, np.integer)):
                raise TypeError(f"Null sentinel for column {name!r} must be an integer")
            info = np.iinfo(spec.dtype)
            if not info.min <= int(null_value) <= info.max:
                raise ValueError(
                    f"Null sentinel for column {name!r}={null_value!r} is outside {spec.dtype} range"
                )
            return
        if isinstance(spec, (float32, float64)):
            if not isinstance(null_value, (int, float, np.integer, np.floating)):
                raise TypeError(f"Null sentinel for column {name!r} must be numeric")
            return
        if isinstance(spec, b2_bool):
            if null_value != 255:
                raise ValueError(f"Null sentinel for nullable bool column {name!r} must be 255")
            return
        if isinstance(spec, (string, Utf8Spec)):
            if not isinstance(null_value, str):
                raise TypeError(f"Null sentinel for string column {name!r} must be str")
            return
        if isinstance(spec, b2_bytes) and not isinstance(null_value, bytes):
            raise TypeError(f"Null sentinel for bytes column {name!r} must be bytes")

    @classmethod
    def _resolve_nullable_specs(
        cls, schema: CompiledSchema, *, validate_column_null_values: bool = True
    ) -> None:
        policy = get_null_policy()
        schema_names = {col.name for col in schema.columns}
        unknown_null_values = set(policy.column_null_values) - schema_names
        if validate_column_null_values and unknown_null_values:
            names = ", ".join(sorted(unknown_null_values))
            raise KeyError(f"column_null_values contains unknown columns: {names}")
        for col in schema.columns:
            spec = col.spec
            if isinstance(spec, NDArraySpec) and getattr(spec, "null_value", None) is not None:
                cls._validate_null_value_for_spec(col.name, spec, spec.null_value)
                if spec.dtype == np.dtype(np.bool_):
                    spec.dtype = np.dtype(np.uint8)
                    spec.itemsize = spec.dtype.itemsize
                    spec.kind = spec.dtype.kind
                    spec.type = spec.dtype.type
                    spec.str = spec.dtype.str
                    spec.name = spec.dtype.name
                col.dtype = getattr(spec, "dtype", None)
                col.display_width = compute_display_width(spec)
                continue
            if (
                isinstance(
                    spec,
                    (ListSpec, VLStringSpec, VLBytesSpec, StructSpec, ObjectSpec, DictionarySpec),
                )
                or getattr(spec, "null_value", None) is not None
            ):
                continue
            if not getattr(spec, "nullable", False):
                continue
            null_value = policy.column_null_values.get(col.name)
            if null_value is None:
                null_value = cls._policy_null_value_for_spec(spec, policy)
            if null_value is None:
                raise TypeError(f"Column {col.name!r} is nullable, but no null policy sentinel is available")
            cls._validate_null_value_for_spec(col.name, spec, null_value)
            spec.null_value = null_value
            if isinstance(spec, string):
                spec.max_length = max(spec.max_length, len(null_value), 1)
                spec.dtype = np.dtype(f"U{spec.max_length}")
            elif isinstance(spec, b2_bytes):
                spec.max_length = max(spec.max_length, len(null_value), 1)
                spec.dtype = np.dtype(f"S{spec.max_length}")
            elif isinstance(spec, b2_bool):
                spec.dtype = np.dtype(np.uint8)
            elif isinstance(spec, NDArraySpec) and spec.dtype == np.dtype(np.bool_):
                spec.dtype = np.dtype(np.uint8)
                spec.itemsize = spec.dtype.itemsize
                spec.kind = spec.dtype.kind
                spec.type = spec.dtype.type
                spec.str = spec.dtype.str
                spec.name = spec.dtype.name
            col.dtype = getattr(spec, "dtype", None)
            col.display_width = compute_display_width(spec)

    def _flush_varlen_columns(self) -> None:
        for col in self._schema.columns:
            if (
                self._is_list_column(col)
                or self._is_varlen_scalar_column(col)
                or self._is_dictionary_column(col)
            ):
                self._cols[col.name].flush()

    # Common itemsizes we snap the representative (median) itemsize to when
    # computing the table-wide shared chunk/block grid.
    _COMMON_ITEMSIZES = (1, 2, 4, 8, 16)

    # Fixed-width string/bytes columns up to this many bytes join the shared
    # grid (so string filters stay on the fast path); ``U32`` is 128 bytes
    # under NumPy's 4-bytes-per-char Unicode dtype.  Larger string columns are
    # better stored as dictionary columns.
    _MAX_ALIGNED_STR_ITEMSIZE = 128

    @staticmethod
    def _snap_itemsize(median: float) -> int:
        """Snap *median* to the nearest value in :attr:`_COMMON_ITEMSIZES`.

        Ties round down (the strict ``<`` comparison keeps the first, smaller
        candidate), so a median of 6 snaps to 4 rather than 8.
        """
        best = CTable._COMMON_ITEMSIZES[0]
        best_dist = abs(median - best)
        for value in CTable._COMMON_ITEMSIZES[1:]:
            dist = abs(median - value)
            if dist < best_dist:
                best, best_dist = value, dist
        return best

    @classmethod
    def _compute_aligned_grid(cls, columns, capacity: int):
        """Compute a single chunk/block grid shared by fixed-size columns.

        All 1-D fixed-size scalar columns (no user-pinned grid) are sized to one
        common ``(chunks, blocks)`` so that lazy expressions over them take the
        ``fast_eval`` path (which requires identical element-unit chunk/block
        grids across operands).  The same grid is also applied to the
        ``_valid_rows`` mask so that ``where()`` keeps the fast path when it
        combines the condition with the mask.

        The grid is sized for the *median* itemsize of the eligible numeric
        columns (the operands that dominate fused arithmetic), snapped to the
        nearest common itemsize.  Numeric columns join the aligned set only if
        the shared grid does not blow their chunk size past ~4x what they would
        pick on their own; this keeps wide numeric columns on per-dtype sizing.

        Fixed-width string/bytes columns are kept *out* of the median (so they
        don't coarsen the numeric grid) but still join the aligned set when
        their itemsize is at most :attr:`_MAX_ALIGNED_STR_ITEMSIZE`, so string
        filters stay on the fast path.  Wider string columns (e.g. ``U183642``)
        keep per-dtype sizing instead of producing multi-GB chunks.

        ``columns`` is an iterable of :class:`CompiledColumn`.  Returns a
        ``(shared_chunks, shared_blocks, included_names)`` tuple, where
        ``included_names`` is the set of column names that should use the shared
        grid.  Returns ``(None, None, set())`` when there is nothing to align.
        """
        numeric, strings = [], []
        for col in columns:
            if (
                cls._is_list_column(col)
                or cls._is_varlen_scalar_column(col)
                or cls._is_dictionary_column(col)
            ):
                continue
            if col.dtype is None:
                continue
            if col.config.chunks is not None or col.config.blocks is not None:
                continue
            if len(cls._column_physical_shape(col, capacity)) != 1:
                continue
            if np.dtype(col.dtype).kind in ("U", "S"):
                strings.append(col)
            else:
                numeric.append(col)

        if not numeric and not strings:
            return None, None, set()

        # Size the grid from the numeric columns; if the table is all strings,
        # fall back to sizing it from the string columns instead.
        basis = numeric or strings
        itemsizes = [np.dtype(col.dtype).itemsize for col in basis]
        snapped = cls._snap_itemsize(float(np.median(itemsizes)))
        shared_chunks, shared_blocks = compute_chunks_blocks((capacity,), dtype=np.dtype(f"V{snapped}"))

        included = set()
        for col in numeric:
            natural_chunks, _ = cls._column_chunks_blocks(col, cls._column_physical_shape(col, capacity))
            # Only align numeric columns whose shared-grid chunk stays within
            # ~4x the chunk they would choose on their own (in rows; itemsize
            # cancels).
            if shared_chunks[0] <= 4 * natural_chunks[0]:
                included.add(col.name)
        for col in strings:
            # Fixed strings join via an absolute byte ceiling, not the relative
            # cap: they fast-path equality filters and a few-MB block is fine.
            if np.dtype(col.dtype).itemsize <= cls._MAX_ALIGNED_STR_ITEMSIZE:
                included.add(col.name)

        if not included:
            return None, None, set()
        return shared_chunks, shared_blocks, included

    def _init_columns(
        self,
        expected_size: int,
        default_chunks,
        default_blocks,
        storage: TableStorage,
        aligned_grid: tuple | None = None,
    ) -> None:
        """Create one physical column per compiled schema column."""
        if aligned_grid is None:
            aligned_grid = self._compute_aligned_grid(self._schema.columns, expected_size)
        shared_chunks, shared_blocks, aligned_names = aligned_grid
        for col in self._schema.columns:
            self.col_names.append(col.name)
            self._col_widths[col.name] = max(len(col.name), col.display_width)
            col_storage = self._resolve_column_storage(col, default_chunks, default_blocks)
            if self._is_list_column(col):
                self._cols[col.name] = storage.create_list_column(
                    col.name,
                    spec=col.spec,
                    cparams=col_storage.get("cparams"),
                    dparams=col_storage.get("dparams"),
                )
                continue
            if self._is_varlen_scalar_column(col):
                self._cols[col.name] = storage.create_varlen_scalar_column(
                    col.name,
                    spec=col.spec,
                    cparams=col_storage.get("cparams"),
                    dparams=col_storage.get("dparams"),
                )
                continue
            if self._is_dictionary_column(col):
                dict_col = storage.create_dictionary_column(
                    col.name,
                    spec=col.spec,
                    cparams=col_storage.get("cparams"),
                    dparams=col_storage.get("dparams"),
                )
                if len(dict_col.codes) < expected_size:
                    dict_col.resize((expected_size,))
                self._cols[col.name] = dict_col
                continue
            # Recompute chunks/blocks using the actual dtype so that wide
            # string columns (e.g. U183642) don't produce multi-GB chunks.
            chunks = col_storage["chunks"]
            blocks = col_storage["blocks"]
            shape = self._column_physical_shape(col, expected_size)
            if col.config.chunks is None and col.config.blocks is None:
                if col.name in aligned_names:
                    # Use the table-wide shared grid so lazy expressions over
                    # this column take the fast_eval path.
                    chunks, blocks = shared_chunks, shared_blocks
                else:
                    chunks, blocks = self._column_chunks_blocks(col, shape)
            self._cols[col.name] = storage.create_column(
                col.name,
                dtype=col.dtype,
                shape=shape,
                chunks=chunks,
                blocks=blocks,
                cparams=col_storage.get("cparams"),
                dparams=col_storage.get("dparams"),
            )

    def _resolve_column_storage(
        self,
        col: CompiledColumn,
        default_chunks,
        default_blocks,
    ) -> dict[str, Any]:
        """Merge table-level and column-level storage settings.

        Column-level settings (from ``b2.field(...)``) take precedence over
        table-level defaults passed to ``CTable.__init__``.
        """
        result: dict[str, Any] = {
            "chunks": col.config.chunks if col.config.chunks is not None else default_chunks,
            "blocks": col.config.blocks if col.config.blocks is not None else default_blocks,
        }
        cparams = col.config.cparams if col.config.cparams is not None else self._table_cparams
        dparams = col.config.dparams if col.config.dparams is not None else self._table_dparams
        if cparams is not None:
            result["cparams"] = cparams
        if dparams is not None:
            result["dparams"] = dparams
        return result

    @staticmethod
    def _flatten_nested_dict(d: dict, prefix: str = "") -> dict:
        """Recursively flatten a nested dict into a dotted-key flat dict.

        Works for both single-row dicts ``{field: value}`` and column-batch
        dicts ``{field: array}``.  Leaves non-dict values unchanged.

        Example::

            {"trip": {"begin": {"lon": 1.0}}} -> {"trip.begin.lon": 1.0}
        """
        result = {}
        for k, v in d.items():
            full_key = join_field_path((*split_field_path(prefix), k)) if prefix else join_field_path((k,))
            if isinstance(v, dict):
                result.update(CTable._flatten_nested_dict(v, full_key))
            else:
                result[full_key] = v
        return result

    def _normalize_row_input(self, data: Any) -> dict[str, Any]:
        """Normalize a row input to a ``{col_name: value}`` dict.

        Accepted shapes:
        - list / tuple  → positional, zipped with stored column names (computed columns skipped)
        - dict          → used as-is (nested dicts are flattened to dotted keys)
        - dataclass     → ``dataclasses.asdict`` (nested fields flattened)
        - np.void / structured scalar → field-name access
        """
        stored = self._append_input_col_names
        if isinstance(data, dict):
            if any(isinstance(v, dict) for v in data.values()):
                return self._flatten_nested_dict(data)
            return data
        if isinstance(data, (list, tuple)):
            return dict(zip(stored, data, strict=False))
        if dataclasses.is_dataclass(data) and not isinstance(data, type):
            d = dataclasses.asdict(data)
            if any(isinstance(v, dict) for v in d.values()):
                return self._flatten_nested_dict(d)
            return d
        if isinstance(data, (np.void, np.record)):
            return {name: data[name] for name in stored}
        # Fallback: try positional indexing
        return {name: data[i] for i, name in enumerate(stored)}

    def _coerce_row_to_storage(self, row: dict[str, Any]) -> dict[str, Any]:
        """Coerce each value in *row* to the column's storage representation."""
        result = {}
        for col in self._schema.columns:
            val = row[col.name]
            if self._is_list_column(col):
                result[col.name] = coerce_list_cell(col.spec, val)
            elif self._is_varlen_scalar_column(col):
                # Coercion is handled inside _ScalarVarLenArray.append.
                result[col.name] = val
            elif self._is_dictionary_column(col):
                # Pass str/None through; DictionaryColumn.__setitem__ encodes.
                result[col.name] = val
            elif self._is_ndarray_column(col):
                result[col.name] = self._coerce_ndarray_value(col.name, col.spec, val)
            elif isinstance(col.spec, timestamp):
                if val is None:
                    result[col.name] = col.spec.null_value
                elif isinstance(val, (np.datetime64, str)) or hasattr(val, "isoformat"):
                    result[col.name] = (
                        np.datetime64(val).astype(f"datetime64[{col.spec.unit}]").astype(np.int64).item()
                    )
                else:
                    result[col.name] = np.array(val, dtype=col.dtype).item()
            else:
                result[col.name] = np.array(val, dtype=col.dtype).item()
        return result

    def _open_column_from_storage(self, storage: TableStorage, name: str):
        """Open one stored column from *storage*."""
        cc = self._schema.columns_by_name[name]
        if self._is_list_column(cc):
            return storage.open_list_column(name)
        if self._is_varlen_scalar_column(cc):
            return storage.open_varlen_scalar_column(name, cc.spec)
        if self._is_dictionary_column(cc):
            return storage.open_dictionary_column(name, cc.spec)
        return storage.open_column(name)

    def _resolve_last_pos(self) -> int:
        """Return the physical index of the next write slot.

        Returns the cached ``_last_pos`` when available.  After a deletion
        ``_last_pos`` is ``None``; this method then walks chunk metadata of
        ``_valid_rows`` from the end (no full decompression) to find the last
        ``True`` position, caches the result, and returns it.
        """
        if self._last_pos is not None:
            return self._last_pos

        arr = self._valid_rows
        chunk_size = arr.chunks[0]
        last_true_pos = -1

        for info in reversed(list(arr.iterchunks_info())):
            actual_size = min(chunk_size, arr.shape[0] - info.nchunk * chunk_size)
            chunk_start = info.nchunk * chunk_size

            if info.special == blosc2.SpecialValue.ZERO:
                continue
            if info.special == blosc2.SpecialValue.VALUE:
                val = np.frombuffer(info.repeated_value, dtype=arr.dtype)[0]
                if not val:
                    continue
                last_true_pos = chunk_start + actual_size - 1
                break

            chunk_data = arr[chunk_start : chunk_start + actual_size]
            nonzero = np.flatnonzero(chunk_data)
            if len(nonzero) == 0:
                continue
            last_true_pos = chunk_start + int(nonzero[-1])
            break

        self._last_pos = last_true_pos + 1
        return self._last_pos

    def trim_capacity(self) -> None:
        """Shrink fixed-width physical storage to the last live row position.

        This removes unused append capacity while preserving holes left by deletes
        before the last live row.  List and variable-length scalar columns already
        grow to their logical length and are left untouched.
        """
        if self._read_only:
            raise ValueError("Table is read-only (opened with mode='r').")
        if self.base is not None:
            raise ValueError("Cannot trim capacity of a view.")

        target = self._resolve_last_pos()
        if target <= 0 or target >= len(self._valid_rows):
            return

        for name, col_arr in self._cols.items():
            cc = self._schema.columns_by_name[name]
            if self._is_list_column(cc) or self._is_varlen_scalar_column(cc):
                continue
            if self._is_dictionary_column(cc):
                col_arr.resize((target,))
                continue
            col_arr.resize(self._column_physical_shape(cc, target))
        self._valid_rows.resize((target,))
        self._last_pos = target

    def _grow(self) -> None:
        """Grow scalar-column capacity and the valid_rows mask by one table chunk."""
        c = len(self._valid_rows)
        growth_rows = min(c, _MAX_GROWTH_ROWS)
        new_capacity = c + growth_rows
        for name, col_arr in self._cols.items():
            cc = self._schema.columns_by_name[name]
            if self._is_list_column(cc) or self._is_varlen_scalar_column(cc):
                continue
            if self._is_dictionary_column(cc):
                col_arr.resize((new_capacity,))
                continue
            col_arr.resize(self._column_physical_shape(cc, new_capacity))
        self._valid_rows.resize((new_capacity,))

    # ------------------------------------------------------------------
    # Display
    # ------------------------------------------------------------------

    def _display_positions(self, display_rows: int | None = None):
        nrows = self._n_rows
        display_rows = _CTABLE_PRINT_OPTIONS["display_rows"] if display_rows is None else display_rows
        if display_rows == 0:
            return np.empty(0, dtype=np.intp), np.empty(0, dtype=np.intp), nrows
        _slp = getattr(self, "_cached_live_positions", None)
        if _slp is not None and self.base is not None:
            all_pos = _slp
        else:
            valid_np = self._valid_rows[:]
            all_pos = np.where(valid_np)[0]
        if display_rows < 0 or nrows <= display_rows:  # -1 (or any negative) shows all rows
            return all_pos, np.array([], dtype=all_pos.dtype), 0

        preview_rows = min(10, display_rows)
        head_rows = (preview_rows + 1) // 2
        tail_rows = preview_rows // 2
        hidden = max(0, nrows - head_rows - tail_rows)
        tail_pos = all_pos[-tail_rows:] if tail_rows else np.array([], dtype=all_pos.dtype)
        return all_pos[:head_rows], tail_pos, hidden

    def _display_widths(self, col_names: list[str] | None = None) -> dict[str, int]:
        widths: dict[str, int] = {}
        col_names = self.col_names if col_names is None else col_names
        single_col = len(col_names) == 1
        for name in col_names:
            if name == "...":
                widths[name] = 3
                continue
            spec = self._schema.columns_by_name.get(name)
            dtype_label = self._dtype_info_label(self._col_dtype(name), spec.spec if spec else None)
            widths[name] = max(self._col_widths[name], len(dtype_label))
            if single_col:
                widths[name] = max(widths[name], 80)
        return widths

    def _display_columns(
        self, *, display_index: bool = False, index_width: int = 0, max_width: int | None = None
    ) -> tuple[list[str], int]:
        """Return width-friendly display columns and hidden count.

        *max_width* is the character budget for column fitting: ``None`` or a
        negative value shows all columns (no truncation); a positive int caps it.
        """
        col_names = list(self.col_names)
        if max_width is None or max_width < 0:  # unlimited: show every column
            return col_names, 0
        widths = self._display_widths(col_names)
        widths["..."] = 3
        total_width = sum(widths[n] + 2 for n in col_names) + 2 * max(0, len(col_names) - 1)
        if display_index:
            total_width += index_width + 2 + 2
        term_width = max_width
        if total_width <= term_width or len(col_names) <= 2:
            return col_names, 0

        selected: list[str] = []
        left = 0
        right = len(col_names) - 1
        used = index_width + 2 + 2 if display_index else 0

        def extra_width(name: str, n_existing: int) -> int:
            return widths[name] + 2 + (2 if n_existing else 0)

        # Account for an ellipsis column between left and right blocks.
        used += widths["..."] + 2
        while left <= right:
            left_name = col_names[left]
            need = extra_width(left_name, len(selected) + 1)
            if used + need > term_width:
                break
            selected.append(left_name)
            used += need
            left += 1
            if left > right:
                break

            right_name = col_names[right]
            need = extra_width(right_name, len(selected) + 1)
            if used + need > term_width:
                break
            selected.append(right_name)
            used += need
            right -= 1

        left_cols = [n for n in col_names if n in selected and col_names.index(n) < left]
        right_cols = [n for n in col_names if n in selected and col_names.index(n) > right]
        display_cols = left_cols + ["..."] + right_cols
        hidden = len(col_names) - len(left_cols) - len(right_cols)
        return display_cols, hidden

    @staticmethod
    def _cell_text(value, float_precision: int | None = None) -> str:
        if isinstance(value, np.datetime64):
            s = str(value).replace("T", " ")
            if s.endswith(".000"):
                s = s[:-4]
            return s
        if isinstance(value, np.ndarray):
            if value.ndim == 1 and value.size <= 6:
                return np.array2string(value, separator=", ", max_line_width=10_000)
            return f"ndarray(shape={value.shape}, dtype={value.dtype})"
        if isinstance(value, (float, np.floating)):
            precision = (
                _CTABLE_PRINT_OPTIONS["display_precision"] if float_precision is None else float_precision
            )
            if _CTABLE_PRINT_OPTIONS["fancy"]:
                return np.format_float_positional(float(value), precision=precision, trim="-")
            return f"{float(value):.{precision}f}"
        return str(value)

    @staticmethod
    def _format_cell(value, width: int, float_precision: int | None = None) -> str:
        s = CTable._cell_text(value, float_precision)
        if len(s) > width:
            s = s[: width - 1] + "…"
        if _CTABLE_PRINT_OPTIONS["fancy"]:
            return f" {s:>{width}} "
        return f"{s:>{width}}"

    @staticmethod
    def _format_index_cell(value, width: int) -> str:
        s = "" if value is None else str(value)
        if len(s) > width:
            s = s[: width - 1] + "…"
        if _CTABLE_PRINT_OPTIONS["fancy"]:
            return f" {s:<{width}} "
        return f"{s:<{width}}"

    @staticmethod
    def _display_index_width(nrows: int, hidden: int, index_name: str) -> int:
        width = max(len(index_name), len(str(max(nrows - 1, 0))))
        if hidden > 0:
            width = max(width, 3)
        return width

    def _format_display_row(
        self,
        values: dict,
        widths: dict[str, int],
        col_names: list[str],
        float_precisions: dict[str, int] | None = None,
    ) -> str:
        float_precisions = {} if float_precisions is None else float_precisions
        return "  ".join(self._format_cell(values[n], widths[n], float_precisions.get(n)) for n in col_names)

    def _format_display_row_with_index(
        self,
        values: dict,
        widths: dict[str, int],
        col_names: list[str],
        index_value,
        index_width: int,
        float_precisions: dict[str, int] | None = None,
    ) -> str:
        return (
            self._format_index_cell(index_value, index_width)
            + "  "
            + self._format_display_row(values, widths, col_names, float_precisions)
        )

    def _prewarm_display_cache(self, display_cols: list[str], head_pos, tail_pos) -> None:
        """Pre-populate the render cache with one combined head+tail gather/col.

        Each storage column is then sparse-read a single time (for head ∪ tail)
        instead of once per slice; the values are split back so the
        ``(col, id(head_pos))`` and ``(col, id(tail_pos))`` lookups made by
        precision detection, width sizing and row rendering all hit the cache.
        Only pays off when both slices are non-empty.
        """
        cache = getattr(self, "_display_fetch_cache", None)
        if cache is None or len(head_pos) == 0 or len(tail_pos) == 0:
            return
        real_cols = [n for n in display_cols if n != "..." and (n in self._cols or n in self._computed_cols)]
        if not real_cols:
            return
        nh = len(head_pos)
        combined = np.concatenate([head_pos, tail_pos])
        for name in real_cols:
            vals = self._fetch_col_at_positions_uncached(name, combined)
            cache[(name, id(head_pos))] = vals[:nh]
            cache[(name, id(tail_pos))] = vals[nh:]

    def _rows_to_dicts(self, positions, col_names: list[str] | None = None) -> list[dict]:
        if len(positions) == 0:
            return []
        col_names = self.col_names if col_names is None else col_names
        real_cols = [n for n in col_names if n != "..."]
        col_data = {n: self._fetch_col_at_positions(n, positions) for n in real_cols}
        rows = []
        for i in range(len(positions)):
            row = {}
            for n in col_names:
                # Keep NumPy scalar types for display so their compact string
                # formatting is preserved (notably float32, e.g. 224.97
                # instead of Python float's 224.97000122070312).
                row[n] = "..." if n == "..." else col_data[n][i]
            rows.append(row)
        return rows

    def _display_separator(
        self,
        widths: dict[str, int],
        display_cols: list[str],
        display_index: bool,
        index_width: int,
        fancy: bool,
    ) -> str | None:
        if not fancy:
            return None
        sep_parts = ["─" * (widths[n] + 2) for n in display_cols]
        if display_index:
            sep_parts.insert(0, "─" * (index_width + 2))
        return "  ".join(sep_parts)

    def _display_dtype_row(self, display_cols: list[str]) -> dict:
        dtype_row = {}
        for n in display_cols:
            if n == "...":
                dtype_row[n] = "..."
            else:
                dtype_row[n] = self._dtype_info_label(
                    self._col_dtype(n),
                    self._schema.columns_by_name[n].spec if n in self._schema.columns_by_name else None,
                )
        return dtype_row

    def _compact_float_precisions(self, display_cols: list[str], head_pos, tail_pos) -> dict[str, int]:
        default_precision = _CTABLE_PRINT_OPTIONS["display_precision"]
        precisions: dict[str, int] = {}
        for n in display_cols:
            finite_float_seen = False
            integer_valued = True
            for positions in (head_pos, tail_pos):
                for row in self._rows_to_dicts(positions, [n]):
                    value = row[n]
                    if not isinstance(value, (float, np.floating)):
                        continue
                    value = float(value)
                    if not np.isfinite(value):
                        continue
                    finite_float_seen = True
                    if not value.is_integer():
                        integer_valued = False
                        break
                if not integer_valued:
                    break
            if finite_float_seen and integer_valued:
                precisions[n] = 1
            else:
                precisions[n] = default_precision
        return precisions

    def _compact_display_widths(
        self,
        display_cols: list[str],
        head_pos,
        tail_pos,
        hidden: int,
        float_precisions: dict[str, int],
    ) -> dict[str, int]:
        widths = {n: len(n) for n in display_cols}
        if hidden > 0:
            for n in display_cols:
                widths[n] = max(widths[n], 3)
        for positions in (head_pos, tail_pos):
            for row in self._rows_to_dicts(positions, display_cols):
                for n, value in row.items():
                    widths[n] = max(widths[n], len(self._cell_text(value, float_precisions.get(n))))
        return widths

    @staticmethod
    def _display_footer(nrows: int, ncols: int, hidden: int, hidden_cols: int, fancy: bool) -> list[str]:
        if not fancy:
            return ["", f"[{nrows} rows x {ncols} columns]"]
        footer = f"{nrows:,} rows × {ncols} columns"
        notes = []
        if hidden > 0:
            notes.append(f"{hidden:,} rows hidden")
        if hidden_cols > 0:
            notes.append(f"{hidden_cols:,} columns hidden")
        if notes:
            footer += f"  ({', '.join(notes)})"
        return [footer]

    def _display_lines_with_index(
        self,
        *,
        display_cols: list[str],
        widths: dict[str, int],
        index_name: str,
        index_width: int,
        head_pos,
        tail_pos,
        hidden: int,
        sep: str | None,
        fancy: bool,
        float_precisions: dict[str, int] | None = None,
    ) -> list[str]:
        header_row = {n: n for n in display_cols}
        lines = [
            self._format_display_row_with_index(
                header_row, widths, display_cols, index_name, index_width, float_precisions
            )
        ]
        if fancy:
            dtype_row = self._display_dtype_row(display_cols)
            lines.append(
                self._format_display_row_with_index(
                    dtype_row, widths, display_cols, None, index_width, float_precisions
                )
            )
        if sep is not None:
            lines.append(sep)
        lines.extend(
            self._format_display_row_with_index(row, widths, display_cols, i, index_width, float_precisions)
            for i, row in enumerate(self._rows_to_dicts(head_pos, display_cols))
        )
        if hidden > 0:
            lines.append(
                self._format_display_row_with_index(
                    dict.fromkeys(display_cols, "..."),
                    widths,
                    display_cols,
                    "...",
                    index_width,
                    float_precisions,
                )
            )
        tail_start = self._n_rows - len(tail_pos)
        lines.extend(
            self._format_display_row_with_index(
                row, widths, display_cols, tail_start + i, index_width, float_precisions
            )
            for i, row in enumerate(self._rows_to_dicts(tail_pos, display_cols))
        )
        return lines

    def _display_lines_without_index(
        self,
        *,
        display_cols: list[str],
        widths: dict[str, int],
        head_pos,
        tail_pos,
        hidden: int,
        sep: str | None,
        fancy: bool,
        float_precisions: dict[str, int] | None = None,
    ) -> list[str]:
        header_row = {n: n for n in display_cols}
        lines = [self._format_display_row(header_row, widths, display_cols, float_precisions)]
        if fancy:
            lines.append(
                self._format_display_row(
                    self._display_dtype_row(display_cols), widths, display_cols, float_precisions
                )
            )
        if sep is not None:
            lines.append(sep)
        lines.extend(
            self._format_display_row(row, widths, display_cols, float_precisions)
            for row in self._rows_to_dicts(head_pos, display_cols)
        )
        if hidden > 0:
            lines.append(
                self._format_display_row(
                    dict.fromkeys(display_cols, "..."), widths, display_cols, float_precisions
                )
            )
        lines.extend(
            self._format_display_row(row, widths, display_cols, float_precisions)
            for row in self._rows_to_dicts(tail_pos, display_cols)
        )
        return lines

    def to_string(
        self,
        *,
        max_rows: int | None = None,
        max_width: int | None = None,
        show_dimensions: bool | str = False,
        display_index: bool | None = None,
        index_name: str = "",
    ) -> str:
        """Return a tabular string representation of the table.

        By default (``max_rows=None``, ``max_width=None``) this renders the
        *whole* table — every row and every column — like ``pandas``'
        ``DataFrame.to_string()``.  This is independent of the global
        :func:`blosc2.set_printoptions`; those only affect the truncated
        ``str``/``repr``/``print`` view.

        Parameters
        ----------
        max_rows:
            Maximum number of rows before truncating to a compact head/tail
            view.  ``None`` (default) shows all rows; ``-1`` also means all,
            ``0`` shows none, a positive int caps it.
        max_width:
            Character budget for column fitting.  ``None`` (default) or ``-1``
            shows all columns; a positive int truncates the middle ones with
            ``...`` to fit.
        show_dimensions:
            Whether to append a ``[N rows x M columns]`` footer.  ``False``
            (default) omits it, matching ``pandas``' ``to_string()``; ``True``
            always shows it; ``"truncate"`` shows it only when the view is
            truncated (the behaviour of ``str``/``repr``).
        display_index:
            Whether to include a pandas-like logical row index column.  If
            ``None`` (default), use the global value configured with
            :func:`blosc2.set_printoptions`.
        index_name:
            Optional label for the displayed index column.
        """
        if display_index is None:
            display_index = _CTABLE_PRINT_OPTIONS["display_index"]
        if not isinstance(display_index, bool):
            raise TypeError("display_index must be a bool or None")
        if not isinstance(index_name, str):
            raise TypeError("index_name must be a str")

        nrows = self._n_rows
        ncols = len(self.col_names)
        rows_arg = -1 if max_rows is None else max_rows  # None ⇒ all rows
        head_pos, tail_pos, hidden = self._display_positions(rows_arg)
        # Memoise per-column sparse gathers for the duration of this render so
        # the repeated (column, head_pos/tail_pos) lookups across precision,
        # width and row formatting only touch storage once.  head_pos/tail_pos
        # stay referenced below, so keying the cache on their id() is safe.
        self._display_fetch_cache = {}
        try:
            return self._to_string_body(
                display_index,
                index_name,
                nrows,
                ncols,
                head_pos,
                tail_pos,
                hidden,
                max_width,
                show_dimensions,
            )
        finally:
            self._display_fetch_cache = None

    def _to_string_body(
        self,
        display_index,
        index_name,
        nrows,
        ncols,
        head_pos,
        tail_pos,
        hidden,
        max_width=None,
        show_dimensions=False,
    ) -> str:
        index_width = self._display_index_width(nrows, hidden, index_name) if display_index else 0
        display_cols, hidden_cols = self._display_columns(
            display_index=display_index, index_width=index_width, max_width=max_width
        )
        # Warm the fetch cache with a single combined head+tail gather per column.
        # head_pos/tail_pos land in different blocks, but folding them into one
        # sparse read still halves the per-call gather overhead vs reading head
        # and tail separately (and every downstream consumer hits the cache).
        self._prewarm_display_cache(display_cols, head_pos, tail_pos)
        fancy = _CTABLE_PRINT_OPTIONS["fancy"]
        float_precisions = {} if fancy else self._compact_float_precisions(display_cols, head_pos, tail_pos)
        widths = (
            self._display_widths(display_cols)
            if fancy
            else self._compact_display_widths(display_cols, head_pos, tail_pos, hidden, float_precisions)
        )
        sep = self._display_separator(widths, display_cols, display_index, index_width, fancy)

        if display_index:
            lines = self._display_lines_with_index(
                display_cols=display_cols,
                widths=widths,
                index_name=index_name,
                index_width=index_width,
                head_pos=head_pos,
                tail_pos=tail_pos,
                hidden=hidden,
                sep=sep,
                fancy=fancy,
                float_precisions=float_precisions,
            )
        else:
            lines = self._display_lines_without_index(
                display_cols=display_cols,
                widths=widths,
                head_pos=head_pos,
                tail_pos=tail_pos,
                hidden=hidden,
                sep=sep,
                fancy=fancy,
                float_precisions=float_precisions,
            )
        if sep is not None:
            lines.append(sep)
        # pandas convention: to_string() omits the dimensions footer
        # (show_dimensions=False); str/repr show it only when truncated.
        truncated = hidden > 0 or hidden_cols > 0
        if show_dimensions is True or (show_dimensions == "truncate" and truncated):
            lines.extend(self._display_footer(nrows, ncols, hidden, hidden_cols, fancy))
        return "\n".join(lines)

    def __str__(self) -> str:
        """Pandas-style tabular display, truncated per :func:`blosc2.set_printoptions`."""
        opts = _CTABLE_PRINT_OPTIONS
        width = opts["display_width"]
        if width is None:  # auto: fit to the current terminal
            width = shutil.get_terminal_size((120, 20)).columns
        return self.to_string(max_rows=opts["display_rows"], max_width=width, show_dimensions="truncate")

    def __repr__(self) -> str:
        """Same truncated table as ``str`` (pandas/polars convention).

        The compact ``CTable<cols>(N rows, …)`` summary is available via
        :attr:`info`.
        """
        return self.__str__()

    def __len__(self):
        """Return the number of live (non-deleted) rows."""
        return self._n_rows

    def __iter__(self):
        """Iterate over live rows in insertion order, yielding namedtuple-like row objects."""
        for i in range(self.nrows):
            yield self._materialize_row(i)

    def _row_namedtuple_type(self):
        visible = tuple(self.col_names)
        if getattr(self, "_row_namedtuple_type_cache_cols", None) != visible:
            self._row_namedtuple_type_cache = _make_namedtuple_row_type(visible)
            self._row_namedtuple_type_cache_cols = visible
        return self._row_namedtuple_type_cache

    def _row_namedtuple_type_for_fields(self, fields: tuple[str, ...]):
        cache = getattr(self, "_row_namedtuple_type_cache_by_fields", None)
        if cache is None:
            cache = {}
            self._row_namedtuple_type_cache_by_fields = cache
        row_type = cache.get(fields)
        if row_type is None:
            row_type = _make_namedtuple_row_type(fields)
            cache[fields] = row_type
        return row_type

    @staticmethod
    def _normalize_scalar_value(value):
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, np.ndarray) and value.ndim == 0:
            return value.item()
        return value

    def _physical_row_value(self, col_name: str, pos: int):
        cc = self._computed_cols.get(col_name)
        if cc is not None:
            return self._normalize_scalar_value(np.asarray(self._build_computed_lazy(cc)[pos]).ravel()[0])
        value = self._normalize_scalar_value(self._cols[col_name][pos])
        spec = self._schema.columns_by_name[col_name].spec
        if isinstance(spec, timestamp):
            return np.datetime64(int(value), spec.unit)
        return value

    def _materialize_row(self, index: int):
        n_rows = self.nrows
        if index < 0:
            index += n_rows
        if not (0 <= index < n_rows):
            raise IndexError(f"row index {index} is out of bounds for table with {n_rows} rows")
        _slp = getattr(self, "_cached_live_positions", None)
        if _slp is not None and self.base is not None:
            pos = int(_slp[index])
        else:
            pos = _find_physical_index(self._valid_rows, index)

        nested_meta = self._schema.metadata.get("nested") if self._schema.metadata else None
        reconstruct = isinstance(nested_meta, dict) and bool(nested_meta.get("reconstruct_rows", False))
        if not reconstruct:
            row_type = self._row_namedtuple_type()
            return row_type(*(self._physical_row_value(name, int(pos)) for name in self.col_names))

        row_dict: dict[str, Any] = {}
        for name in self.col_names:
            value = self._physical_row_value(name, int(pos))
            parts = split_field_path(name)
            if len(parts) <= 1:
                row_dict[name] = value
                continue
            node = row_dict
            for part in parts[:-1]:
                child = node.get(part)
                if not isinstance(child, dict):
                    child = {}
                    node[part] = child
                node = child
            node[parts[-1]] = value

        fields = tuple(row_dict.keys())
        row_type = self._row_namedtuple_type_for_fields(fields)
        return row_type(*(row_dict[f] for f in fields))

    def iter_sorted(
        self,
        cols: str | list[str],
        ascending: bool | list[bool] = True,
        *,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
        batch_size: int = 4096,
    ):
        """Iterate rows in sorted order without materializing a full copy.

        Uses a FULL index when available (no sort needed); otherwise falls
        back to ``np.lexsort`` on live physical positions.  Yields namedtuple-like
        row objects in the same way as ``__iter__``.

        The sorted positions array is stored as a compressed ``blosc2.NDArray``
        to keep RAM usage low for large tables.  ``batch_size`` positions are
        decompressed at a time during iteration.

        Parameters
        ----------
        cols:
            Column name or list of column names to sort by.
        ascending:
            Sort direction.  A single bool applies to all keys; a list must
            have the same length as *cols*.
        start, stop, step:
            Optional slice applied to the sorted sequence before iteration.
            E.g. ``stop=10`` yields only the top-10 rows; ``step=2`` yields
            every other row in sorted order.
        batch_size:
            Number of positions decompressed per iteration step.  Larger
            values reduce decompression overhead; smaller values use less
            transient RAM.  Default is 4096.
        """
        cols, ascending = self._normalise_sort_keys(cols, ascending)

        valid_np = self._valid_rows[:]
        live_pos = np.where(valid_np)[0]
        n = len(live_pos)

        if n == 0:
            return

        sorted_pos = None
        if len(cols) == 1:
            sorted_pos = self._sorted_positions_from_full_index(cols[0], ascending[0])
            if sorted_pos is not None and len(sorted_pos) != n:
                sorted_pos = None

        if sorted_pos is None:
            order = np.lexsort(self._build_lex_keys(cols, ascending, live_pos, n))
            sorted_pos = live_pos[order]

        if start is not None or stop is not None or step is not None:
            sorted_pos = sorted_pos[start:stop:step]

        # Compress positions into an NDArray to reduce RAM usage for large tables.
        # The uncompressed numpy array is released immediately after.
        sorted_pos_nd = blosc2.asarray(np.asarray(sorted_pos, dtype=np.int64))
        del sorted_pos

        # physical → logical index mapping
        phys_to_logical = np.empty(valid_np.shape[0], dtype=np.intp)
        phys_to_logical[live_pos] = np.arange(n, dtype=np.intp)

        total = len(sorted_pos_nd)
        for i in range(0, total, batch_size):
            chunk = sorted_pos_nd[i : i + batch_size]
            for phys in chunk:
                yield self._materialize_row(int(phys_to_logical[phys]))

    # ------------------------------------------------------------------
    # Open existing table (classmethod)
    # ------------------------------------------------------------------

    @classmethod
    def _open_from_existing_filestore(cls, urlpath: str, *, mode: str, store: blosc2.TreeStore) -> CTable:
        """Open a root CTable reusing an already-opened TreeStore."""
        storage = FileTableStorage(urlpath, mode, store=store)
        return cls._open_from_storage(storage)

    @classmethod
    def open(cls, urlpath: str, *, mode: str = "r", mmap_mode: str | None = None) -> CTable:
        """Open a persistent CTable from *urlpath*.

        Parameters
        ----------
        urlpath:
            Path to the table root directory (created by passing ``urlpath``
            to :class:`CTable`).
        mode:
            ``'r'`` (default) — read-only.
            ``'a'`` — read/write.
        mmap_mode:
            ``'r'`` to memory-map the backing store instead of using regular
            file I/O (requires ``mode='r'``).  For a ``.b2z`` container this
            maps the single file once and reads every member — columns and
            index sidecars alike — in place at its offset, with mapped pages
            shared across reader processes.

        Raises
        ------
        FileNotFoundError
            If *urlpath* does not contain a CTable.
        ValueError
            If the metadata at *urlpath* does not identify a CTable, or if
            ``mmap_mode`` is used with a writable ``mode``.
        """
        storage = FileTableStorage(urlpath, mode, mmap_mode=mmap_mode)
        if not storage.table_exists():
            raise FileNotFoundError(f"No CTable found at {urlpath!r}")
        return cls._open_from_storage(storage)

    def to_b2z(self, urlpath: str, *, overwrite: bool = False, compact: bool = False) -> str:
        """Write this table to a compact ``.b2z`` container.

        ``.b2z`` is the compact zip-backed CTable format.  For persistent,
        non-view directory-backed tables and ``compact=False``, this uses a
        fast physical-pack path: the backing :class:`TreeStore` directory is
        zipped with already-compressed leaves stored as-is. This preserves the
        physical layout, including deleted rows and spare capacity, and does
        not recompress columns.  A ``.b2d`` suffix is recommended for
        directory-backed stores, but not required.

        For in-memory tables, views, existing ``.b2z`` tables, or
        ``compact=True``, this falls back to the logical :meth:`save` path,
        materializing only visible/live rows into a new ``.b2z`` store.

        Examples
        --------
        Fast-pack an existing directory-backed table into a compact zip store::

            table = blosc2.CTable.open("data.b2d", mode="r")
            table.to_b2z("data.b2z", overwrite=True)
            table.close()

        Materialize a filtered view into a new compact store::

            view = table.where(table["score"] > 10)
            view.to_b2z("high-score.b2z", overwrite=True)

        Force a logical compacted copy, even for a persistent ``.b2d`` table::

            table.to_b2z("data-compact.b2z", overwrite=True, compact=True)
        """
        if not str(urlpath).endswith(".b2z"):
            raise ValueError("urlpath must have a .b2z extension")

        storage = getattr(self, "_storage", None)
        can_physical_pack = (
            not compact
            and self.base is None
            and isinstance(storage, FileTableStorage)
            and not str(storage._root).endswith(".b2z")
        )
        if can_physical_pack:
            self._flush_varlen_columns()
            store = blosc2.TreeStore(storage._root, mode="r")
            try:
                return store.to_b2z(filename=urlpath, overwrite=overwrite)
            finally:
                store.close()

        if self.base is not None:
            materialized = self.copy(compact=True)
            materialized.save(urlpath, overwrite=overwrite)
        else:
            self.save(urlpath, overwrite=overwrite)
        return os.path.abspath(urlpath)

    def to_b2d(self, urlpath: str, *, overwrite: bool = False, compact: bool = False) -> str:
        """Write this table to a directory-backed store.

        Directory-backed CTable stores may use any path that does not end in
        ``.b2z``; using a ``.b2d`` suffix is recommended for clarity.  For
        persistent, non-view ``.b2z`` tables opened read-only and
        ``compact=False``, this uses a fast physical-unpack path: the zip
        members are extracted as already-compressed leaves. This preserves the
        physical layout, including deleted rows and spare capacity, and does
        not recompress columns.

        For in-memory tables, views, writable ``.b2z`` tables, existing
        directory-backed tables, or ``compact=True``, this falls back to the
        logical :meth:`save` path, materializing only visible/live rows into a
        new directory-backed store.

        Examples
        --------
        Fast-unpack an existing compact zip store into a directory-backed table::

            table = blosc2.CTable.open("data.b2z", mode="r")
            table.to_b2d("data.b2d", overwrite=True)
            table.close()

        Materialize a filtered view into a directory-backed store::

            view = table.where(table["score"] > 10)
            view.to_b2d("high-score.b2d", overwrite=True)

        Force a logical compacted copy, even for a persistent ``.b2z`` table::

            table.to_b2d("data-compact.b2d", overwrite=True, compact=True)
        """
        urlpath = os.fspath(urlpath)
        storage = getattr(self, "_storage", None)
        can_physical_unpack = (
            not compact
            and self.base is None
            and isinstance(storage, FileTableStorage)
            and str(storage._root).endswith(".b2z")
            and storage.open_mode() == "r"
        )
        if can_physical_unpack:
            store = blosc2.TreeStore(storage._root, mode="r")
            try:
                return store.to_b2d(urlpath, overwrite=overwrite)
            finally:
                store.close()

        if self.base is not None:
            materialized = self.copy(compact=True)
            materialized.save(urlpath, overwrite=overwrite)
        else:
            self.save(urlpath, overwrite=overwrite)
        return os.path.abspath(urlpath)

    def to_cframe(self) -> bytes:
        """Serialize this table to a bytes buffer (a CFrame).

        This is the Blosc2-bytes counterpart of :meth:`to_b2z`, mirroring
        :meth:`blosc2.NDArray.to_cframe`.  The table is packed into an
        in-memory :class:`blosc2.EmbedStore` (one entry per column, plus the
        schema, the ``_valid_rows`` mask and any user vlmeta) and serialized to
        a single ``bytes`` object.

        Only live rows are serialized: for a view or slice the result is
        materialized first.  The result is fully self-contained and can be
        rebuilt with :func:`blosc2.ctable_from_cframe` without any temp file,
        which makes it suitable as a transport format (e.g. Caterva2 table
        slice fetch) including under Pyodide.

        Returns
        -------
        bytes
            The serialized table.

        See Also
        --------
        :func:`blosc2.ctable_from_cframe`
        """
        # Materialize live rows for views/slices; for base tables use the live
        # columns directly.  copy() is the canonical materialization path.
        if self.base is not None:
            src = self.copy(compact=True)
        else:
            src = self
            src._flush_varlen_columns()

        from blosc2.ctable_storage import (
            _DICT_SUFFIX,
            _UTF8_DATA_SUFFIX,
            _column_name_to_relpath,
        )

        estore = blosc2.EmbedStore(urlpath=None, mode="w")

        # Manifest: schema + kind/version markers, mirroring FileTableStorage.save_schema
        meta = blosc2.SChunk()
        meta.vlmeta["kind"] = "ctable"
        meta.vlmeta["version"] = 1
        meta.vlmeta["schema"] = json.dumps(src._schema_dict_with_computed())
        estore["/_meta"] = meta
        estore["/_valid_rows"] = src._valid_rows

        # User vlmeta (if any) — best-effort, mirroring FileTableStorage._open_vlmeta
        vlmeta_schunk = getattr(src._storage, "_vlmeta", None)
        if vlmeta_schunk is None:
            vlmeta_schunk = getattr(src._storage, "_vlmeta_schunk", None)
        if isinstance(vlmeta_schunk, blosc2.SChunk):
            estore["/_vlmeta"] = vlmeta_schunk

        for col in src._schema.columns:
            name = col.name
            key = f"/_cols/{_column_name_to_relpath(name)}"
            arr = src._cols[name]
            if self._is_dictionary_column(col):
                estore[key] = arr.codes
                estore[key + _DICT_SUFFIX] = arr._dict_store._backend
            elif self._is_utf8_column(col):
                estore[key] = arr.offsets
                estore[key + _UTF8_DATA_SUFFIX] = arr.data
            elif self._is_varlen_scalar_column(col):
                estore[key] = arr._backend
            else:
                # Scalar NDArray or ListArray — both serialize via to_cframe().
                estore[key] = arr

        return estore.to_cframe()

    def _save_to_storage(  # noqa: C901
        self,
        storage: TableStorage,
        *,
        chunks_override: tuple[int, ...] | None = None,
        blocks_override: tuple[int, ...] | None = None,
        cparams_override: dict[str, Any] | None = None,
    ) -> None:
        """Write all live rows and columns into *storage*.

        The caller is responsible for calling ``storage.close()`` when done.
        This method does **not** close *storage*.
        """
        self._flush_varlen_columns()

        # Collect live physical positions
        valid_np = self._valid_rows[:]
        live_pos = np.where(valid_np)[0]
        n_live = len(live_pos)
        capacity = max(n_live, 1)
        # True when all live positions are [0, 1, ..., n_live-1] — no gaps or tombstones.
        no_deletions = n_live > 0 and int(live_pos[0]) == 0 and int(live_pos[-1]) == n_live - 1

        default_chunks, default_blocks = compute_chunks_blocks((capacity,))
        # Align fixed-size scalar columns (and the _valid_rows mask) on one
        # shared grid so lazy expressions over the saved table take the
        # fast_eval path on reopen.
        shared_chunks, shared_blocks, aligned_names = self._compute_aligned_grid(
            self._schema.columns, capacity
        )

        # --- valid_rows (all True, compacted) ---
        disk_valid = storage.create_valid_rows(
            shape=(capacity,),
            chunks=shared_chunks if shared_chunks is not None else default_chunks,
            blocks=shared_blocks if shared_blocks is not None else default_blocks,
        )
        if n_live > 0:
            disk_valid[:n_live] = True

        # --- columns ---
        for col in self._schema.columns:
            name = col.name
            col_cparams = col.config.cparams if col.config.cparams is not None else self._table_cparams
            eff_cparams = cparams_override if cparams_override is not None else col_cparams
            if self._is_list_column(col):
                src_la = self._cols[name]
                if no_deletions and cparams_override is None and not src_la._pending_cells:
                    # Fast path: C-level chunk transfer, no Python decompression.
                    storage.install_list_column(name, src_la)
                else:
                    disk_col = storage.create_list_column(
                        name,
                        spec=col.spec,
                        cparams=eff_cparams,
                        dparams=col.config.dparams
                        if col.config.dparams is not None
                        else self._table_dparams,
                    )
                    if n_live > 0:
                        items = src_la[:n_live] if no_deletions else (src_la[int(pos)] for pos in live_pos)
                        disk_col.extend(items, validate=False)
                        disk_col.flush()
                continue
            if self._is_varlen_scalar_column(col):
                disk_col = storage.create_varlen_scalar_column(
                    name,
                    spec=col.spec,
                    cparams=eff_cparams,
                    dparams=col.config.dparams if col.config.dparams is not None else self._table_dparams,
                )
                if n_live > 0:
                    src_vl = self._cols[name]
                    items = src_vl[:n_live] if no_deletions else (src_vl[int(pos)] for pos in live_pos)
                    disk_col.extend(items)
                    disk_col.flush()
                continue
            if self._is_dictionary_column(col):
                src_dc = self._cols[name]
                disk_dc = storage.create_dictionary_column(
                    name,
                    spec=col.spec,
                    cparams=eff_cparams,
                    dparams=col.config.dparams if col.config.dparams is not None else self._table_dparams,
                )
                # Copy dictionary values first
                for v in src_dc.dictionary:
                    disk_dc.encode(v)
                disk_dc.flush()
                # Resize codes to the target capacity before writing — the default
                # codes_shape=(4096,) in create_dictionary_column is too small.
                if len(disk_dc.codes) < capacity:
                    disk_dc.codes.resize((capacity,))
                # Copy live codes
                if n_live > 0:
                    pos_slice = np.arange(n_live, dtype=np.int64) if no_deletions else live_pos
                    disk_dc.codes[:n_live] = src_dc.codes[pos_slice]
                continue
            shape = self._column_physical_shape(col, capacity)
            if col.name in aligned_names:
                dtype_chunks, dtype_blocks = shared_chunks, shared_blocks
            else:
                dtype_chunks, dtype_blocks = self._column_chunks_blocks(col, shape)
            col_storage = self._resolve_column_storage(col, dtype_chunks, dtype_blocks)
            src_arr = self._cols[name]
            # Fast reblock path: use NDArray.copy() for a single C-level
            # decompress+recompress pass.  Only safe when the source array has
            # no spare capacity (src shape == n_live); otherwise copy() would
            # include uninitialised rows beyond the live watermark.
            if (
                no_deletions
                and src_arr.shape[0] == n_live
                and (
                    chunks_override is not None
                    or blocks_override is not None
                    or cparams_override is not None
                )
            ):
                copy_kwargs: dict[str, Any] = {}
                if chunks_override is not None:
                    copy_kwargs["chunks"] = chunks_override
                if blocks_override is not None:
                    copy_kwargs["blocks"] = blocks_override
                if cparams_override is not None:
                    copy_kwargs["cparams"] = cparams_override
                new_arr = src_arr.copy(**copy_kwargs)
                storage.install_column(name, new_arr)
            else:
                eff_chunks = chunks_override if chunks_override is not None else col_storage["chunks"]
                if chunks_override is not None and blocks_override is None:
                    _sb = shared_blocks if shared_blocks is not None else default_blocks
                    if _sb is not None and all(b <= c for b, c in zip(_sb, chunks_override, strict=False)):
                        eff_blocks = _sb
                    else:
                        eff_blocks = None
                else:
                    eff_blocks = blocks_override if blocks_override is not None else col_storage["blocks"]
                disk_col = storage.create_column(
                    name,
                    dtype=col.dtype,
                    shape=shape,
                    chunks=eff_chunks,
                    blocks=eff_blocks,
                    cparams=cparams_override if cparams_override is not None else col_storage.get("cparams"),
                    dparams=col_storage.get("dparams"),
                )
                if n_live > 0:
                    # Slice is ~30x faster than fancy-index for sequential no-deletion access.
                    disk_col[:n_live] = src_arr[:n_live] if no_deletions else src_arr[live_pos]

        storage.save_schema(self._schema_dict_with_computed())

    def save(self, urlpath: str, *, overwrite: bool = False) -> None:
        """Persist this table to disk at *urlpath*.

        This writes a standalone copy and returns ``None``; use :meth:`copy`
        directly when the copied :class:`CTable` object is needed.

        Only live rows are written — the on-disk table is always compacted.
        A ``.b2z`` suffix selects the compact zip-backed format; any other
        suffix creates a directory-backed store.  Use a ``.b2d`` suffix for
        directory-backed stores when possible so the format is clear.

        Parameters
        ----------
        urlpath:
            Destination path.  Use a ``.b2z`` suffix for a compact zip-backed
            store; any other suffix creates a directory-backed store.  A
            ``.b2d`` suffix is recommended for directory-backed stores.
        overwrite:
            If ``False`` (default), raise :exc:`ValueError` when *urlpath*
            already exists.  Set to ``True`` to replace an existing table.

        Raises
        ------
        ValueError
            If *urlpath* already exists and ``overwrite=False``.
        """
        if self.base is not None:
            materialized = self.copy(compact=True)
            materialized.save(urlpath, overwrite=overwrite)
            return

        file_storage = FileTableStorage(urlpath, "w")
        target_path = file_storage._root
        if os.path.exists(target_path):
            if not overwrite:
                raise ValueError(f"Path {target_path!r} already exists. Use overwrite=True to replace.")
            if os.path.isdir(target_path):
                shutil.rmtree(target_path)
            else:
                os.remove(target_path)

        self._save_to_storage(file_storage)
        file_storage.close()

    @classmethod
    def _open_from_storage(cls, storage: TableStorage) -> CTable:
        """Construct a :class:`CTable` from an already-configured *storage* backend.

        The caller must have already verified that the storage target exists.
        This is the common open path shared by :meth:`open` and
        :meth:`_open_from_treestore`.
        """
        storage.check_kind()
        schema_dict = storage.load_schema()
        schema = schema_from_dict(schema_dict)
        col_names = [c["name"] for c in schema_dict["columns"]]

        obj = cls.__new__(cls)
        obj._row_type = None
        obj._validate = True
        obj._table_cparams = None
        obj._table_dparams = None
        obj._storage = storage
        obj._read_only = storage.is_read_only()
        obj._schema = schema
        obj._cols = {}
        obj._col_widths = {}
        obj.col_names = col_names
        obj.auto_compact = False
        obj._create_summary_index = schema_dict.get("create_summary_index", True)
        obj._summary_indexes_built = schema_dict.get("summary_indexes_built", False)
        obj.base = None

        obj._valid_rows = storage.open_valid_rows()
        obj._cols = _LazyColumnDict(obj, storage, col_names)
        for name in col_names:
            cc = schema.columns_by_name[name]
            obj._col_widths[name] = max(len(name), cc.display_width)

        obj._n_rows = None
        # Restore cached row count from saved metadata so that
        # where() can skip the _valid_rows intersection for all-valid tables.
        if "n_rows" in schema_dict:
            obj._n_rows_cached = schema_dict["n_rows"]
        obj._last_pos = None
        obj._computed_cols = {}
        obj._materialized_cols = {}
        obj._expr_index_arrays = {}
        obj._load_computed_cols_from_schema(schema_dict)
        obj._load_materialized_cols_from_schema(schema_dict)
        return obj

    def _save_to_treestore(self, store: blosc2.TreeStore, full_key: str) -> None:
        """Save this CTable inline into *store* under *full_key*.

        *full_key* must be the absolute (fully-translated) key within the
        backing DictStore (not a subtree-relative key).
        Internal use only — called by :class:`blosc2.TreeStore`.
        """
        if self.base is not None:
            materialized = self.copy(compact=True)
            materialized._save_to_treestore(store, full_key)
            return
        storage = TreeStoreTableStorage(store, full_key, mode="a", owns_store=False)
        self._save_to_storage(storage)
        # storage is non-owning; outer store handles persistence

    @classmethod
    def _open_from_treestore(cls, store: blosc2.TreeStore, full_key: str) -> CTable:
        """Open an inline CTable from *store* at *full_key*.

        *full_key* must be the absolute key within the backing DictStore.
        Internal use only — called by :class:`blosc2.TreeStore`.
        """
        storage = TreeStoreTableStorage(store, full_key, mode=store.mode, owns_store=False)
        if not storage.table_exists():
            raise FileNotFoundError(f"No inline CTable found at key {full_key!r} in {store.localpath!r}")
        return cls._open_from_storage(storage)

    @classmethod
    def load(cls, urlpath: str) -> CTable:  # noqa: C901
        """Load a persistent table from *urlpath* into RAM.

        The schema is read from the table's metadata — the original Python
        dataclass is not required.  The returned table is fully in-memory and
        read/write.

        Parameters
        ----------
        urlpath:
            Path to the table root directory.

        Raises
        ------
        FileNotFoundError
            If *urlpath* does not contain a CTable.
        ValueError
            If the metadata at *urlpath* does not identify a CTable.
        """
        file_storage = FileTableStorage(urlpath, "r")
        if not file_storage.table_exists():
            raise FileNotFoundError(f"No CTable found at {urlpath!r}")
        file_storage.check_kind()
        schema_dict = file_storage.load_schema()
        schema = schema_from_dict(schema_dict)
        col_names = [c["name"] for c in schema_dict["columns"]]

        disk_valid = file_storage.open_valid_rows()
        disk_cols = {}
        for col in schema.columns:
            if cls._is_list_column(col):
                disk_cols[col.name] = file_storage.open_list_column(col.name)
            elif cls._is_varlen_scalar_column(col):
                disk_cols[col.name] = file_storage.open_varlen_scalar_column(col.name, col.spec)
            elif cls._is_dictionary_column(col):
                disk_cols[col.name] = file_storage.open_dictionary_column(col.name, col.spec)
            else:
                disk_cols[col.name] = file_storage.open_column(col.name)
        phys_size = len(disk_valid)
        n_live = int(blosc2.count_nonzero(disk_valid))
        capacity = max(phys_size, 1)

        mem_storage = InMemoryTableStorage()
        bool_chunks, bool_blocks = compute_chunks_blocks((capacity,), dtype=np.dtype(np.bool_))
        # Align fixed-size scalar columns (and the _valid_rows mask) on one
        # shared grid so lazy expressions over them take the fast_eval path.
        shared_chunks, shared_blocks, aligned_names = cls._compute_aligned_grid(schema.columns, capacity)

        mem_valid = mem_storage.create_valid_rows(
            shape=(capacity,),
            chunks=shared_chunks if shared_chunks is not None else bool_chunks,
            blocks=shared_blocks if shared_blocks is not None else bool_blocks,
        )
        if phys_size > 0:
            mem_valid[:phys_size] = disk_valid[:]

        mem_cols: dict[str, blosc2.NDArray | ListArray | _ScalarVarLenArray] = {}
        for col in schema.columns:
            name = col.name
            if cls._is_list_column(col):
                mem_col = mem_storage.create_list_column(name, spec=col.spec, cparams=None, dparams=None)
                mem_col.extend(disk_cols[name][:])
                mem_col.flush()
                mem_cols[name] = mem_col
                continue
            if cls._is_varlen_scalar_column(col):
                mem_col = mem_storage.create_varlen_scalar_column(name, spec=col.spec)
                mem_col.extend(iter(disk_cols[name]))
                mem_col.flush()
                mem_cols[name] = mem_col
                continue
            if cls._is_dictionary_column(col):
                mem_col = mem_storage.create_dictionary_column(name, spec=col.spec)
                disk_dc = disk_cols[name]
                # Copy dictionary values
                for v in disk_dc.dictionary:
                    mem_col.encode(v)
                # Copy codes
                if phys_size > 0:
                    mem_col.codes[:phys_size] = disk_dc.codes[:phys_size]
                mem_cols[name] = mem_col
                continue
            shape = cls._column_physical_shape(col, capacity)
            if name in aligned_names:
                col_chunks, col_blocks = shared_chunks, shared_blocks
            else:
                col_chunks, col_blocks = cls._column_chunks_blocks(col, shape)
            mem_col = mem_storage.create_column(
                name,
                dtype=col.dtype,
                shape=shape,
                chunks=col_chunks,
                blocks=col_blocks,
                cparams=None,
                dparams=None,
            )
            if phys_size > 0:
                mem_col[:phys_size] = disk_cols[name][:]
            mem_cols[name] = mem_col

        file_storage.close()

        obj = cls.__new__(cls)
        obj._row_type = None
        obj._validate = True
        obj._table_cparams = None
        obj._table_dparams = None
        obj._storage = mem_storage
        obj._read_only = False
        obj._schema = schema
        obj._cols = mem_cols
        obj._col_widths = {col.name: max(len(col.name), col.display_width) for col in schema.columns}
        obj.col_names = col_names
        obj.auto_compact = False
        obj._create_summary_index = schema_dict.get("create_summary_index", True)
        obj._summary_indexes_built = schema_dict.get("summary_indexes_built", False)
        obj.base = None
        obj._valid_rows = mem_valid
        obj._n_rows = n_live
        obj._last_pos = None  # resolve lazily on first write
        obj._computed_cols = {}
        obj._materialized_cols = {}
        obj._expr_index_arrays = {}
        obj._load_computed_cols_from_schema(schema_dict)
        obj._load_materialized_cols_from_schema(schema_dict)
        return obj

    @classmethod
    def _make_view(cls, parent: CTable, new_valid_rows: blosc2.NDArray) -> CTable:
        """Construct a read-only view sharing *parent*'s columns."""
        obj = cls.__new__(cls)
        obj._row_type = parent._row_type
        obj._validate = parent._validate
        obj._table_cparams = parent._table_cparams
        obj._table_dparams = parent._table_dparams
        obj._storage = None
        obj._read_only = parent._read_only  # inherit: only True for mode="r" disk tables
        obj._schema = parent._schema
        obj._cols = parent._cols  # shared — views cannot change row structure
        obj._computed_cols = parent._computed_cols  # shared — LazyExpr refs remain valid
        obj._materialized_cols = parent._materialized_cols
        obj._expr_index_arrays = parent._expr_index_arrays
        obj._cached_index_catalog = None
        obj._cached_live_positions = None
        obj._col_widths = parent._col_widths
        obj.col_names = parent.col_names
        obj.auto_compact = parent.auto_compact
        obj._create_summary_index = parent._create_summary_index
        obj._summary_indexes_built = True  # views never build indexes themselves
        obj.base = parent
        obj._valid_rows = new_valid_rows
        # Keep row counts lazy for views.  Many pipelines (e.g. where(...).sort_by(...))
        # immediately scan the mask for positions, so counting here would duplicate work.
        obj._n_rows = None
        obj._last_pos = None
        return obj

    def _view_from_positions(self, positions: np.ndarray) -> CTable:
        """Return a row-filter view from physical row positions."""
        positions = np.asarray(positions, dtype=np.intp)
        total = len(self._valid_rows)
        if len(positions):
            positions = positions[(positions >= 0) & (positions < total)]
        if len(positions) and self._known_n_rows() != total:
            keep = np.asarray(self._valid_rows[positions], dtype=bool)
            positions = positions[keep]
        result = CTable._make_view(self, self._bool_mask_from_positions(positions, total))
        result._cached_live_positions = positions
        result._n_rows = len(positions)
        return result

    # Rows per chunk when building a positions view mask.  The mask is written
    # one chunk at a time, so this trades peak memory against build time: each
    # touched chunk costs ~one compress (time) and ~2x this many bytes (bool)
    # held transiently (peak).  A full numpy mask is a single fast compress but
    # materializes the whole ``total``-element array (~1 byte/row).  ~4M keeps
    # peak to a few MB while staying few enough chunks that the per-chunk
    # compress cost does not dominate; lower it for less memory at more chunks.
    _MASK_BUILD_CHUNK_ROWS = 4_000_000

    def _bool_mask_from_positions(self, positions: np.ndarray, total: int) -> blosc2.NDArray:
        """Build the compressed boolean row mask for a positions view.

        The mask is written chunk-by-chunk (only chunks that contain a position
        are touched), so it never materializes the full ``total``-element
        uncompressed array — a sparse selection out of millions of rows costs a
        few chunk buffers instead of the whole mask.  The chunk size is capped
        (not the column grid) so peak memory stays bounded without paying a
        compress per column-chunk.
        """
        if total <= 0:
            return blosc2.zeros(max(total, 0), dtype=np.bool_)
        chunk_len = min(total, self._MASK_BUILD_CHUNK_ROWS)
        mask = blosc2.zeros(total, dtype=np.bool_, chunks=(chunk_len,))
        if len(positions) == 0:
            return mask
        chunk_len = int(mask.chunks[0])
        chunk_of = positions // chunk_len
        for c in np.unique(chunk_of):
            lo = int(c) * chunk_len
            hi = min(lo + chunk_len, total)
            buf = np.zeros(hi - lo, dtype=np.bool_)
            buf[positions[chunk_of == c] - lo] = True
            mask[lo:hi] = buf
        return mask

    def view(self, new_valid_rows):
        """Return a row-filter view backed by a boolean mask array without copying data."""
        if isinstance(new_valid_rows, np.ndarray) and new_valid_rows.dtype == np.bool_:
            new_valid_rows = blosc2.asarray(new_valid_rows)
        if not (
            isinstance(new_valid_rows, (blosc2.NDArray, blosc2.LazyExpr))
            and (getattr(new_valid_rows, "dtype", None) == np.bool_)
        ):
            raise TypeError(
                f"Expected boolean blosc2.NDArray or LazyExpr, got {type(new_valid_rows).__name__}"
            )

        new_valid_rows = (
            new_valid_rows.compute() if isinstance(new_valid_rows, blosc2.LazyExpr) else new_valid_rows
        )

        if len(self._valid_rows) != len(new_valid_rows):
            raise ValueError()

        return CTable._make_view(self, new_valid_rows)

    @staticmethod
    def _normalize_row_take_indices(indices, size: int) -> np.ndarray:
        if isinstance(indices, blosc2.NDArray):
            indices = indices[()]
        indices = np.asarray(indices)
        if indices.ndim == 0:
            indices = indices.reshape(1)
        if indices.ndim != 1:
            raise ValueError("CTable.take indices must be a 1-D integer array")
        if indices.size == 0:
            return np.ascontiguousarray(indices, dtype=np.int64)
        if not np.issubdtype(indices.dtype, np.integer):
            raise TypeError("CTable.take indices must be integers")
        normalized = np.ascontiguousarray(indices, dtype=np.int64)
        negative = normalized < 0
        if np.any(negative):
            normalized = normalized.copy()
            normalized[negative] += size
        if np.any((normalized < 0) | (normalized >= size)):
            raise IndexError("CTable.take index out of bounds")
        return normalized

    def take(self, indices, /) -> CTable:
        """Return a compact table containing rows at the requested positions.

        Indices are interpreted as logical row positions among live rows.  The
        returned table preserves the order of ``indices`` and any duplicates,
        unlike mask-based views.
        """
        logical_pos = self._normalize_row_take_indices(indices, self.nrows)
        physical_pos = self._live_positions_from_valid_rows_chunks()[logical_pos]
        n = len(physical_pos)

        result = self._empty_copy(capacity=n)
        for col in self._schema.columns:
            col_name = col.name
            arr = self._cols[col_name]
            if self._is_list_column(col):
                result._cols[col_name].extend((arr[int(pos)] for pos in physical_pos), validate=False)
                result._cols[col_name].flush()
            elif self._is_varlen_scalar_column(col):
                result._cols[col_name].extend(arr[int(pos)] for pos in physical_pos)
                result._cols[col_name].flush()
            elif self._is_dictionary_column(col):
                for v in arr.dictionary:
                    result._cols[col_name].encode(v)
                result._cols[col_name].codes[:n] = arr.codes._take_numpy(physical_pos, axis=0)
            else:
                result._cols[col_name][:n] = arr._take_numpy(physical_pos, axis=0)

        result._valid_rows[:n] = True
        result._valid_rows[n:] = False
        result._n_rows = n
        result._last_pos = n - 1 if n > 0 else None
        return result

    def slice(self, start, stop=None, /, *, copy: bool = True) -> CTable:
        """Return a contiguous range of live (non-deleted) rows.

        The range is given the way :func:`range` takes its bounds, either as a
        single stop (``table.slice(stop)``), as start/stop integers
        (``table.slice(start, stop)``), or as a Python ``slice``
        (``table.slice(slice(start, stop))``).  Negative bounds count from the
        end; ``step`` is not supported.

        Parameters
        ----------
        start, stop:
            Range bounds, interpreted as logical positions among the live rows.
        copy:
            When ``True`` (the default, mirroring :meth:`NDArray.slice`) a compact
            copy of the range is returned.  When ``False`` a zero-copy view is
            returned instead, sharing the parent's column data (read-only, like
            :meth:`head`/:meth:`tail`).

        Returns
        -------
        out: :ref:`CTable`
            The requested rows, re-indexed from 0.
        """
        if isinstance(start, slice):
            if stop is not None:
                raise TypeError("pass either a slice or start/stop integers, not both")
            key = start
        else:
            key = slice(0, start) if stop is None else slice(start, stop)
        if key.step not in (None, 1):
            raise ValueError("CTable.slice does not support a step")
        lo, hi, _ = key.indices(self.nrows)
        hi = max(lo, hi)
        if copy:
            return self.take(np.arange(lo, hi, dtype=np.int64))
        positions = self._live_positions_from_valid_rows_chunks()[lo:hi]
        return self._view_from_positions(np.asarray(positions))

    def head(self, N: int = 5) -> CTable:
        """Return a view of the first *N* live rows (default 5)."""
        if N <= 0:
            return self.view(blosc2.zeros(shape=len(self._valid_rows), dtype=np.bool_))
        _slp = getattr(self, "_cached_live_positions", None)
        if _slp is not None and self.base is not None:
            # A lazily-sorted view: physical row order is not row order, so the
            # first N *logical* rows are the first N entries of the stored
            # permutation, not the first N physical positions.
            return self._view_from_positions(_slp[:N])
        if self._n_rows <= N:
            return self.view(self._valid_rows)

        # Reuse _find_physical_index: physical position of the (N-1)-th live row
        arr = self._valid_rows
        pos_N_true = _find_physical_index(arr, N - 1)

        if pos_N_true < len(arr) // 2:
            mask_arr = blosc2.zeros(shape=len(arr), dtype=np.bool_)
            mask_arr[: pos_N_true + 1] = True
        else:
            mask_arr = blosc2.ones(shape=len(arr), dtype=np.bool_)
            mask_arr[pos_N_true + 1 :] = False

        mask_arr = (mask_arr & self._valid_rows).compute()
        return self.view(mask_arr)

    def tail(self, N: int = 5) -> CTable:
        """Return a view of the last *N* live rows (default 5)."""
        if N <= 0:
            return self.view(blosc2.zeros(shape=len(self._valid_rows), dtype=np.bool_))
        _slp = getattr(self, "_cached_live_positions", None)
        if _slp is not None and self.base is not None:
            # See head(): physical order is not row order for a sorted view.
            return self._view_from_positions(_slp[-N:] if len(_slp) > N else _slp)
        if self._n_rows <= N:
            return self.view(self._valid_rows)

        # Physical position of the first row we want = logical index (nrows - N)
        arr = self._valid_rows
        pos_start = _find_physical_index(arr, self._n_rows - N)

        if pos_start > len(arr) // 2:
            mask_arr = blosc2.zeros(shape=len(arr), dtype=np.bool_)
            mask_arr[pos_start:] = True
        else:
            mask_arr = blosc2.ones(shape=len(arr), dtype=np.bool_)
            if pos_start > 0:
                mask_arr[:pos_start] = False

        mask_arr = (mask_arr & self._valid_rows).compute()
        return self.view(mask_arr)

    def sample(self, n: int, *, seed: int | None = None) -> CTable:
        """Return a read-only view of *n* randomly chosen live rows.

        Parameters
        ----------
        n:
            Number of rows to sample.  If *n* >= number of live rows,
            returns a view of the whole table.
        seed:
            Optional random seed for reproducibility.

        Returns
        -------
        CTable
            A read-only view sharing columns with this table.
        """
        if n <= 0:
            return self.view(blosc2.zeros(shape=len(self._valid_rows), dtype=np.bool_))
        if n >= self._n_rows:
            return self.view(self._valid_rows)

        rng = np.random.default_rng(seed)
        all_pos = self._live_positions_from_valid_rows_chunks()
        chosen = rng.choice(all_pos, size=n, replace=False)

        mask = np.zeros(len(self._valid_rows), dtype=np.bool_)
        mask[chosen] = True
        return self.view(blosc2.asarray(mask))

    def select(self, cols: list[str]) -> CTable:
        """Return a column-projection view exposing only *cols*.

        The returned object shares the underlying NDArrays with this table
        (no data is copied).  Row filtering and value writes work as usual;
        structural mutations (add/drop/rename column, append, …) are blocked.

        Parameters
        ----------
        cols:
            Ordered list of column names to keep.  For tables with **nested
            (dotted) column names**, a struct-prefix name automatically expands
            to all descendant leaves::

                t.select(["trip.begin"])   # expands to trip.begin.lon, trip.begin.lat
                t.select(["trip"])          # expands to all trip.* leaves

        Raises
        ------
        KeyError
            If any name in *cols* is not a column of this table (and does not
            match any struct prefix).
        ValueError
            If *cols* is empty.
        """
        if not cols:
            raise ValueError("select() requires at least one column name.")
        expanded_cols = []
        for name in cols:
            expanded_cols.extend(self._expand_logical_column_selector(name))
        cols = expanded_cols
        for name in cols:
            if name not in self._cols and name not in self._computed_cols:
                raise KeyError(f"No column named {name!r}. Available: {self.col_names}")

        obj = CTable.__new__(CTable)
        obj._row_type = self._row_type
        obj._validate = self._validate
        obj._table_cparams = self._table_cparams
        obj._table_dparams = self._table_dparams
        obj._storage = None
        obj._read_only = self._read_only
        obj._valid_rows = self._valid_rows
        obj._n_rows = self._known_n_rows()
        obj._last_pos = self._last_pos
        obj.auto_compact = self.auto_compact
        obj._create_summary_index = self._create_summary_index
        obj._summary_indexes_built = True  # views never build indexes
        obj.base = self

        # Stored columns — same NDArray objects, no copy.  Project lazily so a
        # column is only opened when the view actually reads it: selecting then
        # touching a subset (or aggregating one column) no longer opens every
        # projected column up front.
        stored_names = [name for name in cols if name in self._cols]
        obj._cols = _LazyColumnDict(obj, self._storage, stored_names, source_cols=self._cols)
        obj.col_names = list(cols)
        obj._materialized_cols = {
            name: dict(self._materialized_cols[name]) for name in cols if name in self._materialized_cols
        }
        obj._expr_index_arrays = self._expr_index_arrays
        obj._cached_index_catalog = None
        obj._cached_live_positions = getattr(self, "_cached_live_positions", None)

        # Computed columns — share the same definitions (LazyExpr refs remain valid)
        obj._computed_cols = {
            name: self._computed_cols[name] for name in cols if name in self._computed_cols
        }

        # Rebuild schema for the selected stored columns only
        stored_sel = [n for n in cols if n in self._cols]
        sel_set = set(stored_sel)
        sel_compiled = [c for c in self._schema.columns if c.name in sel_set]
        # Preserve caller-specified order
        order = {name: i for i, name in enumerate(stored_sel)}
        sel_compiled.sort(key=lambda c: order[c.name])
        obj._schema = CompiledSchema(
            columns=sel_compiled,
            columns_by_name={c.name: c for c in sel_compiled},
            row_cls=self._schema.row_cls,
        )
        obj._col_widths = {name: self._col_widths[name] for name in cols if name in self._col_widths}
        return obj

    def group_by(
        self,
        keys: str | Sequence[str],
        *,
        sort: bool | None = None,
        dropna: bool = True,
        engine: str = "auto",
        chunk_size: int | None = None,
    ):
        """Return a deferred group-by object for this table.

        Parameters
        ----------
        keys:
            Column name or sequence of column names to group by.
        sort:
            Controls the ordering of the output groups:

            * ``True`` -- always return groups sorted by key.
            * ``False`` -- do not sort; groups come out in a deterministic but
              unspecified order (integer/dense keys are still ascending, as that
              order is free).
            * ``None`` (default) -- *auto*: sort only when the path can do so
              cheaply.  Integer and dictionary (string) keys are sorted (free or
              vectorized); float and multi-key results, whose only ordering is a
              Python sort over every distinct group, are left unsorted.  This
              avoids paying an O(G log G) Python sort that can rival the grouping
              cost itself on high-cardinality data.

            .. list-table::
               :header-rows: 1

               * - key kind
                 - ``sort=False``
                 - ``sort=None`` (auto)
                 - ``sort=True``
               * - int / dense
                 - ascending\\ :sup:`*`
                 - ascending
                 - ascending
               * - dictionary / string
                 - first-seen code order
                 - sorted by string
                 - sorted by string
               * - float
                 - unspecified\\ :sup:`**`
                 - unspecified\\ :sup:`**`
                 - sorted
               * - multi-key / generic
                 - unspecified\\ :sup:`**`
                 - unspecified\\ :sup:`**`
                 - sorted

            :sup:`*` integer dense is always ascending, so ``False`` means "no
            ordering promise", not "guaranteed unsorted".

            :sup:`**` deterministic for a given table, but order is an
            implementation detail (hash-bucket or first-appearance depending on
            the path); do not rely on it -- pass ``sort=True`` if you need order.

            .. note::
               This differs from pandas, whose ``groupby`` defaults to
               ``sort=True``.  blosc2 targets large, potentially on-disk data
               where an unconditional group sort is not always cheap relative to
               the scan, hence the ``None`` auto-default.
        dropna:
            If ``True`` (default), rows with null/NaN group keys are skipped.
            If ``False``, null/NaN keys form their own group.
        engine:
            Execution engine for built-in aggregations (``size``, ``count``,
            ``sum``, ``mean``, ``min``, ``max``, ``argmin``, ``argmax``):

            * ``"auto"`` (default) -- currently always the NumPy/Cython
              chunked implementation; may choose automatically in the future.
            * ``"numpy"`` -- explicitly request the NumPy/Cython chunked
              implementation.
            * ``"jit"`` -- reserved for a miniexpr-JIT execution path; not
              implemented yet (raises :class:`NotImplementedError`).
        chunk_size:
            Optional number of physical rows processed per chunk.

        Returns
        -------
        CTableGroupBy
            A lightweight deferred operation builder.  Call methods such as
            ``.size()``, ``.count(column)`` or ``.agg({...})`` to materialize a
            grouped result as a new :class:`CTable`.
        """
        if engine not in ("auto", "numpy", "jit"):
            raise ValueError(f"engine must be 'auto', 'numpy', or 'jit', got {engine!r}")
        if engine == "jit":
            raise NotImplementedError(
                "engine='jit' is reserved for a future miniexpr-JIT execution path; "
                "use 'auto' or 'numpy' for now."
            )
        from blosc2.groupby import CTableGroupBy

        return CTableGroupBy(self, keys, sort=sort, dropna=dropna, engine=engine, chunk_size=chunk_size)

    def describe(self) -> None:
        """Print a per-column statistical summary.

        Numeric columns (int, float): count, mean, std, min, max.
        Bool columns: count, true-count, true-%.
        String columns: count, min (lex), max (lex), n-unique.
        """
        n = self._n_rows
        lines = []
        lines.append(f"CTable  {n:,} rows × {self.ncols} cols")
        lines.append("")

        for name in self.col_names:
            col = self[name]
            dtype = col.dtype
            spec = self._schema.columns_by_name.get(name)
            label = self._dtype_info_label(dtype, spec.spec if spec else None)
            lines.append(f"  {name}  [{label}]")

            if n == 0:
                lines.append("    (empty)")
                lines.append("")
                continue

            nc = col.null_count()
            n_nonnull = n - nc

            if isinstance(spec.spec, NDArraySpec) if spec is not None else False:
                lines.append(f"    count      : {n:,}")
                lines.append(f"    item_shape : {spec.spec.item_shape}")
                lines.append(
                    "    (scalar stats not available for ndarray columns; use column reductions with axis=)"
                )
            elif isinstance(spec.spec, ListSpec) if spec is not None else False:
                lines.append(f"    count : {n:,}")
                lines.append("    (stats not available for list columns)")
            elif dtype.kind in "biufc" and dtype.kind != "c":
                # numeric + bool
                if dtype.kind == "b":
                    arr = col[:]
                    # Exclude null sentinels from true/false counts
                    if col.null_value is not None:
                        arr = arr[col.notnull()]
                    true_n = int(arr.sum())
                    lines.append(f"    count : {n:,}")
                    if nc > 0:
                        lines.append(f"    null  : {nc:,}  ({nc / n * 100:.1f} %)")
                    lines.append(f"    true  : {true_n:,}  ({true_n / n * 100:.1f} %)")
                    lines.append(f"    false : {n - true_n - nc:,}  ({(n - true_n - nc) / n * 100:.1f} %)")
                else:
                    fmt = ".4g"
                    lines.append(f"    count : {n:,}")
                    if nc > 0:
                        lines.append(f"    null  : {nc:,}  ({nc / n * 100:.1f} %)")
                    if n_nonnull > 0:
                        mn = col.min()
                        mx = col.max()
                        avg = col.mean()
                        sd = col.std()
                        lines.append(f"    mean  : {avg:{fmt}}")
                        lines.append(f"    std   : {sd:{fmt}}")
                        lines.append(f"    min   : {mn:{fmt}}")
                        lines.append(f"    max   : {mx:{fmt}}")
                    else:
                        lines.append("    (all values are null)")
            elif dtype.kind in "US":
                nu = len(col.unique())
                lines.append(f"    count   : {n:,}")
                if nc > 0:
                    lines.append(f"    null    : {nc:,}  ({nc / n * 100:.1f} %)")
                lines.append(f"    unique  : {nu:,}")
                if n_nonnull > 0:
                    mn = col.min()
                    mx = col.max()
                    lines.append(f"    min     : {str(mn)!r}")
                    lines.append(f"    max     : {str(mx)!r}")
                else:
                    lines.append("    (all values are null)")
            else:
                lines.append(f"    count : {n:,}")
                lines.append(f"    (stats not available for dtype {dtype})")

            lines.append("")

        print("\n".join(lines))

    def cov(self) -> np.ndarray:
        """Return the covariance matrix as a numpy array.

        Only int, float, and bool columns are supported.  Bool columns are
        cast to int (0/1) before computation.  Complex columns raise
        :exc:`TypeError`.

        Returns
        -------
        numpy.ndarray
            Shape ``(ncols, ncols)``.  Column order matches
            :attr:`col_names`.

        Raises
        ------
        TypeError
            If any column has an unsupported dtype (complex, string, …).
        ValueError
            If the table has fewer than 2 live rows (covariance undefined).
        """
        for name in self.col_names:
            col_info = self._schema.columns_by_name.get(name)
            if col_info is not None and self._is_ndarray_column(col_info):
                raise TypeError(
                    f"Column {name!r} is a fixed-shape ndarray column and is not supported by cov(). "
                    "Materialize scalar generated columns first."
                )
            dtype = self._col_dtype(name)
            if dtype is None or not (
                np.issubdtype(dtype, np.integer) or np.issubdtype(dtype, np.floating) or dtype == np.bool_
            ):
                raise TypeError(
                    f"Column {name!r} has dtype {dtype} which is not supported by cov(). "
                    "Only int, float, and bool columns are allowed."
                )

        if self._n_rows < 2:
            raise ValueError(f"cov() requires at least 2 live rows, got {self._n_rows}.")

        # Build (n_cols, n_rows) matrix — one row per column.
        # Compute a combined null mask: any row that is null in *any* column
        # is excluded from all columns (listwise deletion).
        raw_arrays = []
        null_union = None
        for name in self.col_names:
            col = self[name]
            arr = col[:]
            nm = col._null_mask_for(arr)
            if nm.any():
                null_union = nm if null_union is None else (null_union | nm)
            raw_arrays.append(arr)

        arrays = []
        for arr in raw_arrays:
            if null_union is not None:
                arr = arr[~null_union]
            if arr.dtype == np.bool_:
                arr = arr.astype(np.int8)
            arrays.append(arr.astype(np.float64))

        n_valid = len(arrays[0]) if arrays else 0
        if n_valid < 2:
            raise ValueError(
                f"cov() requires at least 2 non-null rows, got {n_valid} after excluding nulls."
            )

        data = np.stack(arrays, axis=0)  # shape (ncols, n_valid)
        return np.atleast_2d(np.cov(data))

    # ------------------------------------------------------------------
    # Arrow interop
    # ------------------------------------------------------------------

    @staticmethod
    def _require_pyarrow(context: str):
        try:
            import pyarrow as pa
        except ImportError:
            raise ImportError(
                f"pyarrow is required for {context}. Install it with: pip install pyarrow"
            ) from None
        return pa

    @staticmethod
    def _require_pyarrow_parquet(context: str):
        try:
            import pyarrow.parquet as pq
        except ImportError:
            raise ImportError(
                f"pyarrow is required for {context}. Install it with: pip install pyarrow"
            ) from None
        return pq

    @staticmethod
    def _validate_arrow_batch_size(batch_size: int) -> None:
        if batch_size <= 0:
            raise ValueError("batch_size must be greater than 0")

    def _resolve_arrow_columns(self, columns, include_computed: bool = True) -> list[str]:
        if columns is None:
            names = list(self.col_names)
            if not include_computed:
                names = [name for name in names if name not in self._computed_cols]

            # If top-level struct aliases are present in schema metadata (virtual
            # entries not physically stored), prefer exporting them instead of
            # their descendant dotted leaves.
            virtual_structs = [
                n
                for n, cc in self._schema.columns_by_name.items()
                if n not in self.col_names and isinstance(cc.spec, StructSpec)
            ]
            for alias in sorted(virtual_structs, key=len, reverse=True):
                alias_parts = split_field_path(alias)
                children = [
                    n
                    for n in names
                    if split_field_path(n)[: len(alias_parts)] == alias_parts
                    and len(split_field_path(n)) > len(alias_parts)
                ]
                if not children:
                    continue
                first = min(names.index(c) for c in children)
                child_set = set(children)
                names = [n for n in names if n not in child_set]
                names.insert(first, alias)
        else:
            names = []
            for name in columns:
                names.extend(self._expand_logical_column_selector(name))
        if len(set(names)) != len(names):
            raise ValueError("columns must be unique")
        for name in names:
            if name not in self.col_names and name not in self._schema.columns_by_name:
                raise KeyError(f"No column named {name!r}. Available: {self.col_names}")
        return names

    @staticmethod
    def _pa_type_from_spec(pa, spec):
        if isinstance(spec, DictionarySpec):
            return pa.dictionary(pa.int32(), pa.string(), ordered=spec.ordered)
        if isinstance(spec, Utf8Spec):
            # Always large_string: 64-bit offsets match the int64 offsets array,
            # so multi-GB string columns export without int32-offset overflow.
            return pa.large_string()
        if isinstance(spec, VLStringSpec):
            return pa.string()
        if isinstance(spec, VLBytesSpec):
            return pa.large_binary()
        if isinstance(spec, ListSpec):
            return pa.list_(CTable._pa_type_from_spec(pa, spec.item_spec))
        if isinstance(spec, NDArraySpec):
            return pa.list_(pa.from_numpy_dtype(spec.dtype), list_size=int(np.prod(spec.item_shape)))
        if isinstance(spec, timestamp):
            return pa.timestamp(spec.unit, tz=spec.timezone)
        if isinstance(spec, StructSpec):
            return pa.struct(
                [pa.field(name, CTable._pa_type_from_spec(pa, child)) for name, child in spec.fields.items()]
            )
        if isinstance(spec, ObjectSpec):
            raise TypeError(
                "ObjectSpec columns do not have a fixed Arrow type; materialize values explicitly"
            )
        if spec.to_metadata_dict().get("kind") == "bool":
            return pa.bool_()
        dtype = getattr(spec, "dtype", None)
        if dtype is None:
            raise TypeError(f"No Arrow type for blosc2 spec {spec!r}")
        kind = dtype.kind
        if kind == "U":
            return pa.string()
        if kind == "S":
            return pa.large_binary()
        return pa.from_numpy_dtype(dtype)

    def _export_arrow_names(self, names: list[str]) -> list[str]:
        nested = self._schema.metadata.get("nested") if self._schema.metadata else None
        exported = list(names)
        if isinstance(nested, dict):
            root_meta = nested.get("root")
            if isinstance(root_meta, dict):
                physical = root_meta.get("physical")
                if isinstance(physical, str) and physical:
                    exported = ["" if n == physical else n for n in exported]
        for i, n in enumerate(names):
            cc = self._schema.columns_by_name.get(n)
            if n not in self.col_names and cc is not None and isinstance(cc.spec, StructSpec):
                parts = split_field_path(n)
                if len(parts) == 1:
                    exported[i] = parts[0]
        return exported

    def _arrow_schema_for_columns(self, columns=None, *, include_computed: bool = True):
        pa = self._require_pyarrow("to_arrow()/to_parquet()")
        names = self._resolve_arrow_columns(columns, include_computed=include_computed)
        arrow_names = self._export_arrow_names(names)
        fields = []
        for name, arrow_name in zip(names, arrow_names, strict=True):
            cc = self._schema.columns_by_name.get(name)
            metadata = None
            if cc is not None:
                pa_type = self._pa_type_from_spec(pa, cc.spec)
                if isinstance(cc.spec, NDArraySpec):
                    metadata = {b"blosc2:ndarray_shape": json.dumps(list(cc.spec.item_shape)).encode()}
            else:
                pa_type = pa.from_numpy_dtype(np.asarray(self[name][:0]).dtype)
            fields.append(pa.field(arrow_name, pa_type, metadata=metadata))
        return pa.schema(fields)

    def iter_arrow_batches(  # noqa: C901
        self,
        *,
        columns: list[str] | None = None,
        batch_size: int = _BATCH_SIZE_DEFAULT,
        include_computed: bool = True,
    ):
        """Yield live rows as bounded-size :class:`pyarrow.RecordBatch` objects."""
        pa = self._require_pyarrow("iter_arrow_batches()")
        self._validate_arrow_batch_size(batch_size)
        self._flush_varlen_columns()
        names = self._resolve_arrow_columns(columns, include_computed=include_computed)
        arrow_names = self._export_arrow_names(names)

        # Dictionary columns need the physical positions of their live rows.
        # This depends only on self._valid_rows (fixed for the whole call), so
        # it is computed once here instead of once per batch per dictionary
        # column — the previous per-batch recompute was an O(n_rows) scan
        # repeated O(n_rows / batch_size) times.
        dict_real_pos = None
        if any(name in self.col_names and self[name].is_dictionary for name in names):
            dict_real_pos = blosc2.where(self._valid_rows, _arange(len(self._valid_rows))).compute()

        for start in range(0, self._n_rows, batch_size):
            stop = min(start + batch_size, self._n_rows)
            arrays = []
            for name in names:
                cc = self._schema.columns_by_name.get(name)
                if name not in self.col_names and cc is not None and isinstance(cc.spec, StructSpec):
                    values = self[name][start:stop]
                    arrays.append(pa.array(values, type=self._pa_type_from_spec(pa, cc.spec)))
                    continue
                col = self[name]
                if col.is_list:
                    spec = self._schema.columns_by_name[name].spec
                    arrays.append(pa.array(col[start:stop], type=self._pa_type_from_spec(pa, spec)))
                    continue
                if col.is_utf8:
                    spec = self._schema.columns_by_name[name].spec
                    arr8 = self._cols[name]
                    nv = col.null_value
                    if self.base is None and self._last_pos == self._n_rows and stop <= arr8._persisted_rows:
                        # Dense root table: logical rows == persisted rows, so
                        # export straight from the offsets/bytes buffers with
                        # no per-row decode (storage is already Arrow layout).
                        arrays.append(arr8.arrow_slice(pa, start, stop, nv))
                        continue
                    values = col[start:stop]  # StringDType array with sentinel nulls
                    null_mask = col._null_mask_for(values) if nv is not None else None
                    arrays.append(
                        pa.array(
                            values.astype(object),
                            type=self._pa_type_from_spec(pa, spec),
                            mask=null_mask if null_mask is not None and null_mask.any() else None,
                        )
                    )
                    continue
                if col.is_varlen_scalar:
                    spec = self._schema.columns_by_name[name].spec
                    values = col[start:stop]  # list of str/bytes/None
                    arrays.append(pa.array(values, type=self._pa_type_from_spec(pa, spec)))
                    continue
                if col.is_dictionary:
                    dc = self._cols[name]  # DictionaryColumn
                    spec = self._schema.columns_by_name[name].spec
                    # Physical positions for live rows in [start, stop),
                    # precomputed once above (not per batch/column).
                    batch_real_pos = dict_real_pos[start:stop]
                    if len(batch_real_pos) == 0:
                        pa_dict = pa.array(dc.dictionary, type=pa.string())
                        pa_indices = pa.array([], type=pa.int32())
                        arrays.append(
                            pa.DictionaryArray.from_arrays(pa_indices, pa_dict, ordered=spec.ordered)
                        )
                    else:
                        raw_codes = dc.codes[batch_real_pos]
                        null_mask = raw_codes == np.int32(spec.null_code)
                        safe_codes = raw_codes.copy()
                        safe_codes[null_mask] = 0
                        pa_dict = pa.array(dc.dictionary, type=pa.string())
                        pa_indices = pa.array(
                            safe_codes,
                            type=pa.int32(),
                            mask=null_mask if null_mask.any() else None,
                        )
                        arrays.append(
                            pa.DictionaryArray.from_arrays(pa_indices, pa_dict, ordered=spec.ordered)
                        )
                    continue
                if col.is_ndarray:
                    spec = self._schema.columns_by_name[name].spec
                    values = np.asarray(col[start:stop])
                    null_mask = col._null_mask_for(values) if col.null_value is not None else None
                    pa_type = self._pa_type_from_spec(pa, spec)
                    flat_values = np.ascontiguousarray(values.reshape(-1))
                    pa_values = pa.array(flat_values, type=pa_type.value_type)
                    arrays.append(
                        pa.FixedSizeListArray.from_arrays(
                            pa_values,
                            type=pa_type,
                            mask=(
                                pa.array(null_mask, type=pa.bool_())
                                if null_mask is not None and null_mask.any()
                                else None
                            ),
                        )
                    )
                    continue
                arr = np.asarray(col[start:stop])
                nv = col.null_value
                null_mask = col._null_mask_for(arr) if nv is not None else None
                has_nulls = null_mask is not None and bool(null_mask.any())
                if arr.dtype.kind == "U":
                    values = arr.tolist()
                    if has_nulls:
                        values = [None if null_mask[i] else v for i, v in enumerate(values)]
                    arrays.append(pa.array(values, type=pa.string()))
                elif arr.dtype.kind == "S":
                    values = arr.tolist()
                    if has_nulls:
                        values = [None if null_mask[i] else v for i, v in enumerate(values)]
                    arrays.append(pa.array(values, type=pa.large_binary()))
                elif (
                    self._schema.columns_by_name.get(name) is not None
                    and self._schema.columns_by_name[name].spec.to_metadata_dict().get("kind") == "bool"
                ):
                    arrays.append(pa.array(arr == 1, mask=null_mask if has_nulls else None, type=pa.bool_()))
                elif self._schema.columns_by_name.get(name) is not None and isinstance(
                    self._schema.columns_by_name[name].spec, timestamp
                ):
                    spec = self._schema.columns_by_name[name].spec
                    values = arr.astype(f"datetime64[{spec.unit}]")
                    arrays.append(
                        pa.array(
                            values,
                            mask=null_mask if has_nulls else None,
                            type=pa.timestamp(spec.unit, tz=spec.timezone),
                        )
                    )
                else:
                    arrays.append(pa.array(arr, mask=null_mask if has_nulls else None))
            yield pa.RecordBatch.from_arrays(arrays, names=arrow_names)

    def to_arrow(self):
        """Convert all live rows to a :class:`pyarrow.Table`."""
        pa = self._require_pyarrow("to_arrow()")
        batches = list(self.iter_arrow_batches())
        schema = self._arrow_schema_for_columns()
        return pa.Table.from_batches(batches, schema=schema)

    def __arrow_c_stream__(self, requested_schema=None):
        """Arrow PyCapsule protocol: export live rows as a stream of record batches.

        Lets Arrow-native consumers (pyarrow, DuckDB, Polars, pandas) read this
        table directly, pulling decompressed batches lazily with bounded memory.
        Strict zero-copy is impossible here (the data is compressed on disk;
        decompression *is* the copy), but there is zero intermediate
        materialization: only one batch is decompressed at a time.
        """
        pa = self._require_pyarrow("__arrow_c_stream__")
        reader = pa.RecordBatchReader.from_batches(
            self._arrow_schema_for_columns(), self.iter_arrow_batches()
        )
        return reader.__arrow_c_stream__(requested_schema)

    @staticmethod
    def _auto_null_sentinel(pa, pa_type, *, null_policy: NullPolicy):
        return null_policy.sentinel_for_arrow_type(pa, pa_type)

    @staticmethod
    def _arrow_type_needs_object_fallback(pa, pa_type) -> bool:
        """True when *pa_type* has no typed CTable mapping."""
        if pa.types.is_dictionary(pa_type):
            vt = pa_type.value_type
            return not _is_arrow_string_type(pa, vt)
        if _is_arrow_string_type(pa, pa_type):
            return False
        if pa_type in (
            pa.int8(),
            pa.int16(),
            pa.int32(),
            pa.int64(),
            pa.uint8(),
            pa.uint16(),
            pa.uint32(),
            pa.uint64(),
            pa.float32(),
            pa.float64(),
            pa.bool_(),
        ):
            return False
        if _is_arrow_binary_type(pa, pa_type):
            return False
        if pa.types.is_timestamp(pa_type):
            return False
        return not (
            pa.types.is_list(pa_type)
            or pa.types.is_large_list(pa_type)
            or pa.types.is_fixed_size_list(pa_type)
            or pa.types.is_struct(pa_type)
        )

    @staticmethod
    def _arrow_type_to_spec(  # noqa: C901
        pa,
        pa_type,
        arrow_col=None,
        *,
        field_metadata=None,
        string_max_length=None,
        null_value=None,
        nullable=False,
        object_fallback: bool = False,
    ):
        import blosc2.schema as b2s

        # Handle Arrow dictionary types (dict-encoded strings)
        if pa.types.is_fixed_size_list(pa_type):
            shape = None
            if field_metadata:
                encoded = field_metadata.get(b"blosc2:ndarray_shape") or field_metadata.get(
                    "blosc2:ndarray_shape"
                )
                if encoded is not None:
                    if isinstance(encoded, bytes):
                        encoded = encoded.decode()
                    shape = tuple(int(x) for x in json.loads(encoded))
            if shape is None:
                shape = (int(pa_type.list_size),)
            value_type = pa_type.value_type
            value_spec = CTable._arrow_type_to_spec(pa, value_type, object_fallback=object_fallback)
            value_dtype = getattr(value_spec, "dtype", None)
            if value_dtype is None:
                raise TypeError(f"FixedSizeList values must have a fixed NumPy dtype, got {value_type!r}")
            if int(np.prod(shape)) != int(pa_type.list_size):
                raise ValueError(
                    f"Arrow fixed-size-list metadata shape {shape} has size {int(np.prod(shape))}, "
                    f"but the Arrow list size is {pa_type.list_size}."
                )
            return b2s.ndarray(shape, dtype=value_dtype, nullable=nullable, null_value=null_value)

        if pa.types.is_dictionary(pa_type):
            vt = pa_type.value_type
            if _is_arrow_string_type(pa, vt):
                index_type = pa_type.index_type
                # Accept signed and unsigned integer index types; validate fit in int32.
                if not (pa.types.is_integer(index_type) or pa.types.is_unsigned_integer(index_type)):
                    raise TypeError(
                        f"Dictionary column has unsupported index type {index_type!r}; "
                        "expected an integer type."
                    )
                if arrow_col is not None:
                    # Validate all indices fit in signed int32.
                    if pa.types.is_unsigned_integer(index_type):
                        max_idx = arrow_col.combine_chunks().indices.to_pandas().max(skipna=True)
                        if max_idx is not None and max_idx > np.iinfo(np.int32).max:
                            raise ValueError(
                                f"Arrow dictionary column has unsigned indices exceeding int32.max "
                                f"(max={max_idx})."
                            )
                    combined = (
                        arrow_col.combine_chunks() if hasattr(arrow_col, "combine_chunks") else arrow_col
                    )
                    n_cats = len(combined.dictionary)
                    if n_cats > np.iinfo(np.int32).max:
                        raise OverflowError(
                            f"Arrow dictionary has {n_cats} categories, exceeding int32 capacity."
                        )
                return b2s.dictionary(
                    index_type=b2s.int32(),
                    value_type=b2s.vlstring(),
                    ordered=bool(pa_type.ordered),
                    nullable=nullable,
                )
            if object_fallback:
                return b2s.object(nullable=nullable)
            raise TypeError(
                f"No blosc2 spec for Arrow dictionary type {pa_type!r} with "
                f"value type {pa_type.value_type!r}. Only string dictionary values are supported in v1."
            )

        mapping = [
            (pa.int8(), b2s.int8),
            (pa.int16(), b2s.int16),
            (pa.int32(), b2s.int32),
            (pa.int64(), b2s.int64),
            (pa.uint8(), b2s.uint8),
            (pa.uint16(), b2s.uint16),
            (pa.uint32(), b2s.uint32),
            (pa.uint64(), b2s.uint64),
            (pa.float32(), b2s.float32),
            (pa.float64(), b2s.float64),
            (pa.bool_(), b2s.bool),
        ]
        if pa.types.is_timestamp(pa_type):
            return b2s.timestamp(
                unit=pa_type.unit, timezone=pa_type.tz, nullable=nullable, null_value=null_value
            )

        for arrow_t, spec_cls in mapping:
            if pa_type == arrow_t:
                if null_value is not None and hasattr(spec_cls(), "null_value"):
                    return spec_cls(null_value=null_value)
                if null_value is not None and spec_cls is b2s.bool:
                    return spec_cls(null_value=null_value)
                return spec_cls()

        if pa.types.is_list(pa_type) or pa.types.is_large_list(pa_type):
            if arrow_col is not None:
                combined = arrow_col.combine_chunks() if hasattr(arrow_col, "combine_chunks") else arrow_col
                item_arrow_col = combined.values
                nullable = nullable or combined.null_count > 0
            else:
                item_arrow_col = None
                nullable = True
            item_string_max_length = string_max_length
            if _is_arrow_string_type(pa, pa_type.value_type):
                item_string_max_length = max(string_max_length or 1, 1_000_000)
            item_spec = CTable._arrow_type_to_spec(
                pa,
                pa_type.value_type,
                item_arrow_col,
                string_max_length=item_string_max_length,
                object_fallback=object_fallback,
            )
            return b2s.list(item_spec, nullable=nullable, storage="batch", serializer="msgpack")

        if pa.types.is_struct(pa_type):
            fields = {}
            for field in pa_type:
                child_col = None
                if arrow_col is not None:
                    combined = (
                        arrow_col.combine_chunks() if hasattr(arrow_col, "combine_chunks") else arrow_col
                    )
                    child_col = combined.field(field.name)
                child_string_max_length = string_max_length
                if _is_arrow_string_type(pa, field.type):
                    child_string_max_length = max(string_max_length or 1, 1_000_000)
                fields[field.name] = CTable._arrow_type_to_spec(
                    pa,
                    field.type,
                    child_col,
                    string_max_length=child_string_max_length,
                    nullable=field.nullable,
                    object_fallback=object_fallback,
                )
            return b2s.struct(fields, nullable=nullable)

        if _is_arrow_string_type(pa, pa_type):
            if string_max_length is None:
                from blosc2.utf8_array import have_string_dtype

                if not have_string_dtype():
                    # utf8 columns need numpy.dtypes.StringDType (NumPy >= 2.0).
                    # On older NumPy, keep the historical import behavior:
                    # variable-length msgpack strings with native-None nulls.
                    return b2s.vlstring(nullable=nullable)
                # No fixed-width threshold given: store as a variable-length
                # utf8 column (offsets + bytes, StringDType reads).
                return b2s.utf8(nullable=nullable, null_value=null_value)
            max_length = max(string_max_length, len(null_value) if null_value is not None else 1, 1)
            return b2s.string(max_length=max_length, null_value=null_value)

        if _is_arrow_binary_type(pa, pa_type):
            if string_max_length is None:
                # No fixed-width threshold given: store as variable-length scalar bytes.
                return b2s.vlbytes(nullable=nullable)
            max_length = max(string_max_length, len(null_value) if null_value is not None else 1, 1)
            return b2s.bytes(max_length=max_length, null_value=null_value)

        if object_fallback:
            return b2s.object(nullable=nullable)

        raise TypeError(
            f"No blosc2 spec for Arrow type {pa_type!r}. Supported: int8/16/32/64, "
            "uint8/16/32/64, float32/64, bool, string, binary, list, and struct. "
            "Pass object_fallback=True to CTable.from_arrow() to import unsupported Arrow types "
            "as schema-less object columns."
        )

    @staticmethod
    def _string_max_length_for_column(string_max_length, name: str):
        if isinstance(string_max_length, Mapping):
            return string_max_length.get(name)
        return string_max_length

    @classmethod
    def _compiled_columns_from_arrow(
        cls,
        pa,
        schema,
        table_for_inference,
        string_max_length,
        *,
        auto_null_sentinels: bool,
        object_fallback: bool = False,
    ):
        null_policy = get_null_policy()
        column_null_values = null_policy.column_null_values
        schema_names = set(schema.names)
        unknown_null_values = set(column_null_values) - schema_names
        if unknown_null_values:
            names = ", ".join(sorted(unknown_null_values))
            raise KeyError(f"column_null_values contains unknown columns: {names}")
        columns: list[CompiledColumn] = []
        for field in schema:
            name = field.name
            _validate_column_name(name)
            arrow_col = table_for_inference.column(name) if table_for_inference is not None else None
            field_is_ndarray = pa.types.is_fixed_size_list(field.type)
            field_is_list = (
                pa.types.is_list(field.type) or pa.types.is_large_list(field.type)
            ) and not field_is_ndarray
            field_is_struct = pa.types.is_struct(field.type)
            field_is_dictionary = pa.types.is_dictionary(field.type)
            column_string_max_length = cls._string_max_length_for_column(string_max_length, name)
            # Scalar strings without a fixed-width threshold import as utf8
            # columns, which use null sentinels like other scalar columns;
            # only binary columns keep the native-None varlen treatment.
            # On NumPy < 2.0 (no StringDType) utf8 columns are unavailable and
            # scalar strings keep the historical vlstring treatment instead.
            from blosc2.utf8_array import have_string_dtype

            field_is_varlen_scalar = (
                not field_is_list
                and not field_is_struct
                and not field_is_dictionary
                and column_string_max_length is None
                and (
                    _is_arrow_binary_type(pa, field.type)
                    or (not have_string_dtype() and _is_arrow_string_type(pa, field.type))
                )
            )
            field_needs_object_fallback = cls._arrow_type_needs_object_fallback(pa, field.type)
            if field_needs_object_fallback and not object_fallback:
                cls._arrow_type_to_spec(pa, field.type, arrow_col, object_fallback=False)
            field_is_object_fallback = object_fallback and field_needs_object_fallback
            null_value = None
            has_null_value_override = name in column_null_values
            if has_null_value_override and (
                field_is_list or field_is_struct or field_is_dictionary or field_is_object_fallback
            ):
                raise TypeError(
                    f"column_null_values only supports scalar columns and ndarray columns; {name!r} is not scalar"
                )
            if has_null_value_override and field_is_varlen_scalar:
                raise TypeError(
                    f"column_null_values is not supported for vlbytes/vlstring column {name!r}; "
                    "these columns represent nulls as native None."
                )
            if has_null_value_override:
                null_value = column_null_values[name]
            elif (
                auto_null_sentinels
                and field.nullable
                and not (
                    field_is_list
                    or field_is_struct
                    or field_is_dictionary
                    or field_is_varlen_scalar
                    or field_is_object_fallback
                )
            ):
                arrow_type_for_null = field.type.value_type if field_is_ndarray else field.type
                null_value = cls._auto_null_sentinel(pa, arrow_type_for_null, null_policy=null_policy)
            if (
                arrow_col is not None
                and arrow_col.null_count
                and not (
                    field_is_list
                    or field_is_struct
                    or field_is_dictionary
                    or field_is_varlen_scalar
                    or field_is_object_fallback
                )
                and null_value is None
            ):
                raise TypeError(
                    f"Column {name!r} contains Parquet nulls. Provide a CTable schema with a "
                    "null_value sentinel for this column."
                )
            spec = cls._arrow_type_to_spec(
                pa,
                field.type,
                arrow_col,
                field_metadata=field.metadata,
                string_max_length=column_string_max_length,
                null_value=null_value,
                nullable=field.nullable,
                object_fallback=object_fallback,
            )
            if null_value is not None and not (
                field_is_list
                or field_is_struct
                or field_is_dictionary
                or field_is_varlen_scalar
                or field_is_object_fallback
            ):
                cls._validate_null_value_for_spec(name, spec, null_value)
            columns.append(cls._compiled_column_from_spec(name, spec))
        return columns

    @classmethod
    def _compiled_column_from_spec(cls, name: str, spec: SchemaSpec) -> CompiledColumn:
        col_config = ColumnConfig(cparams=None, dparams=None, chunks=None, blocks=None)
        return CompiledColumn(
            name=name,
            py_type=spec.python_type,
            spec=spec,
            dtype=getattr(spec, "dtype", None),
            default=MISSING,
            config=col_config,
            display_width=compute_display_width(spec),
        )

    @classmethod
    def _apply_arrow_column_cparams(
        cls, columns: list[CompiledColumn], column_cparams: Mapping[str, dict[str, Any]] | None
    ) -> None:
        if column_cparams is None:
            return
        unknown = set(column_cparams) - {col.name for col in columns}
        if unknown:
            names = ", ".join(sorted(unknown))
            raise KeyError(f"column_cparams contains unknown columns: {names}")
        for col in columns:
            if col.name in column_cparams:
                if cls._is_list_column(col) or cls._is_varlen_scalar_column(col):
                    raise TypeError(
                        f"column_cparams only supports fixed-width columns; {col.name!r} is not fixed-width"
                    )
                col.config.cparams = dict(column_cparams[col.name])

    @staticmethod
    def _storage_for_arrow_import(urlpath: str | None, mode: str) -> TableStorage:
        if urlpath is None:
            return InMemoryTableStorage()
        if mode == "w" and os.path.exists(urlpath):
            if os.path.isdir(urlpath):
                shutil.rmtree(urlpath)
            else:
                os.remove(urlpath)
        elif mode == "a" and not os.path.exists(urlpath):
            raise FileNotFoundError(
                f"No CTable found at {urlpath!r}: mode='a' opens an existing table; "
                "use mode='w' to create a new one."
            )
        return FileTableStorage(urlpath, mode)

    @classmethod
    def _create_arrow_import_columns(
        cls,
        storage: TableStorage,
        columns: list[CompiledColumn],
        capacity: int,
        cparams,
        dparams,
        chunks_override: tuple[int, ...] | None = None,
        blocks_override: tuple[int, ...] | None = None,
    ):
        default_chunks, default_blocks = compute_chunks_blocks((capacity,))
        # Align fixed-size scalar columns (and the _valid_rows mask) on one
        # shared grid so lazy expressions over them take the fast_eval path.
        shared_chunks, shared_blocks, aligned_names = cls._compute_aligned_grid(columns, capacity)
        valid_chunks = shared_chunks if shared_chunks is not None else default_chunks
        valid_blocks = shared_blocks if shared_blocks is not None else default_blocks
        new_valid = storage.create_valid_rows(shape=(capacity,), chunks=valid_chunks, blocks=valid_blocks)
        new_cols: dict[str, blosc2.NDArray | ListArray | _ScalarVarLenArray | DictionaryColumn] = {}
        for col in columns:
            if cls._is_list_column(col):
                new_cols[col.name] = storage.create_list_column(
                    col.name, spec=col.spec, cparams=cparams, dparams=dparams
                )
            elif cls._is_varlen_scalar_column(col):
                new_cols[col.name] = storage.create_varlen_scalar_column(
                    col.name, spec=col.spec, cparams=cparams, dparams=dparams
                )
            elif cls._is_dictionary_column(col):
                # Create the int32 codes array at full capacity with the aligned
                # grid (codes are int32, so they match the shared numeric grid).
                # This avoids a create-then-resize and the catastrophic 4096-row
                # default chunking, and lets dict-column filters use the fast path.
                if shared_chunks is not None:
                    codes_chunks, codes_blocks = shared_chunks, shared_blocks
                else:
                    codes_chunks, codes_blocks = compute_chunks_blocks((capacity,), dtype=np.dtype(np.int32))
                new_cols[col.name] = storage.create_dictionary_column(
                    col.name,
                    spec=col.spec,
                    cparams=cparams,
                    dparams=dparams,
                    codes_shape=(capacity,),
                    codes_chunks=codes_chunks,
                    codes_blocks=codes_blocks,
                )
            else:
                shape = cls._column_physical_shape(col, capacity)
                chunks, blocks = default_chunks, default_blocks
                if col.name in aligned_names:
                    chunks, blocks = shared_chunks, shared_blocks
                elif col.dtype is not None:
                    chunks, blocks = cls._column_chunks_blocks(col, shape)
                if chunks_override is not None:
                    chunks = chunks_override
                    if blocks_override is None:
                        # Use the shared grid's blocks so all columns stay on the
                        # same (chunks, blocks) pair — required for fast_eval.
                        # Letting each dtype auto-pick its own blocks produces
                        # different values (e.g. 37376 for float32 vs 32768 for
                        # float64), which breaks alignment across columns.
                        # Guard: shared_blocks must fit within chunks_override.
                        eff_blocks = shared_blocks if shared_blocks is not None else default_blocks
                        if eff_blocks is not None and all(
                            b <= c for b, c in zip(eff_blocks, chunks_override, strict=False)
                        ):
                            blocks = eff_blocks
                        else:
                            blocks = None  # chunks too small; let blosc2 auto-pick
                if blocks_override is not None:
                    blocks = blocks_override
                new_cols[col.name] = storage.create_column(
                    col.name,
                    dtype=col.dtype,
                    shape=shape,
                    chunks=chunks,
                    blocks=blocks,
                    cparams=col.config.cparams if col.config.cparams is not None else cparams,
                    dparams=col.config.dparams if col.config.dparams is not None else dparams,
                )
        return new_cols, new_valid

    @classmethod
    def _new_arrow_import_ctable(
        cls,
        compiled,
        storage,
        new_cols,
        new_valid,
        columns,
        *,
        cparams,
        dparams,
        validate,
        create_summary_index=True,
    ):
        obj = cls.__new__(cls)
        obj._row_type = None
        obj._validate = validate
        obj._table_cparams = cparams
        obj._table_dparams = dparams
        obj._storage = storage
        obj._read_only = storage.is_read_only()
        obj._schema = compiled
        obj._cols = new_cols
        obj._col_widths = {col.name: max(len(col.name), col.display_width) for col in columns}
        obj.col_names = [col.name for col in columns]
        obj.auto_compact = False
        obj._create_summary_index = create_summary_index
        obj._summary_indexes_built = False
        obj.base = None
        obj._computed_cols = {}
        obj._materialized_cols = {}
        obj._expr_index_arrays = {}
        obj._valid_rows = new_valid
        obj._n_rows = 0
        obj._last_pos = 0
        return obj

    @staticmethod
    def _timestamp_normalizer_for_spec(spec: SchemaSpec):  # noqa: C901
        """Build a trusted Arrow-import normalizer for timestamp leaves.

        Arrow already validates list/struct values during import, so list columns
        normally skip Python-level coercion.  The exception is nested timestamps:
        ``to_pylist()`` yields ``datetime``/``numpy.datetime64`` objects, while
        msgpack-backed ListArray storage expects integer epoch offsets.  Return a
        small normalizer that descends only into branches containing timestamps,
        or ``None`` when no normalization is needed.
        """
        if isinstance(spec, timestamp):

            def normalize_timestamp(value, unit=spec.unit):
                if value is None:
                    return None
                if isinstance(value, (int, np.integer)):
                    return int(value)
                return np.datetime64(value).astype(f"datetime64[{unit}]").astype(np.int64).item()

            return normalize_timestamp

        if isinstance(spec, ListSpec):
            item_normalizer = CTable._timestamp_normalizer_for_spec(spec.item_spec)
            if item_normalizer is None:
                return None

            def normalize_list(value, item_normalizer=item_normalizer):
                if value is None:
                    return None
                for i, item in enumerate(value):
                    value[i] = item_normalizer(item)
                return value

            return normalize_list

        if isinstance(spec, StructSpec):
            field_normalizers = {
                name: normalizer
                for name, child in spec.fields.items()
                if (normalizer := CTable._timestamp_normalizer_for_spec(child)) is not None
            }
            if not field_normalizers:
                return None

            def normalize_struct(value, field_normalizers=field_normalizers):
                if value is None:
                    return None
                for name, normalizer in field_normalizers.items():
                    if name in value:
                        value[name] = normalizer(value[name])
                return value

            return normalize_struct

        return None

    @classmethod
    def _trim_arrow_import_capacity(cls, obj, n_rows: int) -> None:
        """Shrink append-only Arrow-import columns from capacity to actual row count."""
        obj._last_pos = n_rows
        obj.trim_capacity()

    @classmethod
    def _write_arrow_batches(cls, obj, batches, columns, new_cols, new_valid) -> None:
        pos = 0
        list_normalizers = {
            col.name: cls._timestamp_normalizer_for_spec(col.spec)
            for col in columns
            if cls._is_list_column(col)
        }
        # Buffer fixed-size column writes and flush them chunk-aligned so each
        # chunk is compressed once instead of being merged on every batch.
        writers = {
            col.name: _ChunkAlignedWriter(
                new_cols[col.name],
                new_cols[col.name].chunks[0],
                on_write=obj._summary_feeder(col.name),
            )
            for col in columns
            if not (
                cls._is_list_column(col)
                or cls._is_varlen_scalar_column(col)
                or cls._is_dictionary_column(col)
            )
        }
        for batch in batches:
            end = pos + len(batch)
            while end > len(new_valid):
                obj._grow()
                new_valid = obj._valid_rows
            pos = cls._write_arrow_batch(batch, columns, new_cols, new_valid, pos, list_normalizers, writers)
        for writer in writers.values():
            writer.flush()
        # All imported rows are valid; mark them in a single aligned write.
        if pos:
            new_valid[:pos] = True
        for col in columns:
            if (
                cls._is_list_column(col)
                or cls._is_varlen_scalar_column(col)
                or cls._is_dictionary_column(col)
            ):
                new_cols[col.name].flush()
        cls._trim_arrow_import_capacity(obj, pos)
        obj._n_rows = pos
        obj._last_pos = pos

    @classmethod
    def _write_arrow_batch(
        cls, batch, columns, new_cols, new_valid, pos: int, list_normalizers, writers
    ) -> int:
        m = len(batch)
        if m == 0:
            return pos
        for col in columns:
            arrow_col = batch.column(batch.schema.get_field_index(col.name))
            if cls._is_list_column(col):
                if getattr(col.spec, "serializer", None) == "arrow":
                    new_cols[col.name].extend_arrow(arrow_col)
                    continue
                # Trusted Arrow-import fast path: schema has already been inferred,
                # so avoid Python-level per-item coercion.  If nested timestamps
                # are present, normalize only those leaves before storing.
                values = arrow_col.to_pylist()
                normalizer = list_normalizers[col.name]
                if normalizer is not None:
                    values = [normalizer(value) for value in values]
                new_cols[col.name].extend(values, validate=False)
            elif cls._is_varlen_scalar_column(col):
                new_cols[col.name].extend(arrow_col.to_pylist())
            elif cls._is_dictionary_column(col):
                import pyarrow as _pa

                if _pa.types.is_dictionary(arrow_col.type):
                    # Arrow dictionary array: use unification algorithm.
                    new_cols[col.name].extend_from_arrow(_pa, arrow_col, pos, m, ordered=col.spec.ordered)
                else:
                    # Plain string array: encode values into the dictionary.
                    new_cols[col.name][pos : pos + m] = arrow_col.to_pylist()
            else:
                writers[col.name].append(cls._arrow_column_to_numpy(arrow_col, col))
        return pos + m

    @staticmethod
    def _arrow_column_to_numpy(arrow_col, col: CompiledColumn) -> np.ndarray:
        nv = getattr(col.spec, "null_value", None)
        if col.spec.to_metadata_dict().get("kind") == "bool" and col.dtype == np.dtype(np.uint8):
            return np.array([nv if v is None else int(v) for v in arrow_col.to_pylist()], dtype=np.uint8)
        if isinstance(col.spec, NDArraySpec):
            values = arrow_col.to_pylist()
            arr = CTable._coerce_ndarray_batch(col.name, col.spec, values, len(values))
            return arr.reshape((len(values), *col.spec.item_shape))
        if isinstance(col.spec, timestamp):
            arr = (
                arrow_col.to_numpy(zero_copy_only=False)
                .astype(f"datetime64[{col.spec.unit}]")
                .astype(np.int64)
            )
            if arrow_col.null_count and nv is not None and int(nv) != int(np.iinfo(np.int64).min):
                arr[arr == np.iinfo(np.int64).min] = int(nv)
            return arr.astype(col.dtype, copy=False)
        if col.dtype.kind in "US":
            values = arrow_col.to_pylist()
            if nv is not None:
                values = [nv if v is None else v for v in values]
            max_len = col.spec.max_length
            too_long = [v for v in values if v is not None and len(v) > max_len]
            if too_long:
                raise ValueError(f"Column {col.name!r} contains values longer than max_length={max_len}.")
            return np.array(values, dtype=col.dtype)
        if arrow_col.null_count:
            if nv is None:
                raise TypeError(
                    f"Column {col.name!r} contains Arrow/Parquet nulls. Provide a CTable schema "
                    "with a null_value sentinel for this column."
                )
            arrow_col = arrow_col.fill_null(nv)
        return arrow_col.to_numpy(zero_copy_only=False).astype(col.dtype)

    @staticmethod
    def _arrow_schema_metadata(schema) -> dict[str, Any]:
        import base64

        try:
            schema_ipc = schema.serialize().to_pybytes()
            schema_ipc_base64 = base64.b64encode(schema_ipc).decode("ascii")
        except Exception:
            schema_ipc_base64 = None
        arrow_meta = {}
        if schema_ipc_base64 is not None:
            arrow_meta["schema_ipc_base64"] = schema_ipc_base64
        return {"arrow": arrow_meta}

    @staticmethod
    def _nested_metadata_from_column_names(
        column_names: list[str], *, empty_root_physical: str | None = None
    ) -> dict:
        logical_to_physical = {}
        for name in column_names:
            logical_to_physical[name] = name
        nested = {
            "version": 1,
            "logical_to_physical": logical_to_physical,
        }
        if empty_root_physical:
            logical_to_physical[""] = empty_root_physical
            nested["root"] = {"logical": "", "physical": empty_root_physical}
        return nested

    # ------------------------------------------------------------------
    # Unnamed-root list<struct<...>> detection and flattening helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _detect_unnamed_root_list_struct(pa, schema) -> bool:
        """Return True iff *schema* qualifies for unnamed-root list<struct<...>> flattening.

        Conditions (all must hold):
        * exactly one top-level field;
        * field name is ``""`` (the canonical unnamed Arrow root);
        * field type is ``list<struct<...>>`` or ``large_list<struct<...>>``.
        """
        if len(schema) != 1:
            return False
        field = schema[0]
        if field.name != "":
            return False
        t = field.type
        if not (pa.types.is_list(t) or pa.types.is_large_list(t)):
            return False
        return pa.types.is_struct(t.value_type)

    @staticmethod
    def _inner_schema_for_unnamed_root(pa, schema):
        """Extract the inner struct schema from a single unnamed root list<struct<...>> schema.

        Returns a new Arrow schema whose top-level fields are the struct fields
        of the list value type.  The nullable flag of the original unnamed field
        is not propagated — individual struct child nullability applies.
        """
        field = schema[0]  # the unnamed "" field
        struct_type = field.type.value_type  # struct type inside the list
        return pa.schema(list(struct_type))

    @staticmethod
    def _flatten_root_list_struct_batches(pa, inner_schema, batches, max_rows: int | None = None):
        """Yield flattened :class:`pyarrow.RecordBatch` objects from an unnamed root stream.

        For each incoming batch (which has a single list<struct<...>> column),
        flatten the outer list using ``ListArray.flatten()`` — which skips null
        outer list rows — and convert the resulting struct array into a
        :class:`~pyarrow.RecordBatch` whose columns correspond to the struct fields.

        Parameters
        ----------
        pa:
            The ``pyarrow`` module.
        inner_schema:
            Arrow schema for the inner struct (output of
            :meth:`_inner_schema_for_unnamed_root`).
        batches:
            Iterable of incoming :class:`~pyarrow.RecordBatch` objects from the
            unnamed-root Parquet file.
        max_rows:
            Optional maximum number of flattened element rows to yield.
        """
        rows_seen = 0
        for batch in batches:
            if max_rows is not None and rows_seen >= max_rows:
                break
            list_array = batch.column(0)
            # flatten() skips null outer list rows and concatenates element values
            struct_values = list_array.flatten()
            if len(struct_values) == 0:
                # Emit an empty record batch that still carries the inner schema
                empty_arrays = [pa.array([], type=f.type) for f in inner_schema]
                yield pa.record_batch(empty_arrays, schema=inner_schema)
                continue
            rb = pa.RecordBatch.from_struct_array(struct_values)
            if max_rows is not None:
                # Slice the *record batch* rather than the flattened struct array:
                # ListArray.flatten() can return an offset view, and some pyarrow
                # versions don't honor that offset in from_struct_array(), which
                # would silently import the untruncated batch.  RecordBatch.slice
                # is always respected downstream.
                remaining = max_rows - rows_seen
                if len(rb) > remaining:
                    rb = rb.slice(0, remaining)
            rows_seen += len(rb)
            yield rb

    @staticmethod
    def _flatten_arrow_struct_schema(pa, schema):
        """Flatten top-level struct fields into dotted leaf fields recursively."""

        out_fields = []

        def _walk(field, prefix: tuple[str, ...] = (), parent_nullable: bool = False):
            parts = (*prefix, field.name)
            name = join_field_path(parts)
            nullable = bool(parent_nullable or field.nullable)
            if pa.types.is_struct(field.type):
                for child in field.type:
                    _walk(pa.field(child.name, child.type, nullable=child.nullable), parts, nullable)
            else:
                out_fields.append(pa.field(name, field.type, nullable=nullable))

        for f in schema:
            _walk(f)
        return pa.schema(out_fields, metadata=schema.metadata)

    @staticmethod
    def _flatten_arrow_struct_batch(pa, batch, flat_schema):
        arrays = []

        def _extract(array, arr_type, parts):
            if not parts:
                return array
            head = parts[0]
            if pa.types.is_struct(arr_type):
                return _extract(array.field(head), arr_type[head].type, parts[1:])
            raise KeyError("Invalid flattened path")

        for field in flat_schema:
            parts = split_field_path(field.name)
            col = batch.column(batch.schema.get_field_index(parts[0]))
            arr = _extract(col, col.type, parts[1:])
            arrays.append(arr)
        return pa.RecordBatch.from_arrays(arrays, schema=flat_schema)

    @classmethod
    def _flatten_arrow_struct_input(cls, pa, schema, batches):
        """Return flattened (schema, batches, flattened) for struct-containing Arrow inputs."""
        if not any(pa.types.is_struct(f.type) for f in schema):
            return schema, batches, False
        flat_schema = cls._flatten_arrow_struct_schema(pa, schema)

        def _gen():
            for b in batches:
                yield cls._flatten_arrow_struct_batch(pa, b, flat_schema)

        return flat_schema, _gen(), True

    @classmethod
    def from_arrow(  # noqa: C901
        cls,
        schema,
        batches=None,
        *,
        urlpath: str | None = None,
        mode: str = "w",
        cparams=None,
        dparams=None,
        validate: bool = False,
        capacity_hint: int | None = None,
        string_max_length: int | Mapping[str, int] | None = None,
        auto_null_sentinels: bool = True,
        blosc2_batch_size: int | None = _BATCH_SIZE_DEFAULT,
        blosc2_items_per_block: int | None = None,
        list_serializer: Literal["msgpack", "arrow"] = "msgpack",
        object_fallback: bool = False,
        column_cparams: Mapping[str, dict[str, Any]] | None = None,
        separate_nested_cols: bool = False,
        create_summary_index: bool = True,
        chunks: int | tuple[int, ...] | None = None,
        blocks: int | tuple[int, ...] | None = None,
    ) -> CTable:
        """Build a :class:`CTable` from an Arrow schema and iterable of record batches.

        **Nested struct flattening**: top-level Arrow ``struct<…>`` fields are
        automatically and recursively flattened into dotted leaf columns.  For
        example, a field ``trip: struct<begin: struct<lon: float64, lat: float64>>``
        becomes two CTable columns ``trip.begin.lon`` and ``trip.begin.lat``.
        Each leaf is stored as an independent compressed :class:`~blosc2.NDArray`.
        Row reads via ``t[i]`` reconstruct the original nested dict shape.  Use
        ``t["trip.begin.lon"]`` or ``t.trip.begin.lon`` to access a leaf::

            import pyarrow as pa, blosc2
            trip_type = pa.struct([("begin", pa.struct([("lon", pa.float64())]))])
            schema = pa.schema([pa.field("trip", trip_type)])
            t = blosc2.CTable.from_arrow(schema, batches)
            t.col_names          # ['trip.begin.lon']
            t["trip.begin.lon"].mean()
            t.trip.begin.lon.max()

        When *string_max_length* is ``None`` (the default), scalar Arrow
        ``string`` / ``large_string`` columns are imported as variable-length
        :func:`~blosc2.utf8` columns (offsets + bytes storage; nullable
        columns get a null sentinel string from the active
        :class:`NullPolicy`) and ``binary`` / ``large_binary`` columns are
        imported as :func:`~blosc2.vlbytes` columns.  Non-struct ``struct``
        columns (not containing only scalar leaves) are imported as
        :func:`~blosc2.struct` columns backed by batched variable-length
        storage.  Null values for ``vlbytes``/``struct`` columns are
        represented as native ``None`` with no sentinel needed.

        When *string_max_length* is set to a positive integer, scalar string
        and binary columns are imported as fixed-width
        :func:`~blosc2.string` / :func:`~blosc2.bytes` columns whose dtype is
        sized to *string_max_length* characters/bytes. It may also be a mapping
        from column name to max length; omitted string/binary columns remain
        :func:`~blosc2.utf8` / :func:`~blosc2.vlbytes` columns.

        ``blosc2_batch_size`` controls how many rows are buffered before
        BatchArray-backed imported columns (list columns and variable-length
        scalar columns such as ``vlstring``, ``vlbytes``, ``struct``, and
        schema-less ``object`` columns) are flushed to their backend.  Set it to
        ``None`` to keep those columns pending until the final flush.

        ``list_serializer`` selects the backend serializer for imported list
        columns. ``"msgpack"`` is the default; ``"arrow"`` stores Arrow list
        batches directly and can be much faster for deeply nested list columns.

        Unsupported Arrow types raise by default.  Pass ``object_fallback=True``
        to import such columns as schema-less :func:`~blosc2.object` columns.
        This fallback is intentionally not used by :meth:`from_parquet`.

        ``column_cparams`` optionally maps column names to per-column compression
        parameters. These override the table-level ``cparams`` for fixed-width
        columns imported from Arrow.

        *schema* may also be any object implementing the Arrow PyCapsule
        interchange protocol (``__arrow_c_stream__``) — a pyarrow
        ``Table``/``RecordBatchReader``, a Polars ``DataFrame``, a DuckDB
        result, or another :class:`CTable` — in which case *batches* must be
        omitted and the schema/batches are pulled from the stream::

            t = blosc2.CTable.from_arrow(polars_df)
            t = blosc2.CTable.from_arrow(duckdb_result, urlpath="out.b2z")
        """
        pa = cls._require_pyarrow("from_arrow()")
        if hasattr(schema, "__arrow_c_stream__"):
            if batches is not None:
                raise TypeError(
                    "from_arrow() takes either a single Arrow-stream object "
                    "(implementing __arrow_c_stream__) or a (schema, batches) pair, not both."
                )
            if not hasattr(pa.RecordBatchReader, "from_stream"):
                raise RuntimeError("Importing from an Arrow PyCapsule stream requires pyarrow >= 14.0.")
            reader = pa.RecordBatchReader.from_stream(schema)
            schema, batches = reader.schema, reader
        elif batches is None:
            raise TypeError(
                "from_arrow() requires batches unless the first argument is an "
                "Arrow-stream object implementing __arrow_c_stream__."
            )
        if blosc2_batch_size is not None and blosc2_batch_size <= 0:
            raise ValueError("blosc2_batch_size must be a positive integer or None")
        if blosc2_items_per_block is not None and blosc2_items_per_block <= 0:
            raise ValueError("blosc2_items_per_block must be a positive integer or None")
        if list_serializer not in {"msgpack", "arrow"}:
            raise ValueError("list_serializer must be 'msgpack' or 'arrow'")

        # ------------------------------------------------------------------
        # Unnamed-root list<struct<...>> flattening (opt-in)
        # ------------------------------------------------------------------
        # When the source schema is a single unnamed "" field of type
        # list<struct<...>>, the outer list is a physical Parquet/Awkward
        # chunking artifact, not a semantic column.  Flatten it so that each
        # element becomes a CTable row.  The struct fields become ordinary
        # top-level columns and are further flattened by the struct-leaf
        # machinery below.
        original_root_metadata: dict | None = None
        if separate_nested_cols and cls._detect_unnamed_root_list_struct(pa, schema):
            inner_schema = cls._inner_schema_for_unnamed_root(pa, schema)
            batches = cls._flatten_root_list_struct_batches(pa, inner_schema, batches)
            schema = inner_schema
            original_root_metadata = {
                "kind": "unnamed_list_struct",
                "field_name": "",
                "preserve_grouping": False,
            }

        batches = iter(batches)
        first_batch = None
        table_for_inference = None
        original_top_level_struct_specs: dict[str, SchemaSpec] = {}
        for f in schema:
            if pa.types.is_struct(f.type):
                original_top_level_struct_specs[join_field_path((f.name,))] = cls._arrow_type_to_spec(
                    pa, f.type, nullable=f.nullable, object_fallback=object_fallback
                )
        if string_max_length is None or isinstance(string_max_length, Mapping):
            first_batch = next(batches, None)

        # Flatten top-level Arrow structs into dotted leaf columns so CTable can
        # persist nested scalar leaves as physical columns.
        flattened_structs = False
        if first_batch is not None:
            import itertools as _it

            schema, flat_batches, flattened_structs = cls._flatten_arrow_struct_input(
                pa, schema, _it.chain([first_batch], batches)
            )
            batches = iter(flat_batches)
            first_batch = next(batches, None)
        else:
            schema, batches, flattened_structs = cls._flatten_arrow_struct_input(pa, schema, batches)

        if first_batch is not None:
            table_for_inference = pa.Table.from_batches([first_batch], schema=schema)
        columns = cls._compiled_columns_from_arrow(
            pa,
            schema,
            table_for_inference,
            string_max_length,
            auto_null_sentinels=auto_null_sentinels,
            object_fallback=object_fallback,
        )
        cls._apply_arrow_column_cparams(columns, column_cparams)
        for col in columns:
            if cls._is_list_column(col):
                if getattr(col.spec, "storage", None) == "batch":
                    col.spec.serializer = list_serializer
                    if blosc2_batch_size is not None:
                        col.spec.batch_rows = blosc2_batch_size
                    if blosc2_items_per_block is not None:
                        col.spec.items_per_block = blosc2_items_per_block
            elif cls._is_varlen_scalar_column(col):
                if blosc2_batch_size is not None:
                    col.spec.batch_rows = blosc2_batch_size
                if blosc2_items_per_block is not None:
                    col.spec.items_per_block = blosc2_items_per_block
        metadata = cls._arrow_schema_metadata(schema)
        empty_root_physical = None
        schema_meta = getattr(schema, "metadata", None) or {}
        root_key = b"blosc2_empty_root_physical"
        if root_key in schema_meta:
            raw = schema_meta[root_key]
            empty_root_physical = raw.decode() if isinstance(raw, bytes) else str(raw)
        metadata["nested"] = cls._nested_metadata_from_column_names(
            [col.name for col in columns], empty_root_physical=empty_root_physical
        )
        if flattened_structs:
            metadata["nested"]["reconstruct_rows"] = True
        compiled_columns_by_name = {col.name: col for col in columns}
        for name, spec in original_top_level_struct_specs.items():
            if name in compiled_columns_by_name:
                continue
            compiled_columns_by_name[name] = CompiledColumn(
                name=name,
                py_type=spec.python_type,
                spec=spec,
                dtype=getattr(spec, "dtype", None),
                default=MISSING,
                config=ColumnConfig(cparams=None, dparams=None, chunks=None, blocks=None),
                display_width=compute_display_width(spec),
            )

        compiled = CompiledSchema(
            row_cls=None,
            columns=columns,
            columns_by_name=compiled_columns_by_name,
            metadata=metadata,
        )
        if first_batch is not None:
            import itertools as _it

            batches = _it.chain([first_batch], batches)
        # Use capacity_hint to size initial NDArray chunks/blocks correctly.
        # When capacity_hint is None and we are in the unnamed-root flatten path,
        # fall back to _EXPECTED_SIZE_DEFAULT (1 M) so that compute_chunks_blocks
        # produces a reasonable block size instead of (1,) which causes catastrophic
        # storage fragmentation.  For non-unnamed-root imports capacity_hint is
        # always supplied by from_parquet (pf.metadata.num_rows), so the fallback
        # only matters for direct from_arrow() calls without a hint.
        if capacity_hint is None and original_root_metadata is not None:
            capacity = _EXPECTED_SIZE_DEFAULT
        else:
            capacity = max(capacity_hint or 1, 1)
        _chunks = (chunks,) if isinstance(chunks, int) else chunks
        _blocks = (blocks,) if isinstance(blocks, int) else blocks
        storage = cls._storage_for_arrow_import(urlpath, mode)
        new_cols, new_valid = cls._create_arrow_import_columns(
            storage, columns, capacity, cparams, dparams, _chunks, _blocks
        )
        storage.save_schema(schema_to_dict(compiled))
        obj = cls._new_arrow_import_ctable(
            compiled,
            storage,
            new_cols,
            new_valid,
            columns,
            cparams=cparams,
            dparams=dparams,
            validate=validate,
            create_summary_index=create_summary_index,
        )
        cls._write_arrow_batches(obj, batches, columns, new_cols, new_valid)
        return obj

    def to_parquet(
        self,
        path,
        *,
        columns: list[str] | None = None,
        batch_size: int = _BATCH_SIZE_DEFAULT,
        compression: str | None = "zstd",
        row_group_size: int | None = None,
        include_computed: bool = True,
        **kwargs,
    ) -> None:
        """Write this table to a Parquet file batch-wise using pyarrow."""
        pq = self._require_pyarrow_parquet("to_parquet()")
        pa = self._require_pyarrow("to_parquet()")
        self._validate_arrow_batch_size(batch_size)
        schema = self._arrow_schema_for_columns(columns, include_computed=include_computed)
        with pq.ParquetWriter(path, schema, compression=compression, **kwargs) as writer:
            for batch in self.iter_arrow_batches(
                columns=columns, batch_size=batch_size, include_computed=include_computed
            ):
                table = pa.Table.from_batches([batch], schema=batch.schema)
                writer.write_table(table, row_group_size=row_group_size or len(batch))

    @classmethod
    def from_parquet(  # noqa: C901
        cls,
        path,
        *,
        columns: list[str] | None = None,
        batch_size: int = _BATCH_SIZE_DEFAULT,
        urlpath: str | None = None,
        mode: str = "w",
        cparams=None,
        dparams=None,
        validate: bool = False,
        auto_null_sentinels: bool = True,
        blosc2_batch_size: int | None = _BATCH_SIZE_DEFAULT,
        blosc2_items_per_block: int | None = None,
        list_serializer: Literal["msgpack", "arrow"] = "arrow",
        separate_nested_cols: bool = True,
        max_rows: int | None = None,
        **kwargs,
    ) -> CTable:
        """Read a Parquet file into a :class:`CTable`.

        The Parquet file is streamed batch by batch through :mod:`pyarrow` and then
        converted into a typed :class:`CTable`. By default, the result is created in
        memory, but you can also persist it on disk via ``urlpath``.

        This method delegates the actual table construction to
        :meth:`CTable.from_arrow`, so Arrow schema handling, nullable-column support,
        and Blosc2 write tuning follow the same rules as that method.

        **Nested struct flattening**: top-level Parquet ``struct<…>`` fields are
        automatically and recursively flattened into dotted leaf columns — the same
        as in :meth:`from_arrow`.  For example, a Parquet file that contains a column
        ``trip: struct<begin: struct<lon: double, lat: double>>`` produces two CTable
        columns ``trip.begin.lon`` and ``trip.begin.lat``.  Row reads reconstruct the
        original nested dict shape; individual leaves are accessed via dotted names or
        attribute-chain proxies::

            t = blosc2.CTable.from_parquet("trips.parquet")
            t.col_names               # e.g. ['trip.begin.lon', 'trip.begin.lat', ...]
            t["trip.begin.lon"].mean()
            t.trip.begin.lon.max()

        Unsupported Parquet types are not silently imported as schema-less
        :func:`~blosc2.object` columns; they raise so callers can decide how to
        handle them explicitly.

        Parameters
        ----------
        path : str or path-like
            Path to the source Parquet file.

        columns : list[str] or None, optional
            Subset of columns to read from the Parquet file. If provided, only these
            columns are loaded and their order in the resulting table matches the
            order in this list. Column names must be unique.

        batch_size : int, optional
            Number of rows per Arrow batch read from the Parquet file. This controls
            how much data is pulled from the file at a time before being handed off
            to the CTable builder. Must be greater than 0.

        urlpath : str or None, optional
            Destination storage path for the resulting CTable. If ``None`` (the
            default), the table is created in memory. If provided, the table is backed
            by persistent on-disk storage.

        mode : str, optional
            Storage open mode for ``urlpath``. Defaults to ``"w"``. This is passed
            through to :meth:`CTable.from_arrow`.

        cparams : object, optional
            Compression parameters for the created Blosc2 containers. Passed through
            to :meth:`CTable.from_arrow`.

        dparams : object, optional
            Decompression parameters for the created Blosc2 containers. Passed through
            to :meth:`CTable.from_arrow`.

        validate : bool, optional
            Whether to enable extra internal validation while building the table.
            Defaults to ``False``.

        auto_null_sentinels : bool, optional
            If ``True`` (default), nullable scalar columns imported from Parquet may
            automatically receive per-column null sentinel values when needed. Sentinel
            selection follows the current null-policy rules used by CTable schema
            handling.

        blosc2_batch_size : int or None, optional
            Number of items written to Blosc2 containers per internal write batch.
            Passed through to :meth:`CTable.from_arrow`.

        blosc2_items_per_block : int or None, optional
            Target number of items per internal Blosc2 block. Passed through to
            :meth:`CTable.from_arrow`.  In general, larger number of items
            favors compression ratios but make random access slower.

        list_serializer : {"msgpack", "arrow"}, optional
            Serializer used for imported list columns. The default, ``"arrow"``,
            stores Arrow list batches directly and is much faster for deeply nested
            or ``list<struct<...>>`` columns. The tradeoff is that accessing those
            list columns later requires PyArrow. Use ``"msgpack"`` to keep
            list-column stores independent of PyArrow at read time; it can be
            smaller for simple lists but is much slower and more memory-intensive
            for deeply nested data.

        separate_nested_cols : bool, optional
            Whether to separate qualifying nested columns during import. Defaults to
            ``True``. In particular, a single unnamed top-level
            ``list<struct<...>>`` field is treated as a root record stream: each list
            element becomes a CTable row and struct leaves become ordinary nested
            CTable columns. Use ``separate_nested_cols=False`` when closer fidelity to
            the original Parquet row/schema shape is more important than the separated
            column layout.

        max_rows : int or None, optional
            Maximum number of rows to import. For ordinary Parquet files this limits
            Parquet/CTable rows. For unnamed-root ``list<struct<...>>`` files imported
            with ``separate_nested_cols=True``, this limits flattened element rows.

        **kwargs
            Additional keyword arguments forwarded to ``pyarrow.parquet.ParquetFile``.
            Use these for Parquet-reader-specific options supported by PyArrow.

        Returns
        -------
        CTable
            A new :class:`CTable` populated from the Parquet file. The table contains
            all selected columns and all rows from the file. If ``urlpath`` is
            provided, the returned table is disk-backed; otherwise it is in-memory.

        Raises
        ------
        ImportError
            If :mod:`pyarrow` is not installed.
        ValueError
            If ``batch_size`` is not greater than 0.
        ValueError
            If ``max_rows`` is negative.
        ValueError
            If ``columns`` contains duplicate names.
        Exception
            Any exception raised by :mod:`pyarrow` while opening or reading the Parquet
            file, or by :meth:`CTable.from_arrow` while converting Arrow data into a
            CTable.

        Examples
        --------
        Load an entire Parquet file into an in-memory table:

        >>> import blosc2
        >>> t = blosc2.CTable.from_parquet("data.parquet")

        Load only a subset of columns:

        >>> t = blosc2.CTable.from_parquet(
        ...     "data.parquet",
        ...     columns=["user_id", "amount", "country"],
        ... )

        Create a disk-backed table while reading in batches:

        >>> t = blosc2.CTable.from_parquet(
        ...     "data.parquet",
        ...     batch_size=50_000,
        ...     urlpath="data.ctable",
        ... )

        Pass additional options through to PyArrow's Parquet reader:

        >>> t = blosc2.CTable.from_parquet(
        ...     "data.parquet",
        ...     memory_map=True,
        ... )
        """
        pq = cls._require_pyarrow_parquet("from_parquet()")
        pa = cls._require_pyarrow("from_parquet()")
        cls._validate_arrow_batch_size(batch_size)
        if max_rows is not None and max_rows < 0:
            raise ValueError("max_rows must be non-negative")
        string_max_length = kwargs.pop("string_max_length", None)
        pf = pq.ParquetFile(path, **kwargs)
        arrow_schema = pf.schema_arrow
        if columns is not None:
            if len(set(columns)) != len(columns):
                raise ValueError("columns must be unique")
            fields = [arrow_schema.field(name) for name in columns]
            arrow_schema = pa.schema(fields)
        batches = pf.iter_batches(batch_size=batch_size, columns=columns)

        # Parquet files generated by Awkward-style pipelines may contain an
        # unnamed top-level field (""). When separate_nested_cols=True and the
        # schema qualifies as an unnamed-root list<struct<...>>, skip the
        # rename-to-root logic and pass the original schema directly to
        # from_arrow, which will perform the element-level flattening.
        # Otherwise, normalize empty column names to non-empty names as before.
        _is_unnamed_root_flatten = separate_nested_cols and cls._detect_unnamed_root_list_struct(
            pa, arrow_schema
        )
        if not _is_unnamed_root_flatten and any(name == "" for name in arrow_schema.names):
            used = {n for n in arrow_schema.names if n}

            def _fresh_root_name() -> str:
                base = "root"
                if base not in used:
                    used.add(base)
                    return base
                i = 1
                while True:
                    candidate = f"{base}_{i}"
                    if candidate not in used:
                        used.add(candidate)
                        return candidate
                    i += 1

            original_names = list(arrow_schema.names)
            renamed = [_fresh_root_name() if n == "" else n for n in original_names]
            arrow_schema = pa.schema(
                [arrow_schema.field(i).with_name(renamed[i]) for i in range(len(renamed))]
            )
            # Preserve canonical unnamed-root intent in schema metadata.
            try:
                first_root = next(renamed[i] for i, old in enumerate(original_names) if old == "")
            except StopIteration:
                first_root = renamed[0] if renamed else "root"
            current_meta = dict(arrow_schema.metadata or {})
            current_meta[b"blosc2_empty_root_physical"] = first_root.encode()
            arrow_schema = arrow_schema.with_metadata(current_meta)

            def _renamed_batches(batch_iter, names):
                for b in batch_iter:
                    yield b.rename_columns(names)

            batches = _renamed_batches(batches, renamed)

        def _limited_batches(batch_iter, limit: int):
            rows_seen = 0
            for batch in batch_iter:
                if rows_seen >= limit:
                    break
                remaining = limit - rows_seen
                if len(batch) > remaining:
                    batch = batch.slice(0, remaining)
                rows_seen += len(batch)
                yield batch

        # For unnamed-root flattening, max_rows applies to flattened element rows,
        # not to the outer Parquet rows.  Pre-flatten here when a limit is requested
        # so the limit can be enforced precisely before handing batches to from_arrow.
        if _is_unnamed_root_flatten and max_rows is not None:
            inner_schema = cls._inner_schema_for_unnamed_root(pa, arrow_schema)
            limited_flat_batches = cls._flatten_root_list_struct_batches(
                pa, inner_schema, batches, max_rows=max_rows
            )
            ct = cls.from_arrow(
                inner_schema,
                limited_flat_batches,
                urlpath=urlpath,
                mode=mode,
                cparams=cparams,
                dparams=dparams,
                validate=validate,
                capacity_hint=max_rows,
                string_max_length=string_max_length,
                auto_null_sentinels=auto_null_sentinels,
                blosc2_batch_size=blosc2_batch_size,
                blosc2_items_per_block=blosc2_items_per_block,
                list_serializer=list_serializer,
                separate_nested_cols=False,
            )
            ct._storage.save_schema(schema_to_dict(ct._schema))
            return ct

        if max_rows is not None:
            batches = _limited_batches(batches, max_rows)

        # When flattening a root list<struct<...>>, the actual element count is not
        # known ahead of time.  Pass capacity_hint=None so that from_arrow falls back
        # to _EXPECTED_SIZE_DEFAULT (1 M), which gives compute_chunks_blocks() a
        # reasonable block size instead of the catastrophic (1, 1) produced by
        # capacity=1.  The CLI path computes a better estimate by sampling.
        if _is_unnamed_root_flatten:
            _capacity_hint = None
        elif pf.metadata is not None:
            _capacity_hint = (
                pf.metadata.num_rows if max_rows is None else min(max_rows, pf.metadata.num_rows)
            )
        else:
            _capacity_hint = max_rows

        return cls.from_arrow(
            arrow_schema,
            batches,
            urlpath=urlpath,
            mode=mode,
            cparams=cparams,
            dparams=dparams,
            validate=validate,
            capacity_hint=_capacity_hint,
            string_max_length=string_max_length,
            auto_null_sentinels=auto_null_sentinels,
            blosc2_batch_size=blosc2_batch_size,
            blosc2_items_per_block=blosc2_items_per_block,
            list_serializer=list_serializer,
            separate_nested_cols=separate_nested_cols,
        )

    # ------------------------------------------------------------------
    # CSV interop
    # ------------------------------------------------------------------

    def to_csv(self, path: str | None = None, *, header: bool = True, sep: str = ",") -> str | None:
        """Write all live rows to CSV.

        Uses Python's stdlib ``csv`` module — no extra dependency required.
        Fixed-shape ndarray column cells are serialised as JSON arrays for
        readability and shape safety (e.g. ``"[1.0, 2.0, 3.0]"``).

        Parameters
        ----------
        path:
            Destination file path (created or overwritten).  If ``None`` (the
            default), nothing is written and the CSV is returned as a string,
            like ``pandas``' ``DataFrame.to_csv()``.
        header:
            If ``True`` (default), write column names as the first row.
        sep:
            Field delimiter.  Defaults to ``","``; use ``"\\t"`` for TSV.

        Returns
        -------
        str or None
            The CSV text when *path* is ``None``, otherwise ``None``.
        """
        import csv
        import io

        n = len(self)
        arrays: list = []
        for name in self.col_names:
            col = self[name]
            if col.is_ndarray:
                arr = col[:]
                null_mask = col._null_mask_for(arr)
                json_strings: list[str] = []
                for i in range(n):
                    if null_mask[i]:
                        json_strings.append("")
                    else:
                        json_strings.append(json.dumps(arr[i].tolist()))
                arrays.append(json_strings)
            else:
                arrays.append(col[:])

        def _write(f) -> None:
            writer = csv.writer(f, delimiter=sep)
            if header:
                writer.writerow(self.col_names)
            for row in zip(*arrays, strict=True):
                writer.writerow(row)

        if path is None:
            buf = io.StringIO(newline="")
            _write(buf)
            return buf.getvalue()
        with open(path, "w", newline="") as f:
            _write(f)
        return None

    @staticmethod
    def _csv_ndarray_col_to_array(raw: list[str], col) -> np.ndarray:
        """Convert a list of JSON-array CSV strings to a stacked ndarray for an ndarray column."""
        spec = col.spec
        null_value = getattr(spec, "null_value", None)
        item_shape = spec.item_shape
        dtype = spec.dtype

        rows = []
        for val in raw:
            stripped = val.strip()
            if stripped == "":
                if null_value is not None:
                    rows.append(np.full(item_shape, null_value, dtype=dtype))
                    continue
                raise ValueError(f"Column {col.name!r}: non-nullable column got empty cell")

            try:
                arr = np.array(json.loads(stripped), dtype=dtype)
            except json.JSONDecodeError as exc:
                raise ValueError(f"Column {col.name!r}: invalid JSON array cell {val!r}") from exc

            if arr.shape != item_shape:
                raise ValueError(f"Column {col.name!r}: expected item shape {item_shape}, got {arr.shape}")
            rows.append(arr)

        return np.ascontiguousarray(rows, dtype=dtype)

    @staticmethod
    def _csv_col_to_array(raw: list[str], col, nv) -> np.ndarray:
        """Convert a list of raw CSV strings to a numpy array for *col*."""
        if col.dtype == np.bool_:

            def _parse(v, _nv=nv):
                stripped = v.strip()
                if stripped == "" and _nv is not None:
                    return _nv
                return stripped in ("True", "true", "1")

            return np.array([_parse(v) for v in raw], dtype=np.bool_)
        if col.dtype.kind == "S":
            prepared: list = [nv if (v.strip() == "" and nv is not None) else v.encode() for v in raw]
            return np.array(prepared, dtype=col.dtype)
        prepared2 = [nv if (v.strip() == "" and nv is not None) else v for v in raw]
        return np.array(prepared2, dtype=col.dtype)

    @classmethod
    def from_csv(
        cls,
        path: str,
        row_cls,
        *,
        header: bool = True,
        sep: str = ",",
    ) -> CTable:
        """Build a :class:`CTable` from a CSV file.

        Schema comes from *row_cls* (a dataclass) — CTable is always typed.
        All rows are read in a single pass into per-column Python lists, then
        each column is bulk-written into a pre-allocated NDArray (one slice
        assignment per column, no ``extend()``).

        Parameters
        ----------
        path:
            Source CSV file path.
        row_cls:
            A dataclass whose fields define the column names and types.
        header:
            If ``True`` (default), the first row is treated as a header and
            skipped.  Column order in the file must match *row_cls* field
            order regardless.
        sep:
            Field delimiter.  Defaults to ``","``; use ``"\\t"`` for TSV.

        Returns
        -------
        CTable
            A new in-memory CTable containing all rows from the CSV file.

        Raises
        ------
        TypeError
            If *row_cls* is not a dataclass.
        ValueError
            If a row has a different number of fields than the schema.
        """
        import csv

        schema = compile_schema(row_cls)
        cls._resolve_nullable_specs(schema)
        ncols = len(schema.columns)

        # Accumulate values per column as Python lists (one pass through file)
        col_data: list[list] = [[] for _ in range(ncols)]

        with open(path, newline="") as f:
            reader = csv.reader(f, delimiter=sep)
            if header:
                next(reader)
            for lineno, row in enumerate(reader, start=2 if header else 1):
                if len(row) != ncols:
                    raise ValueError(f"Line {lineno}: expected {ncols} fields, got {len(row)}.")
                for i, val in enumerate(row):
                    col_data[i].append(val)

        n = len(col_data[0]) if ncols > 0 else 0
        capacity = max(n, 1)
        default_chunks, default_blocks = compute_chunks_blocks((capacity,))
        # Align fixed-size scalar columns (and the _valid_rows mask) on one
        # shared grid so lazy expressions over them take the fast_eval path.
        shared_chunks, shared_blocks, aligned_names = cls._compute_aligned_grid(schema.columns, capacity)
        mem_storage = InMemoryTableStorage()

        new_valid = mem_storage.create_valid_rows(
            shape=(capacity,),
            chunks=shared_chunks if shared_chunks is not None else default_chunks,
            blocks=shared_blocks if shared_blocks is not None else default_blocks,
        )
        new_cols: dict[str, blosc2.NDArray] = {}
        for col in schema.columns:
            shape = cls._column_physical_shape(col, capacity)
            if col.name in aligned_names:
                chunks, blocks = shared_chunks, shared_blocks
            else:
                chunks, blocks = cls._column_chunks_blocks(col, shape)
            new_cols[col.name] = mem_storage.create_column(
                col.name,
                dtype=col.dtype,
                shape=shape,
                chunks=chunks,
                blocks=blocks,
                cparams=None,
                dparams=None,
            )

        obj = cls.__new__(cls)
        obj._row_type = row_cls
        obj._validate = True
        obj._table_cparams = None
        obj._table_dparams = None
        obj._storage = mem_storage
        obj._read_only = False
        obj._schema = schema
        obj._cols = new_cols
        obj._col_widths = {col.name: max(len(col.name), col.display_width) for col in schema.columns}
        obj.col_names = [col.name for col in schema.columns]
        obj.auto_compact = False
        obj._create_summary_index = True
        obj._summary_indexes_built = False
        obj.base = None
        obj._computed_cols = {}  # from_csv creates no computed columns
        obj._materialized_cols = {}
        obj._expr_index_arrays = {}
        obj._valid_rows = new_valid
        obj._n_rows = 0
        obj._last_pos = 0

        if n > 0:
            for i, col in enumerate(schema.columns):
                if isinstance(col.spec, NDArraySpec):
                    arr = cls._csv_ndarray_col_to_array(col_data[i], col)
                else:
                    nv = getattr(col.spec, "null_value", None)
                    arr = cls._csv_col_to_array(col_data[i], col, nv)
                new_cols[col.name][:n] = arr
            new_valid[:n] = True
            obj._n_rows = n
            obj._last_pos = n

        return obj

    # ------------------------------------------------------------------
    # Pandas / DataFrame interop
    # ------------------------------------------------------------------

    def to_pandas(self):
        """Convert to a `pandas <https://pandas.pydata.org>`_ DataFrame.

        Scalar columns become regular DataFrame columns.  Fixed-shape ndarray
        columns become ``object``-dtype columns whose cells hold NumPy arrays
        of per-row shape *item_shape*.

        Returns
        -------
        pandas.DataFrame

        Examples
        --------
        >>> import blosc2
        >>> from dataclasses import dataclass
        >>> import numpy as np
        >>> @dataclass
        ... class Row:
        ...     id: int = blosc2.field(blosc2.int64())
        ...     embedding: object = blosc2.field(blosc2.ndarray((3,), dtype=blosc2.float32()))
        >>> t = blosc2.CTable(Row, new_data=[
        ...     (1, np.array([1, 2, 3], dtype=np.float32)),
        ...     (2, np.array([4, 5, 6], dtype=np.float32)),
        ... ])
        >>> df = t.to_pandas()
        >>> df["id"].tolist()
        [1, 2]
        >>> df["embedding"].dtype
        dtype('O')
        >>> np.testing.assert_array_equal(df["embedding"][0], np.array([1, 2, 3], dtype=np.float32))
        """
        import pandas as pd

        data = {}
        for name in self.col_names:
            col = self[name]
            if col.is_ndarray:
                data[name] = list(col)
            else:
                data[name] = col[:]

        return pd.DataFrame(data)

    @classmethod
    def from_pandas(cls, df, row_cls) -> CTable:  # noqa: C901
        """Build a :class:`CTable` from a pandas DataFrame.

        Schema comes from *row_cls* (a dataclass) — CTable is always typed.
        Object-dtype DataFrame columns are **not** automatically inferred as
        ndarray columns; the *row_cls* must explicitly declare
        :func:`blosc2.ndarray` fields.

        Parameters
        ----------
        df:
            Source pandas DataFrame.
        row_cls:
            A dataclass whose fields define the column names and types.

        Returns
        -------
        CTable
            A new CTable containing all DataFrame rows.

        Raises
        ------
        TypeError
            If *row_cls* is not a dataclass.
        ValueError
            If DataFrame columns do not match the *row_cls* schema.
        """
        schema = compile_schema(row_cls)
        cls._resolve_nullable_specs(schema)

        # Validate column names
        schema_names = [col.name for col in schema.columns]
        missing = [name for name in schema_names if name not in df.columns]
        if missing:
            raise ValueError(f"DataFrame missing columns declared in row_cls: {missing}")
        extra = [name for name in df.columns if name not in schema_names]
        if extra:
            raise ValueError(f"DataFrame has extra columns not in row_cls: {extra}")

        n = len(df)
        capacity = max(n, 1)
        default_chunks, default_blocks = compute_chunks_blocks((capacity,))
        # Align fixed-size scalar columns (and the _valid_rows mask) on one
        # shared grid so lazy expressions over them take the fast_eval path.
        shared_chunks, shared_blocks, aligned_names = cls._compute_aligned_grid(schema.columns, capacity)
        mem_storage = InMemoryTableStorage()

        new_valid = mem_storage.create_valid_rows(
            shape=(capacity,),
            chunks=shared_chunks if shared_chunks is not None else default_chunks,
            blocks=shared_blocks if shared_blocks is not None else default_blocks,
        )
        new_cols: dict[str, Any] = {}
        for col in schema.columns:
            if cls._is_list_column(col):
                new_cols[col.name] = mem_storage.create_list_column(
                    col.name,
                    spec=col.spec,
                    cparams=None,
                    dparams=None,
                )
                continue
            if cls._is_varlen_scalar_column(col):
                new_cols[col.name] = mem_storage.create_varlen_scalar_column(
                    col.name,
                    spec=col.spec,
                    cparams=None,
                    dparams=None,
                )
                continue
            if cls._is_dictionary_column(col):
                dict_col = mem_storage.create_dictionary_column(
                    col.name,
                    spec=col.spec,
                    cparams=None,
                    dparams=None,
                )
                if len(dict_col.codes) < capacity:
                    dict_col.resize((capacity,))
                new_cols[col.name] = dict_col
                continue
            shape = cls._column_physical_shape(col, capacity)
            if col.name in aligned_names:
                chunks, blocks = shared_chunks, shared_blocks
            else:
                chunks, blocks = cls._column_chunks_blocks(col, shape)
            new_cols[col.name] = mem_storage.create_column(
                col.name,
                dtype=col.dtype,
                shape=shape,
                chunks=chunks,
                blocks=blocks,
                cparams=None,
                dparams=None,
            )

        obj = cls.__new__(cls)
        obj._row_type = row_cls
        obj._validate = True
        obj._table_cparams = None
        obj._table_dparams = None
        obj._storage = mem_storage
        obj._read_only = False
        obj._schema = schema
        obj._cols = new_cols
        obj._col_widths = {col.name: max(len(col.name), col.display_width) for col in schema.columns}
        obj.col_names = [col.name for col in schema.columns]
        obj.auto_compact = False
        obj._create_summary_index = True
        obj._summary_indexes_built = False
        obj.base = None
        obj._computed_cols = {}
        obj._materialized_cols = {}
        obj._expr_index_arrays = {}
        obj._valid_rows = new_valid
        obj._n_rows = 0
        obj._last_pos = 0

        if n > 0:

            def normalize_pandas_missing(value):
                if value is None:
                    return None
                if isinstance(value, float) and np.isnan(value):
                    return None
                # pandas.NA cannot be compared/coerced reliably; detect it by type name
                # without importing pandas here.
                if type(value).__name__ == "NAType":
                    return None
                return value

            raw_columns = {}
            for col in schema.columns:
                series = df[col.name]
                if isinstance(col.spec, NDArraySpec) and series.values.dtype != object:
                    raise ValueError(
                        f"Column {col.name!r}: expected object dtype in DataFrame "
                        f"for ndarray column, got {series.values.dtype}"
                    )
                if (
                    cls._is_list_column(col)
                    or cls._is_varlen_scalar_column(col)
                    or cls._is_dictionary_column(col)
                    or isinstance(col.spec, NDArraySpec)
                ):
                    raw_columns[col.name] = [normalize_pandas_missing(value) for value in series.tolist()]
                else:
                    raw_columns[col.name] = series.to_numpy(dtype=col.dtype)
            obj.extend(raw_columns, validate=True)

        return obj

    # ------------------------------------------------------------------
    # Schema mutations: add / drop / rename columns
    # ------------------------------------------------------------------

    @staticmethod
    def _column_spec_default_and_config(
        spec_or_field: SchemaSpec | dataclasses.Field,
    ) -> tuple[SchemaSpec, Any, ColumnConfig]:
        """Extract the schema spec, default and storage config for ``add_column()``."""
        if isinstance(spec_or_field, dataclasses.Field):
            meta = get_blosc2_field_metadata(spec_or_field)
            if meta is None:
                raise TypeError("add_column() field descriptors must be created with blosc2.field().")
            spec = copy.deepcopy(meta["spec"])
            if spec_or_field.default is not MISSING:
                default = spec_or_field.default
            elif spec_or_field.default_factory is not MISSING:  # type: ignore[misc]
                default = spec_or_field.default_factory()
            else:
                default = MISSING
            config = ColumnConfig(
                cparams=meta.get("cparams"),
                dparams=meta.get("dparams"),
                chunks=meta.get("chunks"),
                blocks=meta.get("blocks"),
            )
        else:
            spec = spec_or_field
            default = MISSING
            config = ColumnConfig(cparams=None, dparams=None, chunks=None, blocks=None)

        if not isinstance(spec, SchemaSpec):
            raise TypeError(f"add_column() requires a SchemaSpec, got {type(spec)!r}.")
        return spec, default, config

    def add_column(  # noqa: C901
        self,
        name: str,
        spec: SchemaSpec | dataclasses.Field,
    ) -> None:
        """Add a new column filled from the default declared in *spec*.

        Parameters
        ----------
        name:
            Column name.  Must follow the same naming rules as schema fields.
        spec:
            A schema descriptor such as ``b2.int64(ge=0)`` or a field
            descriptor such as ``b2.field(b2.int64(ge=0), default=0)``.
            When the table already has live rows, use ``blosc2.field(...)``
            with a default declared so those rows can be backfilled.

        Raises
        ------
        ValueError
            If the table is read-only, is a view, the column already exists,
            or a non-empty table is given a column with no default declared.
        TypeError
            If a declared default cannot be coerced to *spec*'s dtype.
        """
        if self._read_only:
            raise ValueError("Table is read-only (opened with mode='r').")
        if self.base is not None:
            raise ValueError("Cannot add a column to a view.")
        _validate_column_name(name)
        if name in self._cols:
            raise ValueError(f"Column {name!r} already exists.")
        if name in self._computed_cols:
            raise ValueError(f"A computed column named {name!r} already exists.")

        spec, default, column_config = self._column_spec_default_and_config(spec)

        n_live = self.nrows
        if default is MISSING and n_live > 0:
            raise ValueError(
                "add_column() requires a default declared as blosc2.field(..., default=...) "
                "when the table has live rows."
            )

        compiled_col = self._compiled_column_from_spec(name, spec)
        compiled_col.config = column_config
        self._resolve_nullable_specs(
            CompiledSchema(row_cls=None, columns=[compiled_col], columns_by_name={name: compiled_col}),
            validate_column_null_values=False,
        )
        spec = compiled_col.spec

        if self._is_varlen_scalar_column(compiled_col):
            # Varlen scalar columns don't use fixed-width NDArray storage.
            col_storage = self._resolve_column_storage(compiled_col, None, None)
            new_col = self._storage.create_varlen_scalar_column(
                name,
                spec=spec,
                cparams=col_storage.get("cparams"),
                dparams=col_storage.get("dparams"),
            )
            for _ in range(n_live):
                new_col.append(default)
            new_col.flush()
        elif self._is_list_column(compiled_col):
            raise TypeError(
                "add_column() does not support list columns; use the constructor with a full schema."
            )
        else:
            if default is not MISSING:
                try:
                    if self._is_ndarray_column(compiled_col):
                        default_val = self._coerce_ndarray_value(name, spec, default)
                    else:
                        default_val = spec.dtype.type(default)
                except (ValueError, OverflowError) as exc:
                    raise TypeError(
                        f"Cannot coerce default {default!r} to dtype {spec.dtype!r}: {exc}"
                    ) from exc
            else:
                default_val = None

            capacity = len(self._valid_rows)
            shape = self._column_physical_shape(compiled_col, capacity)
            default_chunks, default_blocks = self._column_chunks_blocks(compiled_col, shape)
            col_storage = self._resolve_column_storage(compiled_col, default_chunks, default_blocks)
            new_col = self._storage.create_column(
                name,
                dtype=spec.dtype,
                shape=shape,
                chunks=col_storage["chunks"],
                blocks=col_storage["blocks"],
                cparams=col_storage.get("cparams"),
                dparams=col_storage.get("dparams"),
            )
            if n_live > 0:
                if self._is_ndarray_column(compiled_col):
                    new_col[self._valid_rows] = np.broadcast_to(default_val, (n_live, *spec.item_shape))
                else:
                    new_col[self._valid_rows] = default_val

        compiled_col.default = default
        self._cols[name] = new_col
        self.col_names.append(name)
        self._col_widths[name] = max(len(name), compiled_col.display_width)

        new_columns = self._schema.columns + [compiled_col]
        self._schema = CompiledSchema(
            row_cls=self._schema.row_cls,
            columns=new_columns,
            columns_by_name={**self._schema.columns_by_name, name: compiled_col},
        )
        if isinstance(self._storage, FileTableStorage):
            self._storage.save_schema(self._schema_dict_with_computed())

    def drop_column(self, name: str) -> None:
        """Remove a column from the table.

        On disk tables the corresponding persisted column leaf is deleted.

        Raises
        ------
        ValueError
            If the table is read-only, is a view, or *name* is the last column.
        KeyError
            If *name* does not exist.
        """
        if self._read_only:
            raise ValueError("Table is read-only (opened with mode='r').")
        if self.base is not None:
            raise ValueError("Cannot drop a column from a view.")
        if name not in self._cols:
            raise KeyError(f"No column named {name!r}. Available: {self.col_names}")
        if len(self._stored_col_names) == 1:
            raise ValueError("Cannot drop the last column.")
        # Guard: refuse if any computed column depends on this column
        dependents = [cc_name for cc_name, cc in self._computed_cols.items() if name in cc["col_deps"]]
        dependents.extend(
            mat_name
            for mat_name, meta in self._materialized_cols.items()
            if name in meta.get("col_deps", ())
        )
        if dependents:
            raise ValueError(
                f"Cannot drop column {name!r}: it is used by computed/generated column(s) "
                + ", ".join(repr(d) for d in dependents)
                + ". Drop those columns first."
            )

        catalog = self._get_index_catalog()
        if name in catalog:
            descriptor = catalog.pop(name)
            self._validate_index_descriptor(name, descriptor)
            self._drop_index_descriptor(name, descriptor)
            self._storage.save_index_catalog(catalog)
            self._invalidate_index_catalog_cache()

        if isinstance(self._storage, FileTableStorage):
            self._storage.delete_column(name)

        self._materialized_cols.pop(name, None)
        del self._cols[name]
        del self._col_widths[name]
        self.col_names.remove(name)

        new_columns = [c for c in self._schema.columns if c.name != name]
        self._schema = CompiledSchema(
            row_cls=self._schema.row_cls,
            columns=new_columns,
            columns_by_name={c.name: c for c in new_columns},
        )
        if isinstance(self._storage, FileTableStorage):
            self._storage.save_schema(self._schema_dict_with_computed())

    def rename_column(self, old: str, new: str) -> None:
        """Rename a column.

        On disk tables the corresponding persisted column leaf is renamed.

        Renaming a flat column to a dotted name (e.g. ``"trip.begin.lon"``)
        promotes it to a nested leaf column: it will be stored under the
        hierarchical path ``/_cols/trip/begin/lon`` on disk and can be
        accessed via ``t["trip.begin.lon"]`` or the attribute-chain proxy
        ``t.trip.begin.lon``.  This is the primary way to define nested
        columns when importing from non-Arrow sources::

            t.rename_column("trip_begin_lon", "trip.begin.lon")
            t["trip.begin.lon"].mean()   # works as a regular Column

        Raises
        ------
        ValueError
            If the table is read-only, is a view, or *new* already exists.
        KeyError
            If *old* does not exist.
        """
        if self._read_only:
            raise ValueError("Table is read-only (opened with mode='r').")
        if self.base is not None:
            raise ValueError("Cannot rename a column in a view.")
        if old not in self._cols and old not in self._computed_cols:
            raise KeyError(f"No column named {old!r}. Available: {self.col_names}")
        if new in self._cols or new in self._computed_cols:
            raise ValueError(f"Column {new!r} already exists.")
        _validate_column_name(new)

        # Computed columns have no physical storage or schema entry.  Renaming
        # them only updates the computed-column registry and visible names.
        if old in self._computed_cols:
            self._computed_cols[new] = self._computed_cols.pop(old)
            idx = self.col_names.index(old)
            self.col_names[idx] = new
            self._col_widths[new] = max(len(new), self._col_widths.pop(old))
            if isinstance(self._storage, FileTableStorage):
                self._storage.save_schema(self._schema_dict_with_computed())
            return

        # Guard: refuse rename if any computed column depends on this stored column
        dependents = [cc_name for cc_name, cc in self._computed_cols.items() if old in cc["col_deps"]]
        if dependents:
            raise ValueError(
                f"Cannot rename column {old!r}: it is used by computed column(s) "
                + ", ".join(repr(d) for d in dependents)
                + ". Drop those computed columns first."
            )

        catalog = self._get_index_catalog()
        rebuild_kwargs = None
        if old in catalog:
            descriptor = catalog.pop(old)
            self._validate_index_descriptor(old, descriptor)
            rebuild_kwargs = self._index_create_kwargs_from_descriptor(descriptor)
            self._drop_index_descriptor(old, descriptor)
            self._storage.save_index_catalog(catalog)
            self._invalidate_index_catalog_cache()

        if isinstance(self._storage, FileTableStorage):
            self._cols[new] = self._rename_stored_column(old, new)
        else:
            self._cols[new] = self._cols[old]
        del self._cols[old]

        idx = self.col_names.index(old)
        self.col_names[idx] = new
        old_compiled = self._schema.columns_by_name[old]
        self._col_widths.pop(old)
        self._col_widths[new] = max(len(new), old_compiled.display_width)

        renamed = CompiledColumn(
            name=new,
            py_type=old_compiled.py_type,
            spec=old_compiled.spec,
            dtype=old_compiled.dtype,
            default=old_compiled.default,
            config=old_compiled.config,
            display_width=old_compiled.display_width,
        )
        new_columns = [renamed if c.name == old else c for c in self._schema.columns]
        self._schema = CompiledSchema(
            row_cls=self._schema.row_cls,
            columns=new_columns,
            columns_by_name={c.name: c for c in new_columns},
        )
        if old in self._materialized_cols:
            self._materialized_cols[new] = self._materialized_cols.pop(old)
        if isinstance(self._storage, FileTableStorage):
            self._storage.save_schema(self._schema_dict_with_computed())
        if rebuild_kwargs is not None:
            self.create_index(new, **rebuild_kwargs)

    def _rename_stored_column(self, old: str, new: str):
        """Rename a stored column's persistent leaves and return the reopened column."""
        old_compiled_col = self._schema.columns_by_name[old]
        if hasattr(self._cols[old], "flush"):
            self._cols[old].flush()
        renamed_col = self._storage.rename_column(old, new)
        if self._is_utf8_column(old_compiled_col):
            # rename_column returns the bare offsets NDArray; reopen the
            # offsets + bytes pair as a proper utf8 column object.
            renamed_col = self._storage.open_varlen_scalar_column(new, old_compiled_col.spec)
        return renamed_col

    # ------------------------------------------------------------------
    # Computed / virtual columns
    # ------------------------------------------------------------------

    @property
    def _stored_col_names(self) -> list[str]:
        """Column names backed by physical NDArrays (excludes computed columns)."""
        return [n for n in self.col_names if n not in self._computed_cols]

    @property
    def _append_input_col_names(self) -> list[str]:
        """Stored columns that callers must normally provide on insert."""
        return [n for n in self._stored_col_names if n not in self._materialized_cols]

    @property
    def computed_columns(self) -> dict[str, dict]:
        """Read-only view of the computed-column definitions.

        Each value is a dict with keys ``expression``, ``col_deps``,
        ``lazy`` (:class:`blosc2.LazyExpr`), and ``dtype``.
        """
        return dict(self._computed_cols)  # shallow copy so callers can't mutate

    def _aligned_grid(self) -> tuple | None:
        """Return ``(chunks, blocks)`` shared by the aligned fixed-size columns.

        Inspects the actual stored 1-D NDArray columns and returns the grid
        used by the largest aligned subset (the fast-path set), or ``None`` when
        there are no such columns.  Works for both freshly created and reopened
        tables since it reads the columns' real chunk/block shapes.
        """
        from collections import Counter

        grids = Counter()
        for name in self._stored_col_names:
            col = self._cols.get(name) if hasattr(self._cols, "get") else self._cols[name]
            chunks = getattr(col, "chunks", None)
            blocks = getattr(col, "blocks", None)
            if chunks is None or blocks is None or len(chunks) != 1:
                continue
            grids[(tuple(chunks), tuple(blocks))] += 1
        if not grids:
            return None
        (chunks, blocks), _ = grids.most_common(1)[0]
        return chunks, blocks

    @property
    def chunks(self) -> tuple | None:
        """Chunk shape shared by the table's aligned fixed-size columns.

        ``None`` if the table has no fixed-size scalar columns.  See
        :attr:`blocks` for the matching block shape.
        """
        grid = self._aligned_grid()
        return grid[0] if grid is not None else None

    @property
    def blocks(self) -> tuple | None:
        """Block shape shared by the table's aligned fixed-size columns.

        ``None`` if the table has no fixed-size scalar columns.  See
        :attr:`chunks` for the matching chunk shape.
        """
        grid = self._aligned_grid()
        return grid[1] if grid is not None else None

    def _ensure_generated_column_not_stale(self, name: str) -> None:
        meta = self._root_table._materialized_cols.get(name)
        if meta is not None and meta.get("stale", False):
            raise ValueError(
                f"Generated column {name!r} is stale because one or more source columns were modified. "
                f"Call refresh_generated_column({name!r}) before using it, or use t[{name!r}].read_stale() "
                "to explicitly read the last stored stale values."
            )

    def _col_dtype(self, name: str) -> np.dtype | None:
        """Return the dtype for *name*, routing through computed cols."""
        cc = self._computed_cols.get(name)
        if cc is not None:
            return cc["dtype"]
        return getattr(self._cols[name], "dtype", None)

    @staticmethod
    def _readable_computed_expr(cc: dict) -> str:
        """Return a human-readable description of a computed column.

        For expression columns the stored string has ``o0``, ``o1``, … replaced
        by their actual column names (``"(o0 * o1)"`` with
        ``col_deps=["price", "qty"]`` → ``"(price * qty)"``).  For DSL columns a
        ``kernel(dep0, dep1)`` call label is returned.
        """
        col_deps = cc["col_deps"]
        if cc.get("kind") == "dsl":
            kernel = cc.get("kernel")
            kname = getattr(kernel, "__name__", "dsl_kernel")
            return f"{kname}({', '.join(col_deps)})"

        def _sub(m: re.Match) -> str:
            idx = int(m.group(1))
            return col_deps[idx] if idx < len(col_deps) else m.group(0)

        return re.sub(r"\bo(\d+)\b", _sub, cc["expression"])

    def _fetch_col_at_positions(self, name: str, positions: np.ndarray):
        """Fetch values at *positions* (physical indices) — used for display.

        During a single ``to_string()`` call the same ``(column, positions)``
        pair is requested repeatedly — by float-precision detection, column
        width sizing and the final row rendering — all sharing the same
        ``head_pos``/``tail_pos`` arrays.  When ``to_string`` installs a
        ``_display_fetch_cache`` we memoise the sparse gather (keyed by the
        positions array identity) so each column is read from storage once
        instead of ~6 times.
        """
        cache = getattr(self, "_display_fetch_cache", None)
        if cache is None:
            return self._fetch_col_at_positions_uncached(name, positions)
        key = (name, id(positions))
        if key not in cache:
            cache[key] = self._fetch_col_at_positions_uncached(name, positions)
        return cache[key]

    def _fetch_col_at_positions_uncached(self, name: str, positions: np.ndarray):
        cc = self._computed_cols.get(name)
        if cc is not None:
            if len(positions) == 0:
                return np.array([], dtype=cc["dtype"])
            lazy = self._build_computed_lazy(cc)
            return np.array(
                [np.asarray(lazy[int(p)]).ravel()[0] for p in positions],
                dtype=cc["dtype"],
            )
        self._ensure_generated_column_not_stale(name)
        col = self._cols[name]
        spec = self._schema.columns_by_name[name].spec
        if self._is_list_spec(spec) or isinstance(
            spec, (VLStringSpec, VLBytesSpec, StructSpec, ObjectSpec, DictionarySpec)
        ):
            return col[positions]
        values = col[positions]
        if isinstance(spec, timestamp):
            return np.asarray(values).astype(f"datetime64[{spec.unit}]")
        return values

    def _schema_dict_with_computed(self) -> dict:
        """Return the schema dict extended with computed/materialized metadata."""
        d = schema_to_dict(self._schema)
        n_rows = self._known_n_rows()
        if n_rows is not None:
            d["n_rows"] = n_rows
        d["create_summary_index"] = getattr(self, "_create_summary_index", True)
        d["summary_indexes_built"] = getattr(self, "_summary_indexes_built", False)
        if self._computed_cols:
            computed = []
            for name, cc in self._computed_cols.items():
                if cc.get("kind") == "dsl":
                    entry = {
                        "name": name,
                        "kind": "dsl",
                        "dsl_source": cc["dsl_source"],
                        "col_deps": cc["col_deps"],
                        "dtype": str(cc["dtype"]),
                    }
                    if cc.get("jit_backend") is not None:
                        entry["jit_backend"] = cc["jit_backend"]
                    computed.append(entry)
                else:
                    computed.append(
                        {
                            "name": name,
                            "kind": "expression",
                            "expression": cc["expression"],
                            "col_deps": cc["col_deps"],
                            "dtype": str(cc["dtype"]),
                        }
                    )
            d["computed_columns"] = computed
        if self._materialized_cols:
            materialized = []
            for name, meta in self._materialized_cols.items():
                entry = {
                    "name": name,
                    "computed_column": meta.get("computed_column"),
                    "expression": meta.get("expression"),
                    "col_deps": meta["col_deps"],
                    "dtype": str(meta["dtype"]),
                    "transformer_kind": meta.get("transformer_kind", "expression"),
                    "stale": bool(meta.get("stale", False)),
                }
                if "transformer" in meta:
                    entry["transformer"] = meta["transformer"]
                if meta.get("dsl_source") is not None:
                    entry["dsl_source"] = meta["dsl_source"]
                if meta.get("jit_backend") is not None:
                    entry["jit_backend"] = meta["jit_backend"]
                materialized.append(entry)
            d["materialized_columns"] = materialized
        return d

    def _save_n_rows_to_meta(self) -> None:
        """Persist the cached row count into the _meta SChunk's vlmeta.

        Updates the vlmeta of the existing _meta SChunk directly and writes
        it back to its backing store.  This avoids going through save_schema()
        which can route through the embed store where SChunk slice writes may
        fail when the backing store has chunksize=-1.
        """
        n_rows = self._known_n_rows()
        if n_rows is None:
            return
        storage = self._storage
        if not hasattr(storage, "_open_meta"):
            return
        try:
            meta = storage._open_meta()
            schema_raw = meta.vlmeta.get("schema")
            if schema_raw is None:
                return
            schema_dict = json.loads(schema_raw)
            schema_dict["n_rows"] = n_rows
            schema_dict["summary_indexes_built"] = getattr(self, "_summary_indexes_built", False)
            meta.vlmeta["schema"] = json.dumps(schema_dict)
            # Persist: for FileTableStorage, rewrite the external _meta.b2f file.
            if hasattr(storage, "_meta_path"):
                meta.save(urlpath=storage._meta_path, mode="w")
            elif hasattr(storage, "_write_leaf"):
                # TreeStoreTableStorage
                storage._write_leaf("/_meta", meta, ".b2f")
        except Exception:
            pass  # best-effort; failure must not prevent close()

    def _is_summary_eligible_column(self, col) -> bool:
        """True if *col* can carry an automatic SUMMARY index.

        Eligible columns are stored, scalar, numeric-or-boolean, and not
        list/varlen/dictionary/computed.  Shared by the close-time builder and
        the incremental accumulator so both agree on the column set.
        """
        name = col.name
        if name in self._computed_cols:
            return False
        if self._is_list_column(col) or self._is_varlen_scalar_column(col):
            return False
        if self._is_dictionary_column(col) or self._is_ndarray_column(col):
            return False
        spec = col.spec
        return hasattr(spec, "dtype") and (
            np.issubdtype(np.dtype(spec.dtype), np.number) or np.issubdtype(np.dtype(spec.dtype), np.bool_)
        )

    def _get_summary_accumulator(self, name: str):
        """Return the live per-block summary accumulator for *name*, or None.

        Lazily creates one the first time an eligible column on a writable
        top-level table is fed.  ``None`` (cached) means the column is not a
        candidate for incremental summaries; an *invalid* accumulator means a
        non-append mutation happened and the close-time builder must fall back
        to the out-of-core path.
        """
        if self.base is not None or getattr(self, "_read_only", False):
            return None
        if not getattr(self, "_create_summary_index", True):
            return None
        accs = self.__dict__.get("_summary_accumulators")
        if accs is None:
            accs = {}
            self._summary_accumulators = accs
        if name in accs:
            return accs[name]
        col = self._schema.columns_by_name.get(name)
        arr = self._cols.get(name)
        if col is None or arr is None or not self._is_summary_eligible_column(col):
            accs[name] = None
            return None
        try:
            block_len = int(arr.blocks[0])
            dtype = arr.dtype
        except Exception:
            accs[name] = None
            return None
        acc = _ColumnSummaryAccumulator(dtype, block_len)
        accs[name] = acc
        return acc

    def _summary_feeder(self, name: str):
        """Return a ``callback(start_pos, values)`` bound to *name*, or None."""
        acc = self._get_summary_accumulator(name)
        if acc is None:
            return None
        return acc.feed

    def _feed_summary(self, name: str, start_pos: int, values: np.ndarray) -> None:
        acc = self._get_summary_accumulator(name)
        if acc is not None:
            acc.feed(start_pos, values)

    def _invalidate_summary_accumulator(self, name: str) -> None:
        accs = self.__dict__.get("_summary_accumulators")
        if accs:
            acc = accs.get(name)
            if acc is not None:
                acc.invalidate()

    def _invalidate_all_summary_accumulators(self) -> None:
        accs = self.__dict__.get("_summary_accumulators")
        if accs:
            for acc in accs.values():
                if acc is not None:
                    acc.invalidate()

    def _precomputed_summary_for(self, name: str):
        """Return ``{"block": summaries}`` for *name* if a valid accumulator
        fully covers the column's physical extent, else None."""
        accs = self.__dict__.get("_summary_accumulators")
        if not accs:
            return None
        acc = accs.get(name)
        if acc is None:
            return None
        arr = self._cols.get(name)
        if arr is None:
            return None
        try:
            expected = int(arr.shape[0])
        except Exception:
            return None
        summaries = acc.finalize(expected)
        if summaries is None:
            return None
        return {"block": summaries}

    def _build_summary_indexes(self) -> None:
        """Create SUMMARY indexes for all eligible scalar columns.

        Called once from :meth:`close` when ``create_summary_index=True`` (the default).
        Skips list, varlen, dictionary, and computed columns.  If all eligible
        columns are already indexed (e.g. from a prior close), this is a no-op.
        Failure on any individual column is silently ignored — it must not
        prevent close().

        When an incremental per-block accumulator fully covers a column (the
        common build-from-empty case), its precomputed summaries are handed to
        ``create_index`` so the column is *not* decompressed again just to
        recompute min/max.  Otherwise the index builder falls back to the
        out-of-core decompress-and-scan path transparently.
        """
        catalog = self._get_index_catalog()
        indexed = set(catalog) if catalog else set()
        eligible = [
            col.name
            for col in self._schema.columns
            if col.name not in indexed and self._is_summary_eligible_column(col)
        ]
        if not eligible:
            self._summary_indexes_built = True
            return

        import warnings

        for name in eligible:
            try:
                precomputed = self._precomputed_summary_for(name)
                if precomputed is not None:
                    self.create_index(name, kind=blosc2.IndexKind.SUMMARY, precomputed_summaries=precomputed)
                else:
                    self.create_index(name, kind=blosc2.IndexKind.SUMMARY)
            except Exception:
                warnings.warn(
                    f"Failed to create SUMMARY index for column {name!r}; skipping.",
                    stacklevel=2,
                )
        self._summary_indexes_built = True

    def _load_computed_cols_from_schema(self, schema_dict: dict) -> None:
        """Reconstruct ``_computed_cols`` from persisted metadata.

        Called from ``__init__``, ``open``, and ``load`` after all stored
        columns have been opened into ``self._cols``.
        """
        for cc_meta in schema_dict.get("computed_columns", []):
            name = cc_meta["name"]
            col_deps = cc_meta["col_deps"]
            dtype = np.dtype(cc_meta["dtype"])
            if cc_meta.get("kind") == "dsl":
                from blosc2.dsl_kernel import kernel_from_source

                dsl_source = cc_meta["dsl_source"]
                self._computed_cols[name] = {
                    "kind": "dsl",
                    "dsl_source": dsl_source,
                    "kernel": kernel_from_source(dsl_source),
                    "col_deps": col_deps,
                    "dtype": dtype,
                    "jit_backend": cc_meta.get("jit_backend"),
                }
            else:
                expression = cc_meta["expression"]
                operands = {f"o{i}": self._cols[dep] for i, dep in enumerate(col_deps)}
                lazy = blosc2.lazyexpr(expression, operands)
                self._computed_cols[name] = {
                    "kind": "expression",
                    "expression": expression,
                    "col_deps": col_deps,
                    "lazy": lazy,
                    "dtype": dtype,
                }
            self.col_names.append(name)
            self._col_widths[name] = max(len(name), 15)

    def _load_materialized_cols_from_schema(self, schema_dict: dict) -> None:
        """Reconstruct ``_materialized_cols`` from persisted metadata."""
        for meta in schema_dict.get("materialized_columns", []):
            loaded = {
                "computed_column": meta.get("computed_column"),
                "expression": meta.get("expression"),
                "col_deps": list(meta["col_deps"]),
                "dtype": np.dtype(meta["dtype"]),
                "transformer_kind": meta.get("transformer_kind", "expression"),
                "stale": bool(meta.get("stale", False)),
            }
            if "transformer" in meta:
                loaded["transformer"] = dict(meta["transformer"])
            if meta.get("dsl_source") is not None:
                loaded["dsl_source"] = meta["dsl_source"]
            if meta.get("jit_backend") is not None:
                loaded["jit_backend"] = meta["jit_backend"]
            self._materialized_cols[meta["name"]] = loaded

    def _require_computed_column(self, name: str) -> dict:
        """Return metadata for computed column *name* or raise ``KeyError``."""
        try:
            return self._computed_cols[name]
        except KeyError:
            raise KeyError(
                f"{name!r} is not a computed column. Computed columns: {list(self._computed_cols)}"
            ) from None

    def _autofill_materialized_row_values(self, row: dict[str, Any]) -> dict[str, Any]:
        """Fill omitted materialized-column values for a single inserted row."""
        row = dict(row)
        for name, meta in self._materialized_cols.items():
            if name in row:
                continue
            missing = [dep for dep in meta["col_deps"] if dep not in row]
            if missing:
                raise ValueError(
                    f"Cannot auto-fill materialized column {name!r}: missing dependency columns {missing!r}."
                )
            if meta.get("transformer_kind") == "row_transformer":
                transformer = RowTransformer.from_metadata(meta["transformer"])
                row[name] = np.asarray(transformer.evaluate_row(row), dtype=meta["dtype"])
            elif meta.get("transformer_kind") == "dsl":
                single = {dep: [row[dep]] for dep in meta["col_deps"]}
                row[name] = self._evaluate_dsl_materialized_batch(meta, single)[0]
            else:
                operands = {f"o{i}": np.asarray([row[dep]]) for i, dep in enumerate(meta["col_deps"])}
                values = blosc2.lazyexpr(meta["expression"], operands)[:]
                row[name] = np.asarray(values, dtype=meta["dtype"])[0]
        return row

    def _validate_no_default_columns_present(self, row: dict[str, Any]) -> None:
        """Raise a clear error when a row omits a column with no default declared."""
        for col in self._schema.columns:
            if col.name in row:
                continue
            is_nullable = getattr(col.spec, "null_value", None) is not None or bool(
                getattr(col.spec, "nullable", False)
            )
            if col.default is MISSING and not is_nullable:
                raise ValueError(f"Column {col.name!r} has no default declared; a value must be provided.")

    def _fill_default_batch_columns(self, raw_columns: dict[str, Any], row_count: int) -> dict[str, Any]:
        """Fill omitted batch columns from defaults, or raise if no default is declared."""
        raw_columns = dict(raw_columns)
        for col in self._schema.columns:
            if col.name in raw_columns:
                continue
            if col.default is MISSING:
                raise ValueError(f"Column {col.name!r} has no default declared; values must be provided.")
            raw_columns[col.name] = [col.default] * row_count
        return raw_columns

    def _autofill_materialized_batch_columns(
        self, raw_columns: dict[str, Any], row_count: int, *, provided_names: set[str]
    ) -> dict[str, Any]:
        """Fill omitted materialized-column arrays for batch inserts."""
        raw_columns = dict(raw_columns)
        for name, meta in self._materialized_cols.items():
            if name in provided_names or name in raw_columns:
                continue
            missing = [dep for dep in meta["col_deps"] if dep not in raw_columns]
            if missing:
                raise ValueError(
                    f"Cannot auto-fill materialized column {name!r}: missing dependency columns {missing!r}."
                )
            if meta.get("transformer_kind") == "row_transformer":
                transformer = RowTransformer.from_metadata(meta["transformer"])
                values = transformer.evaluate_batch(raw_columns)
            elif meta.get("transformer_kind") == "dsl":
                values = self._evaluate_dsl_materialized_batch(meta, raw_columns)
            else:
                operands = {
                    f"o{i}": blosc2.asarray(raw_columns[dep], dtype=self._cols[dep].dtype)
                    for i, dep in enumerate(meta["col_deps"])
                }
                values = blosc2.lazyexpr(meta["expression"], operands)[:]
            values = np.asarray(values, dtype=meta["dtype"])
            if len(values) != row_count:
                raise ValueError(
                    f"Materialized column {name!r} produced {len(values)} values, expected {row_count}."
                )
            raw_columns[name] = values
        return raw_columns

    @staticmethod
    def _coerce_generated_spec(dtype_or_spec, sample: np.ndarray | None = None) -> SchemaSpec:
        """Resolve a generated-column dtype/spec, inferring ndarray shape when needed."""
        if isinstance(dtype_or_spec, SchemaSpec):
            return dtype_or_spec
        if dtype_or_spec is None:
            if sample is None:
                raise TypeError("dtype is required when a generated column has no rows to infer from.")
            arr = np.asarray(sample)
            if arr.ndim <= 1:
                return CTable._schema_spec_from_dtype(arr.dtype)
            return NDArraySpec(item_shape=arr.shape[1:], dtype=arr.dtype)
        return CTable._schema_spec_from_dtype(np.dtype(dtype_or_spec))

    @staticmethod
    def _schema_spec_from_dtype(dtype: np.dtype) -> SchemaSpec:
        """Build a minimal schema spec for a stored column with *dtype*."""
        dtype = np.dtype(dtype)
        spec_factory = _DTYPE_SPEC_FACTORIES.get(dtype)
        if spec_factory is not None:
            return spec_factory()
        if dtype.kind == "U":
            max_length = max(1, dtype.itemsize // np.dtype("U1").itemsize)
            return string(max_length=max_length)
        if dtype.kind == "S":
            return b2_bytes(max_length=max(1, dtype.itemsize))
        raise TypeError(f"Cannot materialize a computed column with unsupported dtype {dtype!r}.")

    def _create_empty_stored_column(
        self,
        name: str,
        dtype: np.dtype | None,
        *,
        spec: SchemaSpec | None = None,
        cparams: dict | None = None,
    ) -> None:
        """Create an empty stored column aligned with the table's physical row space."""
        if spec is None:
            if dtype is None:
                raise TypeError("dtype or spec is required")
            spec = self._schema_spec_from_dtype(dtype)
        dtype = np.dtype(spec.dtype)
        if isinstance(spec, NDArraySpec):
            default = np.zeros(spec.item_shape, dtype=dtype)
        else:
            default = np.array(0, dtype=dtype).item() if dtype.kind not in {"U", "S"} else dtype.type()

        capacity = len(self._valid_rows)
        compiled_col = CompiledColumn(
            name=name,
            py_type=spec.python_type,
            spec=spec,
            dtype=dtype,
            default=default,
            config=ColumnConfig(cparams=cparams, dparams=None, chunks=None, blocks=None),
            display_width=compute_display_width(spec),
        )
        shape = self._column_physical_shape(compiled_col, capacity)
        default_chunks, default_blocks = self._column_chunks_blocks(compiled_col, shape)
        new_col = self._storage.create_column(
            name,
            dtype=dtype,
            shape=shape,
            chunks=default_chunks,
            blocks=default_blocks,
            cparams=cparams,
            dparams=None,
        )
        self._cols[name] = new_col
        self.col_names.append(name)
        self._col_widths[name] = max(len(name), compiled_col.display_width)

        new_columns = self._schema.columns + [compiled_col]
        self._schema = CompiledSchema(
            row_cls=self._schema.row_cls,
            columns=new_columns,
            columns_by_name={**self._schema.columns_by_name, name: compiled_col},
        )
        if isinstance(self._storage, FileTableStorage):
            self._storage.save_schema(self._schema_dict_with_computed())

    def _fill_stored_column_from_computed(
        self,
        target_name: str,
        computed_name: str,
        *,
        dtype: np.dtype,
    ) -> None:
        """Evaluate computed column *computed_name* into stored column *target_name*."""
        cc = self._require_computed_column(computed_name)
        # Expression entries yield a LazyExpr (streamed per slice); DSL entries
        # yield a fully materialized NDArray (the miniexpr DSL path cannot do
        # partial-slice getitem).  Both support the slicing used below.
        lazy = self._build_computed_lazy(cc)
        capacity = len(self._valid_rows)
        step = int(self._valid_rows.chunks[0]) if self._valid_rows.chunks else 65536

        for start in range(0, capacity, step):
            stop = min(start + step, capacity)
            values = lazy[start:stop]
            if isinstance(values, blosc2.NDArray):
                values = values[:]
            try:
                values = np.asarray(values, dtype=dtype)
            except (TypeError, ValueError) as exc:
                raise TypeError(f"Cannot coerce computed values to dtype {dtype!r}: {exc}") from exc
            if values.ndim != 1:
                raise TypeError(
                    f"Computed column {computed_name!r} produced {values.ndim}-D values; expected 1-D slices."
                )
            if len(values) != stop - start:
                raise ValueError(
                    f"Computed column {computed_name!r} produced {len(values)} values for slice "
                    f"[{start}:{stop}], expected {stop - start}."
                )
            self._cols[target_name][start:stop] = values

    def materialize_computed_column(
        self,
        name: str,
        *,
        new_name: str | None = None,
        dtype: np.dtype | None = None,
        cparams: dict | blosc2.CParams | None = None,
    ) -> None:
        """Materialize a computed column into a new stored snapshot column.

        Parameters
        ----------
        name:
            Existing computed column to materialize.
        new_name:
            Name of the new stored column. Defaults to ``f"{name}_stored"``.
        dtype:
            Optional target dtype for the stored column. Defaults to the
            computed column dtype.
        cparams:
            Optional compression parameters for the new stored column.

        Raises
        ------
        ValueError
            If called on a view, on a read-only table, or if the target name
            collides with an existing stored or computed column.
        KeyError
            If *name* is not a computed column.
        TypeError
            If *dtype* is incompatible with the computed values.
        """
        if self.base is not None:
            raise ValueError("Cannot materialize a computed column from a view.")
        if self._read_only:
            raise ValueError("Table is read-only (opened with mode='r').")

        cc = self._require_computed_column(name)
        target_name = new_name or f"{name}_stored"
        _validate_column_name(target_name)
        if target_name in self._cols:
            raise ValueError(f"A stored column named {target_name!r} already exists.")
        if target_name in self._computed_cols:
            raise ValueError(f"A computed column named {target_name!r} already exists.")
        target_dtype = np.dtype(dtype) if dtype is not None else np.dtype(cc["dtype"])

        self._create_empty_stored_column(target_name, target_dtype, cparams=cparams)
        if cc.get("kind") == "dsl":
            self._materialized_cols[target_name] = {
                "computed_column": name,
                "expression": None,
                "dsl_source": cc["dsl_source"],
                "col_deps": list(cc["col_deps"]),
                "dtype": target_dtype,
                "transformer_kind": "dsl",
                "stale": False,
                "jit_backend": cc.get("jit_backend"),
            }
        else:
            self._materialized_cols[target_name] = {
                "computed_column": name,
                "expression": cc["expression"],
                "col_deps": list(cc["col_deps"]),
                "dtype": target_dtype,
            }
        if isinstance(self._storage, FileTableStorage):
            self._storage.save_schema(self._schema_dict_with_computed())
        try:
            self._fill_stored_column_from_computed(target_name, name, dtype=target_dtype)
        except Exception:
            with contextlib.suppress(Exception):
                self.drop_column(target_name)
            raise

    def _normalize_expression_transformer(self, expr) -> tuple[blosc2.LazyExpr, list[str]]:
        if isinstance(expr, RowTransformer):
            raise TypeError(
                "RowTransformer instances cannot be used for computed columns; use add_generated_column()."
            )
        if isinstance(expr, blosc2.LazyExpr):
            lazy = expr
        elif callable(expr):
            lazy = expr(self._cols)
        elif isinstance(expr, str):
            self._guard_scalar_expression(expr)
            operands = self._where_expression_operands(expr)
            expr, operands = self._rewrite_nested_expression(expr, operands)
            lazy = blosc2.lazyexpr(expr, operands)
        else:
            raise TypeError(
                f"expr must be a callable or an expression string (or LazyExpr), got {type(expr).__name__!r}."
            )
        if not isinstance(lazy, blosc2.LazyExpr):
            raise TypeError(f"expr must return a blosc2.LazyExpr, got {type(lazy).__name__!r}.")

        owned_ids = {id(arr): cname for cname, arr in self._cols.items()}
        col_deps = []
        for key in sorted(lazy.operands.keys()):
            arr = lazy.operands[key]
            cname = owned_ids.get(id(arr))
            if cname is None:
                raise ValueError(
                    f"Operand {key!r} in the expression does not reference a stored column of this table."
                )
            self._ensure_generated_column_not_stale(cname)
            col_info = self._schema.columns_by_name.get(cname)
            if col_info is not None and self._is_ndarray_column(col_info):
                raise TypeError(
                    f"Column {cname!r} is a fixed-shape ndarray column. Expression transformers only "
                    "support scalar columns; use a RowTransformer for ndarray row reductions/projections."
                )
            col_deps.append(cname)
        return lazy, col_deps

    def _validate_transformer_dep(self, cname: str) -> blosc2.NDArray:
        """Validate that *cname* is a stored scalar column usable as a transformer
        operand and return its backing NDArray."""
        if cname not in self._cols:
            raise ValueError(f"Column {cname!r} is not a stored column of this table.")
        self._ensure_generated_column_not_stale(cname)
        col_info = self._schema.columns_by_name.get(cname)
        if col_info is not None and self._is_ndarray_column(col_info):
            raise TypeError(
                f"Column {cname!r} is a fixed-shape ndarray column. DSL kernels only "
                "support scalar columns as inputs."
            )
        return self._cols[cname]

    def _dsl_deps_from_lazyudf(self, lazyudf) -> list[str]:
        """Return the stored-column names backing a DSL LazyUDF's inputs, in order."""
        owned_ids = {id(arr): cname for cname, arr in self._cols.items()}
        col_deps = []
        for i, arr in enumerate(lazyudf.inputs):
            cname = owned_ids.get(id(arr))
            if cname is None:
                raise ValueError(
                    f"Input {i} of the DSL kernel does not reference a stored column of this table."
                )
            self._validate_transformer_dep(cname)
            col_deps.append(cname)
        return col_deps

    def _resolve_dsl_kernel(self, kernel, inputs) -> tuple[Any, list[str]]:
        """Validate a bare DSL kernel + its ``inputs`` column bindings."""
        if kernel.dsl_error is not None:
            raise blosc2.DSLSyntaxError(f"Invalid DSL kernel: {kernel.dsl_error}")
        if inputs is None:
            raise TypeError(
                "A DSL kernel passed directly requires inputs=[...] naming one source "
                "column per kernel parameter."
            )
        col_deps = list(inputs)
        expected = kernel.input_names
        if expected is not None and len(col_deps) != len(expected):
            raise ValueError(
                f"DSL kernel expects {len(expected)} input(s) {expected}, "
                f"but inputs={col_deps} provides {len(col_deps)}."
            )
        for d in col_deps:
            self._validate_transformer_dep(d)
        return kernel, col_deps

    def _normalize_transformer(self, expr, inputs=None) -> dict:
        """Resolve *expr* into a transformer descriptor.

        Returns one of::

            {"kind": "expression", "lazy": <LazyExpr>, "col_deps": [...]}
            {"kind": "dsl",        "kernel": <DSLKernel>, "col_deps": [...]}

        A ``@blosc2.dsl_kernel`` is accepted either directly (with *inputs*
        naming one source column per kernel parameter) or as a ``LazyUDF``
        returned by a callable (then operands are matched to columns by
        identity and *inputs* is ignored).
        """
        if isinstance(expr, blosc2.DSLKernel):
            kernel, col_deps = self._resolve_dsl_kernel(expr, inputs)
            return {"kind": "dsl", "kernel": kernel, "col_deps": col_deps}
        # Resolve a callable once (a lambda may return a LazyExpr or a LazyUDF).
        obj = expr(self._cols) if (callable(expr) and not isinstance(expr, blosc2.LazyExpr)) else expr
        if isinstance(obj, blosc2.LazyUDF):
            if not isinstance(obj.func, blosc2.DSLKernel):
                raise TypeError(
                    "Only LazyUDFs backed by a @blosc2.dsl_kernel are supported as CTable columns."
                )
            kernel = obj.func
            if kernel.dsl_error is not None:
                raise blosc2.DSLSyntaxError(f"Invalid DSL kernel: {kernel.dsl_error}")
            return {
                "kind": "dsl",
                "kernel": kernel,
                "col_deps": self._dsl_deps_from_lazyudf(obj),
                "jit_backend": obj.kwargs.get("jit_backend"),
            }
        lazy, col_deps = self._normalize_expression_transformer(obj)
        # Guard: verify the expression string round-trips before storing.
        # An empty string means the LazyExpr was not fully constructed, and a
        # malformed string would silently break on reload.  Catching both here
        # gives an early, actionable error instead of a confusing failure later.
        expression = lazy.expression
        if not expression:
            raise ValueError(
                "The computed-column expression serializes to an empty string "
                "and cannot be persisted. Make sure the lambda returns a "
                "blosc2 expression built from table columns (e.g. cols['x'] * 2)."
            )
        try:
            _ops = {f"o{i}": self._cols[dep] for i, dep in enumerate(col_deps)}
            blosc2.lazyexpr(expression, _ops)
        except Exception as exc:
            raise ValueError(
                f"The computed-column expression {expression!r} cannot be safely "
                f"persisted and reloaded: {exc}"
            ) from exc
        return {"kind": "expression", "lazy": lazy, "col_deps": col_deps}

    def _dsl_result_dtype(self, kernel, col_deps, dtype):
        """Resolve the result dtype for a DSL column.

        When *dtype* is omitted it is inferred by NumPy type promotion of the
        dependency column dtypes — correct for elementwise arithmetic kernels.
        Kernels that change the type (comparisons/``where``/explicit casts)
        should pass *dtype* explicitly.
        """
        if dtype is not None:
            return np.dtype(dtype)
        dep_dtypes = [self._cols[d].dtype for d in col_deps]
        if not dep_dtypes:
            raise TypeError(
                f"Cannot infer dtype for DSL kernel {getattr(kernel, '__name__', '?')!r} "
                "with no column inputs; pass dtype=... explicitly."
            )
        return np.result_type(*dep_dtypes)

    def _build_computed_lazy(self, cc: dict):
        """Return the readable array-like for a computed-column entry *cc*.

        Expression entries return their cached :class:`blosc2.LazyExpr` (which
        supports partial-slice evaluation directly).  DSL entries build a fresh
        :class:`blosc2.LazyUDF` from the current column NDArrays and **eagerly
        materialize** it to a concrete :class:`blosc2.NDArray` via
        ``compute()``: the miniexpr DSL path only supports full-array getitem,
        so a partial slice (used by reads and by ``where()`` per-chunk operand
        access) cannot be evaluated lazily.  Materializing also lets a DSL
        computed column participate in ``where()`` as a plain NDArray operand
        (the all-NDArray miniexpr fast path).  The full column is recomputed on
        each access — acceptable for a virtual, unstored column.
        """
        if cc.get("kind") == "dsl":
            operands = tuple(self._cols[d] for d in cc["col_deps"])
            return blosc2.lazyudf(
                cc["kernel"],
                operands,
                dtype=cc["dtype"],
                jit_backend=cc.get("jit_backend"),
            ).compute()
        return cc["lazy"]

    def _evaluate_expression_materialized_batch(
        self, meta: dict, raw_columns: Mapping[str, Any]
    ) -> np.ndarray:
        operands = {
            f"o{i}": blosc2.asarray(raw_columns[dep], dtype=self._cols[dep].dtype)
            for i, dep in enumerate(meta["col_deps"])
        }
        values = blosc2.lazyexpr(meta["expression"], operands)[:]
        return np.asarray(values, dtype=meta["dtype"])

    def _materialized_dsl_kernel(self, meta: dict):
        """Return the (cached) DSLKernel for a ``transformer_kind == "dsl"`` entry."""
        kernel = meta.get("_kernel")
        if kernel is None:
            from blosc2.dsl_kernel import kernel_from_source

            kernel = kernel_from_source(meta["dsl_source"])
            meta["_kernel"] = kernel  # not serialized (schema dump emits known keys only)
        return kernel

    def _evaluate_dsl_materialized_batch(self, meta: dict, raw_columns: Mapping[str, Any]) -> np.ndarray:
        kernel = self._materialized_dsl_kernel(meta)
        arrays = [np.asarray(raw_columns[dep], dtype=self._cols[dep].dtype) for dep in meta["col_deps"]]
        out_dtype = np.dtype(meta["dtype"])
        n = len(arrays[0]) if arrays else 0
        if n == 0:
            return np.asarray([], dtype=out_dtype)
        # The DSL miniexpr path rejects length-1 (shape ``(1,)``) inputs as
        # "scalar-like"; pad to length 2 and slice the result back.
        pad = n == 1
        if pad:
            arrays = [np.concatenate([arr, arr]) for arr in arrays]
        operands = tuple(blosc2.asarray(arr) for arr in arrays)
        result = blosc2.lazyudf(
            kernel,
            operands,
            dtype=out_dtype,
            jit_backend=meta.get("jit_backend"),
        ).compute()[:]
        result = np.asarray(result, dtype=out_dtype)
        return result[:1] if pad else result

    def _generated_dependency_closure(self, source: str) -> set[str]:
        """Return generated columns transitively depending on *source*."""
        affected: set[str] = set()
        queue = deque([source])
        while queue:
            current = queue.popleft()
            for name, meta in self._materialized_cols.items():
                if name in affected:
                    continue
                if current in meta.get("col_deps", ()):
                    affected.add(name)
                    queue.append(name)
        return affected

    def _mark_generated_columns_stale(self, source: str) -> None:
        affected = self._generated_dependency_closure(source)
        changed = False
        for name in affected:
            meta = self._materialized_cols[name]
            if not meta.get("stale", False):
                meta["stale"] = True
                changed = True
        if changed and isinstance(self._storage, FileTableStorage):
            self._storage.save_schema(self._schema_dict_with_computed())

    def refresh_generated_column(self, name: str) -> None:
        """Recompute a stored generated/materialized column from its source columns."""
        if self._read_only:
            raise ValueError("Table is read-only (opened with mode='r').")
        if name not in self._materialized_cols:
            raise KeyError(f"{name!r} is not a generated/materialized column.")
        meta = self._materialized_cols[name]
        n_live = self.nrows
        if meta.get("transformer_kind") == "row_transformer":
            transformer = RowTransformer.from_metadata(meta["transformer"])
            values = np.asarray(transformer.evaluate_existing(self), dtype=meta["dtype"])
        elif meta.get("transformer_kind") == "dsl":
            raw_columns = {dep: self[dep][:] for dep in meta["col_deps"]}
            values = self._evaluate_dsl_materialized_batch(meta, raw_columns)
        else:
            raw_columns = {dep: self[dep][:] for dep in meta["col_deps"]}
            values = self._evaluate_expression_materialized_batch(meta, raw_columns)
        if len(values) != n_live:
            raise ValueError(f"Generated column {name!r} produced {len(values)} values, expected {n_live}.")
        self._cols[name][self._valid_rows] = values
        meta["stale"] = False
        self._mark_all_indexes_stale()
        if isinstance(self._storage, FileTableStorage):
            self._storage.save_schema(self._schema_dict_with_computed())

    def refresh_generated_columns(self, *, source: str | None = None) -> None:
        """Refresh all generated columns, optionally only those depending on *source*."""
        affected = None if source is None else self._generated_dependency_closure(source)
        for name in list(self._materialized_cols):
            if affected is None or name in affected:
                self.refresh_generated_column(name)

    def apply(
        self,
        func,
        *,
        columns: list[str] | None = None,
        dtype=None,
        engine: str = "auto",
    ) -> np.ndarray:
        """Run a row-batch UDF over live column values and materialize the result.

        Sugar over :func:`blosc2.lazyudf` using this table's columns as
        inputs: ``blosc2.lazyudf(func, tuple(t[c] for c in columns),
        dtype=dtype, jit=...).compute()``. There is no separate execution
        path — this reuses exactly the machinery :meth:`add_computed_column`
        and :meth:`add_generated_column` already use for DSL/UDF columns.

        Parameters
        ----------
        func:
            UDF passed straight to :func:`blosc2.lazyudf`; see that function
            for the expected signature (``func(inputs_tuple, output,
            offset)``) and for how a :func:`blosc2.dsl_kernel`-decorated
            kernel is transpiled instead of run as plain Python.
        columns:
            Names of *stored* columns to bind as *func*'s inputs, in order
            (computed/generated columns are not supported). Defaults to
            every stored column, in schema order.
        dtype:
            Result dtype. Required unless *func* is a
            :func:`blosc2.dsl_kernel` kernel, whose dtype can be inferred by
            NumPy type promotion of the input dtypes -- same rule as
            :func:`blosc2.lazyudf`.
        engine:
            Forwarded to :func:`blosc2.lazyudf` as its ``jit`` policy:
            ``"auto"`` (default) lets it choose, ``"jit"`` forces JIT
            (only effective for a transpilable :func:`blosc2.dsl_kernel`),
            ``"numpy"`` disables JIT.

        Returns
        -------
        numpy.ndarray
            The UDF's result for the live rows only (on a view, the rows
            visible through the view).

        Examples
        --------
        >>> import blosc2
        >>> from dataclasses import dataclass
        >>> @dataclass
        ... class Row:
        ...     price: float = blosc2.field(blosc2.float64())
        ...     qty: float = blosc2.field(blosc2.float64())
        >>> t = blosc2.CTable(Row, new_data=[(10.0, 2.0), (5.0, 3.0)])
        >>> def revenue(inputs, output, offset):
        ...     price, qty = inputs
        ...     output[:] = price * qty
        >>> t.apply(revenue, columns=["price", "qty"], dtype=blosc2.float64().dtype)[:]
        array([20., 15.])
        """
        if engine not in ("auto", "numpy", "jit"):
            raise ValueError(f"engine must be 'auto', 'numpy', or 'jit', got {engine!r}")
        jit = {"auto": None, "numpy": False, "jit": True}[engine]
        names = columns if columns is not None else list(self._stored_col_names)
        missing = [n for n in names if self._logical_to_physical_name(n) not in self._cols]
        if missing:
            raise ValueError(
                f"apply() only accepts stored columns, got {missing!r}. "
                f"Stored columns: {list(self._stored_col_names)!r}."
            )
        # Operands are the raw (full-capacity) storage arrays -- the same
        # inputs add_computed_column()/add_generated_column() pass to
        # lazyudf() for DSL/UDF columns -- so the live-row mask is applied
        # once, here, to the result rather than to every operand.
        operands = tuple(self._cols[self._logical_to_physical_name(name)] for name in names)
        result = blosc2.lazyudf(func, operands, dtype=dtype, jit=jit).compute()
        return result[self._valid_rows]

    def add_generated_column(  # noqa: C901
        self,
        name: str,
        *,
        values: str
        | blosc2.LazyExpr
        | blosc2.DSLKernel
        | Callable[[dict[str, Any]], blosc2.LazyExpr]
        | RowTransformer,
        dtype=None,
        create_index: bool = False,
        inputs: list[str] | None = None,
    ) -> None:
        """Add a stored generated column maintained by the table.

        A generated column is physical storage, not a virtual expression.  The
        initial values are computed for all current live rows, and later
        ``append()`` / ``extend()`` calls automatically compute values for newly
        inserted rows when source columns are provided.  If a source column is
        modified in-place, dependent generated columns are marked stale; call
        :meth:`refresh_generated_column` or :meth:`refresh_generated_columns` to
        recompute them.

        Supported signatures are::

            add_generated_column(name, *, values="price * qty", dtype=..., create_index=False)
            add_generated_column(name, *, values=lazy_expr, dtype=...)
            add_generated_column(name, *, values=dsl_kernel, inputs=["price", "qty"], dtype=...)
            add_generated_column(name, *, values=blosc2.lazyudf(dsl_kernel, (t.price, t.qty)))
            add_generated_column(name, *, values=lambda cols: cols["price"] * 1.21, dtype=...)
            add_generated_column(name, *, values=t.embedding.row_transformer.norm(axis=0), dtype=...)
            add_generated_column(name, *, values=t.image.row_transformer.mean(axis=(0, 1)),
                                 dtype=blosc2.ndarray((3,), dtype=...))

        Parameters
        ----------
        name:
            Name of the generated column to create.  It must be a valid column
            name and must not collide with an existing stored or computed
            column.
        values:
            Definition used to compute the generated values.  Accepted forms:

            * ``str``: scalar expression over stored scalar columns, e.g.
              ``"price * qty"``.  The expression must produce one scalar value
              per row.
            * :class:`blosc2.LazyExpr`: scalar lazy expression over stored
              columns of this table.  It must produce a 1-D scalar stream.
            * :func:`blosc2.dsl_kernel`-decorated kernel passed directly with
              ``inputs=[...]`` — one stored scalar column name per kernel
              parameter, bound positionally.  Produces one scalar per row.
              The kernel source is persisted and recompiled on open; appended
              rows are auto-filled and :meth:`refresh_generated_column`
              recomputes after in-place edits.
            * :class:`blosc2.LazyUDF` built from a :func:`blosc2.dsl_kernel` via
              :func:`blosc2.lazyudf` — column bindings are inferred by identity
              from the operands, so ``inputs=`` is not needed.  Accepts
              :class:`Column` accessors (``t.col1``) or raw NDArrays as
              operands.  Same persistence and auto-fill behaviour as above.
            * callable: called as ``values(self._cols)`` and must return a
              :class:`blosc2.LazyExpr` or a :class:`blosc2.LazyUDF` backed by a
              :func:`blosc2.dsl_kernel`.
            * :class:`RowTransformer`: row-wise projection/reduction bound to a
              fixed-shape ndarray column, e.g.
              ``t.embedding.row_transformer.norm(axis=0)`` or
              ``t.image.row_transformer.mean(axis=(0, 1))``.  Row transformers
              may produce either one scalar per row or one fixed-shape ndarray
              item per row.

            Expression and DSL forms currently cannot depend on computed columns
            and cannot directly consume fixed-shape ndarray columns; use a
            row-transformer for ndarray row projections/reductions.
        dtype:
            Output schema or dtype.  Scalar outputs may pass a NumPy dtype or a
            Blosc2 scalar spec such as ``blosc2.float64()``.  Fixed-shape
            ndarray outputs must pass an ndarray spec such as
            ``blosc2.ndarray((3,), dtype=blosc2.float32())`` unless the table has
            existing rows from which the output shape can be inferred.  When
            omitted, dtype and fixed-shape output shape are inferred from the
            current generated values; this is not possible for an empty table.
        create_index:
            If ``True``, create an index on the generated column immediately.
            Only scalar generated columns can be indexed; fixed-shape ndarray
            generated columns raise :class:`ValueError` when indexing is
            requested.
        inputs:
            Only used when *values* is a bare :func:`blosc2.dsl_kernel`: a list
            of stored scalar column names, one per kernel parameter, bound
            positionally.  Not needed when passing a :class:`blosc2.LazyUDF` or
            a callable — bindings are inferred from the operands in those cases.

        Examples
        --------
        Create and index a scalar generated column from a string expression::

            t.add_generated_column(
                "total",
                values="price * qty",
                dtype=blosc2.float64(),
                create_index=True,
            )

        Use a callable when normal Python composition is more convenient::

            t.add_generated_column(
                "price_with_tax",
                values=lambda cols: cols["price"] * 1.21,
                dtype=blosc2.float64(),
            )

        Generate a scalar from each fixed-shape ndarray row.  For row
        transformers, axes refer to the per-row item shape, so ``axis=0`` is the
        embedding-coordinate axis for ``item_shape=(dim,)``::

            t.add_generated_column(
                "embedding_norm",
                values=t.embedding.row_transformer.norm(axis=0, ord=2),
                dtype=blosc2.float64(),
                create_index=True,
            )

        Generate a fixed-shape ndarray value per row.  Here an image column has
        ``item_shape=(height, width, 3)`` and the generated column stores one RGB
        vector per row::

            t.add_generated_column(
                "image_mean_rgb",
                values=t.image.row_transformer.mean(axis=(0, 1)),
                dtype=blosc2.ndarray((3,), dtype=blosc2.float32()),
            )

        Generated columns are maintained on append/extend::

            t.append((new_id, new_embedding, new_image))
            assert t.embedding_norm[-1] == np.linalg.norm(new_embedding)

        If source values are changed in place, refresh dependent generated
        columns before relying on them::

            t.embedding[0] = new_embedding
            t.refresh_generated_column("embedding_norm")

        Raises
        ------
        ValueError
            If called on a view or read-only table, if *name* already exists,
            if generated output length/shape is incompatible with the table, or
            if ``create_index=True`` is requested for an ndarray generated
            column.
        TypeError
            If *values* has an unsupported form, references unsupported source
            columns, or cannot be coerced to *dtype*.
        KeyError
            If a :class:`RowTransformer` references a missing source column.
        """
        if self.base is not None:
            raise ValueError("Cannot add a generated column to a view.")
        if self._read_only:
            raise ValueError("Table is read-only (opened with mode='r').")
        _validate_column_name(name)
        if name in self._cols:
            raise ValueError(f"A stored column named {name!r} already exists.")
        if name in self._computed_cols:
            raise ValueError(f"A computed column named {name!r} already exists.")

        n_live = self.nrows
        if isinstance(values, RowTransformer):
            transformer = values
            for dep in transformer.source_columns:
                if dep not in self._cols:
                    raise KeyError(f"No source column named {dep!r}.")
                col_info = self._schema.columns_by_name[dep]
                if not self._is_ndarray_column(col_info):
                    raise TypeError(f"RowTransformer source {dep!r} is not an ndarray column.")
            generated_values = (
                transformer.evaluate_existing(self)
                if n_live
                else transformer.evaluate_batch(
                    {
                        transformer.source: np.zeros(
                            (1, *self._schema.columns_by_name[transformer.source].spec.item_shape),
                            dtype=self._cols[transformer.source].dtype,
                        )
                    }
                )[:0]
            )
            spec = self._coerce_generated_spec(dtype, generated_values)
            metadata = {
                "computed_column": None,
                "expression": None,
                "col_deps": list(transformer.source_columns),
                "dtype": np.dtype(spec.dtype),
                "transformer_kind": "row_transformer",
                "transformer": transformer.to_metadata(),
                "stale": False,
            }
        elif (desc := self._normalize_transformer(values, inputs))["kind"] == "dsl":
            kernel = desc["kernel"]
            col_deps = desc["col_deps"]
            compute_dtype = (
                np.dtype(getattr(dtype, "dtype", dtype))
                if dtype is not None
                else self._dsl_result_dtype(kernel, col_deps, None)
            )
            jit_backend = desc.get("jit_backend")
            operands = tuple(self._cols[d] for d in col_deps)
            generated_values = np.asarray(
                blosc2.lazyudf(kernel, operands, dtype=compute_dtype, jit_backend=jit_backend).compute()[:]
            )
            if generated_values.ndim != 1:
                raise TypeError("DSL generated columns must produce a 1-D scalar result.")
            generated_values = (
                generated_values[self._valid_rows[:]]
                if len(generated_values) == len(self._valid_rows)
                else generated_values
            )
            spec = self._coerce_generated_spec(dtype, generated_values)
            metadata = {
                "computed_column": None,
                "expression": None,
                "dsl_source": kernel.dsl_source,
                "col_deps": col_deps,
                "dtype": np.dtype(spec.dtype),
                "transformer_kind": "dsl",
                "stale": False,
                "jit_backend": jit_backend,
            }
        else:
            lazy, col_deps = desc["lazy"], desc["col_deps"]
            generated_values = np.asarray(lazy[:])
            if generated_values.ndim != 1:
                raise TypeError("Expression generated columns must produce a 1-D scalar result.")
            generated_values = (
                generated_values[self._valid_rows[:]]
                if len(generated_values) == len(self._valid_rows)
                else generated_values
            )
            spec = self._coerce_generated_spec(dtype, generated_values)
            metadata = {
                "computed_column": None,
                "expression": lazy.expression,
                "col_deps": col_deps,
                "dtype": np.dtype(spec.dtype),
                "transformer_kind": "expression",
                "stale": False,
            }
        if create_index and isinstance(spec, NDArraySpec):
            raise ValueError("Generated columns intended for indexing must be 1-D scalar columns.")
        generated_values = np.asarray(generated_values, dtype=spec.dtype)
        if len(generated_values) != n_live:
            raise ValueError(
                f"Generated column {name!r} produced {len(generated_values)} values, expected {n_live}."
            )
        if isinstance(spec, NDArraySpec) and generated_values.shape != (n_live, *spec.item_shape):
            raise ValueError(
                f"Generated column {name!r} expected shape {(n_live, *spec.item_shape)}, got {generated_values.shape}."
            )

        self._create_empty_stored_column(name, np.dtype(spec.dtype), spec=spec)
        self._materialized_cols[name] = metadata
        try:
            if n_live:
                self._cols[name][self._valid_rows] = generated_values
            if create_index:
                self.create_index(name)
        except Exception:
            with contextlib.suppress(Exception):
                self.drop_column(name)
            raise
        if isinstance(self._storage, FileTableStorage):
            self._storage.save_schema(self._schema_dict_with_computed())

    def add_computed_column(
        self,
        name: str,
        expr: str | blosc2.LazyExpr | blosc2.DSLKernel | Callable[[dict[str, Any]], blosc2.LazyExpr],
        *,
        dtype: np.dtype | None = None,
        inputs: list[str] | None = None,
    ) -> None:
        """Add a read-only virtual column computed from stored columns.

        A computed column has no physical storage.  It is backed by a
        :class:`blosc2.LazyExpr` and is evaluated when values are read, filtered,
        displayed, exported, or aggregated.  Because it is virtual, it is
        read-only, cannot be indexed directly, and is not supplied in
        ``append()`` / ``extend()`` inputs.  To store and optionally index a
        computed result, use :meth:`add_generated_column` or materialize an
        existing computed column with :meth:`materialize_computed_column`.

        Supported signatures are::

            add_computed_column(name, "price * qty")
            add_computed_column(name, lazy_expr)
            add_computed_column(name, dsl_kernel, inputs=["price", "qty"])
            add_computed_column(name, blosc2.lazyudf(dsl_kernel, (t.price, t.qty)))
            add_computed_column(name, lambda cols: cols["price"] * cols["qty"])

        Parameters
        ----------
        name:
            Name of the virtual computed column.  It must be a valid column name
            and must not collide with an existing stored or computed column.
        expr:
            Definition of the virtual column.  Accepted forms:

            * ``str``: scalar expression over stored scalar columns, e.g.
              ``"price * qty"``.
            * :class:`blosc2.LazyExpr`: lazy expression over stored columns of
              this table.
            * :func:`blosc2.dsl_kernel`-decorated kernel passed directly with
              ``inputs=[...]`` — one stored scalar column name per kernel
              parameter, bound positionally.  The kernel may use loops,
              ``if``/``else`` and ``where(...)``.  Its source is persisted and
              recompiled on open; the column stays virtual/unstored.
            * :class:`blosc2.LazyUDF` built from a :func:`blosc2.dsl_kernel` via
              :func:`blosc2.lazyudf` — column bindings are inferred by identity
              from the operands, so ``inputs=`` is not needed.  Accepted forms
              include ``blosc2.lazyudf(kernel, (t.col1, t.col2))`` (using
              :class:`Column` accessors) or the raw NDArray equivalents.
            * callable: called as ``expr(self._cols)`` and must return a
              :class:`blosc2.LazyExpr` or a :class:`blosc2.LazyUDF` backed by a
              :func:`blosc2.dsl_kernel`.

            DSL columns (last three forms) are persisted — their source is stored
            and recompiled on open — and may be referenced inside :meth:`where`
            predicates.

            Expressions must depend only on stored columns of this table;
            computed columns cannot depend on other computed columns in this
            version.  Fixed-shape ndarray columns are not accepted in computed
            column expressions yet.  For row-wise ndarray projections or
            reductions, use :meth:`add_generated_column` with
            ``values=t.ndarray_col.row_transformer...``.
        dtype:
            Optional dtype override for the computed values.  For expression
            forms it is inferred from the resulting :class:`blosc2.LazyExpr`
            when omitted.  For DSL forms, an omitted dtype is inferred by NumPy
            type promotion of the input column dtypes (correct for elementwise
            arithmetic kernels); pass *dtype* explicitly for kernels that change
            the type (comparisons/``where``/casts) or when the kernel has no
            column inputs.  This changes the dtype reported by the CTable column
            wrapper; it does not create physical storage.
        inputs:
            Only used when *expr* is a bare :func:`blosc2.dsl_kernel`: a list of
            stored scalar column names, one per kernel parameter, bound
            positionally (kernel parameter ``i`` ← ``inputs[i]``).  Not needed
            when passing a :class:`blosc2.LazyUDF` or a callable — bindings are
            inferred from the operands in those cases.

        Examples
        --------
        Add a computed column from a string expression and use it like a normal
        read-only column::

            t.add_computed_column("total", "price * qty")
            assert t.total[:].shape == (t.nrows,)

        Add a computed column from a callable.  The callable receives the table's
        stored column mapping::

            t.add_computed_column(
                "price_with_tax",
                lambda cols: cols["price"] * 1.21,
                dtype=np.float64,
            )

        Callable expressions can use normal Python logic while still returning a
        lazy expression::

            def total_expr(cols):
                base = cols["price"] * cols["qty"]
                return base * 1.21 if include_tax else base

            t.add_computed_column("total", total_expr)

        They are also convenient for reusable, parameterized helpers::

            def ratio(num, den):
                return lambda cols: cols[num] / cols[den]

            t.add_computed_column("margin", ratio("profit", "revenue"))

        Computed columns participate in filters and aggregates::

            expensive = t.where(t.total > 100)
            total_revenue = t.total.sum()

        Computed columns are virtual and read-only and cannot be indexed.  If
        you need to filter or sort by this value frequently, use a generated
        column instead — it is physically stored and can be indexed::

            t.add_generated_column(
                "total_stored",
                values="price * qty",
                dtype=blosc2.float64(),
                create_index=True,
            )

        Or convert an existing computed column to a stored snapshot::

            t.materialize_computed_column("total", new_name="total_stored")
            t.create_index("total_stored")

        Raises
        ------
        ValueError
            If called on a view or read-only table, if *name* already exists,
            or if an expression operand does not reference a stored column of
            this table.
        TypeError
            If *expr* has an unsupported form, does not produce a
            :class:`blosc2.LazyExpr`, references unsupported source columns, or
            if a :class:`RowTransformer` is passed.  Row transformers are only
            accepted by :meth:`add_generated_column`.
        """
        if self.base is not None:
            raise ValueError("Cannot add a computed column to a view.")
        if self._read_only:
            raise ValueError("Table is read-only (opened with mode='r').")
        _validate_column_name(name)
        if name in self._cols:
            raise ValueError(f"A stored column named {name!r} already exists.")
        if name in self._computed_cols:
            raise ValueError(f"A computed column named {name!r} already exists.")

        desc = self._normalize_transformer(expr, inputs)
        if desc["kind"] == "dsl":
            kernel = desc["kernel"]
            col_deps = desc["col_deps"]
            self._computed_cols[name] = {
                "kind": "dsl",
                "dsl_source": kernel.dsl_source,
                "kernel": kernel,
                "col_deps": col_deps,
                "dtype": self._dsl_result_dtype(kernel, col_deps, dtype),
                "jit_backend": desc.get("jit_backend"),
            }
        else:
            lazy = desc["lazy"]
            self._computed_cols[name] = {
                "kind": "expression",
                "expression": lazy.expression,
                "col_deps": desc["col_deps"],
                "lazy": lazy,
                "dtype": np.dtype(dtype) if dtype is not None else lazy.dtype,
            }
        self.col_names.append(name)
        self._col_widths[name] = max(len(name), 15)

        # Persist metadata if backed by a file store
        if isinstance(self._storage, FileTableStorage):
            self._storage.save_schema(self._schema_dict_with_computed())

    def drop_computed_column(self, name: str) -> None:
        """Remove a computed column from the table.

        Parameters
        ----------
        name:
            Name of the computed column to remove.

        Raises
        ------
        KeyError
            If *name* is not a computed column.
        ValueError
            If called on a view.
        """
        if self.base is not None:
            raise ValueError("Cannot drop a computed column from a view.")
        if name not in self._computed_cols:
            raise KeyError(
                f"{name!r} is not a computed column. Computed columns: {list(self._computed_cols)}"
            )
        del self._computed_cols[name]
        self.col_names.remove(name)
        self._col_widths.pop(name, None)

        if isinstance(self._storage, FileTableStorage):
            self._storage.save_schema(self._schema_dict_with_computed())

    @staticmethod
    def _coerce_assign_operand(value):
        """Reduce an assign() value to a form add_computed_column's transformer
        machinery accepts: a LazyExpr, DSLKernel, callable, or string."""
        if isinstance(value, NullableExpr):
            return value._expr
        if isinstance(value, Column):
            value._ensure_queryable()
            raw = value._raw_col
            return raw if isinstance(raw, blosc2.LazyExpr) else blosc2.lazyexpr(raw)
        return value

    def assign(self, **named_exprs) -> CTable:
        """Return a view with additional computed columns, without copying data.

        Each keyword argument names a new computed column; the value defines
        it, in any of the forms :meth:`add_computed_column` accepts (a string
        expression, a :class:`blosc2.LazyExpr`), plus a :class:`Column`,
        :class:`NullableExpr`, or an unbound :func:`col` expression.

        Unlike :meth:`add_computed_column`, which mutates the table in place
        and cannot be called on a view, ``assign()`` never mutates ``self``:
        it returns a new view sharing this table's column storage, with its
        own computed-column metadata. This makes it composable in a chain::

            result = (
                t.assign(profit=col("revenue") - col("cost"))[col("profit") > 0]
                .sort_by("profit", ascending=False)
                .head(10)
            )

        Values are bound against ``self``, before any of this call's new
        columns exist — so a later keyword cannot reference an earlier one
        from the same ``assign()`` call (that raises the usual unknown-column
        error). Chain two ``assign()`` calls for that::

            t2 = t.assign(profit=col("revenue") - col("cost"))
            t3 = t2.assign(margin=col("profit") / col("revenue"))

        Parameters
        ----------
        **named_exprs:
            One computed-column definition per keyword.

        Returns
        -------
        CTable
            A read-only view (see :meth:`select`) with the additional
            computed columns.

        Raises
        ------
        ValueError
            If a name collides with an existing stored or computed column.

        Examples
        --------
        >>> import blosc2
        >>> from blosc2 import col
        >>> from dataclasses import dataclass
        >>> @dataclass
        ... class Row:
        ...     revenue: float = blosc2.field(blosc2.float64())
        ...     cost: float = blosc2.field(blosc2.float64())
        >>> t = blosc2.CTable(Row, new_data=[(100.0, 40.0), (50.0, 60.0)])
        >>> t2 = t.assign(profit=col("revenue") - col("cost"))
        >>> t2.profit[:]
        array([ 60., -10.])
        """
        bound = {}
        for name, expr in named_exprs.items():
            _validate_column_name(name)
            if name in self._cols:
                raise ValueError(f"A stored column named {name!r} already exists.")
            if name in self._computed_cols or name in bound:
                raise ValueError(f"A computed column named {name!r} already exists.")
            value = expr._bind(self) if isinstance(expr, ColExpr) else expr
            bound[name] = self._coerce_assign_operand(value)

        view = CTable._make_view(self, self._valid_rows)
        view._computed_cols = dict(self._computed_cols)
        view.col_names = list(self.col_names)
        view._col_widths = dict(self._col_widths)
        for name, value in bound.items():
            desc = view._normalize_transformer(value)
            if desc["kind"] == "dsl":
                kernel = desc["kernel"]
                col_deps = desc["col_deps"]
                view._computed_cols[name] = {
                    "kind": "dsl",
                    "dsl_source": kernel.dsl_source,
                    "kernel": kernel,
                    "col_deps": col_deps,
                    "dtype": view._dsl_result_dtype(kernel, col_deps, None),
                    "jit_backend": desc.get("jit_backend"),
                }
            else:
                lazy = desc["lazy"]
                view._computed_cols[name] = {
                    "kind": "expression",
                    "expression": lazy.expression,
                    "col_deps": desc["col_deps"],
                    "lazy": lazy,
                    "dtype": lazy.dtype,
                }
            view.col_names.append(name)
            view._col_widths[name] = max(len(name), 15)
        return view

    # ------------------------------------------------------------------
    # Column / row access
    # ------------------------------------------------------------------

    @staticmethod
    def _all_strings(seq) -> bool:
        return all(isinstance(v, str) for v in seq)

    def _getitem_arraylike(self, key):
        if len(key) == 0:
            return self._run_row_logic(key)
        if getattr(key, "dtype", None) is not None:
            if key.dtype == np.bool_:
                return self._run_row_logic(key)
            if np.issubdtype(key.dtype, np.integer):
                return self._run_row_logic(key)
            if key.dtype.kind in {"U", "S"}:
                return self.select(key.tolist())
        values = key.tolist() if hasattr(key, "tolist") else list(key)
        if self._all_strings(values):
            return self.select(values)
        return self._run_row_logic(key)

    def _getitem_row_selector(self, key):
        if isinstance(key, (int, np.integer)) and not isinstance(key, (bool, np.bool_)):
            return self._materialize_row(int(key))
        if isinstance(key, slice):
            return self._run_row_logic(key)
        if isinstance(key, np.ndarray):
            return self._getitem_arraylike(key)
        if isinstance(key, list):
            if key and self._all_strings(key):
                return self.select(key)
            return self._run_row_logic(key)
        if isinstance(key, Iterable) and not isinstance(key, (str, bytes, tuple)):
            key = list(key)
            if key and self._all_strings(key):
                return self.select(key)
            return self._run_row_logic(key)
        raise TypeError(
            "Row selectors must be an int, slice, integer array/list, or boolean mask; "
            f"got {type(key).__name__}"
        )

    def _structured_array_dtype(self) -> np.dtype:
        fields = []
        for name in self.col_names:
            col_info = self._schema.columns_by_name.get(name)
            if col_info is None:
                dtype = np.asarray(self[name][:0]).dtype
            elif (
                self._is_list_column(col_info)
                or self._is_varlen_scalar_column(col_info)
                or self._is_dictionary_column(col_info)
            ):
                dtype = np.dtype(object)
            elif self._is_ndarray_column(col_info):
                fields.append((name, col_info.dtype, col_info.spec.item_shape))
                continue
            else:
                dtype = col_info.dtype if col_info.dtype is not None else np.dtype(object)
            fields.append((name, dtype))
        return np.dtype(fields)

    def __array__(self, dtype=None, copy=None):
        arr = np.empty(self.nrows, dtype=self._structured_array_dtype())
        for name in self.col_names:
            values = self[name][:]
            target_dtype = arr.dtype.fields[name][0]
            if target_dtype == np.dtype(object) and isinstance(values, np.ndarray):
                values = values.tolist()
            arr[name] = values
        if dtype is not None:
            arr = arr.astype(dtype, copy=True if copy is None else copy)
        return arr.copy() if copy else arr

    def _logical_to_physical_name(self, name: str) -> str:
        """Resolve a user/logical column path to a stored physical column name."""
        if name in self._cols or name in self._computed_cols:
            return name
        nested = self._schema.metadata.get("nested") if self._schema.metadata else None
        if isinstance(nested, dict):
            mapping = nested.get("logical_to_physical")
            if isinstance(mapping, dict):
                physical = mapping.get(name)
                if isinstance(physical, str) and (physical in self._cols or physical in self._computed_cols):
                    return physical
        return name

    def _expand_logical_column_selector(self, name: str) -> list[str]:
        """Resolve one logical selector to one or more physical column names.

        If *name* points to a scalar leaf, returns ``[leaf]``. If it points to
        a struct-like prefix (e.g. ``"trip"``), expands to descendant leaves.
        """
        physical = self._logical_to_physical_name(name)
        if physical in self._cols or physical in self._computed_cols:
            return [physical]
        prefix_parts = split_field_path(physical)
        expanded = [
            col for col in self.col_names if split_field_path(col)[: len(prefix_parts)] == prefix_parts
        ]
        if expanded:
            return expanded
        return [physical]

    def __getitem__(self, key):
        """Type-driven indexing for columns, rows, projections, and filters.

        Supported keys are:

        - ``str``: return a :class:`Column` when it matches a stored or computed
          column name; otherwise evaluate it as a boolean expression via
          :meth:`where`.  Dotted names (e.g. ``"trip.begin.lon"``) select
          nested leaf columns directly; a struct-prefix name
          (e.g. ``"trip.begin"``) that matches multiple descendant leaves returns
          a :class:`_StructPathColumn` view.  This item-access form is the
          canonical way to access columns and works for every column name,
          including names that are not valid Python identifiers or that collide
          with existing :class:`CTable` attributes or methods.
        - boolean :class:`blosc2.LazyExpr` or :class:`blosc2.NDArray`: return the
          same filtered view as :meth:`where`, e.g. ``t[t.temperature_f > 70]``.
        - ``int``: return one live row as a namedtuple-like object.
        - ``slice``: return a row-range view.
        - integer array/list: return a gathered-row view.
        - boolean NumPy array/list: return a boolean-mask filtered view.
        - string list: return a column-projection view, equivalent to
          :meth:`select`.

        Examples
        --------
        Access columns and rows::

            temps = t["temperature"]
            first = t[0]
            view = t[10:20]

        Filter rows with a string expression, a stored-column expression, or a
        computed-column expression::

            warm = t["temperature > 20"]
            warm_active = t[(t.temperature > 20) & t.active]
            hot_fahrenheit = t[t.temperature_f > 70]

        Project columns::

            slim = t[["sensor_id", "temperature_f"]]

        Access a nested leaf column with a dotted name or an attribute chain::

            lons = t["trip.begin.lon"]   # Column for the nested leaf
            lons = t.trip.begin.lon      # equivalent attribute-chain form

        Attribute access is only a convenience fallback.  If a column name is
        not a valid identifier, or if it conflicts with an existing table
        attribute or method such as ``nrows``, ``where`` or ``sort_by``, use item
        access instead::

            col = t["where"]             # column named "where"
            method = t.where             # CTable.where method

        Row-range, gathered-row, boolean-mask, sorted (:meth:`sort_by`), and
        column-projection results are all lightweight **views**: they share
        physical storage with the base table instead of copying it. Views
        are read-only — assigning into a value returned by indexing a view
        raises ``ValueError``; use :meth:`take` or :meth:`copy` to obtain an
        independent, writable table. Mutating the base table while a view
        exists leaves the view's row mask frozen at the time the view was
        created, so the view may go stale (it will not see rows appended to
        the base afterwards, and may still reference rows later deleted from
        the base).
        """
        if isinstance(key, ColExpr):
            key = key._bind(self)
        if isinstance(key, str):
            physical = self._logical_to_physical_name(key)
            if physical in self._cols or physical in self._computed_cols:
                return Column(self, physical)
            expanded = self._expand_logical_column_selector(key)
            cc = self._schema.columns_by_name.get(physical)
            if len(expanded) > 1 or (expanded and cc is not None and isinstance(cc.spec, StructSpec)):
                return _StructPathColumn(self, physical, expanded)
            return self.where(key)
        if isinstance(key, (blosc2.NDArray, blosc2.LazyExpr)) and getattr(key, "dtype", None) == np.bool_:
            return self.where(key)
        if isinstance(key, tuple):
            raise TypeError("Tuple indexing is not supported for CTable in V1")
        return self._getitem_row_selector(key)

    def __setitem__(self, key: str, value) -> None:
        """Overwrite all live rows of a stored column.

        ``t["col"] = arr`` is equivalent to ``t["col"][:] = arr``.  *value*
        may be any array-like accepted by :meth:`Column.__setitem__`, including
        a :class:`blosc2.NDArray` (written chunk-by-chunk without full
        decompression when there are no deleted rows).

        Raises ``KeyError`` if *key* is not a stored column name, and
        ``ValueError`` if the table is read-only or a view.

        Examples
        --------
        >>> import blosc2
        >>> from dataclasses import dataclass
        >>> import numpy as np
        >>> @dataclass
        ... class Row:
        ...     price: float = blosc2.field(blosc2.float64())
        ...     embedding: object = blosc2.field(blosc2.ndarray((4,), dtype=blosc2.float32()))
        >>> t = blosc2.CTable(Row, new_data=[
        ...     (0.0, np.zeros(4, dtype=np.float32)),
        ...     (0.0, np.zeros(4, dtype=np.float32)),
        ...     (0.0, np.zeros(4, dtype=np.float32)),
        ... ])

        Overwrite a scalar column from a NumPy array:

        >>> t["price"] = np.array([1.1, 2.2, 3.3])
        >>> t["price"][:]
        array([1.1, 2.2, 3.3])

        Overwrite from a compressed :class:`blosc2.NDArray` without loading the
        full array into memory:

        >>> prices = blosc2.array([1.1, 2.2, 3.3])
        >>> t["price"] = prices
        >>> t["price"][:]
        array([1.1, 2.2, 3.3])

        Overwrite a fixed-shape ndarray column (e.g. embeddings):

        >>> data = blosc2.arange(12, dtype=np.float32, shape=(3, 4))
        >>> t["embedding"] = data
        >>> t["embedding"][:]
        array([[ 0.,  1.,  2.,  3.],
               [ 4.,  5.,  6.,  7.],
               [ 8.,  9., 10., 11.]], dtype=float32)
        """
        if not isinstance(key, str):
            raise TypeError(
                f"CTable.__setitem__ only accepts a column name string, got {type(key).__name__!r}"
            )
        if self.base is not None:
            raise ValueError("Table is a view and cannot be modified.")
        if self._read_only:
            raise ValueError("Table is read-only (opened with mode='r').")
        physical = self._logical_to_physical_name(key)
        if physical not in self._cols:
            raise KeyError(f"Column {key!r} does not exist; use add_column() to add new columns")
        Column(self, physical)[:] = value

    def _nested_namespace(self, prefix: str):
        prefix_parts = split_field_path(prefix)
        for name in self.col_names:
            parts = split_field_path(name)
            if parts[: len(prefix_parts)] == prefix_parts and len(parts) > len(prefix_parts):
                return NestedColumn(self, prefix)
        return None

    def __getattr__(self, s: str):
        """Convenience fallback for attribute-style column access.

        This is called only after normal Python attribute lookup fails.  Thus
        ``t.name`` can return a column only for non-conflicting identifier-like
        column names.  For columns whose names conflict with existing CTable
        attributes/methods, or are not valid identifiers, use the canonical item
        access form ``t["name"]``.
        """
        physical = self._logical_to_physical_name(s)
        if physical in self._cols or physical in self._computed_cols:
            return Column(self, physical)
        ns = self._nested_namespace(s)
        if ns is not None:
            return ns
        return super().__getattribute__(s)

    # ------------------------------------------------------------------
    # Compaction
    # ------------------------------------------------------------------

    def compact(self):
        """Physically rewrite every column array keeping only live rows.

        Closes the gaps left by prior :meth:`delete` calls by shuffling live
        data to the front of each column array.  The underlying NDArray
        allocations are **not resized** — each column retains its original
        capacity.  To actually reclaim memory, use :meth:`copy` with
        ``compact=True`` instead, which allocates fresh arrays sized to the
        live row count.  All existing indexes are dropped and must be
        recreated afterwards.  Raises ``ValueError`` if the table is
        read-only or a view.
        """
        if self._read_only:
            raise ValueError("Table is read-only (opened with mode='r').")
        if self.base is not None:
            raise ValueError("Cannot compact a view.")
        if self._last_pos is not None and self._last_pos == self._n_rows:
            return
        # Compaction rewrites physical column layout; incremental summaries no
        # longer correspond to it.
        self._invalidate_all_summary_accumulators()
        self._flush_varlen_columns()
        real_poss = blosc2.where(self._valid_rows, np.array(range(len(self._valid_rows)))).compute()
        for col in self._schema.columns:
            name = col.name
            v = self._cols[name]
            if self._is_list_column(col):
                compacted = [v[int(pos)] for pos in real_poss[: self._n_rows]]
                replacement = ListArray(spec=col.spec)
                replacement.extend(compacted)
                replacement.flush()
                self._cols[name] = replacement
                continue
            if self._is_utf8_column(col):
                # Clustered fancy-index gather, then a bulk rewrite through the
                # existing backing arrays so a store-backed column stays
                # persistent on disk.
                v.set_all(v[real_poss[: self._n_rows]])
                continue
            if self._is_varlen_scalar_column(col):
                compacted = [v[int(pos)] for pos in real_poss[: self._n_rows]]
                replacement = _ScalarVarLenArray(col.spec)
                replacement.extend(compacted)
                replacement.flush()
                self._cols[name] = replacement
                continue
            if self._is_dictionary_column(col):
                # Keep dictionary values intact; just compact the codes.
                live_codes = v.codes[real_poss[: self._n_rows]]
                v.codes[: self._n_rows] = live_codes
                continue
            start = 0
            block_size = self._valid_rows.blocks[0]
            end = min(block_size, self._n_rows)
            while start < end:
                v[start:end] = v[real_poss[start:end]]
                start += block_size
                end = min(end + block_size, self._n_rows)

        self._valid_rows[: self._n_rows] = True
        self._valid_rows[self._n_rows :] = False
        self._last_pos = self._n_rows
        self._mark_all_indexes_stale()

    @staticmethod
    def _column_selector_name(value: Any) -> str:
        """Return the column name represented by a string or Column-like selector."""
        name = getattr(value, "_col_name", value)
        if not isinstance(name, str):
            raise TypeError(f"Expected a column name or Column object, got {type(value)!r}")
        return name

    def _normalise_sort_keys(
        self,
        cols: str | list[str],
        ascending: bool | list[bool],
    ) -> tuple[list[str], list[bool]]:
        """Validate and normalise sort key arguments; return (cols, ascending)."""
        if isinstance(cols, str) or isinstance(getattr(cols, "_col_name", None), str):
            cols = [self._column_selector_name(cols)]
        else:
            cols = [self._column_selector_name(col) for col in cols]

        resolved_cols: list[str] = []
        for name in cols:
            expanded = self._expand_logical_column_selector(name)
            if len(expanded) != 1:
                raise ValueError(
                    f"Sort key {name!r} resolves to multiple columns {expanded!r}; please choose a leaf column."
                )
            resolved_cols.append(expanded[0])
        cols = resolved_cols
        if isinstance(ascending, bool):
            ascending = [ascending] * len(cols)
        if len(cols) != len(ascending):
            raise ValueError(
                f"'ascending' must have the same length as 'cols' ({len(cols)}), got {len(ascending)}."
            )
        for name in cols:
            if name not in self._cols and name not in self._computed_cols:
                raise KeyError(f"No column named {name!r}. Available: {self.col_names}")
            self._ensure_generated_column_not_stale(name)
            col_info = self._schema.columns_by_name.get(name)
            if col_info is not None and self._is_ndarray_column(col_info):
                raise TypeError(
                    f"Cannot sort by ndarray column {name!r} with per-row shape {col_info.spec.item_shape}. "
                    "Materialize a scalar generated column first, e.g. embedding_norm or embedding_max."
                )
            dtype = self._col_dtype(name)
            if dtype is None:
                if col_info is not None and isinstance(
                    col_info.spec, (VLStringSpec, VLBytesSpec, StructSpec, ObjectSpec)
                ):
                    raise TypeError(
                        f"Column {name!r} is a varlen scalar column and does not support sort ordering."
                    )
                if col_info is not None and self._is_dictionary_column(col_info):
                    pass  # dictionary columns: sorting supported (decoded strings)
                else:
                    raise TypeError(
                        f"Column {name!r} is a list column and does not support sort ordering in V1."
                    )
            if np.issubdtype(dtype, np.complexfloating):
                raise TypeError(
                    f"Column {name!r} has complex dtype {dtype} which does not support ordering."
                )
        return cols, ascending

    def _sorted_positions_from_full_index(self, name: str, ascending: bool) -> np.ndarray | None:  # noqa: C901
        """Return live physical positions from a matching FULL index, if available.

        Reads the pre-sorted positions sidecar directly rather than going through
        the ordered_indices query machinery, which is optimised for selective range
        queries and is much slower for full-table streaming.
        """
        root = self._root_table
        catalog = root._get_index_catalog()
        descriptor = None

        null_value = None
        null_code = None
        is_dict_rank = False
        if name in root._cols:
            col_info = root._schema.columns_by_name.get(name)
            if col_info is not None:
                null_value = getattr(col_info.spec, "null_value", None)
                if isinstance(col_info.spec, DictionarySpec):
                    null_code = col_info.spec.null_code
            descriptor = catalog.get(name)
            if descriptor is None or descriptor.get("kind") != "full" or descriptor.get("stale", False):
                descriptor = None
            else:
                dict_rank_meta = descriptor.get("full", {}).get("dict_rank")
                if dict_rank_meta is not None:
                    if self._dict_rank_index_stale(name, dict_rank_meta):
                        descriptor = None  # ranks no longer match dictionary → lexsort
                    else:
                        is_dict_rank = True
        elif name in root._computed_cols:
            cc = root._computed_cols[name]
            for _lookup_key, candidate in catalog.items():
                target = candidate.get("target") or {}
                if (
                    target.get("source") == "expression"
                    and candidate.get("kind") == "full"
                    and not candidate.get("stale", False)
                    and target.get("expression_key") == cc.get("expression")
                    and list(target.get("dependencies", [])) == list(cc["col_deps"])
                ):
                    descriptor = candidate
                    break
        if descriptor is None:
            return None

        positions_path = descriptor.get("full", {}).get("positions_path")

        # Read pre-sorted positions directly — bypasses the ordered_indices query
        # machinery which is built for selective range queries and is ~70x slower
        # for full-table streaming.
        if positions_path is not None:
            # Persistent table: positions live in a sidecar .b2nd file.  Use the
            # sidecar opener so .b2z (zip) stores are read at their zip offset —
            # blosc2.open() would look for a standalone file that isn't there.
            from blosc2.indexing import _open_sidecar_file

            positions_nd = _open_sidecar_file(positions_path)
        else:
            # In-memory table: positions live in the sidecar handle cache.
            from blosc2.indexing import _SIDECAR_HANDLE_CACHE, _sidecar_handle_cache_key

            target_arr = root._cols.get(name)
            if target_arr is None:
                return None
            token = descriptor["token"]
            cache_key = _sidecar_handle_cache_key(target_arr, token, "full", "positions")
            positions_nd = _SIDECAR_HANDLE_CACHE.get(cache_key)
            if positions_nd is None:
                return None

        positions = np.asarray(positions_nd[:], dtype=np.int64)
        total = len(root._valid_rows)
        # Index sidecars can carry padding positions beyond the live range, so
        # the bounds clip always runs — but the ``.all()`` check skips the copy
        # (and a 24M-element temporary) when there is nothing to clip.
        in_bounds = (positions >= 0) & (positions < total)
        if not bool(in_bounds.all()):
            positions = positions[in_bounds]
        del in_bounds
        # Validity filtering only matters when the table has gaps (deleted rows);
        # for a compact table every clipped position is already live.
        if root._n_rows is None or root._n_rows != total:
            valid = root._valid_rows[:]
            positions = positions[valid[positions]]
        if self is not root:
            current_valid = self._valid_rows[:]
            positions = positions[current_valid[positions]]

        if is_dict_rank:
            # Dict-rank index: positions sorted by rank (int32), nulls have sentinel null_rank.
            # Partition null rows using codes (int32), not decoded strings.
            codes = np.asarray(root._cols[name].codes[:], dtype=np.int32)
            null_phys = codes == null_code
            del codes
            if null_phys.any():
                is_null = null_phys[positions]
                del null_phys
                nulls = positions[is_null]
                nonnull = positions[~is_null]
                del is_null, positions
                if not ascending:
                    nonnull = nonnull[::-1]
                return np.concatenate([nonnull, nulls])
            # No nulls: fall through to simple reverse
        elif null_value is not None:
            # The index sorts by raw value, but sort_by's contract is nulls-last.
            # Partition explicitly so it holds for any sentinel (NaN sorts last,
            # an integer sentinel like INT64_MIN sorts first) and either order.
            # Free each 24M-element temporary as soon as it is consumed to keep
            # peak memory near the size of the permutation itself.
            raw = np.asarray(root._cols[name][:])
            if isinstance(null_value, float) and np.isnan(null_value):
                null_phys = np.isnan(raw)
            else:
                null_phys = raw == null_value
            del raw
            if null_phys.any():
                is_null = null_phys[positions]
                del null_phys
                nulls = positions[is_null]
                nonnull = positions[~is_null]
                del is_null, positions
                if not ascending:
                    nonnull = nonnull[::-1]
                return np.concatenate([nonnull, nulls])

        if not ascending:
            positions = positions[::-1]
        return positions

    def _build_lex_keys(
        self,
        cols: list[str],
        ascending: list[bool],
        live_pos: np.ndarray,
        n: int,
    ) -> list[np.ndarray]:
        """Build the key list for np.lexsort (innermost = last = primary key).

        For nullable columns a null-indicator key (0=non-null, 1=null) is
        inserted immediately after the value key, making it more significant.
        This ensures nulls sort last regardless of ascending/descending order.
        """
        lex_keys = []
        for name, asc in zip(reversed(cols), reversed(ascending), strict=True):
            cc = self._computed_cols.get(name)
            col_info = self._schema.columns_by_name.get(name)
            is_dict_col = False
            if cc is not None:
                # Materialise computed column values at live positions
                raw = np.asarray(self._build_computed_lazy(cc)[:])[live_pos]
            else:
                is_dict_col = col_info is not None and self._is_dictionary_column(col_info)
                if is_dict_col:
                    # Sort dictionary columns by decoded string values.
                    decoded = self._cols[name][live_pos]
                    raw = np.array(decoded, dtype=object)
                    # Replace None with placeholder so lexsort never compares None.
                    # Null indicator key (below) already places nulls last.
                    raw[raw == None] = ""  # noqa: E711
                else:
                    raw = self._cols[name][live_pos]
            nv = getattr(col_info.spec, "null_value", None) if col_info else None

            # Value key
            if not asc:
                if raw.dtype.kind in "USOT":
                    # strings can't be negated — invert via rank
                    rank = np.argsort(np.argsort(raw, kind="stable"), kind="stable")
                    lex_keys.append((n - 1 - rank).astype(np.intp))
                elif np.issubdtype(raw.dtype, np.unsignedinteger):
                    lex_keys.append(-raw.astype(np.int64))
                else:
                    lex_keys.append(-raw)
            else:
                lex_keys.append(raw)

            # Null indicator key — more significant than the value key above,
            # so nulls always sort last (0 before 1 → non-null before null).
            if is_dict_col and col_info.spec.nullable:
                null_code = col_info.spec.null_code
                codes_at_pos = np.asarray(self._cols[name].codes[live_pos], dtype=np.int32)
                null_ind = (codes_at_pos == null_code).astype(np.intp)
                lex_keys.append(null_ind)
            elif nv is not None:
                if isinstance(nv, float) and np.isnan(nv):
                    null_ind = np.isnan(raw).astype(np.intp)
                else:
                    null_ind = (raw == nv).astype(np.intp)
                lex_keys.append(null_ind)

        return lex_keys

    def sort_by(
        self,
        cols: str | list[str],
        ascending: bool | list[bool] = True,
        *,
        inplace: bool = False,
        view: bool = False,
    ) -> CTable:
        """Return the table sorted by one or more columns.

        By default this materialises a new in-memory copy of the sorted rows.
        Pass ``view=True`` to instead get a lightweight **sorted view** that
        shares the parent's column data and gathers rows on demand in sorted
        order — no whole-table copy.  This is ideal for reading a sorted slice
        of a large persistent table (e.g. ``t.sort_by("col", view=True)[:10]``).

        Parameters
        ----------
        cols:
            Column name or list of column names to sort by.  When multiple
            columns are given, the first is the primary key, the second is
            the tiebreaker, and so on.  For tables with **nested (dotted)
            column names**, pass the dotted leaf name directly::

                t.sort_by("trip.begin.lon")
                t.sort_by(["trip.begin.lon", "payment.fare"], ascending=[True, False])

        ascending:
            Sort direction.  A single bool applies to all keys; a list must
            have the same length as *cols*.
        inplace:
            If ``True``, rewrite the physical data in place and return
            ``self`` (like :meth:`compact` but sorted).  If ``False``
            (default), return a new in-memory CTable leaving this one
            untouched.
        view:
            If ``True``, return a zero-copy sorted **view** over this table
            instead of materialising a copy: it shares the parent's columns and
            stores only the sort permutation, gathering rows on demand in sorted
            order.  Slicing the view (``sv[start:stop:step]``) keeps the sorted
            order and touches only the rows read.  A single-column sort backed by
            a non-stale ``FULL`` index reuses its pre-sorted positions (no sort at
            read time); otherwise only the sort-key column(s) are materialised to
            build the permutation — never the whole table.  Mutually exclusive
            with ``inplace``.  Sorting an existing view is always lazy regardless
            of this flag.

        Raises
        ------
        ValueError
            If called on a view or a read-only table when ``inplace=True``, or if
            both ``inplace`` and ``view`` are ``True``.
        KeyError
            If any column name is not found.
        TypeError
            If a column used as a sort key does not support ordering
            (e.g. complex numbers).
        """
        if inplace and view:
            raise ValueError("inplace=True and view=True are mutually exclusive.")
        if self.base is not None and inplace:
            raise ValueError(
                "Cannot sort a view inplace (would modify shared column data). Use sort_by(inplace=False) to get a sorted copy."
            )
        if inplace and self._read_only:
            raise ValueError("Table is read-only (opened with mode='r').")

        cols, ascending = self._normalise_sort_keys(cols, ascending)

        # Live physical positions.  Scan the validity NDArray chunk-wise to avoid
        # materialising the whole mask as a single NumPy array.
        live_pos = self._live_positions_from_valid_rows_chunks()
        n = len(live_pos)

        if n == 0:
            if inplace:
                return self
            return self._empty_copy()

        sorted_pos = None
        if len(cols) == 1:
            sorted_pos = self._sorted_positions_from_full_index(cols[0], ascending[0])
            if sorted_pos is not None and len(sorted_pos) != n:
                sorted_pos = None

        if sorted_pos is None:
            order = np.lexsort(self._build_lex_keys(cols, ascending, live_pos, n))
            sorted_pos = live_pos[order]

        if inplace:
            self._sort_by_inplace(sorted_pos, n)
            return self

        # When sorting a view, return a new lazy sorted view rather than
        # eagerly materialising all columns.  The new view shares _cols with
        # the base and stores the sorted physical positions in
        # _cached_live_positions.  Column reads and row iteration on the result
        # use those positions directly, so columns are fetched on demand and in
        # the correct sorted order — identical performance to pre-projecting
        # with columns= before calling sort_by.
        if self.base is not None or view:
            result = CTable._make_view(self, self._valid_rows)
            result._cached_live_positions = sorted_pos
            result._n_rows = n
            return result

        return self._sorted_copy_from_positions(sorted_pos, n)

    def sorted_slice(self, col: str, key: slice, *, ascending: bool = True) -> CTable:
        """Return rows ``key`` in ``col``-sorted order, reading only the slice window.

        Like ``sort_by(col, ascending=ascending, view=True)[key]`` but, when ``col``
        has a usable FULL index, it reads just the needed window of the index's
        position sidecar instead of materialising the whole 24M-row permutation —
        ideal for small slices (top/bottom *k*).  Falls back to the full sorted
        view (same result) whenever the window path does not apply.
        """
        if not isinstance(key, slice):
            raise TypeError("sorted_slice expects a slice")
        pos = self._sorted_slice_positions(col, ascending, key)
        if pos is None:
            return self.sort_by(col, ascending=ascending, view=True)[key]
        return self._view_from_positions(pos)

    def _sorted_slice_positions(self, name: str, ascending: bool, key: slice) -> np.ndarray | None:
        """Physical positions for the sorted slice ``key``, reading only the window.

        Returns ``None`` (so the caller falls back to the full path) unless this is
        a base table with a non-stale, persistent FULL index over a compact,
        unpadded column indexed by a numeric (or null) sentinel.
        """
        if self.base is not None:
            return None
        descriptor = self._get_index_catalog().get(name)
        if not descriptor or descriptor.get("kind") != "full" or descriptor.get("stale", False):
            return None
        full = descriptor.get("full") or {}
        positions_path = full.get("positions_path")
        if positions_path is None:  # in-memory sidecar: not worth a partial-read path
            return None

        n = self._n_rows
        total = len(self._valid_rows)
        if n is None or n != total:  # deletions → positions are not a clean permutation
            return None

        col_info = self._schema.columns_by_name.get(name)
        null_value = getattr(col_info.spec, "null_value", None) if col_info is not None else None
        # Dict-rank index: use null_rank (int32) as sentinel for null-block location.
        dict_rank = full.get("dict_rank")
        if dict_rank is not None:
            if self._dict_rank_index_stale(name, dict_rank):
                return None  # ranks no longer match dictionary → lexsort
            null_value = dict_rank["null_rank"]
        # Numeric / NaN / string sentinels keep the null rows in one contiguous block
        # once sorted; other non-numeric sentinels (e.g. object) would need a
        # different locator.
        if null_value is not None and not isinstance(null_value, (int, float, str, bytes)):
            return None

        from blosc2.indexing import _open_sidecar_file

        pnd = _open_sidecar_file(positions_path)
        if len(pnd) != total:  # capacity padding → window read would be wrong
            return None

        result_idx = np.arange(*key.indices(n), dtype=np.int64)
        if result_idx.size == 0:
            return np.empty(0, dtype=np.int64)

        # Locate the (contiguous) null block [null_lo, null_hi) in the sorted order.
        null_lo, null_hi = self._null_block_bounds(full, null_value, n)

        # Map each requested result index to its index in the sorted sidecar.  The
        # nulls-last order is the non-null rows (forward or reversed) followed by
        # the null block, where the non-null rows are everything outside the block.
        sidecar_idx = np.empty_like(result_idx)
        if ascending:
            len_below = null_lo  # non-null rows sorted below the null block
            len_above = n - null_hi  # non-null rows sorted above it
            below = result_idx < len_below
            above = (result_idx >= len_below) & (result_idx < len_below + len_above)
            nulls = result_idx >= len_below + len_above
            sidecar_idx[below] = result_idx[below]
            sidecar_idx[above] = null_hi + (result_idx[above] - len_below)
            sidecar_idx[nulls] = null_lo + (result_idx[nulls] - len_below - len_above)
        else:
            len_above = n - null_hi  # largest non-null rows come first
            len_below = null_lo
            above = result_idx < len_above
            below = (result_idx >= len_above) & (result_idx < len_above + len_below)
            nulls = result_idx >= len_above + len_below
            sidecar_idx[above] = (n - 1) - result_idx[above]
            sidecar_idx[below] = (null_lo - 1) - (result_idx[below] - len_above)
            sidecar_idx[nulls] = null_lo + (result_idx[nulls] - len_above - len_below)

        lo = int(sidecar_idx.min())
        hi = int(sidecar_idx.max()) + 1
        window = np.asarray(pnd[lo:hi], dtype=np.int64)
        return window[sidecar_idx - lo]

    def _null_block_bounds(self, full: dict, null_value, n: int) -> tuple[int, int]:
        """Return ``[null_lo, null_hi)``: the null rows' span in the sorted sidecar.

        Empty (``null_lo == null_hi``) when the column is non-nullable or has no
        null rows.  Reads only a handful of sidecar blocks, never the whole array.
        """
        if null_value is None:
            return n, n
        from blosc2.indexing import _open_sidecar_file

        vnd = _open_sidecar_file(full["values_path"])
        if isinstance(null_value, float) and np.isnan(null_value):
            # NaN sorts last and breaks ordered comparisons, so count the trailing
            # block directly, one chunk at a time (peak memory = a single chunk).
            chunk = int(vnd.chunks[0]) if vnd.chunks else len(vnd)
            count = 0
            hi = len(vnd)
            while hi > 0:
                lo = max(0, hi - chunk)
                block = np.isnan(np.asarray(vnd[lo:hi]))
                count += int(block.sum())
                if not block.all():  # reached the non-null region
                    break
                hi = lo
            return n - count, n
        # Ordinary value: the block is wherever the sentinel sorts.  Bisect the
        # sorted values, reading one block per probe.
        return (
            self._sidecar_bisect(vnd, null_value, "left"),
            self._sidecar_bisect(vnd, null_value, "right"),
        )

    @staticmethod
    def _sidecar_bisect(vnd: blosc2.NDArray, value, side: str) -> int:
        """``np.searchsorted`` over an ascending sidecar, reading one element/probe."""
        lo, hi = 0, len(vnd)
        while lo < hi:
            mid = (lo + hi) // 2
            v = vnd[mid : mid + 1][0]
            if (v < value) if side == "left" else (v <= value):
                lo = mid + 1
            else:
                hi = mid
        return lo

    def _sorted_small_copy_from_live_positions(
        self, cols: list[str], ascending: list[bool], live_pos: np.ndarray, n: int
    ) -> CTable:
        """Materialise and sort a small filtered view, avoiding a second gather of sort keys."""
        gathered = {}
        for col in self._schema.columns:
            arr = self._cols[col.name]
            if self._is_dictionary_column(col):
                gathered[col.name] = arr.codes[live_pos]
            else:
                gathered[col.name] = arr[live_pos]

        lex_keys = []
        for name, asc in zip(reversed(cols), reversed(ascending), strict=True):
            col_info = self._schema.columns_by_name.get(name)
            is_dict_col = col_info is not None and self._is_dictionary_column(col_info)
            if is_dict_col:
                raw = np.array(self._cols[name][live_pos], dtype=object)
                # Replace None with placeholder so lexsort never compares None.
                raw[raw == None] = ""  # noqa: E711
            else:
                raw = gathered[name]

            if not asc:
                if raw.dtype.kind in "USO":
                    rank = np.argsort(np.argsort(raw, kind="stable"), kind="stable")
                    lex_keys.append((n - 1 - rank).astype(np.intp))
                elif np.issubdtype(raw.dtype, np.unsignedinteger):
                    lex_keys.append(-raw.astype(np.int64))
                else:
                    lex_keys.append(-raw)
            else:
                lex_keys.append(raw)

            if is_dict_col and col_info.spec.nullable:
                null_code = col_info.spec.null_code
                codes_at_pos = np.asarray(self._cols[name].codes[live_pos], dtype=np.int32)
                null_ind = (codes_at_pos == null_code).astype(np.intp)
                lex_keys.append(null_ind)
            else:
                nv = getattr(col_info.spec, "null_value", None) if col_info else None
                if nv is not None:
                    if isinstance(nv, float) and np.isnan(nv):
                        null_ind = np.isnan(raw).astype(np.intp)
                    else:
                        null_ind = (raw == nv).astype(np.intp)
                    lex_keys.append(null_ind)

        order = np.lexsort(lex_keys)
        result = self._empty_copy(capacity=n)
        for col in self._schema.columns:
            col_name = col.name
            if self._is_dictionary_column(col):
                for v in self._cols[col_name].dictionary:
                    result._cols[col_name].encode(v)
                result._cols[col_name].codes[:n] = gathered[col_name][order]
            else:
                result._cols[col_name][:n] = gathered[col_name][order]
        result._valid_rows[:n] = True
        result._valid_rows[n:] = False
        result._n_rows = n
        result._last_pos = n
        return result

    def _sort_by_inplace(self, sorted_pos: np.ndarray, n: int) -> None:
        for col in self._schema.columns:
            arr = self._cols[col.name]
            if self._is_list_column(col):
                new_arr = ListArray(spec=col.spec)
                new_arr.extend((arr[int(pos)] for pos in sorted_pos), validate=False)
                new_arr.flush()
                self._cols[col.name] = new_arr
            elif self._is_dictionary_column(col):
                sorted_codes = arr.codes[sorted_pos]
                arr.codes[:n] = sorted_codes
            elif self._is_utf8_column(col):
                # Bulk-rewrite through the existing backing arrays so a
                # store-backed column stays persistent on disk.
                arr.set_all(arr[sorted_pos])
            else:
                arr[:n] = arr[sorted_pos]
        self._valid_rows[:n] = True
        self._valid_rows[n:] = False
        self._n_rows = n
        self._last_pos = n
        self._mark_all_indexes_stale()

    def _sorted_copy_from_positions(self, sorted_pos: np.ndarray, n: int) -> CTable:
        # Build a new in-memory table with the sorted rows
        result = self._empty_copy()
        for col in self._schema.columns:
            col_name = col.name
            arr = self._cols[col_name]
            if self._is_list_column(col):
                result._cols[col_name].extend((arr[int(pos)] for pos in sorted_pos), validate=False)
                result._cols[col_name].flush()
            elif self._is_dictionary_column(col):
                # Copy dictionary values, then sorted codes.
                for v in arr.dictionary:
                    result._cols[col_name].encode(v)
                sorted_codes = arr.codes[sorted_pos]
                result._cols[col_name].codes[:n] = sorted_codes
            elif self._is_utf8_column(col):
                result._cols[col_name].set_all(arr[sorted_pos])
            else:
                result._cols[col_name][:n] = arr[sorted_pos]
        result._valid_rows[:n] = True
        result._valid_rows[n:] = False
        result._n_rows = n
        result._last_pos = n
        return result

    def copy(  # noqa: C901
        self,
        compact: bool = True,
        *,
        urlpath: str | os.PathLike[str] | None = None,
        overwrite: bool = False,
        chunks: int | tuple[int, ...] | None = None,
        blocks: int | tuple[int, ...] | None = None,
        cparams: dict[str, Any] | None = None,
    ) -> CTable:
        """Return a new standalone copy of this table.

        This is the only operation that truly reclaims memory: when
        ``compact=True`` the new table allocates fresh arrays sized exactly
        to the live row count, discarding all deleted-row gaps and unused
        capacity.

        Parameters
        ----------
        compact:
            If ``True`` (default), only live (non-deleted) rows are copied.
            The result is a dense table with no tombstones and no parent
            dependency — ideal for materialising a filtered view or freeing
            memory after heavy deletions.
            If ``False``, all physical slots are copied including deleted gaps,
            preserving the tombstone state exactly for in-memory copies.
        urlpath:
            Destination path for a persistent copy.  The ``.b2z`` extension
            selects a compact zip-backed store; any other path uses a
            directory-backed store.  A ``.b2d`` suffix is recommended for
            directory-backed stores.  If ``None`` (default), return an
            in-memory copy.
        overwrite:
            If ``True``, replace an existing persistent destination.
        chunks:
            Chunk size (in items) to use for all scalar columns in the copy.
            Overrides the chunk size inherited from the source schema.  Pass an
            ``int`` for a 1-D chunk or a ``tuple`` for multi-dimensional arrays.
        blocks:
            Block size (in items) to use for all scalar columns in the copy.
            Overrides the block size inherited from the source schema.  Pass an
            ``int`` for a 1-D block or a ``tuple`` for multi-dimensional arrays.
            Summary indexes are rebuilt with the new block granularity.
        cparams:
            Compression parameters (codec, clevel, …) to apply to all columns
            in the copy.  Overrides per-column and table-level settings from
            the source.
        """
        if urlpath is not None:
            urlpath = os.fspath(urlpath)
            if chunks is not None or blocks is not None or cparams is not None:
                # When storage layout changes we must go through _save_to_storage
                # directly — to_b2z/to_b2d may take the physical-pack fast path
                # which zips existing compressed leaves as-is, silently ignoring
                # any chunk/block/cparams override.
                # For views (base is not None) _save_to_storage already limits
                # iteration to self._schema.columns and self._cols, so no in-memory
                # intermediate is needed.
                _chunks = (chunks,) if isinstance(chunks, int) else chunks
                _blocks = (blocks,) if isinstance(blocks, int) else blocks
                file_storage = FileTableStorage(urlpath, "w")
                target_path = file_storage._root
                if os.path.exists(target_path):
                    if not overwrite:
                        raise ValueError(
                            f"Path {target_path!r} already exists. Use overwrite=True to replace."
                        )
                    if os.path.isdir(target_path):
                        shutil.rmtree(target_path)
                    else:
                        os.remove(target_path)
                self._save_to_storage(
                    file_storage,
                    chunks_override=_chunks,
                    blocks_override=_blocks,
                    cparams_override=cparams,
                )
                file_storage.close()
                # Open with mode="a" so _build_summary_indexes() fires automatically,
                # then re-open read-only for the caller.
                result = CTable.open(urlpath, mode="a")
                result.close()
                return CTable.open(urlpath, mode="r")
            if urlpath.endswith(".b2z"):
                self.to_b2z(urlpath, overwrite=overwrite, compact=compact)
            else:
                self.to_b2d(urlpath, overwrite=overwrite, compact=compact)
            return CTable.open(urlpath, mode="r")

        valid_np = self._valid_rows[:]
        live_pos = np.where(valid_np)[0]
        n_live = len(live_pos)

        if compact:
            n = n_live
        else:
            # High watermark: number of slots ever written.
            # List columns are written sequentially with no gaps — their length
            # is the exact high watermark.  For scalar-only tables fall back to
            # the last live position + 1 (writes are always sequential so no
            # deleted slot can exist beyond the last live one).
            n = 0
            for col in self._schema.columns:
                if self._is_list_column(col):
                    n = len(self._cols[col.name])
                    break
            if n == 0:
                n = int(live_pos[-1]) + 1 if n_live > 0 else 0

        # When all live positions are exactly [0, 1, …, n_live-1] a slice read is
        # ~30× faster than fancy indexing.  Check via O(1) boundary test.
        is_dense = compact and n_live > 0 and int(live_pos[0]) == 0 and int(live_pos[-1]) == n_live - 1

        _chunks = (chunks,) if isinstance(chunks, int) else chunks
        _blocks = (blocks,) if isinstance(blocks, int) else blocks
        result = self._empty_copy(
            capacity=n,
            chunks_override=_chunks,
            blocks_override=_blocks,
            cparams_override=cparams,
        )

        for col in self._schema.columns:
            col_name = col.name
            arr = self._cols[col_name]
            if self._is_list_column(col):
                src = (
                    arr[:n_live]
                    if is_dense
                    else (arr[int(pos)] for pos in live_pos)
                    if compact
                    else (arr[i] for i in range(n))
                )
                result._cols[col_name].extend(src, validate=False)
                result._cols[col_name].flush()
            elif self._is_varlen_scalar_column(col):
                # _ScalarVarLenArray.__setitem__ only accepts a single int index
                # (mirroring row-wise append/extend semantics), so bulk copies
                # must go through extend(), same as list columns above.
                src = (
                    arr[:n_live]
                    if is_dense
                    else (arr[int(pos)] for pos in live_pos)
                    if compact
                    else (arr[i] for i in range(n))
                )
                result._cols[col_name].extend(src)
                result._cols[col_name].flush()
            elif self._is_dictionary_column(col):
                # Copy dictionary values, then copy (live) codes.
                for v in arr.dictionary:
                    result._cols[col_name].encode(v)
                pos_slice = (
                    np.arange(n_live, dtype=np.int64)
                    if is_dense
                    else (live_pos if compact else np.arange(n, dtype=np.int64))
                )
                raw_codes = arr.codes[pos_slice]
                result._cols[col_name].codes[:n] = raw_codes
            else:
                result._cols[col_name][:n] = (
                    arr[:n_live] if is_dense else (arr[live_pos] if compact else arr[:n])
                )

        if compact:
            result._valid_rows[:n] = True
            result._n_rows = n
            result._last_pos = n - 1 if n > 0 else None
        else:
            result._valid_rows[:n] = valid_np[:n]
            result._n_rows = n_live
            result._last_pos = None  # recomputed lazily on next append

        return result

    def _empty_copy(
        self,
        capacity: int | None = None,
        *,
        chunks_override: tuple[int, ...] | None = None,
        blocks_override: tuple[int, ...] | None = None,
        cparams_override: dict[str, Any] | None = None,
    ) -> CTable:
        """Return a new empty in-memory CTable with the same schema and capacity."""
        from blosc2 import compute_chunks_blocks

        capacity = max(capacity if capacity is not None else self._n_rows, 1)
        default_chunks, default_blocks = compute_chunks_blocks((capacity,))
        # Align fixed-size scalar columns (and the _valid_rows mask) on one
        # shared grid so lazy expressions over them take the fast_eval path.
        shared_chunks, shared_blocks, aligned_names = self._compute_aligned_grid(
            self._schema.columns, capacity
        )
        mem_storage = InMemoryTableStorage()

        new_valid = mem_storage.create_valid_rows(
            shape=(capacity,),
            chunks=shared_chunks if shared_chunks is not None else default_chunks,
            blocks=shared_blocks if shared_blocks is not None else default_blocks,
        )
        new_cols = {}
        for col in self._schema.columns:
            col_storage = self._resolve_column_storage(col, default_chunks, default_blocks)
            eff_cparams = cparams_override if cparams_override is not None else col_storage.get("cparams")
            if self._is_list_column(col):
                new_cols[col.name] = mem_storage.create_list_column(
                    col.name,
                    spec=col.spec,
                    cparams=eff_cparams,
                    dparams=col_storage.get("dparams"),
                )
            elif self._is_varlen_scalar_column(col):
                new_cols[col.name] = mem_storage.create_varlen_scalar_column(
                    col.name,
                    spec=col.spec,
                    cparams=eff_cparams,
                    dparams=col_storage.get("dparams"),
                )
            elif self._is_dictionary_column(col):
                dict_col = mem_storage.create_dictionary_column(
                    col.name,
                    spec=col.spec,
                    cparams=eff_cparams,
                    dparams=col_storage.get("dparams"),
                )
                if len(dict_col.codes) < capacity:
                    dict_col.codes.resize((capacity,))
                new_cols[col.name] = dict_col
            else:
                shape = self._column_physical_shape(col, capacity)
                chunks = col_storage["chunks"]
                blocks = col_storage["blocks"]
                if col.config.chunks is None and col.config.blocks is None:
                    if col.name in aligned_names:
                        chunks, blocks = shared_chunks, shared_blocks
                    else:
                        chunks, blocks = self._column_chunks_blocks(col, shape)
                if chunks_override is not None:
                    chunks = chunks_override
                    if blocks_override is None:
                        eff_blocks = shared_blocks if shared_blocks is not None else default_blocks
                        if eff_blocks is not None and all(
                            b <= c for b, c in zip(eff_blocks, chunks_override, strict=False)
                        ):
                            blocks = eff_blocks
                        else:
                            blocks = None
                if blocks_override is not None:
                    blocks = blocks_override
                new_cols[col.name] = mem_storage.create_column(
                    col.name,
                    dtype=col.dtype,
                    shape=shape,
                    chunks=chunks,
                    blocks=blocks,
                    cparams=eff_cparams,
                    dparams=col_storage.get("dparams"),
                )

        obj = CTable.__new__(CTable)
        obj._schema = self._schema
        obj._row_type = self._row_type
        obj._table_cparams = self._table_cparams
        obj._table_dparams = self._table_dparams
        obj._storage = mem_storage
        obj._valid_rows = new_valid
        obj._cols = new_cols
        obj._col_widths = self._col_widths.copy()
        obj.col_names = [col.name for col in self._schema.columns]
        obj.auto_compact = self.auto_compact
        obj._create_summary_index = self._create_summary_index
        obj._summary_indexes_built = False  # compact creates a new copy; indexes not yet built
        obj._materialized_cols = {name: dict(meta) for name, meta in self._materialized_cols.items()}
        obj._expr_index_arrays = dict(self._expr_index_arrays)
        # Rebuild computed columns with the new NDArray objects as operands
        obj._computed_cols = {}
        for cc_name, cc in self._computed_cols.items():
            if cc.get("kind") == "dsl":
                # DSL entries hold the live kernel; the LazyUDF is rebuilt on
                # demand from obj._cols, so no operand rebinding is needed here.
                dsl_entry: dict[str, Any] = {
                    "kind": "dsl",
                    "dsl_source": cc["dsl_source"],
                    "kernel": cc["kernel"],
                    "col_deps": cc["col_deps"],
                    "dtype": cc["dtype"],
                    **({"jit_backend": cc["jit_backend"]} if "jit_backend" in cc else {}),
                }
                obj._computed_cols[cc_name] = dsl_entry
            else:
                operands = {f"o{i}": new_cols[dep] for i, dep in enumerate(cc["col_deps"])}
                new_lazy = blosc2.lazyexpr(cc["expression"], operands)
                obj._computed_cols[cc_name] = {
                    "kind": "expression",
                    "expression": cc["expression"],
                    "col_deps": cc["col_deps"],
                    "lazy": new_lazy,
                    "dtype": cc["dtype"],
                }
            obj.col_names.append(cc_name)
            obj._col_widths.setdefault(cc_name, max(len(cc_name), 15))
        obj._n_rows = 0
        obj._last_pos = None
        obj._read_only = False
        obj.base = None
        obj.auto_compact = self.auto_compact
        obj._validate = self._validate
        return obj

    # ------------------------------------------------------------------
    # Properties / info
    # ------------------------------------------------------------------

    @property
    def nrows(self) -> int:
        return self._n_rows

    @property
    def ncols(self) -> int:
        """Total number of columns, including computed (virtual) columns."""
        return len(self.col_names)

    @property
    def cbytes(self) -> int:
        """Total compressed size in bytes (all columns + valid_rows mask)."""
        return sum(col.cbytes for col in self._cols.values()) + self._valid_rows.cbytes

    @property
    def nbytes(self) -> int:
        """Total uncompressed size in bytes (all columns + valid_rows mask)."""
        return sum(col.nbytes for col in self._cols.values()) + self._valid_rows.nbytes

    @property
    def cratio(self) -> float:
        """Compression ratio for the whole table payload."""
        if self.cbytes == 0:
            return float("inf")
        return self.nbytes / self.cbytes

    @property
    def schema(self) -> CompiledSchema:
        """The compiled schema that drives this table's columns and validation."""
        return self._schema

    @property
    def vlmeta(self):
        """Variable-length metadata attached to this table.

        Returns a mapping-like proxy that supports item access, iteration,
        and the ``[:]`` bulk getter.  Values are serialised via msgpack, so
        all standard types (int, float, str, bool, list, dict) are supported.
        The metadata is stored separately from the internal schema metadata
        and persists through ``close()`` / reopen for disk-backed tables.

        Examples
        --------
        >>> import blosc2
        >>> import dataclasses
        >>> @dataclasses.dataclass
        ... class Row:
        ...     x: int = 0
        >>> t = blosc2.CTable(Row)
        >>> t.vlmeta["author"] = "Alice"
        >>> t.vlmeta["tags"] = ["alpha", "beta"]
        >>> t.vlmeta["count"] = 42
        >>> print(t.vlmeta["author"])
        Alice
        >>> print(t.vlmeta[:])
        {'author': 'Alice', 'tags': ['alpha', 'beta'], 'count': 42}
        >>> del t.vlmeta["count"]
        >>> for name in t.vlmeta:
        ...     print(name, t.vlmeta[name])
        ...
        author Alice
        tags ['alpha', 'beta']
        """
        storage = getattr(self, "_storage", None)
        if storage is None:
            raise AttributeError("CTable has no storage backend")
        if not hasattr(storage, "_open_meta"):
            # In-memory table: create a simple SChunk to hold vlmeta lazily
            _tmp = getattr(storage, "_vlmeta_schunk", None)
            if _tmp is None:
                storage._vlmeta_schunk = blosc2.SChunk()
            return storage._vlmeta_schunk.vlmeta
        # Persistent table: use the dedicated user-vlmeta SChunk
        meta = storage._open_vlmeta()
        if meta is None:
            # First access — create an in-memory SChunk; it will be saved
            # to disk when the table is closed.
            meta = blosc2.SChunk()
            storage._vlmeta = meta
        return meta.vlmeta

    def column_schema(self, name: str) -> CompiledColumn:
        """Return the :class:`CompiledColumn` descriptor for *name*.

        Raises
        ------
        KeyError
            If *name* is not a column in this table.
        """
        try:
            return self._schema.columns_by_name[name]
        except KeyError:
            raise KeyError(f"No column named {name!r}. Available: {self.col_names}") from None

    def schema_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible dict describing this table's schema."""
        return schema_to_dict(self._schema)

    # ------------------------------------------------------------------
    # Info reporting
    # ------------------------------------------------------------------

    @property
    def info_items(self) -> list[tuple[str, object]]:
        """Structured summary items used by :meth:`info`."""
        storage_type = "persistent" if isinstance(self._storage, FileTableStorage) else "in-memory"
        urlpath = self._storage._root if isinstance(self._storage, FileTableStorage) else None
        column_summary = {}
        for name in self.col_names:
            if name in self._computed_cols:
                cc = self._computed_cols[name]
                column_summary[name] = _InfoLiteral(
                    f"{cc['dtype']} (computed: {self._readable_computed_expr(cc)})"
                )
            else:
                col_meta = self._schema.columns_by_name.get(name)
                dtype_label = self._dtype_info_label(
                    getattr(self._cols[name], "dtype", None), col_meta.spec if col_meta else None
                )
                cbytes = getattr(self._cols[name], "cbytes", None)
                if cbytes is not None:
                    nbytes = getattr(self._cols[name], "nbytes", None)
                    detail = f"cbytes: {format_nbytes_human(cbytes)}"
                    if nbytes is not None and cbytes:
                        detail += f", cratio: {nbytes / cbytes:.2f}x"
                    column_summary[name] = _InfoLiteral(f"{dtype_label} ({detail})")
                else:
                    column_summary[name] = _InfoLiteral(dtype_label)

        index_summary = {}
        for idx in self.indexes:
            stale = " stale" if idx.stale else ""
            label = f" name={idx.name!r}" if idx.name and idx.name != "__self__" else ""
            stats = idx.storage_stats()
            if stats is None:
                suffix = "(size=n/a, sidecars not directly addressable)"
            else:
                _, cbytes, _ = stats
                suffix = f"({format_nbytes_human(cbytes)})"
            index_summary[idx.col_name] = f"[{idx.kind}{stale}{label}] {suffix}"

        items = [
            ("type", self.__class__.__name__),
            ("storage", storage_type),
            ("view", self.base is not None),
            ("nrows", self.nrows),
            ("ncols", self.ncols),
            ("chunks", self.chunks if self.chunks is not None else "none (no fixed-size columns)"),
            ("blocks", self.blocks if self.blocks is not None else "none (no fixed-size columns)"),
            ("nbytes", format_nbytes_info(self.nbytes)),
            ("cbytes", format_nbytes_info(self.cbytes)),
            ("cratio", f"{self.cratio:.2f}x"),
            ("columns", column_summary),
            ("indexes", index_summary if index_summary else "none"),
        ]
        # Only surface the validity-mask overhead when the table carries
        # uncompacted tombstones, i.e. the physical extent exceeds the number
        # of live rows (mid-table deletions not yet reclaimed by compact()).
        if self.base is None:
            live_rows = int(blosc2.count_nonzero(self._valid_rows))
            if self._resolve_last_pos() > live_rows:
                items.insert(
                    items.index(("columns", column_summary)), ("valid_rows", self._valid_rows_info_label())
                )
        if urlpath is not None:
            items.insert(2, ("urlpath", urlpath))
            open_mode = self._storage.open_mode()
            if open_mode is not None:
                items.insert(3, ("open_mode", open_mode))
        return items

    def _valid_rows_info_label(self) -> str:
        """Storage label for the validity mask (deletion/selection overhead)."""
        vr = self._valid_rows
        dtype_label = str(getattr(vr, "dtype", "bool"))
        cbytes = getattr(vr, "cbytes", None)
        if cbytes is None:
            return dtype_label
        detail = f"cbytes: {format_nbytes_human(cbytes)}"
        nbytes = getattr(vr, "nbytes", None)
        if nbytes is not None and cbytes:
            detail += f", cratio: {nbytes / cbytes:.2f}x"
        return f"{dtype_label} ({detail})"

    @staticmethod
    def _dtype_info_label(dtype: np.dtype | None, spec: SchemaSpec | None = None) -> str:
        """Return a compact dtype label for info reports."""
        if isinstance(spec, DictionarySpec):
            ordered_tag = ", ordered" if spec.ordered else ""
            return f"dictionary[str{ordered_tag}]"
        if isinstance(spec, Utf8Spec):
            return "utf8"
        if isinstance(spec, VLStringSpec):
            return "vlstring"
        if isinstance(spec, VLBytesSpec):
            return "vlbytes"
        if isinstance(spec, StructSpec):
            return spec.display_label()
        if isinstance(spec, ObjectSpec):
            return spec.display_label()
        if isinstance(spec, ListSpec):
            return spec.display_label()
        if isinstance(spec, NDArraySpec):
            return spec.display_label()
        if isinstance(spec, timestamp):
            return (
                f"timestamp[{spec.unit}]"
                if spec.timezone is None
                else f"timestamp[{spec.unit}, {spec.timezone}]"
            )
        if dtype is None:
            return "None"
        if dtype.kind == "U":
            nchars = dtype.itemsize // 4
            return f"U{nchars} (Unicode)"
        if dtype.kind == "S":
            return f"S{dtype.itemsize}"
        return str(dtype)

    @property
    def info(self) -> _CTableInfoReporter:
        """Get information about this table.

        Examples
        --------
        >>> print(t.info)
        >>> t.info()
        """
        return _CTableInfoReporter(self)

    # ------------------------------------------------------------------
    # Mutation: append / extend / delete
    # ------------------------------------------------------------------

    def _load_initial_data(self, new_data) -> None:
        """Dispatch new_data to append() or extend() as appropriate."""
        is_append = False

        if isinstance(new_data, (np.void, np.record)):
            is_append = True
        elif isinstance(new_data, np.ndarray):
            if new_data.dtype.names is not None and new_data.ndim == 0:
                is_append = True
        elif isinstance(new_data, list) and len(new_data) > 0:
            first_elem = new_data[0]
            if isinstance(first_elem, (str, bytes, int, float, bool, complex)):
                is_append = True

        if is_append:
            self.append(new_data)
        else:
            self.extend(new_data)

    def append(self, data: list | np.void | np.ndarray) -> None:
        """Append a single row to the table.

        *data* may be a list, tuple, ``numpy.void``, or structured
        ``numpy.ndarray`` whose fields match the schema column order.
        Materialized columns whose values are omitted are auto-filled from
        their recorded expression.  Raises ``ValueError`` if the table is
        read-only or a view.

        For tables with **nested (dotted) column names** the row dict may be
        supplied either as a flat mapping of dotted keys or as a nested dict
        that mirrors the original struct shape — both are accepted and
        automatically flattened to the physical dotted leaf names::

            # flat dotted keys
            t.append({"trip.begin.lon": -87.6, "trip.begin.lat": 41.8,
                      "payment.fare": 12.5})

            # original nested dict (auto-flattened)
            t.append({"trip": {"begin": {"lon": -87.6, "lat": 41.8}},
                      "payment": {"fare": 12.5}})
        """
        if self._read_only:
            raise ValueError("Table is read-only (opened with mode='r').")
        if self.base is not None:
            raise TypeError("Cannot extend view.")

        # Normalize → validate → coerce
        row = self._normalize_row_input(data)
        row = self._autofill_materialized_row_values(row)
        self._validate_no_default_columns_present(row)
        if self._validate:
            from blosc2.schema_validation import validate_row

            row = validate_row(self._schema, row)
        row = self._coerce_row_to_storage(row)

        pos = self._resolve_last_pos()
        if pos >= len(self._valid_rows):
            self._grow()

        for col in self._schema.columns:
            name = col.name
            col_array = self._cols[name]
            if self._is_list_column(col) or self._is_varlen_scalar_column(col):
                col_array.append(row[name])
            elif self._is_dictionary_column(col):
                col_array[pos] = row[name]  # DictionaryColumn encodes on __setitem__
            else:
                col_array[pos] = row[name]
                acc = self._get_summary_accumulator(name)
                if acc is not None and acc.valid:
                    acc.feed(pos, np.asarray([row[name]], dtype=col_array.dtype))

        n_rows = self.nrows
        self._valid_rows[pos] = True
        self._last_pos = pos + 1
        self._n_rows = n_rows + 1
        self._mark_all_indexes_stale()

    def delete(self, ind: int | slice | str | Iterable) -> None:
        """Mark one or more rows as deleted (tombstone deletion).

        *ind* may be a logical row index (``int``), a slice, or an iterable of
        logical indices.  Deleted rows are excluded from all subsequent queries
        and aggregates.  Physical storage is not reclaimed until
        :meth:`compact` is called.  Raises ``ValueError`` if the table is
        read-only or a view.
        """
        if self._read_only:
            raise ValueError("Table is read-only (opened with mode='r').")
        if self.base is not None:
            raise ValueError("Cannot delete rows from a view.")
        true_pos = self._live_positions_from_valid_rows_chunks()

        if isinstance(ind, Iterable) and not isinstance(ind, (str, bytes)):
            ind = list(ind)
        elif not isinstance(ind, int) and not isinstance(ind, slice):
            raise TypeError(f"Invalid type '{type(ind)}'")

        false_pos = true_pos[ind]
        n_deleted = len(np.unique(false_pos))
        n_rows = self.nrows

        self._valid_rows[false_pos] = False
        self._n_rows = n_rows - n_deleted
        if self._last_pos is None or np.any(false_pos == self._last_pos - 1):
            self._last_pos = None  # last live row deleted; recalculate on next write
        self._storage.bump_visibility_epoch()

    def extend(self, data: list | CTable | Any, *, validate: bool | None = None) -> None:  # noqa: C901
        """Append multiple rows at once.

        *data* may be:

        * a **dict of arrays** ``{"col": array, ...}`` — all arrays must have
          the same length; omitted columns are filled from their declared default;
          columns with no default declared must be provided;
        * a **list of rows**, each compatible with :meth:`append`;
        * another **CTable** — columns are matched by name.

        Pass ``validate=False`` to skip per-row Pydantic validation on trusted
        bulk imports.  Raises ``ValueError`` if the table is read-only or a view.

        For tables with **nested (dotted) column names** both the dict-of-arrays
        and list-of-dicts forms accept the original nested dict shape and
        auto-flatten it to physical dotted leaf names::

            # nested dict of arrays
            t.extend({
                "trip": {"begin": {"lon": lons, "lat": lats}},
                "payment": {"fare": fares},
            })

            # list of nested dicts
            t.extend([
                {"trip": {"begin": {"lon": -87.6, "lat": 41.8}}, "payment": {"fare": 12.5}},
                {"trip": {"begin": {"lon": -87.5, "lat": 41.7}}, "payment": {"fare": 8.0}},
            ])
        """
        if self._read_only:
            raise ValueError("Table is read-only (opened with mode='r').")
        if self.base is not None:
            raise TypeError("Cannot extend view.")
        if len(data) <= 0:
            if isinstance(data, dict):
                raise ValueError("No columns provided for extend().")
            return

        # Resolve effective validate flag: per-call override takes precedence
        do_validate = self._validate if validate is None else validate

        start_pos = self._resolve_last_pos()

        current_col_names = self._stored_col_names  # skip computed columns
        input_col_names = self._append_input_col_names
        new_nrows = 0
        provided_names: set[str] = set()

        if hasattr(data, "_cols") and hasattr(data, "_n_rows"):
            new_nrows = data._n_rows
            raw_columns = {}
            for name in current_col_names:
                if name in data._cols:
                    raw_columns[name] = data._cols[name][: data._n_rows]
                    provided_names.add(name)
        else:
            if isinstance(data, dict):
                if any(isinstance(v, dict) for v in data.values()):
                    data = self._flatten_nested_dict(data)
                known_names = [name for name in current_col_names if name in data]
                if not known_names:
                    raise ValueError("No known stored columns provided for extend().")
                column_lengths = {}
                for name in known_names:
                    try:
                        column_lengths[name] = len(data[name])
                    except TypeError as exc:
                        raise TypeError(f"Column {name!r} does not have a length.") from exc
                new_nrows = column_lengths[known_names[0]]
                mismatched = {name: n for name, n in column_lengths.items() if n != new_nrows}
                if mismatched:
                    details = ", ".join(f"{name}={n}" for name, n in mismatched.items())
                    raise ValueError(
                        f"All provided columns must have the same length; "
                        f"expected {new_nrows}, got {details}."
                    )
                provided_names = set(known_names)
                raw_columns = {name: data[name] for name in known_names}
            elif isinstance(data, np.ndarray) and data.dtype.names is not None:
                new_nrows = len(data)
                raw_columns = {name: data[name] for name in data.dtype.names if name in current_col_names}
                provided_names = set(raw_columns)
            elif data and isinstance(data[0], dict):
                # List of dicts: flatten any nested dicts and pivot to column arrays.
                flat_rows = [
                    self._flatten_nested_dict(row) if any(isinstance(v, dict) for v in row.values()) else row
                    for row in data
                ]
                new_nrows = len(flat_rows)
                col_set = set(input_col_names)
                raw_columns = {
                    name: [row[name] for row in flat_rows]
                    for name in input_col_names
                    if name in flat_rows[0]
                }
                provided_names = set(raw_columns)
                # Fill any remaining columns from the rows (may include extra keys)
                for row in flat_rows:
                    for key in row:
                        if key in col_set and key not in raw_columns:
                            raw_columns[key] = [r.get(key) for r in flat_rows]
                            provided_names.add(key)
            else:
                new_nrows = len(data)
                batch_columns = list(zip(*data, strict=False))
                raw_columns = {
                    input_col_names[i]: batch_columns[i]
                    for i in range(min(len(input_col_names), len(batch_columns)))
                }
                provided_names = set(raw_columns)

        raw_columns = self._autofill_materialized_batch_columns(
            raw_columns, new_nrows, provided_names=provided_names
        )
        raw_columns = self._fill_default_batch_columns(raw_columns, new_nrows)

        # Validate constraints column-by-column before writing
        if do_validate:
            from blosc2.schema_vectorized import validate_column_batch

            validate_column_batch(self._schema, raw_columns)

        scalar_processed_cols: dict[str, blosc2.NDArray] = {}
        list_processed_cols: dict[str, list] = {}
        varlen_scalar_processed_cols: dict[str, list] = {}
        dict_processed_cols: dict[str, list] = {}
        for name in current_col_names:
            col_meta = self._schema.columns_by_name[name]
            if self._is_list_column(col_meta):
                list_processed_cols[name] = list(raw_columns[name])
            elif self._is_varlen_scalar_column(col_meta):
                varlen_scalar_processed_cols[name] = list(raw_columns[name])
            elif self._is_dictionary_column(col_meta):
                dict_processed_cols[name] = list(raw_columns[name])
            else:
                target_dtype = self._cols[name].dtype
                if isinstance(col_meta.spec, timestamp):
                    values = np.asarray(raw_columns[name])
                    if np.issubdtype(values.dtype, np.datetime64):
                        values = values.astype(f"datetime64[{col_meta.spec.unit}]").astype(np.int64)
                    elif values.dtype.kind in "OUS":
                        values = np.array(
                            [
                                col_meta.spec.null_value
                                if v is None
                                else np.datetime64(v)
                                .astype(f"datetime64[{col_meta.spec.unit}]")
                                .astype(np.int64)
                                if isinstance(v, (np.datetime64, str)) or hasattr(v, "isoformat")
                                else v
                                for v in values
                            ],
                            dtype=target_dtype,
                        )
                    scalar_processed_cols[name] = np.ascontiguousarray(values, dtype=target_dtype)
                elif self._is_ndarray_column(col_meta):
                    scalar_processed_cols[name] = self._coerce_ndarray_batch(
                        name, col_meta.spec, raw_columns[name], new_nrows
                    )
                else:
                    raw = raw_columns[name]
                    if isinstance(raw, blosc2.NDArray):
                        # Keep as-is; written chunk-by-chunk in the write loop below.
                        # validate_column_batch() above also scans NDArray columns
                        # chunk-by-chunk, so validation never fully decompresses them.
                        scalar_processed_cols[name] = raw
                    else:
                        scalar_processed_cols[name] = np.ascontiguousarray(raw, dtype=target_dtype)

        end_pos = start_pos + new_nrows

        if self.auto_compact and end_pos >= len(self._valid_rows):
            self.compact()  # sets _last_pos = _n_rows
            start_pos = self._last_pos
            end_pos = start_pos + new_nrows

        while end_pos > len(self._valid_rows):
            self._grow()

        for name in current_col_names:
            col_meta = self._schema.columns_by_name[name]
            if self._is_list_column(col_meta):
                self._cols[name].extend(list_processed_cols[name], validate=do_validate)
            elif self._is_varlen_scalar_column(col_meta):
                self._cols[name].extend(varlen_scalar_processed_cols[name])
            elif self._is_dictionary_column(col_meta):
                # DictionaryColumn.__setitem__ with a slice encodes all values.
                self._cols[name][start_pos:end_pos] = dict_processed_cols[name]
            else:
                values = scalar_processed_cols[name]
                if isinstance(values, blosc2.NDArray):
                    # Decompress one chunk at a time to bound peak memory usage.
                    tgt = self._cols[name]
                    chunk_size = values.chunks[0] if values.chunks else 65536
                    for c in range(0, new_nrows, chunk_size):
                        c_end = min(c + chunk_size, new_nrows)
                        chunk = np.ascontiguousarray(values[c:c_end], dtype=tgt.dtype)
                        tgt[start_pos + c : start_pos + c_end] = chunk
                        self._feed_summary(name, start_pos + c, chunk)
                else:
                    self._cols[name][start_pos:end_pos] = values[:]
                    self._feed_summary(name, start_pos, values)

        n_rows = self.nrows
        self._valid_rows[start_pos:end_pos] = True
        self._last_pos = end_pos
        self._n_rows = n_rows + new_nrows
        self._mark_all_indexes_stale()

    # ------------------------------------------------------------------
    # Filtering
    # ------------------------------------------------------------------

    def _where_expression_operands(
        self, expr: str | None = None
    ) -> dict[str, blosc2.NDArray | blosc2.LazyExpr]:
        operands = {}
        for name, arr in self._cols.items():
            col = self._schema.columns_by_name.get(name)
            if col is not None and not (
                self._is_list_column(col)
                or self._is_varlen_scalar_column(col)
                or self._is_dictionary_column(col)
                or self._is_ndarray_column(col)
            ):
                operands[name] = arr
        for name, cc in self._computed_cols.items():
            if expr is None or self._expression_references_name(expr, name):
                operands[name] = self._build_computed_lazy(cc)
        return operands

    # Quoted string literal (may contain commas/spaces/escapes), single or double.
    _STR_LITERAL = r"""(\"(?:[^"\\]|\\.)*\"|'(?:[^'\\]|\\.)*')"""

    def _rewrite_dictionary_predicates(
        self, expr: str, operands: dict[str, blosc2.NDArray | blosc2.LazyExpr]
    ) -> tuple[str, dict[str, blosc2.NDArray | blosc2.LazyExpr]]:
        """Rewrite dictionary-column string predicates into integer-code comparisons.

        Dictionary columns are excluded from the plain operand namespace because
        the expression engine would compare raw int32 codes against the string
        literal (never matching).  Two predicate forms are resolved here, against
        the dictionary's codes:

        - ``dictcol == "literal"`` / ``!=`` -> ``dictcol ==/!= <code>``; an absent
          literal maps to a sentinel code no row carries (``==`` matches nothing,
          ``!=`` everything).
        - ``"literal" in dictcol`` -> substring search: an ``OR`` of ``==`` over
          every dictionary value containing *literal* (none -> matches nothing).

        The column's codes array is then supplied as the operand, so the result is
        an ordinary numeric expression that combines with the rest (``and``/``or``,
        precedence) natively.  Other uses of a dictionary column are left untouched
        and still raise ``Unknown symbol``.
        """
        absent_code = int(np.iinfo(np.int32).min)  # a code no live row carries
        rewritten = expr
        new_operands = dict(operands)
        for col in self._schema.columns:
            if not self._is_dictionary_column(col) or not self._expression_references_name(expr, col.name):
                continue
            dc = self._cols[col.name]  # DictionaryColumn
            name = col.name

            def eq_repl(match: re.Match, _dc=dc, _name=name) -> str:
                value = ast.literal_eval(match.group(2))
                try:
                    code = int(_dc.value_to_code(value))
                except KeyError:
                    code = absent_code
                return f"{_name} {match.group(1)} {code}"

            def in_repl(match: re.Match, _dc=dc, _name=name) -> str:
                needle = ast.literal_eval(match.group(1))
                codes = [c for c, value in enumerate(_dc.dictionary) if needle in value]
                if not codes:
                    return f"({_name} == {absent_code})"
                # Per-term parens: bitwise ``|`` binds tighter than ``==``.
                return "(" + " | ".join(f"({_name} == {c})" for c in codes) + ")"

            eq_pattern = r"(?<![\w.])" + re.escape(name) + r"\s*(==|!=)\s*" + self._STR_LITERAL
            in_pattern = self._STR_LITERAL + r"\s+in\s+(?<![\w.])" + re.escape(name) + r"(?![\w.])"
            new_expr = re.sub(eq_pattern, eq_repl, rewritten)
            new_expr = re.sub(in_pattern, in_repl, new_expr)
            if new_expr != rewritten:
                # numexpr needs all operands chunked alike; the codes array uses a
                # different chunkshape than the regular columns, so when mixed with
                # one (e.g. ``dictcol == "x" and other > y``) re-chunk it to match.
                codes = dc.codes
                target = next(
                    (
                        v
                        for v in operands.values()
                        if isinstance(v, blosc2.NDArray)
                        and v.shape == codes.shape
                        and v.chunks != codes.chunks
                    ),
                    None,
                )
                if target is not None:
                    codes = blosc2.asarray(codes[:], chunks=target.chunks, blocks=target.blocks)
                new_operands[name] = codes
                rewritten = new_expr
        return rewritten, new_operands

    def _rewrite_nested_expression(
        self, expr: str, operands: dict[str, blosc2.NDArray | blosc2.LazyExpr]
    ) -> tuple[str, dict[str, blosc2.NDArray | blosc2.LazyExpr]]:
        """Rewrite dotted nested names in *expr* to safe identifiers.

        `blosc2.lazyexpr` does not accept dotted identifiers, but nested leaf
        columns are naturally addressed as dotted paths (e.g. ``trip.begin.lon``).
        This maps them to temporary aliases and returns rewritten expression and
        operand mapping.
        """
        dotted = [name for name in operands if "." in name]
        if not dotted:
            return expr, operands

        rewritten = expr
        new_operands = dict(operands)
        # Longest names first so trip.begin.lon is rewritten before trip.begin.
        for i, name in enumerate(sorted(dotted, key=len, reverse=True)):
            alias = f"__nf{i}"
            pattern = rf"(?<![\w.]){re.escape(name)}(?![\w.])"
            replaced = re.sub(pattern, alias, rewritten)
            if replaced != rewritten:
                rewritten = replaced
                new_operands[alias] = new_operands.pop(name)
        return rewritten, new_operands

    @staticmethod
    def _expression_references_name(expr: str, name: str) -> bool:
        return re.search(rf"(?<![\w.]){re.escape(name)}(?![\w.])", expr) is not None

    def _guard_scalar_expression(self, expr: str) -> None:
        for name, meta in self._root_table._materialized_cols.items():
            if meta.get("stale", False) and self._expression_references_name(expr, name):
                raise ValueError(
                    f"Generated column {name!r} is stale because one or more source columns were modified. "
                    f"Call refresh_generated_column({name!r}) before using it in expressions, or use "
                    f"t[{name!r}].read_stale() to explicitly read the last stored stale values."
                )
        for col in self._schema.columns:
            if self._is_ndarray_column(col) and self._expression_references_name(expr, col.name):
                raise TypeError(
                    f"Column {col.name!r} is a fixed-shape ndarray column. String expressions only "
                    "support scalar columns. Use an element projection or a row-wise reduction first."
                )
            if self._is_utf8_column(col) and self._expression_references_name(expr, col.name):
                raise NotImplementedError(
                    f"Column {col.name!r} is a variable-length utf8 column; "
                    "string expressions on utf8 columns are not supported yet."
                )
            if self._is_varlen_scalar_column(col) and self._expression_references_name(expr, col.name):
                raise NotImplementedError(
                    f"Column {col.name!r} is a variable-length scalar column (vlstring/vlbytes/struct/object); "
                    "lazy expressions are not supported yet."
                )

    def _guard_varlen_scalar_expression(self, expr: str) -> None:
        self._guard_scalar_expression(expr)

    def _is_nullable_column(self, name: str) -> bool:
        col = self[name]
        return col.null_value is not None or col.is_dictionary or col.is_varlen_scalar

    def dropna(self, subset: list[str] | None = None) -> CTable:
        """Return a view excluding rows where any column in *subset* is null.

        Parameters
        ----------
        subset:
            Column names to check for nulls. Defaults to every nullable
            column (sentinel-backed, dictionary, or variable-length scalar).

        Returns
        -------
        CTable
            A read-only view over the live, non-null rows (see :meth:`where`).
        """
        names = subset if subset is not None else [n for n in self.col_names if self._is_nullable_column(n)]
        mask = np.ones(self.nrows, dtype=np.bool_)
        for name in names:
            mask &= self[name].notnull()
        return self.where(mask)

    def where(  # noqa: C901
        self,
        expr_result: str | np.ndarray | blosc2.NDArray | blosc2.LazyExpr | blosc2.LazyUDF | Column | ColExpr,
        *,
        columns: list[str] | tuple[str, ...] | None = None,
    ) -> CTable:
        """Return a row-filtered view matching a boolean predicate.

        Signature::

            where(expr_result) -> CTable

        The predicate can be supplied as a boolean :class:`blosc2.LazyExpr`,
        a boolean :class:`blosc2.NDArray`, a boolean NumPy array, a boolean
        ``Column``, an unbound :func:`col` expression, a :class:`blosc2.LazyUDF`
        (including those backed by a :func:`blosc2.dsl_kernel`), or a string
        expression evaluated against this table's columns.  String expressions
        can reference stored and computed columns directly by name.

        The returned object is a :class:`CTable` view sharing the original
        column data.  The row-selection mask is evaluated immediately and
        intersected with the table's current live rows; selected column data is
        not copied.

        Parameters
        ----------
        expr_result:
            Boolean predicate selecting rows.  Strings are converted to a
            lazy expression with table columns as operands, e.g.
            ``"value * category >= 150"``.  Column objects can also be used in
            Python expressions, e.g. ``(t.value * t.category) >= 150``.

        Returns
        -------
        CTable
            A view over the same columns containing only rows where the
            predicate is true and the source row is live.  When ``columns`` is
            provided, the returned view is additionally projected to that
            ordered subset of columns.

        Raises
        ------
        TypeError
            If *expr_result* does not evaluate to a boolean Blosc2/NumPy
            array or lazy expression.

        Examples
        --------
        Filter using a string expression::

            view = t.where("value * category >= 150")
            slim = t.where("value * category >= 150", columns=["value", "category"])

        Filter using column arithmetic::

            view = t.where((t.value * t.category) >= 150)

        Blosc2 lazy functions can be used in column expressions::

            view = t.where(((t.value + 2) * blosc2.sin(t.category)) >= 10)

        For column names that are not valid Python identifiers, use item
        access::

            view = t.where((t["unit price"] * t["quantity"]) > 100)

        For tables with **nested (dotted) column names**, dotted leaf names and
        attribute-chain proxies work in both string and expression forms::

            view = t.where("trip.begin.lon > -87.7 and payment.fare > 10")
            view = t.where(t.trip.begin.lon > -87.7)

        Notes
        -----
        Use bitwise operators (``&``, ``|``, ``~``) or string expressions for
        element-wise boolean logic.  Python's logical operators ``and``, ``or``
        and ``not`` cannot be overloaded and therefore do not build lazy column
        expressions.

        Use::

            t.where((t.x > 0) & (t.y < 10))
            t.where(~t.returned)
            t.where("not returned")

        not::

            t.where((t.x > 0) and (t.y < 10))
            t.where(not t.returned)
        """
        if isinstance(expr_result, ColExpr):
            expr_result = expr_result._bind(self)
        if isinstance(expr_result, str):
            self._guard_varlen_scalar_expression(expr_result)
            operands = self._where_expression_operands(expr_result)
            expr_result, operands = self._rewrite_dictionary_predicates(expr_result, operands)
            expr_result, operands = self._rewrite_nested_expression(expr_result, operands)
            expr_result = blosc2.lazyexpr(expr_result, operands)
        if isinstance(expr_result, np.ndarray) and expr_result.dtype == np.bool_:
            expr_result = blosc2.asarray(expr_result)
        if isinstance(expr_result, Column):
            expr_result = (
                expr_result._raw_col == 1 if expr_result._is_nullable_bool else expr_result._raw_col
            )
        if isinstance(expr_result, blosc2.LazyUDF):
            # DSL miniexpr only supports full-array getitem, so we cannot stream
            # a LazyUDF chunk-by-chunk the way LazyExpr does.  Materialise the
            # full boolean array upfront and let the NDArray path handle it.
            expr_result = expr_result.compute()

        if not (
            isinstance(expr_result, (blosc2.NDArray, blosc2.LazyExpr))
            and (getattr(expr_result, "dtype", None) == np.bool_)
        ):
            raise TypeError(f"Expected boolean blosc2.NDArray or LazyExpr, got {type(expr_result).__name__}")

        # Attempt index-accelerated filtering before falling back to a full scan.
        if isinstance(expr_result, blosc2.LazyExpr):
            index_result = self._try_index_where(expr_result)
            if index_result is not None:
                if isinstance(index_result, blosc2.NDArray):
                    # Mask-producing index (SUMMARY/BUCKET): keep the compressed
                    # boolean mask and route it through the same downstream path a
                    # plain scan uses (no positions <-> mask round-trip).  The
                    # view extracts live positions lazily if sort_by/gather needs
                    # them.  Position-producing indexes (FULL/PARTIAL/OPSI) still
                    # return positions and go through _view_from_positions.
                    expr_result = index_result
                else:
                    result = self._view_from_positions(index_result)
                    return result if columns is None else result.select(list(columns))

        target_len = len(self._valid_rows)
        known_n_rows = self._known_n_rows()
        all_rows_valid = known_n_rows == target_len
        filter_intersected = False

        # Prefer a compressed boolean mask for LazyExpr filters so temporary
        # mask materialization stays compact even for medium-sized selections.
        if isinstance(expr_result, blosc2.LazyExpr):
            filter = expr_result.compute()
        else:
            filter = expr_result

        if getattr(filter, "ndim", 1) != 1:
            raise ValueError(
                "CTable.where() requires a 1-D row mask. Reduce ndarray-column predicates to one "
                "boolean per row before filtering."
            )

        filter_len = len(filter)
        if filter_len != target_len:
            if filter_len == self.nrows:
                physical = blosc2.zeros(target_len, dtype=np.bool_)
                physical[self._valid_rows] = filter[:]
                filter = physical
                filter_intersected = True
            elif filter_len > target_len:
                filter = filter[:target_len]
                filter_intersected = False
            else:
                padding = blosc2.zeros(target_len, dtype=np.bool_)
                padding[:filter_len] = filter[:]
                filter = padding
                filter_intersected = False

        if not filter_intersected and not all_rows_valid:
            if isinstance(filter, np.ndarray):
                filter &= self._valid_rows[:]
            else:
                filter = (filter & self._valid_rows).compute()

        result = self.view(filter)
        return result if columns is None else result.select(list(columns))

    def _run_row_logic(self, ind: int | slice | str | Iterable) -> CTable:
        true_pos = self._live_positions_from_valid_rows_chunks()

        # A boolean mask is set-like: it selects rows in physical-ascending order
        # and cannot represent duplicates.  Every other selector — a slice (incl.
        # negative-step reverse), an integer array, or a list of positions —
        # carries an explicit order and may repeat rows.  Detect the mask before
        # any list() coercion so we keep its dtype.
        is_bool_mask = isinstance(ind, np.ndarray) and ind.dtype == np.bool_

        if isinstance(ind, Iterable) and not isinstance(ind, (str, bytes)):
            ind = list(ind)

        mant_pos = np.asarray(true_pos[ind])

        # Carry the positions forward whenever order matters: for an ordered view
        # (sorted view or position view), or for any order-carrying selector.  A
        # boolean mask is physical-order and set-like, so going through it would
        # silently drop both the requested order and any duplicates.
        if getattr(self, "_cached_live_positions", None) is not None or not is_bool_mask:
            return self._view_from_positions(mant_pos)

        new_mask_np = np.zeros(len(self._valid_rows), dtype=bool)
        new_mask_np[mant_pos] = True

        new_mask = blosc2.asarray(new_mask_np)
        return self.view(new_mask)


def ctable_from_cframe(cframe: bytes, *, copy: bool = True) -> CTable:
    """Deserialize a CFrame into a :class:`CTable`.

    The counterpart of :meth:`CTable.to_cframe`.  The cframe is decoded into an
    in-memory :class:`blosc2.EmbedStore` and opened through
    :class:`~blosc2.ctable_storage.EmbedStoreTableStorage`, so the result is a
    standalone in-memory table with no file dependency.

    Parameters
    ----------
    cframe : bytes
        The serialized table, as produced by :meth:`CTable.to_cframe`.
    copy : bool, optional
        If ``True``, copy the underlying buffers so the result does not share
        memory with *cframe*.  Default is ``False``.

    Returns
    -------
    CTable
        The deserialized table.

    See Also
    --------
    :meth:`blosc2.CTable.to_cframe`
    """
    from blosc2.ctable_storage import EmbedStoreTableStorage

    # Probe the cframe type with a cheap non-copying open; bail early on
    # non-EmbedStore / non-CTable frames so callers can try-fallback.
    probe = blosc2.schunk_from_cframe(cframe, copy=False)
    if "b2embed" not in probe.meta:
        raise ValueError("Not an EmbedStore cframe (no b2embed marker)")
    estore = blosc2.from_cframe(cframe, copy=copy)
    storage = EmbedStoreTableStorage(estore)
    storage.check_kind()  # raise if not a CTable
    return CTable._open_from_storage(storage)
