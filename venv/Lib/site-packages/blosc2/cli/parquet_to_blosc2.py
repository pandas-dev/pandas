#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################

"""Import/export Parquet datasets through a Blosc2 CTable store.

The installed ``parquet-to-blosc2`` utility supports three modes:

* default import: parquet -> ``.b2z`` / ``.b2d``
* ``--export``: existing ``.b2z`` / ``.b2d`` -> parquet
* ``--roundtrip``: parquet -> ``.b2z`` / ``.b2d`` -> parquet and compare

The output extension selects the storage layout: ``.b2z`` is compact/zip-backed,
while ``.b2d`` is sparse directory-backed.

Scalar string columns are stored as ``utf8`` (variable-length, no length
limit; nullable columns get a sentinel string). Scalar binary columns are
stored as ``vlbytes``, whose nullable columns represent nulls with native
``None`` — no sentinel value is needed.

Struct-valued columns are wrapped as ``list<struct>`` (one-element lists) so
that they round-trip through the list-column machinery. True list columns pass
through unchanged. Timestamp columns are imported as semantic CTable
``timestamp`` columns. Unsupported types (nested lists, durations, etc.) are
skipped.
"""

from __future__ import annotations

import argparse
import base64
import contextlib
import cProfile
import gc
import io
import os
import pstats
import shutil
import sys
import time
from pathlib import Path
from typing import Any

import blosc2
from blosc2.schema_compiler import _validate_column_name, schema_to_dict

DEFAULT_BATCH_SIZE = 2048
MAX_ELEMENT_WRITE_BATCH = 5_000_000  # cap on flattened elements yielded per write
UNNAMED_ROOT_CAPACITY_SAFETY = 1.15  # first-batch estimates are often a little low
# Target in-memory size of one Arrow read batch for the unnamed-root flatten
# path.  Nested list<struct> batches amplify ~10x downstream (flatten + cast +
# write buffers + Arrow pool), so an auto parquet batch size is capped to keep
# this Arrow batch small enough that peak RSS stays well under ~1 GB.
PARQUET_BATCH_ARROW_BUDGET = 48 * 2**20  # 48 MiB


def require_pyarrow():
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as exc:
        raise ImportError(
            "parquet-to-blosc2 requires pyarrow; install it with: pip install 'blosc2[parquet]'"
        ) from exc
    return pa, pq


def _default_import_output(input_path: Path) -> Path:
    return input_path.with_suffix(".b2z")


def _default_export_output(input_path: Path) -> Path:
    return input_path.with_suffix(".parquet")


def _default_roundtrip_output(input_path: Path) -> Path:
    return input_path.with_name(f"{input_path.stem}-roundtrip.parquet")


def _format_bytes(n: int | None) -> str:
    if n is None:
        return "n/a"
    value = float(n)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if abs(value) < 1024 or unit == "TiB":
            return f"{value:.1f} {unit}"
        value /= 1024
    return f"{value:.1f} TiB"


def _peak_rss_bytes() -> int:
    import resource

    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return int(peak)
    return int(peak) * 1024


def _current_rss_bytes() -> int | None:
    try:
        import psutil
    except ImportError:
        return None
    return int(psutil.Process(os.getpid()).memory_info().rss)


def memory_report(label: str, pa=None) -> None:
    arrow_allocated = None
    arrow_pool = None
    if pa is not None:
        try:
            arrow_allocated = int(pa.total_allocated_bytes())
            arrow_pool = int(pa.default_memory_pool().bytes_allocated())
        except Exception:
            pass
    parts = [
        f"[mem] {label}",
        f"rss={_format_bytes(_current_rss_bytes())}",
        f"peak={_format_bytes(_peak_rss_bytes())}",
    ]
    if arrow_allocated is not None:
        parts.append(f"arrow_total={_format_bytes(arrow_allocated)}")
    if arrow_pool is not None:
        parts.append(f"arrow_pool={_format_bytes(arrow_pool)}")
    print("  ".join(parts), flush=True)


def maybe_memory_report(args, label: str, pa=None) -> None:
    if args.mem_report:
        memory_report(label, pa)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Import/export Parquet datasets via Blosc2 CTable (.b2z compact or .b2d sparse).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--export", action="store_true", help="Export input .b2z/.b2d to output parquet.")
    mode.add_argument(
        "--roundtrip", action="store_true", help="Run parquet -> .b2z/.b2d -> parquet and compare."
    )
    parser.add_argument(
        "input_path", type=Path, help="Input parquet file or Blosc2 store, depending on mode."
    )
    parser.add_argument(
        "output_path",
        nargs="?",
        type=Path,
        default=None,
        help="Output path. Defaults depend on the mode and input path.",
    )
    parser.add_argument(
        "--parquet-batch-size",
        type=int,
        default=None,
        help="Rows per Parquet read batch. Defaults to the source Parquet average row-group size.",
    )
    parser.add_argument(
        "--fixed-str-maxlen",
        type=int,
        default=None,
        help=(
            "Pre-scan string columns and import columns whose maximum character length is at most "
            "this value as fixed-width, indexable strings. Other string columns remain utf8."
        ),
    )
    parser.add_argument(
        "--fixed-bytes-maxlen",
        type=int,
        default=None,
        help=(
            "Pre-scan binary columns and import columns whose maximum byte length is at most this value "
            "as fixed-width, indexable bytes. Other binary columns remain vlbytes."
        ),
    )
    parser.add_argument(
        "--max-rows",
        type=int,
        default=None,
        help=(
            "Maximum number of CTable rows to import.  "
            "In normal mode this equals the number of Parquet rows read.  "
            "With separate nested columns enabled for an unnamed-root list<struct<...>> "
            "file, the unit is list elements "
            "(i.e. the number of rows in the resulting CTable), "
            "not outer Parquet rows."
        ),
    )
    parser.add_argument(
        "--batch-size",
        dest="parquet_batch_size",
        type=int,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--blosc2-batch-size",
        type=int,
        default=None,
        help="Internal batch_rows for BatchArray/varlen columns in the imported CTable. "
        "Defaults to the blocks value from blosc2.compute_chunks_blocks() based on "
        "the estimated CTable row count.",
    )
    parser.add_argument(
        "--blosc2-items-per-block",
        type=int,
        default=None,
        help=(
            "Items per internal BatchArray block for imported Blosc2 varlen/list columns. "
            "Defaults to BatchArray's automatic heuristic."
        ),
    )
    parser.add_argument(
        "--list-serializer",
        choices=["msgpack", "arrow"],
        default="arrow",
        help=(
            "Serializer for imported list columns. 'arrow' is the default and stores Arrow list "
            "batches directly, which is much faster for deeply nested lists but requires PyArrow "
            "when reading those columns later. Use 'msgpack' to avoid that read-time dependency."
        ),
    )
    parser.add_argument("--use-dict", action="store_true", help="Enable C-Blosc2 dictionary compression.")
    parser.add_argument(
        "--float-trunc-prec",
        action="append",
        default=[],
        metavar="BITS|COLUMN=BITS",
        help=(
            "Apply the Blosc2 TRUNC_PREC filter to imported float32/float64 columns. "
            "Pass an integer to affect all float columns, or COLUMN=integer for a single column. "
            "May be repeated; column-specific entries override the global value."
        ),
    )
    parser.add_argument(
        "--timestamp-unit",
        choices=["s", "ms", "us", "ns", "auto"],
        default=None,
        help=(
            "Import timestamp columns using this unit. Explicit units use Arrow's safe cast and fail "
            "if conversion would lose precision. 'auto' pre-scans timestamp columns and chooses the "
            "coarsest lossless unit per column."
        ),
    )
    parser.add_argument(
        "--chunks",
        type=int,
        default=None,
        help=(
            "Chunk size (in rows) for all scalar columns in the imported CTable. "
            "Overrides the automatic chunk size chosen by blosc2.compute_chunks_blocks(). "
            "Only affects fixed-width scalar columns; list, varlen, and dictionary columns "
            "use their own internal chunking."
        ),
    )
    parser.add_argument(
        "--blocks",
        type=int,
        default=None,
        help=(
            "Block size (in rows) for all scalar columns in the imported CTable. "
            "Overrides the automatic block size chosen by blosc2.compute_chunks_blocks(). "
            "Must be <= chunks; if omitted when --chunks is given, blosc2 picks a suitable block size."
        ),
    )
    parser.add_argument("--codec", type=str, default="ZSTD", choices=[c.name for c in blosc2.Codec])
    parser.add_argument("--clevel", type=int, default=5)
    parser.add_argument(
        "--reduce-mem",
        action="store_true",
        help=(
            "Shrink an auto-chosen Parquet batch size so a single Arrow read batch fits a "
            "small memory budget, lowering peak RSS at the cost of import speed. "
            "Only affects auto batch sizing for unnamed-root list<struct<...>> flattening; "
            "an explicit --batch-size is always left untouched."
        ),
    )
    parser.add_argument(
        "--mem-report",
        action="store_true",
        help="Print process/Arrow memory diagnostics at import phases and during batch processing.",
    )
    parser.add_argument(
        "--mem-every",
        type=int,
        default=1,
        help="With --mem-report, print batch memory diagnostics every N batches.",
    )
    parser.add_argument(
        "--batch-report-every",
        type=int,
        default=1,
        help="With --progress, print progress every N batches; the final batch is always reported.",
    )
    parser.add_argument(
        "--progress",
        action="store_true",
        help="Print import progress lines. By default, only the import summary is shown.",
    )
    parser.add_argument(
        "--profile",
        action="store_true",
        help="Run the selected operation under cProfile and print cumulative timing stats.",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--decode-dictionaries",
        action="store_true",
        help=(
            "Decode Arrow dictionary-encoded columns to plain utf8 instead of preserving "
            "the dictionary encoding.  By default, supported dictionary columns "
            "(string values with integer indices) are imported as Blosc2 dictionary columns."
        ),
    )
    parser.add_argument(
        "--separate-nested-cols",
        action=argparse.BooleanOptionalAction,
        default=True,
        dest="separate_nested_cols",
        help=(
            "Import nested columns as separate CTable columns where possible. "
            "Top-level struct fields are flattened recursively into dotted leaf columns "
            "(e.g. trip.begin.lon).  For a single unnamed top-level list<struct<…>> "
            "field (the Awkward Array / Chicago-taxi layout), flatten the outer list "
            "so that each element becomes a CTable row. Enabled by default; use "
            "--no-separate-nested-cols when closer Parquet schema fidelity is desired."
        ),
    )
    parser.add_argument(
        "--no-summary-index",
        action="store_false",
        dest="create_summary_index",
        default=True,
        help=(
            "Disable automatic SUMMARY index creation on close. "
            "By default, SUMMARY indexes are built for all eligible scalar columns, "
            "which costs <0.1%% of column size and accelerates WHERE queries."
        ),
    )
    return parser


def prepare_output(path: Path, overwrite: bool) -> None:
    if not path.exists():
        return
    if not overwrite:
        raise FileExistsError(f"Output already exists: {path} (use --overwrite to replace)")
    if path.is_dir():
        shutil.rmtree(path)
    else:
        path.unlink()


def encode_arrow_schema(schema) -> str:
    return base64.b64encode(schema.serialize().to_pybytes()).decode("ascii")


def decode_arrow_schema(pa, encoded: str):
    return pa.ipc.read_schema(pa.BufferReader(base64.b64decode(encoded)))


def _release_arrow_temporaries(pa) -> None:
    gc.collect()
    with contextlib.suppress(Exception):
        pa.default_memory_pool().release_unused()


def ctable_column_name_map(schema) -> dict[str, str]:
    """Return a mapping from Arrow field names to CTable-safe column names.

    Remaps invalid names (empty strings, names starting with '_', names
    containing '/') to safe substitutes like ``column_0``.
    """
    used: set[str] = set()
    result: dict[str, str] = {}
    for i, field in enumerate(schema):
        original = field.name
        try:
            _validate_column_name(original)
            candidate = original
        except ValueError:
            candidate = f"column_{i}"
        if candidate in used:
            base = candidate
            suffix = 1
            while f"{base}_{suffix}" in used:
                suffix += 1
            candidate = f"{base}_{suffix}"
        used.add(candidate)
        result[original] = candidate
    return result


def classify_columns(  # noqa: C901
    pa,
    schema,
    fixed_string_lengths: dict[str, int] | None = None,
    fixed_bytes_lengths: dict[str, int] | None = None,
    *,
    decode_dictionaries: bool = False,
    separate_nested_cols: bool = True,
):
    """Classify Parquet schema columns into importable categories."""
    fixed_cols: dict[str, object] = {}
    struct_wrap_cols: dict[str, object] = {}
    conversions: dict[str, dict[str, Any]] = {}
    nullable_scalars: list[str] = []
    fixed_string_lengths = fixed_string_lengths or {}
    fixed_bytes_lengths = fixed_bytes_lengths or {}

    for field in schema:
        t = field.type
        if pa.types.is_struct(t):
            if separate_nested_cols:
                # Let CTable.from_arrow() apply its normal struct flattening so
                # top-level structs become dotted leaf columns.
                fixed_cols[field.name] = field
                conversions[field.name] = {"conversion": "struct_flattened_to_columns"}
            else:
                struct_wrap_cols[field.name] = pa.list_(t)
                conversions[field.name] = {"conversion": "struct_wrapped_as_singleton_list"}
            continue
        if pa.types.is_list(t) or pa.types.is_large_list(t):
            value_type = t.value_type
            if pa.types.is_list(value_type) or pa.types.is_large_list(value_type):
                conversions[field.name] = {"conversion": "skipped", "reason": f"nested list: {t}"}
            else:
                fixed_cols[field.name] = field
            continue
        if pa.types.is_dictionary(t):
            vt = t.value_type
            if vt in (pa.string(), pa.large_string(), pa.utf8(), pa.large_utf8()):
                if decode_dictionaries:
                    # Decode to plain utf8.
                    fixed_cols[field.name] = pa.field(
                        field.name, pa.string(), nullable=field.nullable, metadata=field.metadata
                    )
                    conversions[field.name] = {
                        "conversion": "dictionary_decoded_to_utf8",
                        "ordered": bool(t.ordered),
                    }
                else:
                    fixed_cols[field.name] = field
                    conversions[field.name] = {
                        "conversion": "dictionary_preserved",
                        "ordered": bool(t.ordered),
                    }
            else:
                conversions[field.name] = {
                    "conversion": "skipped",
                    "reason": f"unsupported dictionary value type: {vt}",
                }
            continue
        if pa.types.is_boolean(t):
            fixed_cols[field.name] = field
            if field.nullable:
                nullable_scalars.append(field.name)
                conversions[field.name] = {"conversion": "nullable_scalar_sentinel"}
            continue
        if pa.types.is_integer(t) or pa.types.is_floating(t):
            fixed_cols[field.name] = field
            if field.nullable:
                nullable_scalars.append(field.name)
                conversions[field.name] = {"conversion": "nullable_scalar_sentinel"}
            continue
        if pa.types.is_timestamp(t):
            fixed_cols[field.name] = field
            if field.nullable:
                nullable_scalars.append(field.name)
                conversions[field.name] = {"conversion": "timestamp_nullable"}
            else:
                conversions[field.name] = {"conversion": "timestamp"}
            continue
        if pa.types.is_string(t) or pa.types.is_large_string(t):
            fixed_cols[field.name] = field
            if field.name in fixed_string_lengths:
                conversions[field.name] = {
                    "conversion": "fixed_string_nullable" if field.nullable else "fixed_string",
                    "max_length": fixed_string_lengths[field.name],
                }
            else:
                conversions[field.name] = {"conversion": "utf8_nullable" if field.nullable else "utf8"}
            continue
        if pa.types.is_binary(t) or pa.types.is_large_binary(t):
            fixed_cols[field.name] = field
            if field.name in fixed_bytes_lengths:
                conversions[field.name] = {
                    "conversion": "fixed_bytes_nullable" if field.nullable else "fixed_bytes",
                    "max_length": fixed_bytes_lengths[field.name],
                }
            else:
                conversions[field.name] = {"conversion": "vlbytes_nullable" if field.nullable else "vlbytes"}
            continue
        conversions[field.name] = {"conversion": "skipped", "reason": f"unsupported: {t}"}

    return fixed_cols, struct_wrap_cols, conversions, nullable_scalars


def build_import_schema(
    pa,
    original_schema,
    fixed_cols: dict,
    struct_wrap_cols: dict,
    timestamp_units: dict[str, str] | None = None,
    column_name_map: dict[str, str] | None = None,
):
    """Build the Arrow schema passed to CTable.from_arrow()."""
    timestamp_units = timestamp_units or {}
    column_name_map = column_name_map or {}
    fields = []
    for field in original_schema:
        ctable_name = column_name_map.get(field.name, field.name)
        if field.name in struct_wrap_cols:
            fields.append(pa.field(ctable_name, struct_wrap_cols[field.name], nullable=True))
        elif field.name in fixed_cols:
            unit = timestamp_units.get(field.name)
            if unit is not None:
                fields.append(
                    pa.field(ctable_name, pa.timestamp(unit, tz=field.type.tz), nullable=field.nullable)
                )
            else:
                # Use the field from fixed_cols in case it was remapped (e.g. dict→string)
                fc = fixed_cols[field.name]
                if hasattr(fc, "type") and fc.type != field.type:
                    # fc has the remapped type; use ctable_name for the field name
                    fields.append(
                        pa.field(
                            ctable_name,
                            fc.type,
                            nullable=fc.nullable,
                            metadata=fc.metadata if fc.metadata else None,
                        )
                    )
                elif ctable_name != field.name:
                    fields.append(
                        pa.field(ctable_name, field.type, nullable=field.nullable, metadata=field.metadata)
                    )
                else:
                    fields.append(field)
    return pa.schema(fields)


def candidate_fixed_scalar_columns(pa, schema, *, scan_strings: bool, scan_bytes: bool) -> list[str]:
    columns = []
    for field in schema:
        if (scan_strings and (pa.types.is_string(field.type) or pa.types.is_large_string(field.type))) or (
            scan_bytes and (pa.types.is_binary(field.type) or pa.types.is_large_binary(field.type))
        ):
            columns.append(field.name)
    return columns


def update_string_and_bytes_max_lengths(pa, pc, batch, max_lengths: dict[str, int]) -> None:
    for field in batch.schema:
        arr = batch.column(field.name)
        if pa.types.is_string(field.type) or pa.types.is_large_string(field.type):
            lengths = pc.utf8_length(arr)
        else:
            lengths = pc.binary_length(arr)
        batch_max = pc.max(lengths).as_py()
        if batch_max is not None:
            max_lengths[field.name] = max(max_lengths[field.name], int(batch_max))


def nullable_sentinel_adjusted_length(pa, field, max_length: int, null_policy) -> int:
    if not field.nullable:
        return max_length
    null_value = null_policy.column_null_values.get(
        field.name, null_policy.sentinel_for_arrow_type(pa, field.type)
    )
    return max(max_length, len(null_value)) if null_value is not None else max_length


def parse_float_trunc_prec_options(args) -> tuple[int | None, dict[str, int]]:
    """Parse --float-trunc-prec entries into (global_bits, per_column_bits)."""
    global_bits = None
    per_column: dict[str, int] = {}
    for raw in args.float_trunc_prec:
        if "=" in raw:
            name, value = raw.split("=", 1)
            name = name.strip()
            if not name:
                raise ValueError("--float-trunc-prec column name cannot be empty")
            try:
                bits = int(value)
            except ValueError as exc:
                raise ValueError(f"Invalid --float-trunc-prec value for column {name!r}: {value!r}") from exc
            if bits < 1 or bits > 64:
                raise ValueError("--float-trunc-prec bits must be in the range 1..64")
            per_column[name] = bits
        else:
            try:
                bits = int(raw)
            except ValueError as exc:
                raise ValueError(f"Invalid --float-trunc-prec value: {raw!r}") from exc
            if bits < 1 or bits > 64:
                raise ValueError("--float-trunc-prec bits must be in the range 1..64")
            global_bits = bits
    args.float_trunc_prec_global = global_bits
    args.float_trunc_prec_columns = per_column
    return global_bits, per_column


def build_float_trunc_column_cparams(pa, schema, args) -> dict[str, dict[str, Any]]:
    """Return per-column cparams for float columns selected by --float-trunc-prec."""
    global_bits = getattr(args, "float_trunc_prec_global", None)
    per_column = getattr(args, "float_trunc_prec_columns", {})
    if global_bits is None and not per_column:
        return {}

    fields_by_name = {field.name: field for field in schema}
    unknown = set(per_column) - set(fields_by_name)
    if unknown:
        names = ", ".join(sorted(unknown))
        raise KeyError(f"--float-trunc-prec references unknown imported columns: {names}")

    result: dict[str, dict[str, Any]] = {}
    for field in schema:
        if not pa.types.is_floating(field.type):
            if field.name in per_column:
                raise TypeError(
                    f"--float-trunc-prec can only be used with float columns; {field.name!r} is {field.type}"
                )
            continue
        bits = per_column.get(field.name, global_bits)
        if bits is None:
            continue
        max_bits = 23 if field.type.bit_width == 32 else 52
        if bits > max_bits:
            raise ValueError(
                f"--float-trunc-prec for column {field.name!r} is {bits}, "
                f"but float{field.type.bit_width} columns support at most {max_bits}"
            )
        result[field.name] = {
            "codec": blosc2.Codec[args.codec].value,
            "clevel": args.clevel,
            "use_dict": args.use_dict,
            "typesize": field.type.bit_width // 8,
            "filters": [blosc2.Filter.TRUNC_PREC.value, blosc2.Filter.SHUFFLE.value],
            "filters_meta": [bits, 0],
        }
    return result


def fixed_string_and_bytes_lengths_from_scan(pa, schema, args, max_lengths: dict[str, int]):
    from blosc2.ctable import get_null_policy

    null_policy = get_null_policy()
    fixed_string_lengths = {}
    fixed_bytes_lengths = {}
    for field in schema:
        max_length = max_lengths.get(field.name)
        if max_length is None:
            continue
        max_length = nullable_sentinel_adjusted_length(pa, field, max_length, null_policy)
        if (
            args.fixed_str_maxlen is not None
            and (pa.types.is_string(field.type) or pa.types.is_large_string(field.type))
            and max_length <= args.fixed_str_maxlen
        ):
            fixed_string_lengths[field.name] = args.fixed_str_maxlen
        elif (
            args.fixed_bytes_maxlen is not None
            and (pa.types.is_binary(field.type) or pa.types.is_large_binary(field.type))
            and max_length <= args.fixed_bytes_maxlen
        ):
            fixed_bytes_lengths[field.name] = args.fixed_bytes_maxlen
    return fixed_string_lengths, fixed_bytes_lengths


_TIMESTAMP_UNIT_NS = {"s": 1_000_000_000, "ms": 1_000_000, "us": 1_000, "ns": 1}
_TIMESTAMP_UNITS_COARSE_TO_FINE = ("s", "ms", "us", "ns")


def timestamp_columns(pa, schema) -> list[str]:
    return [field.name for field in schema if pa.types.is_timestamp(field.type)]


def initial_timestamp_divisibility(units: dict[str, str]) -> dict[str, dict[str, bool]]:
    return {
        name: {
            unit: True
            for unit in _TIMESTAMP_UNITS_COARSE_TO_FINE
            if _TIMESTAMP_UNIT_NS[unit] >= _TIMESTAMP_UNIT_NS[source_unit]
        }
        for name, source_unit in units.items()
    }


def update_timestamp_divisibility(
    batch, units: dict[str, str], divisible: dict[str, dict[str, bool]]
) -> None:
    for name, source_unit in units.items():
        arr = batch.column(batch.schema.get_field_index(name)).drop_null()
        if len(arr) == 0:
            continue
        values = arr.to_numpy(zero_copy_only=False).astype(f"datetime64[{source_unit}]").astype("int64")
        for unit in list(divisible[name]):
            if unit == source_unit:
                continue
            factor = _TIMESTAMP_UNIT_NS[unit] // _TIMESTAMP_UNIT_NS[source_unit]
            if factor > 1 and not bool((values % factor == 0).all()):
                divisible[name][unit] = False


def choose_timestamp_units(units: dict[str, str], divisible: dict[str, dict[str, bool]]) -> dict[str, str]:
    result = {}
    for name, source_unit in units.items():
        result[name] = source_unit
        for unit in _TIMESTAMP_UNITS_COARSE_TO_FINE:
            if unit in divisible[name] and divisible[name][unit]:
                result[name] = unit
                break
    return result


def infer_timestamp_units(pa, pf, args, schema) -> dict[str, str]:
    """Return target timestamp units for import according to --timestamp-unit."""
    columns = timestamp_columns(pa, schema)
    if args.timestamp_unit is None or not columns:
        return {}
    if args.timestamp_unit != "auto":
        return dict.fromkeys(columns, args.timestamp_unit)

    print("Pre-scanning timestamp units...")
    fields_by_name = {field.name: field for field in schema}
    units = {name: fields_by_name[name].type.unit for name in columns}
    divisible = initial_timestamp_divisibility(units)
    rows_done = 0
    total = pf.metadata.num_rows if args.max_rows is None else min(args.max_rows, pf.metadata.num_rows)
    for batch in pf.iter_batches(batch_size=args.parquet_batch_size, columns=columns):
        remaining = total - rows_done
        if remaining <= 0:
            break
        if len(batch) > remaining:
            batch = batch.slice(0, remaining)
        update_timestamp_divisibility(batch, units, divisible)
        rows_done += len(batch)

    result = choose_timestamp_units(units, divisible)
    changed = {name: unit for name, unit in result.items() if unit != units[name]}
    print(f"  timestamp columns: {len(columns):,}; unit changes: {len(changed):,}")
    for name, unit in sorted(changed.items()):
        print(f"    - {name}: {units[name]} -> {unit}")
    return result


def scan_string_and_bytes_lengths(pa, pf, args, schema) -> tuple[dict[str, int], dict[str, int]]:
    if args.fixed_str_maxlen is None and args.fixed_bytes_maxlen is None:
        return {}, {}

    import pyarrow.compute as pc

    columns = candidate_fixed_scalar_columns(
        pa,
        schema,
        scan_strings=args.fixed_str_maxlen is not None,
        scan_bytes=args.fixed_bytes_maxlen is not None,
    )
    if not columns:
        return {}, {}

    print("Pre-scanning string/binary column lengths...")
    rows_done = 0
    total = pf.metadata.num_rows if args.max_rows is None else min(args.max_rows, pf.metadata.num_rows)
    max_lengths = dict.fromkeys(columns, 0)
    for batch in pf.iter_batches(batch_size=args.parquet_batch_size, columns=columns):
        remaining = total - rows_done
        if remaining <= 0:
            break
        if len(batch) > remaining:
            batch = batch.slice(0, remaining)
        update_string_and_bytes_max_lengths(pa, pc, batch, max_lengths)
        rows_done += len(batch)

    fixed_string_lengths, fixed_bytes_lengths = fixed_string_and_bytes_lengths_from_scan(
        pa, schema, args, max_lengths
    )
    print(
        f"  fixed string columns: {len(fixed_string_lengths):,}; "
        f"fixed bytes columns: {len(fixed_bytes_lengths):,}"
    )
    return fixed_string_lengths, fixed_bytes_lengths


def transform_batch(
    pa,
    batch,
    selected_cols: list[str],
    struct_wrap_cols: dict,
    timestamp_units: dict[str, str],
    import_schema=None,
):
    """Apply import-time Arrow conversions; pass everything else through."""
    arrays = list(batch.columns)
    for name, unit in timestamp_units.items():
        idx = batch.schema.get_field_index(name)
        if idx < 0:
            continue
        field = batch.schema.field(idx)
        target_type = pa.timestamp(unit, tz=field.type.tz)
        arrays[idx] = batch.column(idx).cast(target_type, safe=True)
    for name, target_type in struct_wrap_cols.items():
        try:
            idx = batch.schema.get_field_index(name)
        except KeyError:
            continue
        if idx < 0:
            continue
        arr = batch.column(idx)
        arrays[idx] = pa.array([[v] if v is not None else None for v in arr.to_pylist()], type=target_type)
    if import_schema is not None:
        # Cast / rename arrays to match import_schema (e.g. dict→string, renamed columns).
        for i, field in enumerate(import_schema):
            if not arrays[i].type.equals(field.type):
                arrays[i] = arrays[i].cast(field.type, safe=True)
        return pa.record_batch(arrays, schema=import_schema)
    if not struct_wrap_cols and not timestamp_units:
        return batch
    return pa.record_batch(arrays, names=selected_cols)


def store_original_arrow_metadata(
    ct, original_schema, conversions: dict, column_name_map: dict | None = None
) -> None:
    column_name_map = column_name_map or {}
    fields_meta = {}
    for field in original_schema:
        entry = conversions.get(field.name)
        if entry is None:
            continue
        entry = dict(entry)
        ctable_name = column_name_map.get(field.name, field.name)
        if ctable_name != field.name:
            entry["ctable_name"] = ctable_name
        fields_meta[field.name] = entry
    ct._schema.metadata = {
        "arrow": {
            "schema_ipc_base64": encode_arrow_schema(original_schema),
            "fields": fields_meta,
        }
    }
    ct._storage.save_schema(schema_to_dict(ct._schema))


def ctable_store_kind(path: Path) -> str:
    if path.suffix == ".b2d":
        return "sparse directory (.b2d)"
    if path.suffix == ".b2z":
        return "compact zip (.b2z)"
    return f"unknown ({path.suffix or 'no suffix'})"


def print_import_plan(
    args,
    input_path,
    output_path,
    pf,
    parquet_schema,
    fixed_cols,
    struct_wrap_cols,
    conversions,
    nullable_scalars,
):
    utf8_cols = [n for n, e in conversions.items() if e.get("conversion") in {"utf8", "utf8_nullable"}]
    vlbytes_cols = [
        n for n, e in conversions.items() if e.get("conversion") in {"vlbytes", "vlbytes_nullable"}
    ]
    fixed_string_cols = [
        n for n, e in conversions.items() if e.get("conversion") in {"fixed_string", "fixed_string_nullable"}
    ]
    fixed_bytes_cols = [
        n for n, e in conversions.items() if e.get("conversion") in {"fixed_bytes", "fixed_bytes_nullable"}
    ]
    dict_cols = [n for n, e in conversions.items() if e.get("conversion") == "dictionary_preserved"]
    dict_decoded_cols = [
        n for n, e in conversions.items() if e.get("conversion") == "dictionary_decoded_to_utf8"
    ]
    flattened_structs = [
        n for n, e in conversions.items() if e.get("conversion") == "struct_flattened_to_columns"
    ]
    wrapped_structs = list(struct_wrap_cols)
    skipped = {n: e for n, e in conversions.items() if e.get("conversion") == "skipped"}
    print(f"Input:                 {input_path} ({input_path.stat().st_size / 1e6:.1f} MB)")
    print(f"Output:                {output_path}")
    print(f"CTable store:          {ctable_store_kind(output_path)}")
    print(f"Rows:                  {pf.metadata.num_rows:,}")
    if args.max_rows is not None:
        print(f"Rows to import:        {min(args.max_rows, pf.metadata.num_rows):,} (Parquet rows)")
    print(f"Parquet columns:       {len(parquet_schema)}")
    print(f"Imported columns:      {len(fixed_cols) + len(struct_wrap_cols)}")
    n_fixed_non_string = (
        len(fixed_cols) - len(utf8_cols) - len(vlbytes_cols) - len(dict_cols) - len(dict_decoded_cols)
    )
    print(f"  Fixed-width:         {n_fixed_non_string}")
    print(f"  Fixed strings:       {len(fixed_string_cols)}")
    print(f"  Fixed bytes:         {len(fixed_bytes_cols)}")
    print(f"  utf8:                {len(utf8_cols)}")
    print(f"  vlbytes:             {len(vlbytes_cols)}")
    print(f"  Dictionary:          {len(dict_cols)}")
    if dict_decoded_cols:
        print(f"  Dict→utf8:           {len(dict_decoded_cols)}")
    print(f"  Struct→columns:      {len(flattened_structs)}")
    print(f"  Struct→list:         {len(wrapped_structs)}")
    print(f"  Nullable scalars:    {len(nullable_scalars)}")
    print(f"  Skipped unsupported: {len(skipped)}")
    for name, entry in skipped.items():
        print(f"    - {name}: {entry['reason']}")
    if args.fixed_str_maxlen is not None:
        print(f"Fixed string maxlen:   {args.fixed_str_maxlen:,} characters")
    if args.fixed_bytes_maxlen is not None:
        print(f"Fixed bytes maxlen:    {args.fixed_bytes_maxlen:,} bytes")
    print(f"Parquet batch size:    {args.parquet_batch_size:,}")
    print(f"Blosc2 batch size:     {args.blosc2_batch_size:,}")
    if args.blosc2_items_per_block is not None:
        print(f"Blosc2 items/block:    {args.blosc2_items_per_block:,}")
    print(f"List serializer:       {args.list_serializer}")
    print(f"Codec / level:         {args.codec} / {args.clevel}")
    print(f"Use dict:              {args.use_dict}")
    if args.chunks is not None:
        print(f"Chunks:                {args.chunks:,}")
    if args.blocks is not None:
        print(f"Blocks:                {args.blocks:,}")
    trunc_global = getattr(args, "float_trunc_prec_global", None)
    trunc_columns = getattr(args, "float_trunc_prec_columns", {})
    if trunc_global is not None:
        print(f"Float trunc precision: {trunc_global} bits (all float columns)")
    if trunc_columns:
        formatted = ", ".join(f"{name}={bits}" for name, bits in sorted(trunc_columns.items()))
        print(f"Float trunc columns:   {formatted}")
    if args.timestamp_unit is not None:
        print(f"Timestamp unit:        {args.timestamp_unit}")
    print()


def progress_batches(pa, pf, args, selected_cols, struct_wrap_cols, timestamp_units, import_schema=None):
    rows_done = 0
    t0 = time.perf_counter()
    total = pf.metadata.num_rows if args.max_rows is None else min(args.max_rows, pf.metadata.num_rows)
    for batch_n, raw_batch in enumerate(
        pf.iter_batches(batch_size=args.parquet_batch_size, columns=selected_cols), start=1
    ):
        remaining = total - rows_done
        if remaining <= 0:
            break
        if len(raw_batch) > remaining:
            raw_batch = raw_batch.slice(0, remaining)
        report_batch_mem = args.mem_report and batch_n % args.mem_every == 0
        if report_batch_mem:
            memory_report(f"batch {batch_n} after parquet read", pa)
        batch = transform_batch(
            pa, raw_batch, selected_cols, struct_wrap_cols, timestamp_units, import_schema
        )
        if report_batch_mem:
            memory_report(f"batch {batch_n} after transform", pa)
        rows_done += len(batch)
        elapsed = time.perf_counter() - t0
        rate = rows_done / elapsed if elapsed > 0 else 0.0
        eta = (total - rows_done) / rate if rate > 0 else 0.0
        if args.progress and (batch_n % args.batch_report_every == 0 or rows_done >= total):
            print(
                f"  batch {batch_n:4d}  {rows_done:>12,}/{total:,}  "
                f"{elapsed:7.1f}s  {rate / 1e3:7.1f}k rows/s  ETA {eta:6.0f}s",
                flush=True,
            )
        if report_batch_mem:
            memory_report(f"batch {batch_n} before ctable write", pa)
        yield batch
        if report_batch_mem:
            memory_report(f"batch {batch_n} after ctable write", pa)


def _flatten_root_batches_with_progress(
    pa,
    pf,
    inner_schema,
    args,
    capacity_hint=None,
):
    """Yield flattened :class:`pyarrow.RecordBatch` objects from an unnamed-root Parquet file.

    Reads Parquet batches, flattens the outer ``list<struct<...>>`` column via
    ``ListArray.flatten()``, and honours ``args.max_rows`` as an element-level
    row limit.  When ``args.progress`` is enabled, progress is printed per
    Parquet batch according to ``args.batch_report_every``.

    Each flattened Parquet batch is yielded as a single write to CTable so that
    the per-write Python/Arrow overhead is amortised over as many rows as
    possible.  Batches exceeding ``MAX_ELEMENT_WRITE_BATCH`` are split into
    cap-sized chunks to bound memory usage.
    """
    rows_done = 0
    max_rows = args.max_rows
    t0 = time.perf_counter()
    # total_str is the CTable-row (element) limit for the progress display.
    total_str = f"{max_rows:,} CTable rows" if max_rows is not None else "?"
    # Use capacity_hint as the estimated total for ETA when max_rows is not set.
    estimated_total = max_rows if max_rows is not None else capacity_hint

    for parquet_batch_n, raw_batch in enumerate(
        pf.iter_batches(batch_size=args.parquet_batch_size), start=1
    ):
        if max_rows is not None and rows_done >= max_rows:
            break

        report_batch_mem = args.mem_report and parquet_batch_n % args.mem_every == 0
        if report_batch_mem:
            memory_report(f"batch {parquet_batch_n} after parquet read", pa)

        list_array = raw_batch.column(0)
        struct_values = list_array.flatten()  # skips null outer-list rows

        if len(struct_values) == 0:
            continue

        if max_rows is not None:
            remaining = max_rows - rows_done
            if len(struct_values) > remaining:
                struct_values = struct_values.slice(0, remaining)

        # Yield the whole flattened batch as one write; split only when it
        # exceeds MAX_ELEMENT_WRITE_BATCH to bound peak memory.
        n_elems = len(struct_values)

        elapsed = time.perf_counter() - t0
        rate = rows_done / elapsed if elapsed > 0 and rows_done > 0 else 0.0
        eta_str = (
            f"  ETA {(estimated_total - rows_done) / rate:6.0f}s"
            if rate > 0 and estimated_total is not None
            else ""
        )
        report_progress = parquet_batch_n % args.batch_report_every == 0 or (
            max_rows is not None and rows_done + n_elems >= max_rows
        )
        n_writes = (n_elems + MAX_ELEMENT_WRITE_BATCH - 1) // MAX_ELEMENT_WRITE_BATCH
        if args.progress and report_progress:
            print(
                f"  parquet batch {parquet_batch_n:4d}: "
                f"{n_elems:>12,} CTable rows -> {n_writes:,} write(s)  "
                f"done {rows_done:>12,}/{total_str}  "
                f"{elapsed:7.1f}s  {rate / 1e3:7.1f}k rows/s{eta_str}",
                flush=True,
            )

        for offset in range(0, n_elems, MAX_ELEMENT_WRITE_BATCH):
            chunk = struct_values.slice(offset, min(MAX_ELEMENT_WRITE_BATCH, n_elems - offset))
            sub_batch = pa.RecordBatch.from_struct_array(chunk)
            rows_done += len(sub_batch)
            yield sub_batch

        if report_batch_mem:
            memory_report(f"batch {parquet_batch_n} after flatten+write", pa)

        if max_rows is not None and rows_done >= max_rows:
            break


def _apply_parquet_batch_memory_budget(args, sample, n_outer_sampled: int) -> None:
    """Shrink an auto parquet batch size so one Arrow read batch fits the budget.

    Nested list<struct> batches amplify several-fold downstream (flatten + cast
    + write buffers + Arrow pool), so an auto-chosen parquet batch size is capped
    to keep peak RSS well under ~1 GB.  An explicit --parquet-batch-size is left
    untouched.

    Opt-in via --reduce-mem: it trades import speed for lower peak RSS, so the
    default keeps the original (larger) auto batch sizes.
    """
    if not getattr(args, "reduce_mem", False):
        return
    if not getattr(args, "parquet_batch_size_auto", False):
        return
    bytes_per_outer = sample.nbytes / n_outer_sampled
    if bytes_per_outer > 0:
        budget_rows = max(1, int(PARQUET_BATCH_ARROW_BUDGET / bytes_per_outer))
        args.parquet_batch_size = min(args.parquet_batch_size, budget_rows)


def import_unnamed_root_separate_cols(  # noqa: C901
    args,
    input_path: Path,
    output_path: Path,
    pa,
    pf,
    parquet_schema,
) -> list[str]:
    """Import an unnamed-root ``list<struct<...>>`` Parquet file with nested column separation.

    Each element of the unnamed root list becomes a CTable row.  Struct leaves
    are stored as separate physical columns with dotted logical paths such as
    ``trip.begin.lon`` and ``payment.fare``.

    Returns the list of imported CTable column names.
    """

    inner_schema = blosc2.CTable._inner_schema_for_unnamed_root(pa, parquet_schema)
    flat_inner_schema = blosc2.CTable._flatten_arrow_struct_schema(pa, inner_schema)
    float_trunc_column_cparams = build_float_trunc_column_cparams(pa, flat_inner_schema, args)
    total_parquet_rows = pf.metadata.num_rows if pf.metadata is not None else None

    # ------------------------------------------------------------------
    # Estimate total element count by sampling the first Parquet batch.
    # This is used as capacity_hint so that compute_chunks_blocks() picks
    # chunk/block sizes proportional to the actual data volume rather than
    # defaulting to (1, 1) when the element count is unknown.
    # pf.iter_batches() creates a fresh iterator each call, so sampling
    # here does not affect the import iterator created later.
    # ------------------------------------------------------------------
    capacity_hint = None
    estimated_batch_rows = None
    if total_parquet_rows is not None and total_parquet_rows > 0:
        try:
            # Sample only a few outer rows: enough for the per-outer-row ratio
            # and byte estimate, while avoiding a large transient Arrow batch
            # (which the Arrow pool would retain and inflate peak RSS).
            sample_rows = min(args.parquet_batch_size, total_parquet_rows, 64)
            sample = next(pf.iter_batches(batch_size=sample_rows), None)
            if sample is not None and len(sample) > 0:
                n_outer_sampled = len(sample)
                n_elems_sampled = len(sample.column(0).flatten())
                avg_per_outer_row = n_elems_sampled / n_outer_sampled
                _apply_parquet_batch_memory_budget(args, sample, n_outer_sampled)
                estimated_batch_rows = max(1, round(args.parquet_batch_size * avg_per_outer_row))
                estimate = round(total_parquet_rows * avg_per_outer_row)
                if args.max_rows is None:
                    estimate = round(estimate * UNNAMED_ROOT_CAPACITY_SAFETY)
                else:
                    estimate = min(estimate, args.max_rows)
                capacity_hint = max(1, estimate)
        except Exception:
            pass  # sampling failure is non-fatal; from_arrow falls back to _EXPECTED_SIZE_DEFAULT

    if args.blosc2_batch_size is None:
        if args.list_serializer == "arrow":
            # Arrow list storage appends incoming Arrow chunks directly, without
            # materializing Python nested-list objects.  Use the natural flattened
            # Parquet-batch scale (about 1M rows for Chicago taxi), capped only for
            # pathological batches, so the displayed BatchArray size matches the
            # actual write granularity better than the absolute cap would.
            args.blosc2_batch_size = min(
                MAX_ELEMENT_WRITE_BATCH,
                estimated_batch_rows if estimated_batch_rows is not None else MAX_ELEMENT_WRITE_BATCH,
            )
        else:
            # Msgpack list storage materializes nested Arrow list data as Python objects
            # before serializing.  Keep its internal BatchArray batch_rows at Blosc2's
            # cache-tuned block granularity instead of the larger Arrow write scale.
            if capacity_hint is not None:
                _, blocks = blosc2.compute_chunks_blocks((capacity_hint,))
                args.blosc2_batch_size = max(1, blocks[0])
            else:
                args.blosc2_batch_size = DEFAULT_BATCH_SIZE

    print(f"Input:                 {input_path} ({input_path.stat().st_size / 1e6:.1f} MB)")
    print(f"Output:                {output_path}")
    print(f"CTable store:          {ctable_store_kind(output_path)}")
    print("Mode:                  unnamed-root list<struct> flattening")
    print("Nested columns:        separated into dotted CTable columns")
    if total_parquet_rows is not None:
        print(f"Parquet rows:          {total_parquet_rows:,}")
    if capacity_hint is not None:
        print(f"Est. CTable rows:      ~{capacity_hint:,}")
    n_inner = len(inner_schema)
    print(f"Inner struct fields:   {n_inner}")
    for f in inner_schema:
        print(f"  {f.name}: {f.type}")
    if args.max_rows is not None:
        print(f"Max CTable rows:       {args.max_rows:,} (list elements)")
    print(f"Parquet batch size:    {args.parquet_batch_size:,} outer rows")
    blosc2_batch_note = (
        f"auto, max: {MAX_ELEMENT_WRITE_BATCH:,}"
        if getattr(args, "blosc2_batch_size_auto", False)
        else f"max: {MAX_ELEMENT_WRITE_BATCH:,}"
    )
    print(f"Blosc2 batch size:     {args.blosc2_batch_size:,} BatchArray rows ({blosc2_batch_note})")
    if args.blosc2_items_per_block is not None:
        print(f"Blosc2 items/block:    {args.blosc2_items_per_block:,}")
    print(f"List serializer:       {args.list_serializer}")
    print(f"Codec / level:         {args.codec} / {args.clevel}")
    print(f"Use dict:              {args.use_dict}")
    if args.chunks is not None:
        print(f"Chunks:                {args.chunks:,}")
    if args.blocks is not None:
        print(f"Blocks:                {args.blocks:,}")
    trunc_global = getattr(args, "float_trunc_prec_global", None)
    trunc_columns = getattr(args, "float_trunc_prec_columns", {})
    if trunc_global is not None:
        print(f"Float trunc precision: {trunc_global} bits (all float columns)")
    if trunc_columns:
        formatted = ", ".join(f"{name}={bits}" for name, bits in sorted(trunc_columns.items()))
        print(f"Float trunc columns:   {formatted}")
    print()

    cparams = blosc2.CParams(codec=blosc2.Codec[args.codec], clevel=args.clevel, use_dict=args.use_dict)
    t0 = time.perf_counter()
    maybe_memory_report(args, "before CTable import", pa)

    ct = blosc2.CTable.from_arrow(
        inner_schema,
        _flatten_root_batches_with_progress(pa, pf, inner_schema, args, capacity_hint=capacity_hint),
        urlpath=str(output_path),
        mode="w",
        cparams=cparams,
        capacity_hint=capacity_hint,
        auto_null_sentinels=True,
        blosc2_batch_size=args.blosc2_batch_size,
        blosc2_items_per_block=args.blosc2_items_per_block,
        list_serializer=args.list_serializer,
        column_cparams=float_trunc_column_cparams or None,
        create_summary_index=args.create_summary_index,
        chunks=args.chunks,
        blocks=args.blocks,
    )

    maybe_memory_report(args, "after CTable import", pa)

    maybe_memory_report(args, "after metadata save", pa)

    elapsed = time.perf_counter() - t0
    rows = len(ct)
    cols = len(ct.col_names)
    col_names = list(ct.col_names)
    ct.close()

    maybe_memory_report(args, "after CTable close", pa)

    output_size = (
        output_path.stat().st_size
        if output_path.is_file()
        else sum(f.stat().st_size for f in output_path.rglob("*") if f.is_file())
    )
    print(f"Done in {elapsed:.2f}s")
    print(f"Element rows imported: {rows:,}")
    print(f"Columns imported:      {cols}")
    print(f"Output size:           {output_size / 1e6:.1f} MB")
    return col_names


def import_parquet_to_ctable(args, input_path: Path, output_path: Path):  # noqa: C901
    if args.parquet_batch_size <= 0:
        raise ValueError("--parquet-batch-size must be positive")
    if args.blosc2_batch_size is not None and args.blosc2_batch_size <= 0:
        raise ValueError("--blosc2-batch-size must be positive")
    if args.blosc2_items_per_block is not None and args.blosc2_items_per_block <= 0:
        raise ValueError("--blosc2-items-per-block must be positive")
    if args.fixed_str_maxlen is not None and args.fixed_str_maxlen <= 0:
        raise ValueError("--fixed-str-maxlen must be positive")
    if args.fixed_bytes_maxlen is not None and args.fixed_bytes_maxlen <= 0:
        raise ValueError("--fixed-bytes-maxlen must be positive")
    if args.chunks is not None and args.chunks <= 0:
        raise ValueError("--chunks must be positive")
    if args.blocks is not None and args.blocks <= 0:
        raise ValueError("--blocks must be positive")
    if args.chunks is not None and args.blocks is not None and args.blocks > args.chunks:
        raise ValueError("--blocks cannot be greater than --chunks")
    parse_float_trunc_prec_options(args)
    if args.max_rows is not None and args.max_rows < 0:
        raise ValueError("--max-rows must be non-negative")
    if args.mem_every <= 0:
        raise ValueError("--mem-every must be positive")
    if args.batch_report_every <= 0:
        raise ValueError("--batch-report-every must be positive")
    if output_path.suffix not in {".b2z", ".b2d"}:
        raise ValueError("output_path must use the .b2z (compact) or .b2d (sparse) extension")
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    prepare_output(output_path, args.overwrite)

    pa, pq = require_pyarrow()
    maybe_memory_report(args, "after pyarrow import", pa)
    pf = pq.ParquetFile(input_path)
    maybe_memory_report(args, "after ParquetFile open", pa)
    parquet_schema = pf.schema_arrow

    # ------------------------------------------------------------------
    # Early dispatch: --separate-nested-cols for unnamed-root datasets
    # ------------------------------------------------------------------
    if getattr(args, "separate_nested_cols", False) and blosc2.CTable._detect_unnamed_root_list_struct(
        pa, parquet_schema
    ):
        return import_unnamed_root_separate_cols(args, input_path, output_path, pa, pf, parquet_schema)

    fixed_string_lengths, fixed_bytes_lengths = scan_string_and_bytes_lengths(pa, pf, args, parquet_schema)
    maybe_memory_report(args, "after string/binary length scan", pa)

    timestamp_units = infer_timestamp_units(pa, pf, args, parquet_schema)
    maybe_memory_report(args, "after timestamp unit scan", pa)

    fixed_cols, struct_wrap_cols, conversions, nullable_scalars = classify_columns(
        pa,
        parquet_schema,
        fixed_string_lengths,
        fixed_bytes_lengths,
        decode_dictionaries=getattr(args, "decode_dictionaries", False),
        separate_nested_cols=getattr(args, "separate_nested_cols", True),
    )
    maybe_memory_report(args, "after column classification", pa)

    selected_cols = [f.name for f in parquet_schema if f.name in fixed_cols or f.name in struct_wrap_cols]
    column_name_map = ctable_column_name_map(parquet_schema)
    import_schema = build_import_schema(
        pa, parquet_schema, fixed_cols, struct_wrap_cols, timestamp_units, column_name_map
    )
    fixed_scalar_lengths = {
        column_name_map.get(name, name): length
        for name, length in {**fixed_string_lengths, **fixed_bytes_lengths}.items()
    } or None
    float_trunc_column_cparams = build_float_trunc_column_cparams(pa, import_schema, args)
    maybe_memory_report(args, "after import schema build", pa)

    print_import_plan(
        args,
        input_path,
        output_path,
        pf,
        parquet_schema,
        fixed_cols,
        struct_wrap_cols,
        conversions,
        nullable_scalars,
    )

    t0 = time.perf_counter()
    maybe_memory_report(args, "before CTable import", pa)

    ct = blosc2.CTable.from_arrow(
        import_schema,
        progress_batches(pa, pf, args, selected_cols, struct_wrap_cols, timestamp_units, import_schema),
        urlpath=str(output_path),
        mode="w",
        cparams=blosc2.CParams(codec=blosc2.Codec[args.codec], clevel=args.clevel, use_dict=args.use_dict),
        capacity_hint=(
            pf.metadata.num_rows if args.max_rows is None else min(args.max_rows, pf.metadata.num_rows)
        ),
        string_max_length=fixed_scalar_lengths,
        auto_null_sentinels=True,
        blosc2_batch_size=args.blosc2_batch_size,
        blosc2_items_per_block=args.blosc2_items_per_block,
        list_serializer=args.list_serializer,
        column_cparams=float_trunc_column_cparams or None,
        create_summary_index=args.create_summary_index,
        chunks=args.chunks,
        blocks=args.blocks,
    )
    maybe_memory_report(args, "after CTable import", pa)
    store_original_arrow_metadata(ct, parquet_schema, conversions, column_name_map)
    maybe_memory_report(args, "after metadata save", pa)
    elapsed = time.perf_counter() - t0
    rows = len(ct)
    cols = len(ct.col_names)
    ct.close()
    maybe_memory_report(args, "after CTable close", pa)

    output_size = (
        output_path.stat().st_size
        if output_path.is_file()
        else sum(f.stat().st_size for f in output_path.rglob("*") if f.is_file())
    )
    print(f"Done in {elapsed:.2f}s")
    print(f"Rows imported:         {rows:,}")
    print(f"Columns imported:      {cols}")
    print(f"Output size:           {output_size / 1e6:.1f} MB")
    return selected_cols


def original_schema_from_ctable(pa, ct):
    arrow_meta = ct._schema.metadata.get("arrow", {})
    encoded = arrow_meta.get("schema_ipc_base64")
    if encoded:
        return decode_arrow_schema(pa, encoded)
    return None


def unwrap_singleton_list(pa, arr, arrow_type):
    return pa.array(
        [None if cell is None or len(cell) == 0 else cell[0] for cell in arr.to_pylist()], type=arrow_type
    )


def export_ctable_to_parquet(input_path: Path, output_path: Path, *, batch_size: int, overwrite: bool):
    pa, pq = require_pyarrow()
    if batch_size <= 0:
        raise ValueError("--parquet-batch-size must be positive")
    prepare_output(output_path, overwrite)
    ct = blosc2.CTable.open(str(input_path))
    original_schema = original_schema_from_ctable(pa, ct)
    fields_meta = ct._schema.metadata.get("arrow", {}).get("fields", {})
    export_names = [
        name
        for name in (original_schema.names if original_schema is not None else ct.col_names)
        if name in ct.col_names
    ]
    export_schema = (
        pa.schema([original_schema.field(name) for name in export_names])
        if original_schema is not None
        else ct._arrow_schema_for_columns(export_names)
    )

    singleton_list_conversions = {
        "struct_wrapped_as_singleton_list",
        "nullable_scalar_wrapped_as_singleton_list",
        "long_nullable_scalar_wrapped_as_singleton_list",
        "scalar_string_promoted_after_overflow",
    }

    t0 = time.perf_counter()
    with pq.ParquetWriter(output_path, export_schema, compression="zstd") as writer:
        for batch in ct.iter_arrow_batches(columns=export_names, batch_size=batch_size):
            arrays = []
            for name in export_names:
                arr = batch.column(name)
                meta = fields_meta.get(name, {})
                field = export_schema.field(name)
                conversion = meta.get("conversion", "")
                if conversion in singleton_list_conversions:
                    arr = unwrap_singleton_list(pa, arr, field.type)
                elif conversion in {"utf8", "utf8_nullable", "vlbytes", "vlbytes_nullable"}:
                    if str(arr.type) != str(field.type):
                        arr = arr.cast(field.type)
                elif conversion in {"dictionary_preserved"}:
                    # CTable emits dictionary<int32, string>; restore original type if needed.
                    if str(arr.type) != str(field.type):
                        arr = arr.cast(field.type, safe=True)
                elif conversion in {"dictionary_decoded_to_utf8"}:
                    # Was decoded to utf8 on import; restore as dictionary type on export.
                    if pa.types.is_dictionary(field.type):
                        encoded = pa.DictionaryArray.from_arrays(
                            *pa.array(arr.to_pylist())
                            .dictionary_encode()
                            .unify_dictionaries([pa.array(arr.to_pylist()).dictionary_encode()]),
                            ordered=field.type.ordered,
                        )
                        arr = encoded.cast(field.type)
                    elif str(arr.type) != str(field.type):
                        arr = arr.cast(field.type)
                elif str(arr.type) != str(field.type):
                    arr = pa.array(arr.to_pylist(), type=field.type)
                arrays.append(arr)
            out_batch = pa.record_batch(arrays, schema=export_schema)
            writer.write_table(pa.Table.from_batches([out_batch]), row_group_size=len(out_batch))
    elapsed = time.perf_counter() - t0
    rows = len(ct)
    ct.close()
    print(f"Exported {rows:,} rows and {len(export_names)} columns to {output_path} in {elapsed:.2f}s")
    return export_names


def read_parquet_prefix(pa, pq, path: Path, columns: list[str], max_rows: int | None):
    if max_rows is None:
        return pq.read_table(path, columns=columns)
    pf = pq.ParquetFile(path)
    schema = pa.schema([pf.schema_arrow.field(name) for name in columns])
    batches = []
    rows_done = 0
    for batch in pf.iter_batches(batch_size=DEFAULT_BATCH_SIZE, columns=columns):
        remaining = max_rows - rows_done
        if remaining <= 0:
            break
        if len(batch) > remaining:
            batch = batch.slice(0, remaining)
        batches.append(batch)
        rows_done += len(batch)
    return pa.Table.from_batches(batches, schema=schema)


def assess_parquet_difference(
    original_path: Path, roundtrip_path: Path, exported_cols: list[str], max_rows: int | None = None
):
    pa, pq = require_pyarrow()
    orig_pf = pq.ParquetFile(original_path)
    rt_pf = pq.ParquetFile(roundtrip_path)
    original_schema = orig_pf.schema_arrow
    roundtrip_schema = rt_pf.schema_arrow
    common = [
        name for name in exported_cols if name in original_schema.names and name in roundtrip_schema.names
    ]
    missing = [name for name in original_schema.names if name not in roundtrip_schema.names]

    orig = read_parquet_prefix(pa, pq, original_path, common, max_rows)
    rt = pq.read_table(roundtrip_path, columns=common)
    differing = []
    type_diffs = []
    null_diffs = []
    for name in common:
        if str(original_schema.field(name).type) != str(roundtrip_schema.field(name).type):
            type_diffs.append(name)
        if orig[name].null_count != rt[name].null_count:
            null_diffs.append((name, orig[name].null_count, rt[name].null_count))
        if not orig[name].equals(rt[name]):
            differing.append(name)

    print("\nRoundtrip assessment")
    print(f"  Original rows:       {orig_pf.metadata.num_rows:,}")
    if max_rows is not None:
        print(f"  Original rows compared: {orig.num_rows:,}")
    print(f"  Roundtrip rows:      {rt_pf.metadata.num_rows:,}")
    print(f"  Original columns:    {len(original_schema)}")
    print(f"  Roundtrip columns:   {len(roundtrip_schema)}")
    print(f"  Missing columns:     {len(missing)}")
    for name in missing:
        print(f"    - {name}: not imported/exported")
    print(f"  Type differences:    {len(type_diffs)}")
    for name in type_diffs:
        print(f"    - {name}: {original_schema.field(name).type} -> {roundtrip_schema.field(name).type}")
    print(f"  Null-count diffs:    {len(null_diffs)}")
    for name, a, b in null_diffs[:20]:
        print(f"    - {name}: {a} -> {b}")
    if len(null_diffs) > 20:
        print(f"    ... {len(null_diffs) - 20} more")
    print(f"  Value differences:   {len(differing)} of {len(common)} compared columns")
    if differing:
        print("  First value-different columns:")
        for name in differing[:20]:
            print(f"    - {name}")
    print(f"  Original size:       {original_path.stat().st_size / 1e6:.1f} MB")
    print(f"  Roundtrip size:      {roundtrip_path.stat().st_size / 1e6:.1f} MB")


def _run_command(args) -> int:
    if args.export:
        input_path = args.input_path
        output_path = args.output_path or _default_export_output(input_path)
        export_ctable_to_parquet(
            input_path, output_path, batch_size=args.parquet_batch_size, overwrite=args.overwrite
        )
        return 0
    if args.roundtrip:
        input_path = args.input_path
        b2_path = args.output_path or _default_import_output(input_path)
        roundtrip_path = _default_roundtrip_output(input_path)
        selected = import_parquet_to_ctable(args, input_path, b2_path)
        exported = export_ctable_to_parquet(
            b2_path, roundtrip_path, batch_size=args.parquet_batch_size, overwrite=True
        )
        assess_parquet_difference(input_path, roundtrip_path, exported or selected, max_rows=args.max_rows)
        return 0

    output_path = args.output_path or _default_import_output(args.input_path)
    import_parquet_to_ctable(args, args.input_path, output_path)
    return 0


def _run_profiled(args) -> int:
    profiler = cProfile.Profile()
    profiler.enable()
    try:
        return _run_command(args)
    finally:
        profiler.disable()
        stream = io.StringIO()
        stats = pstats.Stats(profiler, stream=stream).sort_stats("cumulative")
        stats.print_stats(50)
        print("\n[cProfile] Top cumulative-time functions\n")
        print(stream.getvalue().rstrip())


def _option_present(argv: list[str], option: str) -> bool:
    return any(arg == option or arg.startswith(option + "=") for arg in argv)


def average_parquet_row_group_size(input_path: Path) -> int | None:
    if input_path.suffix != ".parquet" or not input_path.exists():
        return None
    try:
        _, pq = require_pyarrow()
        pf = pq.ParquetFile(input_path)
    except Exception:
        return None
    metadata = pf.metadata
    if metadata is None or metadata.num_row_groups <= 0 or metadata.num_rows <= 0:
        return None
    return max(1, round(metadata.num_rows / metadata.num_row_groups))


def is_unnamed_root_parquet_input(input_path: Path) -> bool:
    if input_path.suffix != ".parquet" or not input_path.exists():
        return False
    try:
        pa, pq = require_pyarrow()
        pf = pq.ParquetFile(input_path)
        return blosc2.CTable._detect_unnamed_root_list_struct(pa, pf.schema_arrow)
    except Exception:
        return False


def resolve_default_batch_sizes(args, *, parquet_specified: bool, blosc2_specified: bool) -> None:
    if getattr(args, "separate_nested_cols", False) and is_unnamed_root_parquet_input(args.input_path):
        # In separate-nested mode the two batch-size options use different units:
        # Parquet batches are outer rows, while Blosc2 batches are flattened
        # CTable rows.  Keep them independent so a large write batch does not
        # accidentally imply a huge Parquet read batch (and vice versa).
        args.parquet_batch_size_auto = not parquet_specified
        if not parquet_specified:
            args.parquet_batch_size = average_parquet_row_group_size(args.input_path) or DEFAULT_BATCH_SIZE
        if not blosc2_specified:
            # Defer separate-nested defaults until import, where we have a sampled
            # estimate of flattened CTable rows per Parquet batch.  Arrow uses that
            # natural per-Parquet-batch scale; msgpack uses a smaller blocks-based
            # scale because it materializes nested Python objects before serializing.
            args.blosc2_batch_size = None
        return

    if parquet_specified and not blosc2_specified:
        args.blosc2_batch_size = args.parquet_batch_size
    elif blosc2_specified and not parquet_specified:
        args.parquet_batch_size = args.blosc2_batch_size
    elif not parquet_specified and not blosc2_specified:
        default = average_parquet_row_group_size(args.input_path) or DEFAULT_BATCH_SIZE
        args.parquet_batch_size = default
        args.blosc2_batch_size = default


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else list(argv)
    args = build_parser().parse_args(argv)

    parquet_specified = _option_present(argv, "--parquet-batch-size") or _option_present(
        argv, "--batch-size"
    )
    blosc2_specified = _option_present(argv, "--blosc2-batch-size")
    args.blosc2_batch_size_auto = not blosc2_specified
    resolve_default_batch_sizes(args, parquet_specified=parquet_specified, blosc2_specified=blosc2_specified)

    if args.profile:
        return _run_profiled(args)
    return _run_command(args)


if __name__ == "__main__":
    raise SystemExit(main())
