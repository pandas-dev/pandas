from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

from dataclasses import dataclass
import os
from pathlib import Path
import re

_PANDAS_ROOT: Path = Path(__file__).resolve().parents[1]
_PANDAS_ROOT_STR: str = os.path.realpath(str(_PANDAS_ROOT)).replace(os.sep, "/").lower()

DEFAULT_PHASE = "other pandas"
RULES_VERSION = "v3-module-prefix-no-pd-import"


def _norm_path(s: str) -> str:
    return os.path.realpath(s).replace(os.sep, "/").lower()


def _split_location(location: str) -> tuple[str, str]:
    loc = location
    if ":" in loc:
        head, tail = loc.rsplit(":", 1)
        if "/" in head or "\\" in head:
            return head, tail
    return loc, ""


_ABI_TAG_RE = re.compile(r"\.(cpython|pypy|abi3)[0-9a-zA-Z_\-]*$", re.IGNORECASE)


def _strip_ext_and_abi(rel: str) -> str:
    if rel.endswith(".py"):
        rel = rel[:-3]
    else:
        for ext in (".so", ".pyd", ".dll", ".dylib"):
            if rel.endswith(ext):
                rel = rel[: -len(ext)]
                break
    rel = _ABI_TAG_RE.sub("", rel)
    return rel


def _file_to_module(filename: str) -> str:
    fn = _norm_path(filename)
    if fn.startswith(_PANDAS_ROOT_STR):
        rel = fn[len(_PANDAS_ROOT_STR) :].lstrip("/")
    else:
        idx = fn.find("/pandas/")
        if idx == -1:
            return ""
        rel = fn[idx + 1 :].lstrip("/")

    rel = rel.strip("/")
    if not rel:
        return "pandas"

    rel = _strip_ext_and_abi(rel)

    dotted = "pandas." + rel.replace("/", ".")
    if dotted.endswith(".__init__"):
        dotted = dotted[: -len(".__init__")]
    return dotted


@dataclass(frozen=True)
class BucketRule:
    name: str
    prefixes: tuple[str, ...]
    keywords: tuple[str, ...] = ()


DEFAULT_RULES: list[BucketRule] = [
    # Parsers (CSV/FWF/table) live under pandas.io.parsers (readers.py etc.)
    BucketRule(
        "io/parsers",
        prefixes=("pandas.io.parsers.",),
        keywords=("read_csv", "read_table", "read_fwf", "TextFileReader"),
    ),
    # CSV formatting/writing is under pandas.io.formats.csvs
    BucketRule(
        "io/csv-write",
        prefixes=("pandas.io.formats.csvs",),
        keywords=("to_csv", "CSVFormatter"),
    ),
    # Excel I/O is under pandas.io.excel.*
    BucketRule(
        "io/excel",
        prefixes=("pandas.io.excel.",),
        keywords=("read_excel", "to_excel", "ExcelFile", "ExcelWriter"),
    ),
    # Parquet I/O is under pandas.io.parquet
    BucketRule(
        "io/parquet",
        prefixes=("pandas.io.parquet",),
        keywords=("read_parquet", "to_parquet"),
    ),
    # Other binary-ish formats commonly hit perf/mem
    BucketRule(
        "io/feather-orc",
        prefixes=("pandas.io.feather_format", "pandas.io.orc"),
        keywords=("read_feather", "to_feather", "read_orc", "to_orc"),
    ),
    # SQL I/O
    BucketRule(
        "io/sql",
        prefixes=("pandas.io.sql",),
        keywords=("read_sql", "read_sql_table", "read_sql_query", "to_sql"),
    ),
    # JSON/HTML/XML
    BucketRule(
        "io/json-html-xml",
        prefixes=("pandas.io.json.", "pandas.io.html", "pandas.io.xml."),
        keywords=(
            "read_json",
            "to_json",
            "json_normalize",
            "read_html",
            "to_html",
            "read_xml",
            "to_xml",
        ),
    ),
    # Stata/SAS/SPSS
    BucketRule(
        "io/stata-sas-spss",
        prefixes=("pandas.io.stata", "pandas.io.sas", "pandas.io.spss"),
        keywords=("read_stata", "read_sas", "read_spss"),
    ),
    # Pickle/clipboard
    BucketRule(
        "io/misc",
        prefixes=("pandas.io.pickle", "pandas.io.clipboards"),
        keywords=("read_pickle", "to_pickle", "read_clipboard"),
    ),
    # Catch-all for anything else in pandas.io.*
    BucketRule(
        "io/other",
        prefixes=("pandas.io.",),
        keywords=("read_", "to_"),
    ),
    # Merge/join is in pandas.core.reshape.merge
    BucketRule(
        "merge/join",
        prefixes=("pandas.core.reshape.merge",),
        keywords=("merge", "join", "merge_asof", "merge_ordered"),
    ),
    # Concat is in pandas.core.reshape.concat
    BucketRule(
        "concat",
        prefixes=("pandas.core.reshape.concat",),
        keywords=("concat",),
    ),
    # Reshape/pivot/melt/stack/unstack are spread across core.reshape.*
    BucketRule(
        "reshape/pivot-melt-stack",
        prefixes=(
            "pandas.core.reshape.melt",
            "pandas.core.reshape.reshape",
            "pandas.core.reshape.pivot",
            "pandas.core.reshape.tile",
        ),
        keywords=(
            "pivot",
            "pivot_table",
            "melt",
            "stack",
            "unstack",
            "wide_to_long",
            "lreshape",
            "crosstab",
            "get_dummies",
        ),
    ),
    # Groupby lives under pandas.core.groupby.*
    BucketRule(
        "groupby/agg",
        prefixes=("pandas.core.groupby.",),
        keywords=("groupby", "aggregate", "agg", "transform", "apply"),
    ),
    # Rolling/ewm/expanding live under pandas.core.window.*
    BucketRule(
        "window/rolling",
        prefixes=("pandas.core.window.",),
        keywords=("rolling", "ewm", "expanding", "window"),
    ),
    # Resample is currently a single module path pandas.core.resample
    BucketRule(
        "resample",
        prefixes=("pandas.core.resample",),
        keywords=("resample",),
    ),
    # Indexing & alignment is split across indexing/indexes/indexers
    BucketRule(
        "indexing/alignment",
        prefixes=(
            "pandas.core.indexing",
            "pandas.core.indexes.",
            "pandas.core.indexers.",
        ),
        keywords=("reindex", "align", "take", "loc", "iloc", "xs"),
    ),
    # Sorting/ranking
    BucketRule(
        "sort/rank",
        prefixes=("pandas.core.sorting",),
        keywords=("sort", "rank", "lexsort"),
    ),
    # Missing data paths (fill/drop/isna/interpolate)
    BucketRule(
        "missing-data",
        prefixes=("pandas.core.missing", "pandas.core.nanops"),
        keywords=("fillna", "dropna", "isna", "notna", "interpolate"),
    ),
    # Dtype conversion / inference / casting
    BucketRule(
        "astype/convert",
        prefixes=(
            "pandas.core.dtypes.",
            "pandas.core.arrays.",
            "pandas.core.tools.",
            "pandas.arrays.",
        ),
        keywords=(
            "astype",
            "to_datetime",
            "to_timedelta",
            "to_numeric",
            "convert_dtypes",
        ),
    ),
    # String accessor operations (often expensive)
    BucketRule(
        "strings",
        prefixes=("pandas.core.strings.",),
        keywords=("StringMethods", ".str"),
    ),
    # Internals / block manager / array manager
    BucketRule(
        "internals/blocks",
        prefixes=("pandas.core.internals.",),
        keywords=("Block", "BlockManager", "ArrayManager"),
    ),
    # Expression eval/query engine
    BucketRule(
        "computation/eval",
        prefixes=("pandas.core.computation.",),
        keywords=("eval", "query"),
    ),
    # Time series / offsets / frequencies
    BucketRule(
        "time-series",
        prefixes=("pandas.tseries.",),
        keywords=(
            "date_range",
            "bdate_range",
            "period_range",
            "timedelta_range",
            "infer_freq",
        ),
    ),
    # Cython/compiled hot paths
    BucketRule(
        "pandas._libs",
        prefixes=("pandas._libs",),
        keywords=("pandas._libs",),
    ),
]


def bucket_for(location: str, rules: Sequence[BucketRule] | None = None) -> str:
    rules_list = DEFAULT_RULES if rules is None else list(rules)
    filename, func = _split_location(location)

    mod = _file_to_module(filename)
    if not mod.startswith("pandas"):
        lf = (func or "").lower()
        if "pandas._libs" in lf:
            return "pandas._libs"
        return ""

    text = (mod + ":" + (func or "")).lower()

    for r in rules_list:
        if any(mod.startswith(p) for p in r.prefixes):
            return r.name

    for r in rules_list:
        if r.keywords and any(k in text for k in r.keywords):
            return r.name

    return DEFAULT_PHASE


__all__ = ["DEFAULT_PHASE", "RULES_VERSION", "bucket_for"]
