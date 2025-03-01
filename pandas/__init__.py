from __future__ import annotations

__docformat__ = "restructuredtext"

# Let users know if they're missing any of our hard dependencies
_hard_dependencies = ("numpy", "dateutil")
_missing_dependencies = []

for _dependency in _hard_dependencies:
    try:
        __import__(_dependency)
    except ImportError as _e:  # pragma: no cover
        _missing_dependencies.append(f"{_dependency}: {_e}")

if _missing_dependencies:  # pragma: no cover
    raise ImportError(
        "Unable to import required dependencies:\n" + "\n".join(_missing_dependencies)
    )
del _hard_dependencies, _dependency, _missing_dependencies

try:
    # numpy compat
    from pandas.compat import (  # pyright: ignore[reportUnusedImport] # noqa: F401
        is_numpy_dev as _is_numpy_dev,
    )
except ImportError as _err:  # pragma: no cover
    _module = _err.name
    raise ImportError(
        f"C extension: {_module} not built. If you want to import "
        "pandas from the source directory, you may need to run "
        "'python -m pip install -ve . --no-build-isolation -Ceditable-verbose=true' "
        "to build the C extensions first."
    ) from _err

from pandas._config import (
    describe_option,
    get_option,
    option_context,
    options,
    reset_option,
    set_option,
)

from pandas.util._print_versions import show_versions
from pandas.util._tester import test

from pandas.core.dtypes.dtypes import SparseDtype

from pandas import (
    api,
    arrays,
    errors,
    io,
    plotting,
    testing,
    tseries,
)
from pandas.core.api import (  # dtype; missing; indexes; tseries; conversion; misc
    NA,
    ArrowDtype,
    BooleanDtype,
    Categorical,
    CategoricalDtype,
    CategoricalIndex,
    DataFrame,
    DateOffset,
    DatetimeIndex,
    DatetimeTZDtype,
    Flags,
    Float32Dtype,
    Float64Dtype,
    Grouper,
    Index,
    IndexSlice,
    Int8Dtype,
    Int16Dtype,
    Int32Dtype,
    Int64Dtype,
    Interval,
    IntervalDtype,
    IntervalIndex,
    MultiIndex,
    NamedAgg,
    NaT,
    Period,
    PeriodDtype,
    PeriodIndex,
    RangeIndex,
    Series,
    StringDtype,
    Timedelta,
    TimedeltaIndex,
    Timestamp,
    UInt8Dtype,
    UInt16Dtype,
    UInt32Dtype,
    UInt64Dtype,
    array,
    bdate_range,
    date_range,
    factorize,
    interval_range,
    isna,
    isnull,
    notna,
    notnull,
    period_range,
    set_eng_float_format,
    timedelta_range,
    to_datetime,
    to_numeric,
    to_timedelta,
    unique,
)
from pandas.core.computation.api import eval

# let init-time option registration happen
import pandas.core.config_init  # pyright: ignore[reportUnusedImport] # noqa: F401
from pandas.core.reshape.api import (
    concat,
    crosstab,
    cut,
    from_dummies,
    get_dummies,
    lreshape,
    melt,
    merge,
    merge_asof,
    merge_ordered,
    pivot,
    pivot_table,
    qcut,
    wide_to_long,
)

from pandas.io.api import (  # excel; parsers; pickle; pytables; sql; misc
    ExcelFile,
    ExcelWriter,
    HDFStore,
    read_clipboard,
    read_csv,
    read_excel,
    read_feather,
    read_fwf,
    read_hdf,
    read_html,
    read_json,
    read_orc,
    read_parquet,
    read_pickle,
    read_sas,
    read_spss,
    read_sql,
    read_sql_query,
    read_sql_table,
    read_stata,
    read_table,
    read_xml,
    to_pickle,
)
from pandas.io.json._normalize import json_normalize
from pandas.tseries import offsets
from pandas.tseries.api import infer_freq

# use the closest tagged version if possible
_built_with_meson = False
try:
    from pandas._version_meson import (  # pyright: ignore [reportMissingImports]
        __git_version__,
        __version__,
    )

    _built_with_meson = True
except ImportError:
    from pandas._version import get_versions

    v = get_versions()
    __version__ = v.get("closest-tag", v["version"])
    __git_version__ = v.get("full-revisionid")
    del get_versions, v


# module level doc-string
__doc__ = """
pandas - a powerful data analysis and manipulation library for Python
=====================================================================

**pandas** is a Python package providing fast, flexible, and expressive data
structures designed to make working with "relational" or "labeled" data both
easy and intuitive. It aims to be the fundamental high-level building block for
doing practical, **real world** data analysis in Python. Additionally, it has
the broader goal of becoming **the most powerful and flexible open source data
analysis / manipulation tool available in any language**. It is already well on
its way toward this goal.

Main Features
-------------
Here are just a few of the things that pandas does well:

  - Easy handling of missing data in floating point as well as non-floating
    point data.
  - Size mutability: columns can be inserted and deleted from DataFrame and
    higher dimensional objects
  - Automatic and explicit data alignment: objects can be explicitly aligned
    to a set of labels, or the user can simply ignore the labels and let
    `Series`, `DataFrame`, etc. automatically align the data for you in
    computations.
  - Powerful, flexible group by functionality to perform split-apply-combine
    operations on data sets, for both aggregating and transforming data.
  - Make it easy to convert ragged, differently-indexed data in other Python
    and NumPy data structures into DataFrame objects.
  - Intelligent label-based slicing, fancy indexing, and subsetting of large
    data sets.
  - Intuitive merging and joining data sets.
  - Flexible reshaping and pivoting of data sets.
  - Hierarchical labeling of axes (possible to have multiple labels per tick).
  - Robust IO tools for loading data from flat files (CSV and delimited),
    Excel files, databases, and saving/loading data from the ultrafast HDF5
    format.
  - Time series-specific functionality: date range generation and frequency
    conversion, moving window statistics, date shifting and lagging.
"""

# Use __all__ to let type checkers know what is part of the public API.
# Pandas is not (yet) a py.typed library: the public API is determined
# based on the documentation.
__all__ = [
    "NA",
    "ArrowDtype",
    "BooleanDtype",
    "Categorical",
    "CategoricalDtype",
    "CategoricalIndex",
    "DataFrame",
    "DateOffset",
    "DatetimeIndex",
    "DatetimeTZDtype",
    "ExcelFile",
    "ExcelWriter",
    "Flags",
    "Float32Dtype",
    "Float64Dtype",
    "Grouper",
    "HDFStore",
    "Index",
    "IndexSlice",
    "Int8Dtype",
    "Int16Dtype",
    "Int32Dtype",
    "Int64Dtype",
    "Interval",
    "IntervalDtype",
    "IntervalIndex",
    "MultiIndex",
    "NaT",
    "NamedAgg",
    "Period",
    "PeriodDtype",
    "PeriodIndex",
    "RangeIndex",
    "Series",
    "SparseDtype",
    "StringDtype",
    "Timedelta",
    "TimedeltaIndex",
    "Timestamp",
    "UInt8Dtype",
    "UInt16Dtype",
    "UInt32Dtype",
    "UInt64Dtype",
    "api",
    "array",
    "arrays",
    "bdate_range",
    "concat",
    "crosstab",
    "cut",
    "date_range",
    "describe_option",
    "errors",
    "eval",
    "factorize",
    "from_dummies",
    "get_dummies",
    "get_option",
    "infer_freq",
    "interval_range",
    "io",
    "isna",
    "isnull",
    "json_normalize",
    "lreshape",
    "melt",
    "merge",
    "merge_asof",
    "merge_ordered",
    "notna",
    "notnull",
    "offsets",
    "option_context",
    "options",
    "period_range",
    "pivot",
    "pivot_table",
    "plotting",
    "qcut",
    "read_clipboard",
    "read_csv",
    "read_excel",
    "read_feather",
    "read_fwf",
    "read_hdf",
    "read_html",
    "read_json",
    "read_orc",
    "read_parquet",
    "read_pickle",
    "read_sas",
    "read_spss",
    "read_sql",
    "read_sql_query",
    "read_sql_table",
    "read_stata",
    "read_table",
    "read_xml",
    "reset_option",
    "set_eng_float_format",
    "set_option",
    "show_versions",
    "test",
    "testing",
    "timedelta_range",
    "to_datetime",
    "to_numeric",
    "to_pickle",
    "to_timedelta",
    "tseries",
    "unique",
    "wide_to_long",
]
