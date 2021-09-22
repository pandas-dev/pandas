# flake8: noqa

__docformat__ = "restructuredtext"

# Let users know if they're missing any of our hard dependencies
hard_dependencies = ("numpy", "pytz", "dateutil")
missing_dependencies = []

for dependency in hard_dependencies:
    try:
        __import__(dependency)
    except ImportError as e:
        missing_dependencies.append(f"{dependency}: {e}")

if missing_dependencies:
    raise ImportError(
        "Unable to import required dependencies:\n" + "\n".join(missing_dependencies)
    )
del hard_dependencies, dependency, missing_dependencies

# numpy compat
from pandas.compat import is_numpy_dev as _is_numpy_dev

try:
    from pandas._libs import hashtable as _hashtable, lib as _lib, tslib as _tslib
except ImportError as e:  # pragma: no cover
    module = e.name
    raise ImportError(
        f"C extension: {module} not built. If you want to import "
        "pandas from the source directory, you may need to run "
        "'python setup.py build_ext --force' to build the C extensions first."
    ) from e

# Use redundant imports (X as X) for type checkers to know what is part of the
# public API. Pandas is not (yet) a py.typed library: the public API is determined
# based on the documentation.
from pandas._config import (
    get_option as get_option,
    set_option as set_option,
    reset_option as reset_option,
    describe_option as describe_option,
    option_context as option_context,
    options as options,
)

# let init-time option registration happen
import pandas.core.config_init

from pandas.core.api import (
    # dtype
    Int8Dtype as Int8Dtype,
    Int16Dtype as Int16Dtype,
    Int32Dtype as Int32Dtype,
    Int64Dtype as Int64Dtype,
    UInt8Dtype as UInt8Dtype,
    UInt16Dtype as UInt16Dtype,
    UInt32Dtype as UInt32Dtype,
    UInt64Dtype as UInt64Dtype,
    Float32Dtype as Float32Dtype,
    Float64Dtype as Float64Dtype,
    CategoricalDtype as CategoricalDtype,
    PeriodDtype as PeriodDtype,
    IntervalDtype as IntervalDtype,
    DatetimeTZDtype as DatetimeTZDtype,
    StringDtype as StringDtype,
    BooleanDtype as BooleanDtype,
    # missing
    NA as NA,
    isna as isna,
    isnull as isnull,
    notna as notna,
    notnull as notnull,
    # indexes
    Index as Index,
    CategoricalIndex as CategoricalIndex,
    RangeIndex as RangeIndex,
    NumericIndex as NumericIndex,
    MultiIndex as MultiIndex,
    IntervalIndex as IntervalIndex,
    TimedeltaIndex as TimedeltaIndex,
    DatetimeIndex as DatetimeIndex,
    PeriodIndex as PeriodIndex,
    IndexSlice as IndexSlice,
    # tseries
    NaT as NaT,
    Period as Period,
    period_range as period_range,
    Timedelta as Timedelta,
    timedelta_range as timedelta_range,
    Timestamp as Timestamp,
    date_range as date_range,
    bdate_range as bdate_range,
    Interval as Interval,
    interval_range as interval_range,
    DateOffset as DateOffset,
    # conversion
    to_numeric as to_numeric,
    to_datetime as to_datetime,
    to_timedelta as to_timedelta,
    # misc
    Flags as Flags,
    Grouper as Grouper,
    factorize as factorize,
    unique as unique,
    value_counts as value_counts,
    NamedAgg as NamedAgg,
    array as array,
    Categorical as Categorical,
    set_eng_float_format as set_eng_float_format,
    Series as Series,
    DataFrame as DataFrame,
)

from pandas.core.arrays.sparse import SparseDtype as SparseDtype

from pandas.tseries.api import infer_freq as infer_freq
from pandas.tseries import offsets as offsets

from pandas.core.computation.api import eval as eval

from pandas.core.reshape.api import (
    concat as concat,
    lreshape as lreshape,
    melt as melt,
    wide_to_long as wide_to_long,
    merge as merge,
    merge_asof as merge_asof,
    merge_ordered as merge_ordered,
    crosstab as crosstab,
    pivot as pivot,
    pivot_table as pivot_table,
    get_dummies as get_dummies,
    cut as cut,
    qcut as qcut,
)

import pandas.api
from pandas.util._print_versions import show_versions as show_versions

from pandas.io.api import (
    # excel
    ExcelFile as ExcelFile,
    ExcelWriter as ExcelWriter,
    read_excel as read_excel,
    # parsers
    read_csv as read_csv,
    read_fwf as read_fwf,
    read_table as read_table,
    # pickle
    read_pickle as read_pickle,
    to_pickle as to_pickle,
    # pytables
    HDFStore as HDFStore,
    read_hdf as read_hdf,
    # sql
    read_sql as read_sql,
    read_sql_query as read_sql_query,
    read_sql_table as read_sql_table,
    # misc
    read_clipboard as read_clipboard,
    read_parquet as read_parquet,
    read_orc as read_orc,
    read_feather as read_feather,
    read_gbq as read_gbq,
    read_html as read_html,
    read_xml as read_xml,
    read_json as read_json,
    read_stata as read_stata,
    read_sas as read_sas,
    read_spss as read_spss,
)

from pandas.io.json import _json_normalize as json_normalize

from pandas.util._tester import test as test
import pandas.testing
import pandas.arrays

# use the closest tagged version if possible
from pandas._version import get_versions as get_versions

v = get_versions()
__version__ = v.get("closest-tag", v["version"])
__git_version__ = v.get("full-revisionid")
del get_versions, v


# GH 27101
__deprecated_num_index_names = ["Float64Index", "Int64Index", "UInt64Index"]


def __dir__():
    # GH43028
    # Int64Index etc. are deprecated, but we still want them to be available in the dir.
    # Remove in Pandas 2.0, when we remove Int64Index etc. from the code base.
    return list(globals().keys()) + __deprecated_num_index_names


def __getattr__(name):
    import warnings

    if name in __deprecated_num_index_names:
        warnings.warn(
            f"pandas.{name} is deprecated "
            "and will be removed from pandas in a future version. "
            "Use pandas.NumericIndex with the appropriate dtype instead.",
            FutureWarning,
            stacklevel=2,
        )
        from pandas.core.api import Float64Index, Int64Index, UInt64Index

        return {
            "Float64Index": Float64Index,
            "Int64Index": Int64Index,
            "UInt64Index": UInt64Index,
        }[name]
    elif name == "datetime":
        warnings.warn(
            "The pandas.datetime class is deprecated "
            "and will be removed from pandas in a future version. "
            "Import from datetime module instead.",
            FutureWarning,
            stacklevel=2,
        )

        from datetime import datetime as dt

        return dt

    elif name == "np":

        warnings.warn(
            "The pandas.np module is deprecated "
            "and will be removed from pandas in a future version. "
            "Import numpy directly instead.",
            FutureWarning,
            stacklevel=2,
        )
        import numpy as np

        return np

    elif name in {"SparseSeries", "SparseDataFrame"}:
        warnings.warn(
            f"The {name} class is removed from pandas. Accessing it from "
            "the top-level namespace will also be removed in the next version.",
            FutureWarning,
            stacklevel=2,
        )

        return type(name, (), {})

    elif name == "SparseArray":

        warnings.warn(
            "The pandas.SparseArray class is deprecated "
            "and will be removed from pandas in a future version. "
            "Use pandas.arrays.SparseArray instead.",
            FutureWarning,
            stacklevel=2,
        )
        from pandas.core.arrays.sparse import SparseArray as _SparseArray

        return _SparseArray

    raise AttributeError(f"module 'pandas' has no attribute '{name}'")


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
