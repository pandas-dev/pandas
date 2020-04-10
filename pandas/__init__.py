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

from pandas._config import (
    describe_option,
    get_option,
    option_context,
    options,
    reset_option,
    set_option,
)

# numpy compat
from pandas.compat.numpy import (
    _is_numpy_dev,
    _np_version_under1p14,
    _np_version_under1p15,
    _np_version_under1p16,
    _np_version_under1p17,
    _np_version_under1p18,
)
from pandas.util._print_versions import show_versions
from pandas.util._tester import test

# let init-time option registration happen
import pandas.api
import pandas.arrays
from pandas.core.api import (  # dtype; missing; indexes; tseries; conversion; misc
    NA,
    BooleanDtype,
    Categorical,
    CategoricalDtype,
    CategoricalIndex,
    DataFrame,
    DateOffset,
    DatetimeIndex,
    DatetimeTZDtype,
    Float64Index,
    Grouper,
    Index,
    IndexSlice,
    Int8Dtype,
    Int16Dtype,
    Int32Dtype,
    Int64Dtype,
    Int64Index,
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
    UInt64Index,
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
    value_counts,
)
from pandas.core.arrays.sparse import SparseDtype
from pandas.core.computation.api import eval
import pandas.core.config_init
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
import pandas.testing

from pandas.io.api import (  # excel; parsers; pickle; pytables; sql; misc
    ExcelFile,
    ExcelWriter,
    HDFStore,
    read_clipboard,
    read_csv,
    read_excel,
    read_feather,
    read_fwf,
    read_gbq,
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
    to_pickle,
)
from pandas.io.json import _json_normalize as json_normalize
from pandas.tseries import offsets
from pandas.tseries.api import infer_freq

# use the closest tagged version if possible
from ._version import get_versions

try:
    from pandas._libs import hashtable as _hashtable, lib as _lib, tslib as _tslib
except ImportError as e:  # pragma: no cover
    # hack but overkill to use re
    module = str(e).replace("cannot import name ", "")
    raise ImportError(
        f"C extension: {module} not built. If you want to import "
        "pandas from the source directory, you may need to run "
        "'python setup.py build_ext --inplace --force' to build the C extensions first."
    ) from e













v = get_versions()
__version__ = v.get("closest-tag", v["version"])
__git_version__ = v.get("full-revisionid")
del get_versions, v

# GH 27101
# TODO: remove Panel compat in 1.0
if pandas.compat.PY37:

    def __getattr__(name):
        import warnings

        if name == "Panel":

            warnings.warn(
                "The Panel class is removed from pandas. Accessing it "
                "from the top-level namespace will also be removed in the next version",
                FutureWarning,
                stacklevel=2,
            )

            class Panel:
                pass

            return Panel

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
                "Import numpy directly instead",
                FutureWarning,
                stacklevel=2,
            )
            import numpy as np

            return np

        elif name in {"SparseSeries", "SparseDataFrame"}:
            warnings.warn(
                f"The {name} class is removed from pandas. Accessing it from "
                "the top-level namespace will also be removed in the next version",
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


else:

    class Panel:
        pass

    class SparseDataFrame:
        pass

    class SparseSeries:
        pass

    class __numpy:
        def __init__(self):
            import numpy as np
            import warnings

            self.np = np
            self.warnings = warnings

        def __getattr__(self, item):
            self.warnings.warn(
                "The pandas.np module is deprecated "
                "and will be removed from pandas in a future version. "
                "Import numpy directly instead",
                FutureWarning,
                stacklevel=2,
            )

            try:
                return getattr(self.np, item)
            except AttributeError as err:
                raise AttributeError(f"module numpy has no attribute {item}") from err

    np = __numpy()

    class __Datetime(type):

        from datetime import datetime as dt

        datetime = dt

        def __getattr__(cls, item):
            cls.emit_warning()

            try:
                return getattr(cls.datetime, item)
            except AttributeError as err:
                raise AttributeError(
                    f"module datetime has no attribute {item}"
                ) from err

        def __instancecheck__(cls, other):
            return isinstance(other, cls.datetime)

    class __DatetimeSub(metaclass=__Datetime):
        def emit_warning(dummy=0):
            import warnings

            warnings.warn(
                "The pandas.datetime class is deprecated "
                "and will be removed from pandas in a future version. "
                "Import from datetime instead.",
                FutureWarning,
                stacklevel=3,
            )

        def __new__(cls, *args, **kwargs):
            cls.emit_warning()
            from datetime import datetime as dt

            return dt(*args, **kwargs)

    datetime = __DatetimeSub

    class __SparseArray(type):

        from pandas.core.arrays.sparse import SparseArray as sa

        SparseArray = sa

        def __instancecheck__(cls, other):
            return isinstance(other, cls.SparseArray)

    class __SparseArraySub(metaclass=__SparseArray):
        def emit_warning(dummy=0):
            import warnings

            warnings.warn(
                "The pandas.SparseArray class is deprecated "
                "and will be removed from pandas in a future version. "
                "Use pandas.arrays.SparseArray instead.",
                FutureWarning,
                stacklevel=3,
            )

        def __new__(cls, *args, **kwargs):
            cls.emit_warning()
            from pandas.core.arrays.sparse import SparseArray as sa

            return sa(*args, **kwargs)

    SparseArray = __SparseArraySub


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
