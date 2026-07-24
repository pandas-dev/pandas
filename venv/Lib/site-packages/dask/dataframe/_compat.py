from __future__ import annotations

from dask._compatibility import import_optional_dependency

import_optional_dependency("pandas")
import_optional_dependency("numpy")

from packaging.version import Version

from dask._pandas_compat import (
    PANDAS_GE_202,
    PANDAS_GE_210,
    PANDAS_GE_220,
    PANDAS_GE_230,
    PANDAS_GE_300,
    PANDAS_GE_310,
    PANDAS_VERSION,
    IndexingError,
    assert_categorical_equal,
    assert_numpy_array_equal,
    check_apply_dataframe_deprecation,
    check_convert_dtype_deprecation,
    check_groupby_axis_deprecation,
    check_observed_deprecation,
    is_any_real_numeric_dtype,
    is_string_dtype,
    makeDataFrame,
    makeDateIndex,
    makeMissingDataFrame,
    makeMixedDataFrame,
    makeTimeDataFrame,
    makeTimedeltaIndex,
    makeTimeSeries,
    tm,
)

# re-export all the things we need from dask._pandas_compat


try:
    import pyarrow as pa

    PYARROW_VERSION = Version(pa.__version__)
    HAS_PYARROW = PYARROW_VERSION >= Version("16.0.0")
except ImportError:
    PYARROW_VERSION = Version("0.0.0")
    HAS_PYARROW = False


__all__ = [
    "PANDAS_VERSION",
    "PANDAS_GE_202",
    "PANDAS_GE_210",
    "PANDAS_GE_220",
    "PANDAS_GE_230",
    "PANDAS_GE_300",
    "PANDAS_GE_310",
    "assert_categorical_equal",
    "assert_numpy_array_equal",
    "makeDataFrame",
    "makeTimeDataFrame",
    "makeTimeSeries",
    "makeDateIndex",
    "makeTimedeltaIndex",
    "makeMissingDataFrame",
    "makeMixedDataFrame",
    "check_groupby_axis_deprecation",
    "check_observed_deprecation",
    "check_convert_dtype_deprecation",
    "check_apply_dataframe_deprecation",
    "is_any_real_numeric_dtype",
    "is_string_dtype",
    "IndexingError",
    "HAS_PYARROW",
    "PYARROW_VERSION",
    "tm",
]
