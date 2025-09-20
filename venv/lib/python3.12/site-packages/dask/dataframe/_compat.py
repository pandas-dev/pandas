from __future__ import annotations

from dask._compatibility import import_optional_dependency

import_optional_dependency("pandas")
import_optional_dependency("numpy")
# import_optional_dependency("pyarrow")

from packaging.version import Version

from dask._pandas_compat import (
    PANDAS_GE_201,
    PANDAS_GE_202,
    PANDAS_GE_210,
    PANDAS_GE_211,
    PANDAS_GE_220,
    PANDAS_GE_230,
    PANDAS_GE_300,
    PANDAS_VERSION,
    IndexingError,
    assert_categorical_equal,
    assert_numpy_array_equal,
    check_apply_dataframe_deprecation,
    check_convert_dtype_deprecation,
    check_groupby_axis_deprecation,
    check_observed_deprecation,
    check_reductions_runtime_warning,
    is_any_real_numeric_dtype,
    is_string_dtype,
    makeDataFrame,
    makeDateIndex,
    makeMissingDataframe,
    makeMixedDataFrame,
    makeTimeDataFrame,
    makeTimedeltaIndex,
    makeTimeSeries,
    tm,
)

# re-export all the things we need from dask._pandas_compat


try:
    import pyarrow as pa
except ImportError:
    HAS_PYARROW = False
else:
    HAS_PYARROW = True


if HAS_PYARROW:
    PYARROW_VERSION: Version | None = Version(pa.__version__)
    # we know that Version should be non-None when PYARROW_VERSION is True.
    PYARROW_GE_1500: bool | None = PYARROW_VERSION.release >= (15, 0, 0)  # type: ignore[union-attr]
    PYARROW_GE_2101: bool | None = PYARROW_VERSION.release >= (21, 0, 1)  # type: ignore[union-attr]
else:
    PYARROW_VERSION = None
    PYARROW_GE_1500 = None
    PYARROW_GE_2101 = None


__all__ = [
    "PANDAS_VERSION",
    "PANDAS_GE_201",
    "PANDAS_GE_202",
    "PANDAS_GE_210",
    "PANDAS_GE_211",
    "PANDAS_GE_220",
    "PANDAS_GE_230",
    "PANDAS_GE_300",
    "assert_categorical_equal",
    "assert_numpy_array_equal",
    "makeDataFrame",
    "makeTimeDataFrame",
    "makeTimeSeries",
    "makeDateIndex",
    "makeTimedeltaIndex",
    "makeMissingDataframe",
    "makeMixedDataFrame",
    "check_groupby_axis_deprecation",
    "check_observed_deprecation",
    "check_convert_dtype_deprecation",
    "check_apply_dataframe_deprecation",
    "check_reductions_runtime_warning",
    "is_any_real_numeric_dtype",
    "is_string_dtype",
    "IndexingError",
    "HAS_PYARROW",
    "PYARROW_VERSION",
    "PYARROW_GE_1500",
    "PYARROW_GE_2101",
    "tm",
]
