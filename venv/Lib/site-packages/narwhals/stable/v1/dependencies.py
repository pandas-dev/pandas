from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from typing_extensions import TypeIs

    from narwhals.stable.v1.typing import IntoDataFrameT, IntoLazyFrameT, IntoSeriesT


from narwhals.dependencies import (
    _is_cudf_dataframe as is_cudf_dataframe,
    _is_cudf_series as is_cudf_series,
    _is_dask_dataframe as is_dask_dataframe,
    _is_duckdb_relation as is_duckdb_relation,
    _is_ibis_table as is_ibis_table,
    _is_modin_dataframe as is_modin_dataframe,
    _is_modin_series as is_modin_series,
    _is_native_dataframe,
    _is_native_lazyframe,
    _is_native_series,
    _is_pandas_dataframe as is_pandas_dataframe,
    _is_pandas_like_dataframe as is_pandas_like_dataframe,
    _is_pandas_like_series as is_pandas_like_series,
    _is_pandas_series as is_pandas_series,
    _is_polars_dataframe as is_polars_dataframe,
    _is_polars_lazyframe as is_polars_lazyframe,
    _is_polars_series as is_polars_series,
    _is_pyarrow_chunked_array as is_pyarrow_chunked_array,
    _is_pyarrow_table as is_pyarrow_table,
    _is_pyspark_connect_dataframe as is_pyspark_connect_dataframe,
    _is_pyspark_dataframe as is_pyspark_dataframe,
    _is_sqlframe_dataframe as is_sqlframe_dataframe,
    get_cudf,
    get_dask,
    get_dask_dataframe,
    get_duckdb,
    get_ibis,
    get_modin,
    get_numpy,
    get_pandas,
    get_polars,
    get_pyarrow,
    get_pyspark,
    get_pyspark_connect,
    get_pyspark_sql,
    get_sqlframe,
    is_cudf_index,
    is_modin_index,
    is_narwhals_dataframe,
    is_narwhals_lazyframe,
    is_narwhals_series,
    is_numpy_array,
    is_pandas_index,
    is_pandas_like_index,
)


# TODO(Unassigned): For duckdb and ibis backends:
#   * is_into_dataframe(native_frame) returns False
#   * nw_v1.from_native(native_frame) returns a narwhals.stable.v1.DataFrame
#   * Therefore is_into_dataframe(nw_v1.from_native(native_frame)) returns True
# See discussion https://github.com/narwhals-dev/narwhals/pull/3613#discussion_r3288440039
def is_into_dataframe(native_dataframe: Any | IntoDataFrameT) -> TypeIs[IntoDataFrameT]:
    """Check whether `native_dataframe` can be converted to a narwhals.stable.v1.DataFrame."""
    from narwhals.stable.v1 import DataFrame

    return isinstance(native_dataframe, DataFrame) or _is_native_dataframe(
        native_dataframe
    )


def is_into_lazyframe(native_lazyframe: Any | IntoLazyFrameT) -> TypeIs[IntoLazyFrameT]:
    """Check whether `native_lazyframe` can be converted to a narwhals.stable.v1.LazyFrame."""
    from narwhals.stable.v1 import LazyFrame

    return isinstance(native_lazyframe, LazyFrame) or _is_native_lazyframe(
        native_lazyframe
    )


def is_into_series(native_series: Any | IntoSeriesT) -> TypeIs[IntoSeriesT]:
    """Check whether `native_series` can be converted to a narwhals.stable.v1.Series."""
    from narwhals.stable.v1 import Series

    return isinstance(native_series, Series) or _is_native_series(native_series)


__all__ = [
    "get_cudf",
    "get_dask",
    "get_dask_dataframe",
    "get_duckdb",
    "get_ibis",
    "get_modin",
    "get_numpy",
    "get_pandas",
    "get_polars",
    "get_pyarrow",
    "get_pyspark",
    "get_pyspark_connect",
    "get_pyspark_sql",
    "get_sqlframe",
    "is_cudf_dataframe",
    "is_cudf_index",
    "is_cudf_series",
    "is_dask_dataframe",
    "is_duckdb_relation",
    "is_ibis_table",
    "is_into_dataframe",
    "is_into_lazyframe",
    "is_into_series",
    "is_modin_dataframe",
    "is_modin_index",
    "is_modin_series",
    "is_narwhals_dataframe",
    "is_narwhals_lazyframe",
    "is_narwhals_series",
    "is_numpy_array",
    "is_pandas_dataframe",
    "is_pandas_index",
    "is_pandas_like_dataframe",
    "is_pandas_like_index",
    "is_pandas_like_series",
    "is_pandas_series",
    "is_polars_dataframe",
    "is_polars_lazyframe",
    "is_polars_series",
    "is_pyarrow_chunked_array",
    "is_pyarrow_table",
    "is_pyspark_connect_dataframe",
    "is_pyspark_dataframe",
    "is_sqlframe_dataframe",
]
