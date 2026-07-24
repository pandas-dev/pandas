from __future__ import annotations

from typing import TYPE_CHECKING, Any

from narwhals.dependencies import (
    _is_native_dataframe,
    _is_native_lazyframe,
    _is_native_series,
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
    is_cudf_dataframe,
    is_cudf_index,
    is_cudf_series,
    is_dask_dataframe,
    is_duckdb_relation,
    is_ibis_table,
    is_modin_dataframe,
    is_modin_index,
    is_modin_series,
    is_narwhals_dataframe,
    is_narwhals_lazyframe,
    is_narwhals_series,
    is_numpy_array,
    is_pandas_dataframe,
    is_pandas_index,
    is_pandas_like_dataframe,
    is_pandas_like_index,
    is_pandas_like_series,
    is_pandas_series,
    is_polars_dataframe,
    is_polars_lazyframe,
    is_polars_series,
    is_pyarrow_chunked_array,
    is_pyarrow_table,
    is_pyspark_connect_dataframe,
    is_pyspark_dataframe,
    is_sqlframe_dataframe,
)

if TYPE_CHECKING:
    from typing_extensions import TypeIs

    from narwhals.stable.v2.typing import IntoDataFrameT, IntoLazyFrameT, IntoSeriesT


def is_into_dataframe(native_dataframe: Any | IntoDataFrameT) -> TypeIs[IntoDataFrameT]:
    """Check whether `native_dataframe` can be converted to a narwhals.stable.v2.DataFrame."""
    from narwhals.stable.v2 import DataFrame

    return isinstance(native_dataframe, DataFrame) or _is_native_dataframe(
        native_dataframe
    )


def is_into_lazyframe(native_lazyframe: Any | IntoLazyFrameT) -> TypeIs[IntoLazyFrameT]:
    """Check whether `native_lazyframe` can be converted to a narwhals.stable.v2.LazyFrame."""
    from narwhals.stable.v2 import LazyFrame

    return isinstance(native_lazyframe, LazyFrame) or _is_native_lazyframe(
        native_lazyframe
    )


def is_into_series(native_series: Any | IntoSeriesT) -> TypeIs[IntoSeriesT]:
    """Check whether `native_series` can be converted to a narwhals.stable.v2.Series."""
    from narwhals.stable.v2 import Series

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
