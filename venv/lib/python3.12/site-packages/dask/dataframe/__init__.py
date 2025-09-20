from __future__ import annotations


def _dask_expr_enabled() -> bool:
    import dask

    use_dask_expr = dask.config.get("dataframe.query-planning")
    if use_dask_expr is False:
        raise NotImplementedError("The legacy implementation is no longer supported")

    return True


_dask_expr_enabled()


import dask.array._array_expr._backends  # Import this to register array dispatch # noqa: F401

# Ensure that dtypes are registered
import dask.dataframe._dtypes
import dask.dataframe._pyarrow_compat
from dask._dispatch import get_collection_type
from dask.base import compute
from dask.dataframe import backends, dispatch
from dask.dataframe.dask_expr import (
    DataFrame,
    Index,
    Scalar,
    Series,
    concat,
    from_array,
    from_dask_array,
    from_delayed,
    from_dict,
    from_graph,
    from_map,
    from_pandas,
    get_dummies,
    isna,
    map_overlap,
    map_partitions,
    melt,
    merge,
    merge_asof,
    pivot_table,
    read_parquet,
    repartition,
    to_bag,
    to_datetime,
    to_numeric,
    to_parquet,
    to_records,
    to_timedelta,
)
from dask.dataframe.groupby import Aggregation
from dask.dataframe.io import (
    demo,
    read_csv,
    read_fwf,
    read_hdf,
    read_json,
    read_orc,
    read_sql,
    read_sql_query,
    read_sql_table,
    read_table,
    to_csv,
    to_hdf,
    to_json,
    to_orc,
    to_sql,
)
from dask.dataframe.utils import assert_eq

__all__ = [
    "DataFrame",
    "Index",
    "Scalar",
    "Series",
    "Aggregation",
    "backends",
    "compute",
    "concat",
    "dispatch",
    "from_array",
    "from_dask_array",
    "from_delayed",
    "from_dict",
    "from_graph",
    "from_map",
    "from_pandas",
    "get_collection_type",
    "get_dummies",
    "isna",
    "map_overlap",
    "map_partitions",
    "melt",
    "merge",
    "merge_asof",
    "pivot_table",
    "demo",
    "read_csv",
    "read_fwf",
    "read_hdf",
    "read_json",
    "read_orc",
    "read_parquet",
    "read_sql",
    "read_sql_query",
    "read_sql_table",
    "read_table",
    "repartition",
    "to_bag",
    "to_csv",
    "to_datetime",
    "to_hdf",
    "to_json",
    "to_numeric",
    "to_orc",
    "to_parquet",
    "to_records",
    "to_sql",
    "to_timedelta",
    "assert_eq",
]
