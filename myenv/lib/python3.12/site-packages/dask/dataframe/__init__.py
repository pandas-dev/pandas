from __future__ import annotations

import importlib
import warnings

from packaging.version import Version

# The "dataframe.query-planning" config can only be processed once
DASK_EXPR_ENABLED: bool | None = None


def _dask_expr_enabled() -> bool:
    import pandas as pd

    import dask

    global DASK_EXPR_ENABLED

    use_dask_expr = dask.config.get("dataframe.query-planning")
    if DASK_EXPR_ENABLED is not None:
        if (use_dask_expr is True and DASK_EXPR_ENABLED is False) or (
            use_dask_expr is False and DASK_EXPR_ENABLED is True
        ):
            warnings.warn(
                "The 'dataframe.query-planning' config is now set to "
                f"{use_dask_expr}, but query planning is already "
                f"{'enabled' if DASK_EXPR_ENABLED else 'disabled'}. "
                "The query-planning config can only be changed before "
                "`dask.dataframe` is first imported!"
            )
        return DASK_EXPR_ENABLED

    if (
        use_dask_expr is False
        or use_dask_expr is None
        and Version(pd.__version__).major < 2
    ):
        return (DASK_EXPR_ENABLED := False)
    try:
        import dask_expr  # noqa: F401
    except ImportError:
        msg = """
Dask dataframe query planning is disabled because dask-expr is not installed.

You can install it with `pip install dask[dataframe]` or `conda install dask`.
This will raise in a future version.
"""
        if use_dask_expr is None:
            warnings.warn(msg, FutureWarning)
            return (DASK_EXPR_ENABLED := False)
        else:
            raise ImportError(msg)
    return (DASK_EXPR_ENABLED := True)


try:
    import dask.dataframe._pyarrow_compat
    from dask.base import compute
    from dask.dataframe import backends, dispatch, methods, rolling
    from dask.dataframe._testing import test_dataframe
    from dask.dataframe.core import (
        DataFrame,
        Index,
        Series,
        _Frame,
        map_partitions,
        repartition,
        to_datetime,
        to_timedelta,
    )
    from dask.dataframe.groupby import Aggregation
    from dask.dataframe.io import (
        demo,
        from_array,
        from_dask_array,
        from_delayed,
        from_dict,
        from_map,
        from_pandas,
        read_csv,
        read_fwf,
        read_hdf,
        read_json,
        read_sql,
        read_sql_query,
        read_sql_table,
        read_table,
        to_bag,
        to_csv,
        to_hdf,
        to_json,
        to_records,
        to_sql,
    )
    from dask.dataframe.multi import concat, merge, merge_asof
    from dask.dataframe.numeric import to_numeric
    from dask.dataframe.optimize import optimize
    from dask.dataframe.reshape import get_dummies, melt, pivot_table
    from dask.dataframe.rolling import map_overlap
    from dask.dataframe.utils import assert_eq

    try:
        from dask.dataframe.io import read_parquet, to_parquet
    except ImportError:
        pass
    try:
        from dask.dataframe.io import read_orc, to_orc
    except ImportError:
        pass
    try:
        from dask.dataframe.core import isna
    except ImportError:
        pass

    if _dask_expr_enabled():
        import dask_expr as dd

        # trigger loading of dask-expr which will in-turn import dask.dataframe and run remainder
        # of this module's init updating attributes to be dask-expr
        # note: needs reload, in case dask-expr imported before dask.dataframe; works fine otherwise
        dd = importlib.reload(dd)
except ImportError as e:
    msg = (
        "Dask dataframe requirements are not installed.\n\n"
        "Please either conda or pip install as follows:\n\n"
        "  conda install dask                     # either conda install\n"
        '  python -m pip install "dask[dataframe]" --upgrade  # or python -m pip install'
    )
    raise ImportError(msg) from e


if _dask_expr_enabled():
    try:
        from dask_expr import (  # type: ignore
            DataFrame,
            Index,
            Series,
            concat,
            from_array,
            from_dask_array,
            from_dask_dataframe,
            from_delayed,
            from_dict,
            from_graph,
            from_legacy_dataframe,
            from_map,
            from_pandas,
            get_collection_type,
            get_dummies,
            isna,
            map_overlap,
            map_partitions,
            melt,
            merge,
            merge_asof,
            pivot_table,
            read_csv,
            read_fwf,
            read_hdf,
            read_json,
            read_orc,
            read_parquet,
            read_sql,
            read_sql_query,
            read_sql_table,
            read_table,
            repartition,
            to_bag,
            to_csv,
            to_datetime,
            to_hdf,
            to_json,
            to_numeric,
            to_orc,
            to_parquet,
            to_records,
            to_sql,
            to_timedelta,
        )

        import dask.dataframe._pyarrow_compat
        from dask.base import compute
        from dask.dataframe import backends, dispatch
        from dask.dataframe.groupby import Aggregation
        from dask.dataframe.io import demo
        from dask.dataframe.utils import assert_eq

        def raise_not_implemented_error(attr_name):
            def inner_func(*args, **kwargs):
                raise NotImplementedError(
                    f"Function {attr_name} is not implemented for dask-expr."
                )

            return inner_func

        _Frame = raise_not_implemented_error("_Frame")  # type: ignore

    # Due to the natural circular imports caused from dask-expr
    # wanting to import things from dask.dataframe, this module's init
    # can be run multiple times as it walks code trying to import
    # dask-expr while dask-expr is also trying to import from dask.dataframe
    # Each time this happens and hits a circular import, we can reload
    # dask.dataframe to update itself until dask-expr is fully initialized.
    # TODO: This can go away when dask-expr is merged into dask
    except ImportError:
        import dask.dataframe as dd

        dd = importlib.reload(dd)
