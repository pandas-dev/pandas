from __future__ import annotations

import pandas as pd

from dask._dispatch import get_collection_type
from dask.backends import CreationDispatch
from dask.dataframe.backends import DataFrameBackendEntrypoint
from dask.dataframe.dask_expr._expr import ToBackend
from dask.dataframe.dispatch import to_pandas_dispatch

dataframe_creation_dispatch = CreationDispatch(
    module_name="dataframe",
    default="pandas",
    entrypoint_root="dask_expr",
    entrypoint_class=DataFrameBackendEntrypoint,
    name="dataframe_creation_dispatch",
)


class ToPandasBackend(ToBackend):
    @staticmethod
    def operation(df, options):
        return to_pandas_dispatch(df, **options)

    def _simplify_down(self):
        if isinstance(self.frame._meta, (pd.DataFrame, pd.Series, pd.Index)):
            # We already have pandas data
            return self.frame


class PandasBackendEntrypoint(DataFrameBackendEntrypoint):
    """Pandas-Backend Entrypoint Class for Dask-Expressions

    Note that all DataFrame-creation functions are defined
    and registered 'in-place'.
    """

    @classmethod
    def to_backend(cls, data, **kwargs):
        from dask.dataframe.dask_expr._collection import new_collection

        return new_collection(ToPandasBackend(data, kwargs))


dataframe_creation_dispatch.register_backend("pandas", PandasBackendEntrypoint())


@get_collection_type.register(pd.Series)
def get_collection_type_series(_):
    from dask.dataframe.dask_expr._collection import Series

    return Series


@get_collection_type.register(pd.DataFrame)
def get_collection_type_dataframe(_):
    from dask.dataframe.dask_expr._collection import DataFrame

    return DataFrame


@get_collection_type.register(pd.Index)
def get_collection_type_index(_):
    from dask.dataframe.dask_expr._collection import Index

    return Index


######################################
# cuDF: Pandas Dataframes on the GPU #
######################################


@get_collection_type.register_lazy("cudf")
def _register_cudf():
    import dask_cudf  # noqa: F401
