from __future__ import annotations

import functools

from dask.dataframe.accessor import _bind_method, _bind_property, maybe_wrap_pandas
from dask.dataframe.dask_expr._expr import Elemwise, Expr
from dask.dataframe.dispatch import make_meta, meta_nonempty


class Accessor:
    """
    Base class for pandas Accessor objects cat, dt, and str.

    Notes
    -----
    Subclasses should define ``_accessor_name``, ``_accessor_methods``, and
    ``_accessor_properties``.
    """

    def __init__(self, series):
        from dask.dataframe.dask_expr import Series

        if not isinstance(series, Series):
            raise ValueError("Accessor cannot be initialized")

        series_meta = series._meta
        if hasattr(series_meta, "to_series"):  # is index-like
            series_meta = series_meta.to_series()
        meta = getattr(series_meta, self._accessor_name)

        self._meta = meta
        self._series = series

    def __init_subclass__(cls, **kwargs):
        """Bind all auto-generated methods & properties"""
        import pandas as pd

        super().__init_subclass__(**kwargs)
        pd_cls = getattr(pd.Series, cls._accessor_name)
        for item in cls._accessor_methods:
            attr, min_version = item if isinstance(item, tuple) else (item, None)
            if not hasattr(cls, attr):
                _bind_method(cls, pd_cls, attr, min_version)
        for item in cls._accessor_properties:
            attr, min_version = item if isinstance(item, tuple) else (item, None)
            if not hasattr(cls, attr):
                _bind_property(cls, pd_cls, attr, min_version)

    @staticmethod
    def _delegate_property(obj, accessor, attr):
        out = getattr(getattr(obj, accessor, obj), attr)
        return maybe_wrap_pandas(obj, out)

    @staticmethod
    def _delegate_method(obj, accessor, attr, args, kwargs):
        out = getattr(getattr(obj, accessor, obj), attr)(*args, **kwargs)
        return maybe_wrap_pandas(obj, out)

    def _function_map(self, attr, *args, **kwargs):
        from dask.dataframe.dask_expr._collection import Index, new_collection

        if isinstance(self._series, Index):
            return new_collection(
                FunctionMapIndex(self._series, self._accessor_name, attr, args, kwargs)
            )

        return new_collection(
            FunctionMap(self._series, self._accessor_name, attr, args, kwargs)
        )

    def _property_map(self, attr, *args, **kwargs):
        from dask.dataframe.dask_expr._collection import Index, new_collection

        if isinstance(self._series, Index):
            return new_collection(
                PropertyMapIndex(self._series, self._accessor_name, attr)
            )

        return new_collection(PropertyMap(self._series, self._accessor_name, attr))


class PropertyMap(Elemwise):
    _parameters = [
        "frame",
        "accessor",
        "attr",
    ]

    @staticmethod
    def operation(obj, accessor, attr):
        out = getattr(getattr(obj, accessor, obj), attr)
        return maybe_wrap_pandas(obj, out)


class PropertyMapIndex(PropertyMap):
    def _divisions(self):
        # TODO: We can do better here
        return (None,) * (self.frame.npartitions + 1)


class FunctionMap(Elemwise):
    _parameters = ["frame", "accessor", "attr", "args", "kwargs"]

    @functools.cached_property
    def _meta(self):
        args = [
            meta_nonempty(op._meta) if isinstance(op, Expr) else op for op in self._args
        ]
        return make_meta(self.operation(*args, **self._kwargs))

    @staticmethod
    def operation(obj, accessor, attr, args, kwargs):
        out = getattr(getattr(obj, accessor, obj), attr)(*args, **kwargs)
        return maybe_wrap_pandas(obj, out)


class FunctionMapIndex(FunctionMap):
    def _divisions(self):
        # TODO: We can do better here
        return (None,) * (self.frame.npartitions + 1)
