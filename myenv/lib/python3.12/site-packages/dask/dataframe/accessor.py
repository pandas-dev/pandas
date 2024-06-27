from __future__ import annotations

import functools
import warnings

import numpy as np
import pandas as pd

from dask.dataframe._compat import check_to_pydatetime_deprecation
from dask.utils import derived_from


def _bind_method(cls, pd_cls, attr, min_version=None):
    def func(self, *args, **kwargs):
        return self._function_map(attr, *args, **kwargs)

    func.__name__ = attr
    func.__qualname__ = f"{cls.__name__}.{attr}"
    try:
        func.__wrapped__ = getattr(pd_cls, attr)
    except Exception:
        pass
    setattr(cls, attr, derived_from(pd_cls, version=min_version)(func))


def _bind_property(cls, pd_cls, attr, min_version=None):
    def func(self):
        return self._property_map(attr)

    func.__name__ = attr
    func.__qualname__ = f"{cls.__name__}.{attr}"
    try:
        # Attempt to determine the method we are wrapping
        original_prop = getattr(pd_cls, attr)
        if isinstance(original_prop, property):
            method = original_prop.fget
        elif isinstance(original_prop, functools.cached_property):
            method = original_prop.func
        else:
            method = original_prop
            func.__wrapped__ = method
    except Exception:
        # If we can't then no matter, the function still works.
        pass
    setattr(cls, attr, property(derived_from(pd_cls, version=min_version)(func)))


def maybe_wrap_pandas(obj, x):
    if isinstance(x, np.ndarray):
        if isinstance(obj, pd.Series):
            return pd.Series(x, index=obj.index, dtype=x.dtype)
        return pd.Index(x)
    return x


class Accessor:
    """
    Base class for pandas Accessor objects cat, dt, and str.

    Notes
    -----
    Subclasses should define ``_accessor_name``, ``_accessor_methods``, and
    ``_accessor_properties``.
    """

    def __init__(self, series):
        from dask.dataframe.core import Series

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
    def _delegate_method(
        obj, accessor, attr, args, kwargs, catch_deprecation_warnings: bool = False
    ):
        with check_to_pydatetime_deprecation(catch_deprecation_warnings):
            with warnings.catch_warnings():
                # Falling back on a non-pyarrow code path which may decrease performance
                warnings.simplefilter("ignore", pd.errors.PerformanceWarning)
                out = getattr(getattr(obj, accessor, obj), attr)(*args, **kwargs)
                return maybe_wrap_pandas(obj, out)

    def _property_map(self, attr):
        meta = self._delegate_property(self._series._meta, self._accessor_name, attr)
        token = f"{self._accessor_name}-{attr}"
        return self._series.map_partitions(
            self._delegate_property, self._accessor_name, attr, token=token, meta=meta
        )

    def _function_map(self, attr, *args, **kwargs):
        if "meta" in kwargs:
            meta = kwargs.pop("meta")
        else:
            meta = self._delegate_method(
                self._series._meta_nonempty, self._accessor_name, attr, args, kwargs
            )
        token = f"{self._accessor_name}-{attr}"
        return self._series.map_partitions(
            self._delegate_method,
            self._accessor_name,
            attr,
            args,
            kwargs,
            catch_deprecation_warnings=True,
            meta=meta,
            token=token,
        )


class DatetimeAccessor(Accessor):
    """Accessor object for datetimelike properties of the Series values.

    Examples
    --------

    >>> s.dt.microsecond  # doctest: +SKIP
    """

    _accessor_name = "dt"

    _accessor_methods = (
        "asfreq",
        "ceil",
        "day_name",
        "floor",
        "month_name",
        "normalize",
        "round",
        "strftime",
        "to_period",
        "to_pydatetime",
        "to_pytimedelta",
        "to_timestamp",
        "total_seconds",
        "tz_convert",
        "tz_localize",
    )

    _accessor_properties = (
        "components",
        "date",
        "day",
        "day_of_week",
        "day_of_year",
        "dayofweek",
        "dayofyear",
        "days",
        "days_in_month",
        "daysinmonth",
        "end_time",
        "freq",
        "hour",
        "is_leap_year",
        "is_month_end",
        "is_month_start",
        "is_quarter_end",
        "is_quarter_start",
        "is_year_end",
        "is_year_start",
        "microsecond",
        "microseconds",
        "minute",
        "month",
        "nanosecond",
        "nanoseconds",
        "quarter",
        "qyear",
        "second",
        "seconds",
        "start_time",
        "time",
        "timetz",
        "tz",
        "week",
        "weekday",
        "weekofyear",
        "year",
    )

    @derived_from(pd.Series.dt)
    def isocalendar(self):
        # Sphinx can't solve types with dask-expr available so define explicitly, see
        # https://github.com/sphinx-doc/sphinx/issues/4961
        return self._function_map("isocalendar")


class StringAccessor(Accessor):
    """Accessor object for string properties of the Series values.

    Examples
    --------

    >>> s.str.lower()  # doctest: +SKIP
    """

    _accessor_name = "str"

    _accessor_methods = (
        "capitalize",
        "casefold",
        "center",
        "contains",
        "count",
        "decode",
        "encode",
        "find",
        "findall",
        "fullmatch",
        "get",
        "index",
        "isalnum",
        "isalpha",
        "isdecimal",
        "isdigit",
        "islower",
        "isnumeric",
        "isspace",
        "istitle",
        "isupper",
        "join",
        "len",
        "ljust",
        "lower",
        "lstrip",
        "match",
        "normalize",
        "pad",
        "partition",
        ("removeprefix", "1.4"),
        ("removesuffix", "1.4"),
        "repeat",
        "replace",
        "rfind",
        "rindex",
        "rjust",
        "rpartition",
        "rstrip",
        "slice",
        "slice_replace",
        "strip",
        "swapcase",
        "title",
        "translate",
        "upper",
        "wrap",
        "zfill",
    )
    _accessor_properties = ()

    def _split(self, method, pat=None, n=-1, expand=False):
        if expand:
            if n == -1:
                raise NotImplementedError(
                    "To use the expand parameter you must specify the number of "
                    "expected splits with the n= parameter. Usually n splits "
                    "result in n+1 output columns."
                )
            else:
                delimiter = " " if pat is None else pat
                meta = self._series._meta._constructor(
                    [delimiter.join(["a"] * (n + 1))],
                    index=self._series._meta_nonempty.iloc[:1].index,
                )
                meta = getattr(meta.str, method)(n=n, expand=expand, pat=pat)
        else:
            meta = (self._series.name, object)
        return self._function_map(method, pat=pat, n=n, expand=expand, meta=meta)

    @derived_from(
        pd.Series.str,
        inconsistencies="``expand=True`` with unknown ``n`` will raise a ``NotImplementedError``",
    )
    def split(self, pat=None, n=-1, expand=False):
        """Known inconsistencies: ``expand=True`` with unknown ``n`` will raise a ``NotImplementedError``."""
        return self._split("split", pat=pat, n=n, expand=expand)

    @derived_from(pd.Series.str)
    def rsplit(self, pat=None, n=-1, expand=False):
        return self._split("rsplit", pat=pat, n=n, expand=expand)

    @derived_from(pd.Series.str)
    def cat(self, others=None, sep=None, na_rep=None):
        from dask.dataframe.core import Index, Series

        if others is None:

            def str_cat_none(x):
                if isinstance(x, (Series, Index)):
                    x = x.compute()

                return x.str.cat(sep=sep, na_rep=na_rep)

            return self._series.reduction(chunk=str_cat_none, aggregate=str_cat_none)

        valid_types = (Series, Index, pd.Series, pd.Index)
        if isinstance(others, valid_types):
            others = [others]
        elif not all(isinstance(a, valid_types) for a in others):
            raise TypeError("others must be Series/Index")

        return self._series.map_partitions(
            str_cat, *others, sep=sep, na_rep=na_rep, meta=self._series._meta
        )

    @derived_from(pd.Series.str)
    def extractall(self, pat, flags=0):
        return self._series.map_partitions(
            str_extractall, pat, flags, token="str-extractall"
        )

    def __getitem__(self, index):
        return self._series.map_partitions(str_get, index, meta=self._series._meta)

    @derived_from(pd.Series.str)
    def extract(self, *args, **kwargs):
        # Sphinx can't solve types with dask-expr available so define explicitly, see
        # https://github.com/sphinx-doc/sphinx/issues/4961
        return self._function_map("extract", *args, **kwargs)

    @derived_from(pd.Series.str)
    def startswith(self, *args, **kwargs):
        # Sphinx can't solve types with dask-expr available so define explicitly, see
        # https://github.com/sphinx-doc/sphinx/issues/4961
        return self._function_map("startswith", *args, **kwargs)

    @derived_from(pd.Series.str)
    def endswith(self, *args, **kwargs):
        # Sphinx can't solve types with dask-expr available so define explicitly, see
        # https://github.com/sphinx-doc/sphinx/issues/4961
        return self._function_map("endswith", *args, **kwargs)


def str_extractall(series, pat, flags):
    return series.str.extractall(pat, flags=flags)


def str_get(series, index):
    """Implements series.str[index]"""
    return series.str[index]


def str_cat(self, *others, **kwargs):
    return self.str.cat(others=others, **kwargs)


# Ported from pandas
# https://github.com/pandas-dev/pandas/blob/master/pandas/core/accessor.py
class CachedAccessor:
    """
    Custom property-like object (descriptor) for caching accessors.

    Parameters
    ----------
    name : str
        The namespace this will be accessed under, e.g. ``df.foo``
    accessor : cls
        The class with the extension methods. The class' __init__ method
        should expect one of a ``Series``, ``DataFrame`` or ``Index`` as
        the single argument ``data``
    """

    def __init__(self, name, accessor):
        self._name = name
        self._accessor = accessor

    def __get__(self, obj, cls):
        if obj is None:
            # we're accessing the attribute of the class, i.e., Dataset.geo
            return self._accessor
        accessor_obj = self._accessor(obj)
        # Replace the property with the accessor object. Inspired by:
        # http://www.pydanny.com/cached-property.html
        # We need to use object.__setattr__ because we overwrite __setattr__ on
        # NDFrame
        object.__setattr__(obj, self._name, accessor_obj)
        return accessor_obj


def _register_accessor(name, cls):
    def decorator(accessor):
        if hasattr(cls, name):
            warnings.warn(
                "registration of accessor {!r} under name {!r} for type "
                "{!r} is overriding a preexisting attribute with the same "
                "name.".format(accessor, name, cls),
                UserWarning,
                stacklevel=2,
            )
        setattr(cls, name, CachedAccessor(name, accessor))
        cls._accessors.add(name)
        return accessor

    return decorator


def register_dataframe_accessor(name):
    """
    Register a custom accessor on :class:`dask.dataframe.DataFrame`.

    See :func:`pandas.api.extensions.register_dataframe_accessor` for more.
    """
    from dask.dataframe import DataFrame

    return _register_accessor(name, DataFrame)


def register_series_accessor(name):
    """
    Register a custom accessor on :class:`dask.dataframe.Series`.

    See :func:`pandas.api.extensions.register_series_accessor` for more.
    """
    from dask.dataframe import Series

    return _register_accessor(name, Series)


def register_index_accessor(name):
    """
    Register a custom accessor on :class:`dask.dataframe.Index`.

    See :func:`pandas.api.extensions.register_index_accessor` for more.
    """
    from dask.dataframe import Index

    return _register_accessor(name, Index)
