# Only tests that raise an error and have no better location should go here.
# Tests for specific groupby methods should go in their respective
# test file.

import datetime

import pytest

from pandas import DataFrame
from pandas.tests.groupby import get_groupby_method_args


@pytest.mark.parametrize("how", ["method", "agg", "transform"])
def test_groupby_raises_string(how, groupby_func, as_index, sort):
    df = DataFrame(
        {
            "a": [1, 1, 1, 2, 2],
            "b": range(5),
            "c": list("xyzwt"),
        }
    )
    args = get_groupby_method_args(groupby_func, df)
    gb = df.groupby("a", as_index=as_index, sort=sort)

    klass, msg = {
        "all": (None, ""),
        "any": (None, ""),
        "bfill": (None, ""),
        "corrwith": (TypeError, "Could not convert"),
        "count": (None, ""),
        "cumcount": (None, ""),
        "cummax": (NotImplementedError, "function is not implemented for this dtype"),
        "cummin": (NotImplementedError, "function is not implemented for this dtype"),
        "cumprod": (NotImplementedError, "function is not implemented for this dtype"),
        "cumsum": (NotImplementedError, "function is not implemented for this dtype"),
        "diff": (TypeError, "unsupported operand type"),
        "ffill": (None, ""),
        "fillna": (None, ""),
        "first": (None, ""),
        "idxmax": (TypeError, "'argmax' not allowed for this dtype"),
        "idxmin": (TypeError, "'argmin' not allowed for this dtype"),
        "last": (None, ""),
        "max": (None, ""),
        "mean": (TypeError, "Could not convert xyz to numeric"),
        "median": (TypeError, "could not convert string to float"),
        "min": (None, ""),
        "ngroup": (None, ""),
        "nunique": (None, ""),
        "pct_change": (TypeError, "unsupported operand type"),
        "prod": (TypeError, "can't multiply sequence by non-int of type 'str'"),
        "quantile": (TypeError, "cannot be performed against 'object' dtypes!"),
        "rank": (None, ""),
        "sem": (ValueError, "could not convert string to float"),
        "shift": (None, ""),
        "size": (None, ""),
        "skew": (TypeError, "could not convert string to float"),
        "std": (ValueError, "could not convert string to float"),
        "sum": (None, ""),
        "var": (TypeError, "could not convert string to float"),
    }[groupby_func]

    if klass is None:
        if how == "method":
            getattr(gb, groupby_func)(*args)
        elif how == "agg":
            gb.agg(groupby_func, *args)
        else:
            gb.transform(groupby_func, *args)
    else:
        with pytest.raises(klass, match=msg):
            if how == "method":
                getattr(gb, groupby_func)(*args)
            elif how == "agg":
                gb.agg(groupby_func, *args)
            else:
                gb.transform(groupby_func, *args)


@pytest.mark.parametrize("how", ["agg", "transform"])
def test_groupby_raises_string_udf(how):
    df = DataFrame(
        {
            "a": [1, 1, 1, 2, 2],
            "b": range(5),
            "c": list("xyzwt"),
        }
    )
    gb = df.groupby("a")

    def func(x):
        raise TypeError("Test error message")

    with pytest.raises(TypeError, match="Test error message"):
        getattr(gb, how)(func)


@pytest.mark.parametrize("how", ["method", "agg", "transform"])
def test_groupby_raises_datetime(how, groupby_func, as_index, sort):
    df = DataFrame(
        {
            "a": [1, 1, 1, 2, 2],
            "b": range(5),
            "c": datetime.datetime(2005, 1, 1, 10, 30, 23, 540000),
        }
    )
    args = get_groupby_method_args(groupby_func, df)
    gb = df.groupby("a", as_index=as_index, sort=sort)

    klass, msg = {
        "all": (None, ""),
        "any": (None, ""),
        "bfill": (None, ""),
        "corrwith": (TypeError, "cannot perform __mul__ with this index type"),
        "count": (None, ""),
        "cumcount": (None, ""),
        "cummax": (None, ""),
        "cummin": (None, ""),
        "cumprod": (TypeError, "datetime64 type does not support cumprod operations"),
        "cumsum": (TypeError, "datetime64 type does not support cumsum operations"),
        "diff": (None, ""),
        "ffill": (None, ""),
        "fillna": (None, ""),
        "first": (None, ""),
        "idxmax": (None, ""),
        "idxmin": (None, ""),
        "last": (None, ""),
        "max": (None, ""),
        "mean": (None, ""),
        "median": (None, ""),
        "min": (None, ""),
        "ngroup": (None, ""),
        "nunique": (None, ""),
        "pct_change": (TypeError, "cannot perform __truediv__ with this index type"),
        "prod": (TypeError, "datetime64 type does not support prod"),
        "quantile": (None, ""),
        "rank": (None, ""),
        "sem": (TypeError, "Cannot cast DatetimeArray to dtype float64"),
        "shift": (None, ""),
        "size": (None, ""),
        "skew": (TypeError, r"dtype datetime64\[ns\] does not support reduction"),
        "std": (TypeError, "Cannot cast DatetimeArray to dtype float64"),
        "sum": (TypeError, "datetime64 type does not support sum operations"),
        "var": (None, ""),
    }[groupby_func]

    if klass is None:
        if how == "method":
            getattr(gb, groupby_func)(*args)
        elif how == "agg":
            gb.agg(groupby_func, *args)
        else:
            gb.transform(groupby_func, *args)
    else:
        with pytest.raises(klass, match=msg):
            if how == "method":
                getattr(gb, groupby_func)(*args)
            elif how == "agg":
                gb.agg(groupby_func, *args)
            else:
                gb.transform(groupby_func, *args)


@pytest.mark.parametrize("how", ["agg", "transform"])
def test_groupby_raises_datetime_udf(how):
    df = DataFrame(
        {
            "a": [1, 1, 1, 2, 2],
            "b": range(5),
            "c": datetime.datetime(2005, 1, 1, 10, 30, 23, 540000),
        }
    )
    gb = df.groupby("a")

    def func(x):
        raise TypeError("Test error message")

    with pytest.raises(TypeError, match="Test error message"):
        getattr(gb, how)(func)
