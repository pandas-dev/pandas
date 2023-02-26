# Only tests that raise an error and have no better location should go here.
# Tests for specific groupby methods should go in their respective
# test file.

import datetime

import numpy as np
import pytest

from pandas import (
    Categorical,
    DataFrame,
    Grouper,
    Series,
)
from pandas.tests.groupby import get_groupby_method_args


@pytest.fixture(
    params=[
        "a",
        ["a"],
        ["a", "b"],
        Grouper(key="a"),
        lambda x: x % 2,
        [0, 0, 0, 1, 2, 2, 2, 3, 3],
        np.array([0, 0, 0, 1, 2, 2, 2, 3, 3]),
        dict(zip(range(9), [0, 0, 0, 1, 2, 2, 2, 3, 3])),
        Series([1, 1, 1, 1, 1, 2, 2, 2, 2]),
        [Series([1, 1, 1, 1, 1, 2, 2, 2, 2]), Series([3, 3, 4, 4, 4, 4, 4, 3, 3])],
    ]
)
def by(request):
    return request.param


@pytest.fixture(params=[True, False])
def groupby_series(request):
    return request.param


@pytest.mark.parametrize("how", ["method", "agg", "transform"])
def test_groupby_raises_string(how, by, groupby_series, groupby_func):
    df = DataFrame(
        {
            "a": [1, 1, 1, 1, 1, 2, 2, 2, 2],
            "b": [3, 3, 4, 4, 4, 4, 4, 3, 3],
            "c": range(9),
            "d": list("xyzwtyuio"),
        }
    )
    args = get_groupby_method_args(groupby_func, df)
    gb = df.groupby(by=by)

    if groupby_series:
        gb = gb["d"]

        if groupby_func == "corrwith":
            assert not hasattr(gb, "corrwith")
            return

    klass, msg = {
        "all": (None, ""),
        "any": (None, ""),
        "bfill": (None, ""),
        "corrwith": (TypeError, "Could not convert"),
        "count": (None, ""),
        "cumcount": (None, ""),
        "cummax": (
            (NotImplementedError, TypeError),
            "(function|cummax) is not (implemented|supported) for (this|object) dtype",
        ),
        "cummin": (
            (NotImplementedError, TypeError),
            "(function|cummin) is not (implemented|supported) for (this|object) dtype",
        ),
        "cumprod": (
            (NotImplementedError, TypeError),
            "(function|cumprod) is not (implemented|supported) for (this|object) dtype",
        ),
        "cumsum": (
            (NotImplementedError, TypeError),
            "(function|cumsum) is not (implemented|supported) for (this|object) dtype",
        ),
        "diff": (TypeError, "unsupported operand type"),
        "ffill": (None, ""),
        "fillna": (None, ""),
        "first": (None, ""),
        "idxmax": (TypeError, "'argmax' not allowed for this dtype"),
        "idxmin": (TypeError, "'argmin' not allowed for this dtype"),
        "last": (None, ""),
        "max": (None, ""),
        "mean": (TypeError, "Could not convert xy?z?w?t?y?u?i?o? to numeric"),
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
def test_groupby_raises_string_udf(how, by, groupby_series):
    df = DataFrame(
        {
            "a": [1, 1, 1, 1, 1, 2, 2, 2, 2],
            "b": [3, 3, 4, 4, 4, 4, 4, 3, 3],
            "c": range(9),
            "d": list("xyzwtyuio"),
        }
    )
    gb = df.groupby(by=by)

    if groupby_series:
        gb = gb["d"]

    def func(x):
        raise TypeError("Test error message")

    with pytest.raises(TypeError, match="Test error message"):
        getattr(gb, how)(func)


@pytest.mark.parametrize("how", ["agg", "transform"])
@pytest.mark.parametrize("groupby_func_np", [np.sum, np.mean])
def test_groupby_raises_string_np(how, by, groupby_series, groupby_func_np):
    # GH#50749
    df = DataFrame(
        {
            "a": [1, 1, 1, 1, 1, 2, 2, 2, 2],
            "b": [3, 3, 4, 4, 4, 4, 4, 3, 3],
            "c": range(9),
            "d": list("xyzwtyuio"),
        }
    )
    gb = df.groupby(by=by)

    if groupby_series:
        gb = gb["d"]

    klass, msg = {
        np.sum: (None, ""),
        np.mean: (TypeError, "Could not convert xy?z?w?t?y?u?i?o? to numeric"),
    }[groupby_func_np]

    if klass is None:
        getattr(gb, how)(groupby_func_np)
    else:
        with pytest.raises(klass, match=msg):
            getattr(gb, how)(groupby_func_np)


@pytest.mark.parametrize("how", ["method", "agg", "transform"])
def test_groupby_raises_datetime(how, by, groupby_series, groupby_func):
    df = DataFrame(
        {
            "a": [1, 1, 1, 1, 1, 2, 2, 2, 2],
            "b": [3, 3, 4, 4, 4, 4, 4, 3, 3],
            "c": range(9),
            "d": datetime.datetime(2005, 1, 1, 10, 30, 23, 540000),
        }
    )
    args = get_groupby_method_args(groupby_func, df)
    gb = df.groupby(by=by)

    if groupby_series:
        gb = gb["d"]

        if groupby_func == "corrwith":
            assert not hasattr(gb, "corrwith")
            return

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
        "sem": (None, ""),
        "shift": (None, ""),
        "size": (None, ""),
        "skew": (TypeError, r"dtype datetime64\[ns\] does not support reduction"),
        "std": (None, ""),
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
def test_groupby_raises_datetime_udf(how, by, groupby_series):
    df = DataFrame(
        {
            "a": [1, 1, 1, 1, 1, 2, 2, 2, 2],
            "b": [3, 3, 4, 4, 4, 4, 4, 3, 3],
            "c": range(9),
            "d": datetime.datetime(2005, 1, 1, 10, 30, 23, 540000),
        }
    )

    gb = df.groupby(by=by)

    if groupby_series:
        gb = gb["d"]

    def func(x):
        raise TypeError("Test error message")

    with pytest.raises(TypeError, match="Test error message"):
        getattr(gb, how)(func)


@pytest.mark.parametrize("how", ["agg", "transform"])
@pytest.mark.parametrize("groupby_func_np", [np.sum, np.mean])
def test_groupby_raises_datetime_np(how, by, groupby_series, groupby_func_np):
    # GH#50749
    df = DataFrame(
        {
            "a": [1, 1, 1, 1, 1, 2, 2, 2, 2],
            "b": [3, 3, 4, 4, 4, 4, 4, 3, 3],
            "c": range(9),
            "d": datetime.datetime(2005, 1, 1, 10, 30, 23, 540000),
        }
    )
    gb = df.groupby(by=by)

    if groupby_series:
        gb = gb["d"]

    klass, msg = {
        np.sum: (TypeError, "datetime64 type does not support sum operations"),
        np.mean: (None, ""),
    }[groupby_func_np]

    if klass is None:
        getattr(gb, how)(groupby_func_np)
    else:
        with pytest.raises(klass, match=msg):
            getattr(gb, how)(groupby_func_np)


@pytest.mark.parametrize("how", ["method", "agg", "transform"])
def test_groupby_raises_category(
    how, by, groupby_series, groupby_func, using_copy_on_write
):
    # GH#50749
    df = DataFrame(
        {
            "a": [1, 1, 1, 1, 1, 2, 2, 2, 2],
            "b": [3, 3, 4, 4, 4, 4, 4, 3, 3],
            "c": range(9),
            "d": Categorical(
                ["a", "a", "a", "a", "b", "b", "b", "b", "c"],
                categories=["a", "b", "c", "d"],
                ordered=True,
            ),
        }
    )
    args = get_groupby_method_args(groupby_func, df)
    gb = df.groupby(by=by)

    if groupby_series:
        gb = gb["d"]

        if groupby_func == "corrwith":
            assert not hasattr(gb, "corrwith")
            return

    klass, msg = {
        "all": (None, ""),
        "any": (None, ""),
        "bfill": (None, ""),
        "corrwith": (
            TypeError,
            r"unsupported operand type\(s\) for \*: 'Categorical' and 'int'",
        ),
        "count": (None, ""),
        "cumcount": (None, ""),
        "cummax": (
            (NotImplementedError, TypeError),
            "(category type does not support cummax operations|"
            + "category dtype not supported|"
            + "cummax is not supported for category dtype)",
        ),
        "cummin": (
            (NotImplementedError, TypeError),
            "(category type does not support cummin operations|"
            + "category dtype not supported|"
            "cummin is not supported for category dtype)",
        ),
        "cumprod": (
            (NotImplementedError, TypeError),
            "(category type does not support cumprod operations|"
            + "category dtype not supported|"
            "cumprod is not supported for category dtype)",
        ),
        "cumsum": (
            (NotImplementedError, TypeError),
            "(category type does not support cumsum operations|"
            + "category dtype not supported|"
            "cumsum is not supported for category dtype)",
        ),
        "diff": (
            TypeError,
            r"unsupported operand type\(s\) for -: 'Categorical' and 'Categorical'",
        ),
        "ffill": (None, ""),
        "fillna": (
            TypeError,
            r"Cannot setitem on a Categorical with a new category \(0\), "
            + "set the categories first",
        )
        if not using_copy_on_write
        else (None, ""),  # no-op with CoW
        "first": (None, ""),
        "idxmax": (None, ""),
        "idxmin": (None, ""),
        "last": (None, ""),
        "max": (None, ""),
        "mean": (
            TypeError,
            "'Categorical' with dtype category does not support reduction 'mean'",
        ),
        "median": (
            TypeError,
            "'Categorical' with dtype category does not support reduction 'median'",
        ),
        "min": (None, ""),
        "ngroup": (None, ""),
        "nunique": (None, ""),
        "pct_change": (
            TypeError,
            r"unsupported operand type\(s\) for /: 'Categorical' and 'Categorical'",
        ),
        "prod": (TypeError, "category type does not support prod operations"),
        "quantile": (TypeError, "No matching signature found"),
        "rank": (None, ""),
        "sem": (ValueError, "Cannot cast object dtype to float64"),
        "shift": (None, ""),
        "size": (None, ""),
        "skew": (
            TypeError,
            "'Categorical' with dtype category does not support reduction 'skew'",
        ),
        "std": (ValueError, "Cannot cast object dtype to float64"),
        "sum": (TypeError, "category type does not support sum operations"),
        "var": (
            TypeError,
            "'Categorical' with dtype category does not support reduction 'var'",
        ),
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
def test_groupby_raises_category_udf(how, by, groupby_series):
    # GH#50749
    df = DataFrame(
        {
            "a": [1, 1, 1, 1, 1, 2, 2, 2, 2],
            "b": [3, 3, 4, 4, 4, 4, 4, 3, 3],
            "c": range(9),
            "d": Categorical(
                ["a", "a", "a", "a", "b", "b", "b", "b", "c"],
                categories=["a", "b", "c", "d"],
                ordered=True,
            ),
        }
    )
    gb = df.groupby(by=by)

    if groupby_series:
        gb = gb["d"]

    def func(x):
        raise TypeError("Test error message")

    with pytest.raises(TypeError, match="Test error message"):
        getattr(gb, how)(func)


@pytest.mark.parametrize("how", ["agg", "transform"])
@pytest.mark.parametrize("groupby_func_np", [np.sum, np.mean])
def test_groupby_raises_category_np(how, by, groupby_series, groupby_func_np):
    # GH#50749
    df = DataFrame(
        {
            "a": [1, 1, 1, 1, 1, 2, 2, 2, 2],
            "b": [3, 3, 4, 4, 4, 4, 4, 3, 3],
            "c": range(9),
            "d": Categorical(
                ["a", "a", "a", "a", "b", "b", "b", "b", "c"],
                categories=["a", "b", "c", "d"],
                ordered=True,
            ),
        }
    )
    gb = df.groupby(by=by)

    if groupby_series:
        gb = gb["d"]

    klass, msg = {
        np.sum: (TypeError, "category type does not support sum operations"),
        np.mean: (
            TypeError,
            "'Categorical' with dtype category does not support reduction 'mean'",
        ),
    }[groupby_func_np]

    if klass is None:
        getattr(gb, how)(groupby_func_np)
    else:
        with pytest.raises(klass, match=msg):
            getattr(gb, how)(groupby_func_np)


@pytest.mark.parametrize("how", ["method", "agg", "transform"])
def test_groupby_raises_category_on_category(
    how, by, groupby_series, groupby_func, observed, using_copy_on_write
):
    # GH#50749
    df = DataFrame(
        {
            "a": Categorical(
                ["a", "a", "a", "a", "b", "b", "b", "b", "c"],
                categories=["a", "b", "c", "d"],
                ordered=True,
            ),
            "b": [3, 3, 4, 4, 4, 4, 4, 3, 3],
            "c": range(9),
            "d": Categorical(
                ["a", "a", "a", "a", "b", "b", "c", "c", "c"],
                categories=["a", "b", "c", "d"],
                ordered=True,
            ),
        }
    )
    args = get_groupby_method_args(groupby_func, df)
    gb = df.groupby(by=by, observed=observed)

    if groupby_series:
        gb = gb["d"]

        if groupby_func == "corrwith":
            assert not hasattr(gb, "corrwith")
            return

    empty_groups = any(group.empty for group in gb.groups.values())

    klass, msg = {
        "all": (None, ""),
        "any": (None, ""),
        "bfill": (None, ""),
        "corrwith": (
            TypeError,
            r"unsupported operand type\(s\) for \*: 'Categorical' and 'int'",
        ),
        "count": (None, ""),
        "cumcount": (None, ""),
        "cummax": (
            (NotImplementedError, TypeError),
            "(cummax is not supported for category dtype|"
            + "category dtype not supported|"
            + "category type does not support cummax operations)",
        ),
        "cummin": (
            (NotImplementedError, TypeError),
            "(cummin is not supported for category dtype|"
            + "category dtype not supported|"
            "category type does not support cummin operations)",
        ),
        "cumprod": (
            (NotImplementedError, TypeError),
            "(cumprod is not supported for category dtype|"
            + "category dtype not supported|"
            "category type does not support cumprod operations)",
        ),
        "cumsum": (
            (NotImplementedError, TypeError),
            "(cumsum is not supported for category dtype|"
            + "category dtype not supported|"
            + "category type does not support cumsum operations)",
        ),
        "diff": (TypeError, "unsupported operand type"),
        "ffill": (None, ""),
        "fillna": (
            TypeError,
            r"Cannot setitem on a Categorical with a new category \(0\), "
            + "set the categories first",
        )
        if not using_copy_on_write
        else (None, ""),  # no-op with CoW
        "first": (None, ""),
        "idxmax": (ValueError, "attempt to get argmax of an empty sequence")
        if empty_groups
        else (None, ""),
        "idxmin": (ValueError, "attempt to get argmin of an empty sequence")
        if empty_groups
        else (None, ""),
        "last": (None, ""),
        "max": (None, ""),
        "mean": (
            TypeError,
            "'Categorical' with dtype category does not support reduction 'mean'",
        ),
        "median": (
            TypeError,
            "'Categorical' with dtype category does not support reduction 'median'",
        ),
        "min": (None, ""),
        "ngroup": (None, ""),
        "nunique": (None, ""),
        "pct_change": (TypeError, "unsupported operand type"),
        "prod": (TypeError, "category type does not support prod operations"),
        "quantile": (TypeError, ""),
        "rank": (None, ""),
        "sem": (ValueError, "Cannot cast object dtype to float64"),
        "shift": (None, ""),
        "size": (None, ""),
        "skew": (
            TypeError,
            "'Categorical' with dtype category does not support reduction 'skew'",
        ),
        "std": (ValueError, "Cannot cast object dtype to float64"),
        "sum": (TypeError, "category type does not support sum operations"),
        "var": (
            TypeError,
            "'Categorical' with dtype category does not support reduction 'var'",
        ),
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


def test_subsetting_columns_axis_1_raises():
    # GH 35443
    df = DataFrame({"a": [1], "b": [2], "c": [3]})
    gb = df.groupby("a", axis=1)
    with pytest.raises(ValueError, match="Cannot subset columns when using axis=1"):
        gb["b"]
