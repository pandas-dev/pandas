# Only tests that raise an error and have no better location should go here.
# Tests for specific groupby methods should go in their respective
# test file.

import datetime
import re

import numpy as np
import pytest

from pandas.errors import Pandas4Warning

from pandas import (
    Categorical,
    DataFrame,
    Grouper,
    Series,
)
import pandas._testing as tm
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
        dict(zip(range(9), [0, 0, 0, 1, 2, 2, 2, 3, 3], strict=True)),
        Series([1, 1, 1, 1, 1, 2, 2, 2, 2]),
        [Series([1, 1, 1, 1, 1, 2, 2, 2, 2]), Series([3, 3, 4, 4, 4, 4, 4, 3, 3])],
    ]
)
def by(request):
    return request.param


@pytest.fixture(params=[True, False])
def groupby_series(request):
    return request.param


@pytest.fixture
def df_with_string_col():
    df = DataFrame(
        {
            "a": [1, 1, 1, 1, 1, 2, 2, 2, 2],
            "b": [3, 3, 4, 4, 4, 4, 4, 3, 3],
            "c": range(9),
            "d": list("xyzwtyuio"),
        }
    )
    return df


@pytest.fixture
def df_with_datetime_col():
    df = DataFrame(
        {
            "a": [1, 1, 1, 1, 1, 2, 2, 2, 2],
            "b": [3, 3, 4, 4, 4, 4, 4, 3, 3],
            "c": range(9),
            "d": datetime.datetime(2005, 1, 1, 10, 30, 23, 540000),
        }
    )
    return df


@pytest.fixture
def df_with_cat_col():
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
    return df


def _call_and_check(
    klass, msg, how, gb, groupby_func, args, warn_category=None, warn_msg=""
):
    with tm.assert_produces_warning(
        warn_category, match=warn_msg, check_stacklevel=False
    ):
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


@pytest.mark.parametrize("how", ["method", "agg", "transform"])
def test_groupby_raises_string(
    how, by, groupby_series, groupby_func, df_with_string_col, using_infer_string
):
    df = df_with_string_col
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
        "first": (None, ""),
        "idxmax": (None, ""),
        "idxmin": (None, ""),
        "last": (None, ""),
        "max": (None, ""),
        "mean": (
            TypeError,
            re.escape("agg function failed [how->mean,dtype->object]"),
        ),
        "median": (
            TypeError,
            re.escape("agg function failed [how->median,dtype->object]"),
        ),
        "min": (None, ""),
        "ngroup": (None, ""),
        "nunique": (None, ""),
        "pct_change": (TypeError, "unsupported operand type"),
        "prod": (
            TypeError,
            re.escape("agg function failed [how->prod,dtype->object]"),
        ),
        "quantile": (TypeError, "dtype 'object' does not support operation 'quantile'"),
        "rank": (None, ""),
        "sem": (ValueError, "could not convert string to float"),
        "shift": (None, ""),
        "size": (None, ""),
        "skew": (ValueError, "could not convert string to float"),
        "kurt": (ValueError, "could not convert string to float"),
        "std": (ValueError, "could not convert string to float"),
        "sum": (None, ""),
        "var": (
            TypeError,
            re.escape("agg function failed [how->var,dtype->"),
        ),
    }[groupby_func]

    if using_infer_string:
        if groupby_func in [
            "prod",
            "mean",
            "median",
            "cumsum",
            "cumprod",
            "std",
            "sem",
            "var",
            "skew",
            "kurt",
            "quantile",
        ]:
            msg = f"dtype 'str' does not support operation '{groupby_func}'"
            if groupby_func in ["sem", "std", "skew", "kurt"]:
                # The object-dtype raises ValueError when trying to convert to numeric.
                klass = TypeError
        elif groupby_func == "pct_change" and df["d"].dtype.storage == "pyarrow":
            # This doesn't go through EA._groupby_op so the message isn't controlled
            #  there.
            msg = "operation 'truediv' not supported for dtype 'str' with dtype 'str'"
        elif groupby_func == "diff" and df["d"].dtype.storage == "pyarrow":
            # This doesn't go through EA._groupby_op so the message isn't controlled
            #  there.
            msg = "operation 'sub' not supported for dtype 'str' with dtype 'str'"

        elif groupby_func in ["cummin", "cummax"]:
            msg = msg.replace("object", "str")
        elif groupby_func == "corrwith":
            msg = "Cannot perform reduction 'mean' with string dtype"

    if groupby_func == "corrwith":
        warn_category = Pandas4Warning
        warn_msg = "DataFrameGroupBy.corrwith is deprecated"
    else:
        warn_category = None
        warn_msg = ""
    _call_and_check(klass, msg, how, gb, groupby_func, args, warn_category, warn_msg)


@pytest.mark.parametrize("how", ["agg", "transform"])
def test_groupby_raises_string_udf(how, by, groupby_series, df_with_string_col):
    df = df_with_string_col
    gb = df.groupby(by=by)

    if groupby_series:
        gb = gb["d"]

    def func(x):
        raise TypeError("Test error message")

    with pytest.raises(TypeError, match="Test error message"):
        getattr(gb, how)(func)


@pytest.mark.parametrize("how", ["agg", "transform"])
@pytest.mark.parametrize("groupby_func_np", [np.sum, np.mean])
def test_groupby_raises_string_np(
    how,
    by,
    groupby_series,
    groupby_func_np,
    df_with_string_col,
    using_infer_string,
):
    # GH#50749
    df = df_with_string_col
    gb = df.groupby(by=by)

    if groupby_series:
        gb = gb["d"]

    klass, msg = {
        np.sum: (None, ""),
        np.mean: (
            TypeError,
            "Could not convert string .* to numeric|"
            "Cannot perform reduction 'mean' with string dtype",
        ),
    }[groupby_func_np]

    if using_infer_string:
        if groupby_func_np is np.mean:
            klass = TypeError
        msg = f"Cannot perform reduction '{groupby_func_np.__name__}' with string dtype"

    _call_and_check(klass, msg, how, gb, groupby_func_np, ())


@pytest.mark.parametrize("how", ["method", "agg", "transform"])
def test_groupby_raises_datetime(
    how, by, groupby_series, groupby_func, df_with_datetime_col
):
    df = df_with_datetime_col
    args = get_groupby_method_args(groupby_func, df)
    gb = df.groupby(by=by)

    if groupby_series:
        gb = gb["d"]

        if groupby_func == "corrwith":
            assert not hasattr(gb, "corrwith")
            return

    klass, msg = {
        "all": (TypeError, "'all' with datetime64 dtypes is no longer supported"),
        "any": (TypeError, "'any' with datetime64 dtypes is no longer supported"),
        "bfill": (None, ""),
        "corrwith": (TypeError, "cannot perform __mul__ with this index type"),
        "count": (None, ""),
        "cumcount": (None, ""),
        "cummax": (None, ""),
        "cummin": (None, ""),
        "cumprod": (TypeError, "datetime64 type does not support operation 'cumprod'"),
        "cumsum": (TypeError, "datetime64 type does not support operation 'cumsum'"),
        "diff": (None, ""),
        "ffill": (None, ""),
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
        "prod": (TypeError, "datetime64 type does not support operation 'prod'"),
        "quantile": (None, ""),
        "rank": (None, ""),
        "sem": (None, ""),
        "shift": (None, ""),
        "size": (None, ""),
        "skew": (
            TypeError,
            "|".join(
                [
                    r"dtype datetime64\[ns\] does not support operation",
                    "datetime64 type does not support operation 'skew'",
                ]
            ),
        ),
        "kurt": (
            TypeError,
            "|".join(
                [
                    r"dtype datetime64\[ns\] does not support operation",
                    "datetime64 type does not support operation 'kurt'",
                ]
            ),
        ),
        "std": (None, ""),
        "sum": (TypeError, "datetime64 type does not support operation 'sum"),
        "var": (TypeError, "datetime64 type does not support operation 'var'"),
    }[groupby_func]

    if groupby_func == "corrwith":
        warn_category = Pandas4Warning
        warn_msg = "DataFrameGroupBy.corrwith is deprecated"
    else:
        warn_category = None
        warn_msg = ""
    _call_and_check(klass, msg, how, gb, groupby_func, args, warn_category, warn_msg)


@pytest.mark.parametrize("how", ["agg", "transform"])
def test_groupby_raises_datetime_udf(how, by, groupby_series, df_with_datetime_col):
    df = df_with_datetime_col
    gb = df.groupby(by=by)

    if groupby_series:
        gb = gb["d"]

    def func(x):
        raise TypeError("Test error message")

    with pytest.raises(TypeError, match="Test error message"):
        getattr(gb, how)(func)


@pytest.mark.parametrize("how", ["agg", "transform"])
@pytest.mark.parametrize("groupby_func_np", [np.sum, np.mean])
def test_groupby_raises_datetime_np(
    how, by, groupby_series, groupby_func_np, df_with_datetime_col
):
    # GH#50749
    df = df_with_datetime_col
    gb = df.groupby(by=by)

    if groupby_series:
        gb = gb["d"]

    klass, msg = {
        np.sum: (
            TypeError,
            re.escape("datetime64[us] does not support operation 'sum'"),
        ),
        np.mean: (None, ""),
    }[groupby_func_np]
    _call_and_check(klass, msg, how, gb, groupby_func_np, ())


@pytest.mark.parametrize("func", ["prod", "cumprod", "skew", "kurt", "var"])
def test_groupby_raises_timedelta(func):
    df = DataFrame(
        {
            "a": [1, 1, 1, 1, 1, 2, 2, 2, 2],
            "b": [3, 3, 4, 4, 4, 4, 4, 3, 3],
            "c": range(9),
            "d": datetime.timedelta(days=1),
        }
    )
    gb = df.groupby(by="a")

    _call_and_check(
        TypeError,
        "timedelta64 type does not support .* operations",
        "method",
        gb,
        func,
        [],
    )


@pytest.mark.parametrize("how", ["method", "agg", "transform"])
def test_groupby_raises_category(
    how, by, groupby_series, groupby_func, df_with_cat_col
):
    # GH#50749
    df = df_with_cat_col
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
            "category dtype not supported|"
            "cummax is not supported for category dtype)",
        ),
        "cummin": (
            (NotImplementedError, TypeError),
            "(category type does not support cummin operations|"
            "category dtype not supported|"
            "cummin is not supported for category dtype)",
        ),
        "cumprod": (
            (NotImplementedError, TypeError),
            "(category type does not support cumprod operations|"
            "category dtype not supported|"
            "cumprod is not supported for category dtype)",
        ),
        "cumsum": (
            (NotImplementedError, TypeError),
            "(category type does not support cumsum operations|"
            "category dtype not supported|"
            "cumsum is not supported for category dtype)",
        ),
        "diff": (
            TypeError,
            r"unsupported operand type\(s\) for -: 'Categorical' and 'Categorical'",
        ),
        "ffill": (None, ""),
        "first": (None, ""),
        "idxmax": (None, ""),
        "idxmin": (None, ""),
        "last": (None, ""),
        "max": (None, ""),
        "mean": (
            TypeError,
            "|".join(
                [
                    "'Categorical' .* does not support operation 'mean'",
                    "category dtype does not support aggregation 'mean'",
                ]
            ),
        ),
        "median": (
            TypeError,
            "|".join(
                [
                    "'Categorical' .* does not support operation 'median'",
                    "category dtype does not support aggregation 'median'",
                ]
            ),
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
        "sem": (
            TypeError,
            "|".join(
                [
                    "'Categorical' .* does not support operation 'sem'",
                    "category dtype does not support aggregation 'sem'",
                ]
            ),
        ),
        "shift": (None, ""),
        "size": (None, ""),
        "skew": (
            TypeError,
            "|".join(
                [
                    "dtype category does not support operation 'skew'",
                    "category type does not support skew operations",
                ]
            ),
        ),
        "kurt": (
            TypeError,
            "|".join(
                [
                    "dtype category does not support operation 'kurt'",
                    "category type does not support kurt operations",
                ]
            ),
        ),
        "std": (
            TypeError,
            "|".join(
                [
                    "'Categorical' .* does not support operation 'std'",
                    "category dtype does not support aggregation 'std'",
                ]
            ),
        ),
        "sum": (TypeError, "category type does not support sum operations"),
        "var": (
            TypeError,
            "|".join(
                [
                    "'Categorical' .* does not support operation 'var'",
                    "category dtype does not support aggregation 'var'",
                ]
            ),
        ),
    }[groupby_func]

    if groupby_func == "corrwith":
        warn_category = Pandas4Warning
        warn_msg = "DataFrameGroupBy.corrwith is deprecated"
    else:
        warn_category = None
        warn_msg = ""
    _call_and_check(klass, msg, how, gb, groupby_func, args, warn_category, warn_msg)


@pytest.mark.parametrize("how", ["agg", "transform"])
def test_groupby_raises_category_udf(how, by, groupby_series, df_with_cat_col):
    # GH#50749
    df = df_with_cat_col
    gb = df.groupby(by=by)

    if groupby_series:
        gb = gb["d"]

    def func(x):
        raise TypeError("Test error message")

    with pytest.raises(TypeError, match="Test error message"):
        getattr(gb, how)(func)


@pytest.mark.parametrize("how", ["agg", "transform"])
@pytest.mark.parametrize("groupby_func_np", [np.sum, np.mean])
def test_groupby_raises_category_np(
    how, by, groupby_series, groupby_func_np, df_with_cat_col
):
    # GH#50749
    df = df_with_cat_col
    gb = df.groupby(by=by)

    if groupby_series:
        gb = gb["d"]

    klass, msg = {
        np.sum: (TypeError, "dtype category does not support operation 'sum'"),
        np.mean: (
            TypeError,
            "dtype category does not support operation 'mean'",
        ),
    }[groupby_func_np]
    _call_and_check(klass, msg, how, gb, groupby_func_np, ())


@pytest.mark.filterwarnings("ignore:In a future version, the keys")
@pytest.mark.parametrize("how", ["method", "agg", "transform"])
def test_groupby_raises_category_on_category(
    how,
    by,
    groupby_series,
    groupby_func,
    observed,
    df_with_cat_col,
):
    # GH#50749
    df = df_with_cat_col
    df["a"] = Categorical(
        ["a", "a", "a", "a", "b", "b", "b", "b", "c"],
        categories=["a", "b", "c", "d"],
        ordered=True,
    )
    args = get_groupby_method_args(groupby_func, df)
    gb = df.groupby(by=by, observed=observed)

    if groupby_series:
        gb = gb["d"]

        if groupby_func == "corrwith":
            assert not hasattr(gb, "corrwith")
            return

    empty_groups = not observed and any(group.empty for group in gb.groups.values())
    if how == "transform":
        # empty groups will be ignored
        empty_groups = False

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
            "category dtype not supported|"
            "category type does not support cummax operations)",
        ),
        "cummin": (
            (NotImplementedError, TypeError),
            "(cummin is not supported for category dtype|"
            "category dtype not supported|"
            "category type does not support cummin operations)",
        ),
        "cumprod": (
            (NotImplementedError, TypeError),
            "(cumprod is not supported for category dtype|"
            "category dtype not supported|"
            "category type does not support cumprod operations)",
        ),
        "cumsum": (
            (NotImplementedError, TypeError),
            "(cumsum is not supported for category dtype|"
            "category dtype not supported|"
            "category type does not support cumsum operations)",
        ),
        "diff": (TypeError, "unsupported operand type"),
        "ffill": (None, ""),
        "first": (None, ""),
        "idxmax": (ValueError, "empty group due to unobserved categories")
        if empty_groups
        else (None, ""),
        "idxmin": (ValueError, "empty group due to unobserved categories")
        if empty_groups
        else (None, ""),
        "last": (None, ""),
        "max": (None, ""),
        "mean": (TypeError, "category dtype does not support aggregation 'mean'"),
        "median": (TypeError, "category dtype does not support aggregation 'median'"),
        "min": (None, ""),
        "ngroup": (None, ""),
        "nunique": (None, ""),
        "pct_change": (TypeError, "unsupported operand type"),
        "prod": (TypeError, "category type does not support prod operations"),
        "quantile": (TypeError, "No matching signature found"),
        "rank": (None, ""),
        "sem": (
            TypeError,
            "|".join(
                [
                    "'Categorical' .* does not support operation 'sem'",
                    "category dtype does not support aggregation 'sem'",
                ]
            ),
        ),
        "shift": (None, ""),
        "size": (None, ""),
        "skew": (
            TypeError,
            "|".join(
                [
                    "category type does not support skew operations",
                    "dtype category does not support operation 'skew'",
                ]
            ),
        ),
        "kurt": (
            TypeError,
            "|".join(
                [
                    "category type does not support kurt operations",
                    "dtype category does not support operation 'kurt'",
                ]
            ),
        ),
        "std": (
            TypeError,
            "|".join(
                [
                    "'Categorical' .* does not support operation 'std'",
                    "category dtype does not support aggregation 'std'",
                ]
            ),
        ),
        "sum": (TypeError, "category type does not support sum operations"),
        "var": (
            TypeError,
            "|".join(
                [
                    "'Categorical' .* does not support operation 'var'",
                    "category dtype does not support aggregation 'var'",
                ]
            ),
        ),
    }[groupby_func]

    if groupby_func == "corrwith":
        warn_category = Pandas4Warning
        warn_msg = "DataFrameGroupBy.corrwith is deprecated"
    else:
        warn_category = None
        warn_msg = ""
    _call_and_check(klass, msg, how, gb, groupby_func, args, warn_category, warn_msg)
