"""
Copyright (c) Virtus Lab sp. z o.o. (Ltd.)

Distributed under the terms of the MIT license.

The full license is in the STUBS_LICENSE file, distributed with this software.
"""
# flake8: noqa: F841
# TODO: many functions need return types annotations for pyright
# to run with reportGeneralTypeIssues = true
from pathlib import Path
from typing import (
    List,
    cast,
)

import numpy as np

from pandas._typing import (
    ArrayLike,
    Scalar,
)

import pandas as pd
import pandas._testing as tm
from pandas.core.window import ExponentialMovingWindow
from pandas.util import _test_decorators as td


def test_types_init() -> None:
    pd.Series(1)
    pd.Series((1, 2, 3))
    pd.Series(np.array([1, 2, 3]))
    pd.Series(data=[1, 2, 3, 4], name="series")
    pd.Series(data=[1, 2, 3, 4], dtype=np.int8)
    pd.Series(data={"row1": [1, 2], "row2": [3, 4]})
    pd.Series(data=[1, 2, 3, 4], index=[4, 3, 2, 1], copy=True)


def test_types_any() -> None:
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # bool]", variable has type "bool")
    res1: bool = pd.Series([False, False]).any()  # type: ignore[assignment]
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # bool]", variable has type "bool")
    res2: bool = pd.Series([False, False]).any(  # type: ignore[assignment]
        bool_only=False
    )
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # bool]", variable has type "bool")
    res3: bool = pd.Series([np.nan]).any(skipna=False)  # type: ignore[assignment]


def test_types_all() -> None:
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # bool]", variable has type "bool")
    res1: bool = pd.Series([False, False]).all()  # type: ignore[assignment]
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # bool]", variable has type "bool")
    res2: bool = pd.Series([False, False]).all(  # type: ignore[assignment]
        bool_only=False
    )
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # bool]", variable has type "bool")
    res3: bool = pd.Series([np.nan]).all(skipna=False)  # type: ignore[assignment]


def test_types_csv() -> None:
    s = pd.Series(data=[1, 2, 3])
    # error: Incompatible types in assignment (expression has type "Optional[str]",
    # variable has type "str")
    csv_df: str = s.to_csv()  # type: ignore[assignment]

    with tm.ensure_clean() as path:
        s.to_csv(cast(str, path))
        s2: pd.DataFrame = pd.read_csv(cast(str, path))

    with tm.ensure_clean() as path:
        s.to_csv(Path(cast(str, path)))
        s3: pd.DataFrame = pd.read_csv(Path(cast(str, path)))

    # This keyword was added in 1.1.0
    # https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
    with tm.ensure_clean() as path:
        s.to_csv(cast(str, path), errors="replace")
        s4: pd.DataFrame = pd.read_csv(cast(str, path))


def test_types_copy() -> None:
    s = pd.Series(data=[1, 2, 3, 4])
    s2: pd.Series = s.copy()


def test_types_select() -> None:
    s = pd.Series(data={"row1": 1, "row2": 2})
    s[0]
    s[1:]


def test_types_iloc_iat() -> None:
    s = pd.Series(data={"row1": 1, "row2": 2})
    s2 = pd.Series(data=[1, 2])
    s.loc["row1"]
    s.iat[0]
    s2.loc[0]
    s2.iat[0]


def test_types_loc_at() -> None:
    s = pd.Series(data={"row1": 1, "row2": 2})
    s2 = pd.Series(data=[1, 2])
    s.loc["row1"]
    s.at["row1"]
    s2.loc[1]
    s2.at[1]


def test_types_boolean_indexing() -> None:
    s = pd.Series([0, 1, 2])
    s[s > 1]
    s[s]


def test_types_head_tail() -> None:
    s = pd.Series([0, 1, 2])
    s.head(1)
    s.tail(1)


def test_types_sample() -> None:
    s = pd.Series([0, 1, 2])
    s.sample(frac=0.5)
    s.sample(n=1)


def test_types_nlargest_nsmallest() -> None:
    s = pd.Series([0, 1, 2])
    s.nlargest(1)
    s.nlargest(1, "first")
    s.nsmallest(1, "last")
    s.nsmallest(1, "all")


def test_types_filter() -> None:
    s = pd.Series(data=[1, 2, 3, 4], index=["cow", "coal", "coalesce", ""])
    s.filter(items=["cow"])
    s.filter(regex="co.*")
    s.filter(like="al")


def test_types_setting() -> None:
    s = pd.Series([0, 1, 2])
    s[3] = 4
    s[s == 1] = 5
    s[:] = 3


def test_types_drop() -> None:
    s = pd.Series([0, 1, 2])
    res: pd.Series = s.drop(0)
    res2: pd.Series = s.drop([0, 1])
    res3: pd.Series = s.drop(0, axis=0)
    # error: Incompatible types in assignment (expression has type "Series", variable
    # has type "None")
    res4: None = s.drop(  # type: ignore[assignment]
        [0, 1], inplace=True, errors="raise"
    )
    # error: Incompatible types in assignment (expression has type "Series", variable
    # has type "None")
    res5: None = s.drop(  # type: ignore[assignment]
        [0, 1], inplace=True, errors="ignore"
    )


def test_types_drop_multilevel() -> None:
    index = pd.MultiIndex(
        levels=[["top", "bottom"], ["first", "second", "third"]],
        codes=[[0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 1, 2]],
    )
    s = pd.Series(data=[1, 2, 3, 4, 5, 6], index=index)
    res: pd.Series = s.drop(labels="first", level=1)


def test_types_dropna() -> None:
    s = pd.Series([1, np.nan, np.nan])
    res: pd.Series = s.dropna()
    res2: None = s.dropna(axis=0, inplace=True)


def test_types_fillna() -> None:
    s = pd.Series([1, np.nan, np.nan, 3])
    res: pd.Series = s.fillna(0)
    res2: pd.Series = s.fillna(0, axis="index")
    res3: pd.Series = s.fillna(method="backfill", axis=0)
    res4: None = s.fillna(method="bfill", inplace=True)
    res5: pd.Series = s.fillna(method="pad")
    res6: pd.Series = s.fillna(method="ffill", limit=1)


def test_types_sort_index() -> None:
    s = pd.Series([1, 2, 3], index=[2, 3, 1])
    res: pd.Series = s.sort_index()
    res2: None = s.sort_index(ascending=False, inplace=True)
    res3: pd.Series = s.sort_index(kind="mergesort")


# This was added in 1.1.0 https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
def test_types_sort_index_with_key() -> None:
    s = pd.Series([1, 2, 3], index=["a", "B", "c"])
    res: pd.Series = s.sort_index(key=lambda k: k.str.lower())


def test_types_sort_values() -> None:
    s = pd.Series([4, 2, 1, 3])
    res: pd.Series = s.sort_values(axis=0)
    res2: pd.Series = s.sort_values(ascending=False)
    res3: None = s.sort_values(inplace=True, kind="quicksort")
    res4: pd.Series = s.sort_values(na_position="last")
    res5: pd.Series = s.sort_values(ignore_index=True)


# This was added in 1.1.0 https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
def test_types_sort_values_with_key() -> None:
    s = pd.Series([1, 2, 3], index=[2, 3, 1])
    res: pd.Series = s.sort_values(key=lambda k: -k)


def test_types_shift() -> None:
    s = pd.Series([1, 2, 3])
    s.shift()
    s.shift(axis=0, periods=1)
    s.shift(-1, fill_value=0)


def test_types_rank() -> None:
    s = pd.Series([1, 1, 2, 5, 6, np.nan])
    s.rank()
    s.rank(axis=0, na_option="bottom")
    s.rank(method="min", pct=True)
    s.rank(method="dense", ascending=True)
    s.rank(method="first", numeric_only=True)


def test_types_mean() -> None:
    s = pd.Series([1, 2, 3, np.nan])
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # float]", variable has type "float")
    f1: float = s.mean()  # type: ignore[assignment]
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # float]", variable has type "Series")
    s1: pd.Series = s.mean(axis=0)  # type: ignore[assignment]
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # float]", variable has type "float")
    f2: float = s.mean(skipna=False)  # type: ignore[assignment]
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # float]", variable has type "float")
    f3: float = s.mean(numeric_only=False)  # type: ignore[assignment]


def test_types_median() -> None:
    s = pd.Series([1, 2, 3, np.nan])
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # float]", variable has type "float")
    f1: float = s.median()  # type: ignore[assignment]
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # float]", variable has type "Series")
    s1: pd.Series = s.median(axis=0)  # type: ignore[assignment]
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # float]", variable has type "float")
    f2: float = s.median(skipna=False)  # type: ignore[assignment]
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # float]", variable has type "float")
    f3: float = s.median(numeric_only=False)  # type: ignore[assignment]


def test_types_sum() -> None:
    s = pd.Series([1, 2, 3, np.nan])
    s.sum()
    s.sum(axis=0)
    s.sum(skipna=False)
    s.sum(numeric_only=False)
    s.sum(min_count=4)


def test_types_cumsum() -> None:
    s = pd.Series([1, 2, 3, np.nan])
    s.cumsum()
    s.cumsum(axis=0)
    s.cumsum(skipna=False)


def test_types_min() -> None:
    s = pd.Series([1, 2, 3, np.nan])
    s.min()
    s.min(axis=0)
    s.min(skipna=False)


def test_types_max() -> None:
    s = pd.Series([1, 2, 3, np.nan])
    s.max()
    s.max(axis=0)
    s.max(skipna=False)


def test_types_quantile() -> None:
    s = pd.Series([1, 2, 3, 10])
    s.quantile([0.25, 0.5])
    s.quantile(0.75)
    s.quantile()
    s.quantile(interpolation="nearest")


def test_types_clip() -> None:
    s = pd.Series([-10, 2, 3, 10])
    s.clip(lower=0, upper=5)
    s.clip(lower=0, upper=5, inplace=True)


def test_types_abs() -> None:
    s = pd.Series([-10, 2, 3, 10])
    s.abs()


def test_types_var() -> None:
    s = pd.Series([-10, 2, 3, 10])
    s.var()
    s.var(axis=0, ddof=1)
    s.var(skipna=True, numeric_only=False)


def test_types_std() -> None:
    s = pd.Series([-10, 2, 3, 10])
    s.std()
    s.std(axis=0, ddof=1)
    s.std(skipna=True, numeric_only=False)


def test_types_idxmin() -> None:
    s = pd.Series([-10, 2, 3, 10])
    s.idxmin()
    s.idxmin(axis=0)


def test_types_idxmax() -> None:
    s = pd.Series([-10, 2, 3, 10])
    s.idxmax()
    s.idxmax(axis=0)


def test_types_value_counts() -> None:
    s = pd.Series([1, 2])
    s.value_counts()


def test_types_unique() -> None:
    s = pd.Series([-10, 2, 2, 3, 10, 10])
    s.unique()


def test_types_apply() -> None:
    s = pd.Series([-10, 2, 2, 3, 10, 10])
    s.apply(lambda x: x**2)
    s.apply(np.exp)
    s.apply(str)


def test_types_element_wise_arithmetic() -> None:
    s = pd.Series([0, 1, -10])
    s2 = pd.Series([7, -5, 10])

    s + s2
    s.add(s2, fill_value=0)

    s - s2
    s.sub(s2, fill_value=0)

    s * s2
    s.mul(s2, fill_value=0)

    s / s2
    # error: Unexpected keyword argument "fill_value"
    s.div(s2, fill_value=0)  # type: ignore[call-arg]

    s // s2
    s.floordiv(s2, fill_value=0)

    s % s2
    s.mod(s2, fill_value=0)


def test_types_groupby() -> None:
    s = pd.Series([4, 2, 1, 8], index=["a", "b", "a", "b"])
    s.groupby(["a", "b", "a", "b"])
    s.groupby(level=0)
    s.groupby(s > 2)


# This added in 1.1.0 https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
def test_types_group_by_with_dropna_keyword() -> None:
    s = pd.Series([1, 2, 3, 3], index=["col1", "col2", "col3", np.nan])
    s.groupby(level=0, dropna=True).sum()
    s.groupby(level=0, dropna=False).sum()
    s.groupby(level=0).sum()


@td.skip_if_no("matplotlib")  # type: ignore[misc]
def test_types_plot() -> None:
    s = pd.Series([0, 1, 1, 0, -10])
    s.plot.hist()


def test_types_window() -> None:
    s = pd.Series([0, 1, 1, 0, 5, 1, -10])
    s.expanding()
    s.expanding(axis=0)

    s.rolling(2)
    s.rolling(2, axis=0, center=True)


def test_types_cov() -> None:
    s1 = pd.Series([0, 1, 1, 0, 5, 1, -10])
    s2 = pd.Series([0, 2, 12, -4, 7, 9, 2])
    s1.cov(s2)
    s1.cov(s2, min_periods=1)
    # ddof param was added in 1.1.0 https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
    s1.cov(s2, ddof=2)


def test_update() -> None:
    s1 = pd.Series([0, 1, 1, 0, 5, 1, -10])
    s1.update(pd.Series([0, 2, 12]))
    # Series.update() accepting objects that can be coerced to a
    # Series was added in 1.1.0
    # https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
    s1.update([1, 2, -4, 3])
    s1.update([1, "b", "c", "d"])
    s1.update({1: 9, 3: 4})


# error: Untyped decorator makes function "test_to_markdown" untyped
@td.skip_if_no("tabulate")  # type: ignore[misc]
def test_to_markdown() -> None:
    s = pd.Series([0, 1, 1, 0, 5, 1, -10])
    s.to_markdown()
    s.to_markdown(buf=None, mode="wt")
    # index param was added in 1.1.0 https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
    s.to_markdown(index=False)


# compare() method added in 1.1.0 https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
def test_types_compare() -> None:
    s1 = pd.Series([0, 1, 1, 0, 5, 1, -10])
    s2 = pd.Series([0, 2, 12, -4, 7, 9, 2])
    s1.compare(s2)
    s2.compare(s1, align_axis="columns", keep_shape=True, keep_equal=True)


def test_types_agg() -> None:
    s = pd.Series([1, 2, 3], index=["col1", "col2", "col3"])
    s.agg("min")
    s.agg(x=max, y="min", z=np.mean)
    s.agg("mean", axis=0)


def test_types_describe() -> None:
    s = pd.Series([1, 2, 3, np.datetime64("2000-01-01")])
    s.describe()
    s.describe(percentiles=[0.5], include="all")
    s.describe(exclude=np.number)
    # datetime_is_numeric param added in 1.1.0
    # https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
    s.describe(datetime_is_numeric=True)


def test_types_resample() -> None:
    s = pd.Series(range(9), index=pd.date_range("1/1/2000", periods=9, freq="T"))
    s.resample("3T").sum()
    # origin and offset params added in 1.1.0
    # https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
    s.resample("20min", origin="epoch", offset=pd.Timedelta(value=2, unit="minutes"))


# set_flags() method added in 1.2.0 https://pandas.pydata.org/docs/whatsnew/v1.2.0.html
def test_types_set_flags() -> None:
    pd.Series([1, 2], index=["a", "b"]).set_flags(allows_duplicate_labels=False)
    pd.Series([3, 4], index=["a", "a"]).set_flags(allows_duplicate_labels=True)
    pd.Series([5, 2], index=["a", "a"])


def test_types_getitem() -> None:
    s = pd.Series({"key": [0, 1, 2, 3]})
    key: List[int] = s["key"]
    s2 = pd.Series([0, 1, 2, 3])
    value: int = s2[0]
    s3: pd.Series = s[:2]


def test_types_eq() -> None:
    s1 = pd.Series([1, 2, 3])
    res1: pd.Series = s1 == 1
    s2 = pd.Series([1, 2, 4])
    res2: pd.Series = s1 == s2


def test_types_rename_axis() -> None:
    s: pd.Series = pd.Series([1, 2, 3]).rename_axis("A")


def test_types_values() -> None:
    n1: ArrayLike = pd.Series([1, 2, 3]).values
    n2: ArrayLike = pd.Series(list("aabc")).values
    n3: ArrayLike = pd.Series(list("aabc")).astype("category").values
    n4: ArrayLike = pd.Series(
        pd.date_range("20130101", periods=3, tz="US/Eastern")
    ).values


def test_types_rename() -> None:
    # Scalar
    s1 = pd.Series([1, 2, 3]).rename("A")
    # Hashable Sequence
    s2 = pd.Series([1, 2, 3]).rename(("A", "B"))
    # Optional
    s3 = pd.Series([1, 2, 3]).rename(None)

    # Functions
    def add1(x: int) -> int:
        return x + 1

    s4 = pd.Series([1, 2, 3]).rename(add1)

    # Dictionary
    s5 = pd.Series([1, 2, 3]).rename({1: 10})
    # inplace
    # error: Incompatible types in assignment (expression has type "Optional[Series]",
    # variable has type "None")
    s6: None = pd.Series([1, 2, 3]).rename(  # type: ignore[assignment]
        "A", inplace=True
    )


def test_types_ne() -> None:
    s1 = pd.Series([1, 2, 3])
    s2 = pd.Series([1, 2, 4])
    s3: pd.Series = s1 != s2


def test_types_bfill() -> None:
    s1 = pd.Series([1, 2, 3])
    # error: Incompatible types in assignment (expression has type "Optional[Series]",
    # variable has type "Series")
    s2: pd.Series = s1.bfill(inplace=False)  # type: ignore[assignment]
    # error: Incompatible types in assignment (expression has type "Optional[Series]",
    # variable has type "None")
    s3: None = s1.bfill(inplace=True)  # type: ignore[assignment]


def test_types_ewm() -> None:
    s1 = pd.Series([1, 2, 3])
    w1: ExponentialMovingWindow = s1.ewm(
        com=0.3, min_periods=0, adjust=False, ignore_na=True, axis=0
    )
    w2: ExponentialMovingWindow = s1.ewm(alpha=0.4)
    w3: ExponentialMovingWindow = s1.ewm(span=1.6)
    w4: ExponentialMovingWindow = s1.ewm(halflife=0.7)


def test_types_ffill() -> None:
    s1 = pd.Series([1, 2, 3])
    # error: Incompatible types in assignment (expression has type "Optional[Series]",
    # variable has type "Series")
    s2: pd.Series = s1.ffill(inplace=False)  # type: ignore[assignment]
    # error: Incompatible types in assignment (expression has type "Optional[Series]",
    # variable has type "None")
    s3: None = s1.ffill(inplace=True)  # type: ignore[assignment]


def test_types_as_type() -> None:
    s1 = pd.Series([1, 2, 8, 9])
    s2: pd.Series = s1.astype("int32")


def test_types_dot() -> None:
    s1 = pd.Series([0, 1, 2, 3])
    s2 = pd.Series([-1, 2, -3, 4])
    df1 = pd.DataFrame([[0, 1], [-2, 3], [4, -5], [6, 7]])
    n1 = np.array([[0, 1], [1, 2], [-1, -1], [2, 0]])
    sc1: Scalar = s1.dot(s2)
    sc2: Scalar = s1 @ s2
    s3: pd.Series = s1.dot(df1)
    s4: pd.Series = s1 @ df1
    n2: np.ndarray = s1.dot(n1)
    n3: np.ndarray = s1 @ n1
