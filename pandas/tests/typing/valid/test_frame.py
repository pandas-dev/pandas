"""
Copyright (c) Virtus Lab sp. z o.o. (Ltd.)

Distributed under the terms of the MIT license.

The full license is in the STUBS_LICENSE file, distributed with this software.
"""
# flake8: noqa: F841
# TODO: many functions need return types annotations for pyright
# to run with reportGeneralTypeIssues = true
import io
from pathlib import Path
import tempfile
from typing import (
    Any,
    Iterable,
    List,
    Tuple,
)

import numpy as np

import pandas as pd
from pandas.util import _test_decorators as td


def test_types_init() -> None:
    pd.DataFrame(data={"col1": [1, 2], "col2": [3, 4]})
    pd.DataFrame(data={"col1": [1, 2], "col2": [3, 4]}, index=[2, 1])
    pd.DataFrame(data=[1, 2, 3, 4], dtype=np.int8)
    pd.DataFrame(
        np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
        columns=["a", "b", "c"],
        dtype=np.int8,
        copy=True,
    )


def test_types_to_csv() -> None:
    df = pd.DataFrame(data={"col1": [1, 2], "col2": [3, 4]})
    # error: Incompatible types in assignment (expression has type "Optional[str]",
    # variable has type "str")
    csv_df: str = df.to_csv()  # type: ignore[assignment]

    with tempfile.NamedTemporaryFile() as file:
        df.to_csv(file.name)
        df2: pd.DataFrame = pd.read_csv(file.name)

    with tempfile.NamedTemporaryFile() as file:
        df.to_csv(Path(file.name))
        df3: pd.DataFrame = pd.read_csv(Path(file.name))

    # This keyword was added in 1.1.0
    # https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
    with tempfile.NamedTemporaryFile() as file:
        df.to_csv(file.name, errors="replace")
        df4: pd.DataFrame = pd.read_csv(file.name)

    # Testing support for binary file handles, added in 1.2.0
    # https://pandas.pydata.org/docs/whatsnew/v1.2.0.html
    df.to_csv(io.BytesIO(), encoding="utf-8", compression="gzip")


def test_types_to_csv_when_path_passed() -> None:
    df = pd.DataFrame(data={"col1": [1, 2], "col2": [3, 4]})
    path: Path = Path("./dummy_path.txt")
    try:
        assert not path.exists()
        df.to_csv(path)
        df5: pd.DataFrame = pd.read_csv(path)
    finally:
        path.unlink()


def test_types_copy() -> None:
    df = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    df2: pd.DataFrame = df.copy()


def test_types_getitem() -> None:
    df = pd.DataFrame(data={"col1": [1, 2], "col2": [3, 4], 5: [6, 7]})
    i = pd.Index(["col1", "col2"])
    s = pd.Series(["col1", "col2"])
    select_df = pd.DataFrame({"col1": [True, True], "col2": [False, True]})
    a = np.array(["col1", "col2"])
    df["col1"]
    df[5]
    df[["col1", "col2"]]
    df[1:]
    df[s]
    df[a]
    df[select_df]
    df[i]


def test_slice_setitem() -> None:
    # Due to the bug in pandas 1.2.3
    # (https://github.com/pandas-dev/pandas/issues/40440),
    # this is in separate test case
    df = pd.DataFrame(data={"col1": [1, 2], "col2": [3, 4], 5: [6, 7]})
    df[1:] = ["a", "b", "c"]


def test_types_setitem() -> None:
    df = pd.DataFrame(data={"col1": [1, 2], "col2": [3, 4], 5: [6, 7]})
    i = pd.Index(["col1", "col2"])
    s = pd.Series(["col1", "col2"])
    a = np.array(["col1", "col2"])
    df["col1"] = [1, 2]
    df[5] = [5, 6]
    df[["col1", "col2"]] = [[1, 2], [3, 4]]
    df[s] = [5, 6]
    df[a] = [[1, 2], [3, 4]]
    df[i] = [8, 9]


def test_types_setitem_mask() -> None:
    df = pd.DataFrame(data={"col1": [1, 2], "col2": [3, 4], 5: [6, 7]})
    select_df = pd.DataFrame({"col1": [True, True], "col2": [False, True]})
    df[select_df] = [1, 2, 3]


def test_types_iloc_iat() -> None:
    df = pd.DataFrame(data={"col1": [1, 2], "col2": [3, 4]})
    df.iloc[1, 1]
    df.iloc[[1], [1]]
    df.iat[0, 0]


def test_types_loc_at() -> None:
    df = pd.DataFrame(data={"col1": [1, 2], "col2": [3, 4]})
    df.loc[[0], "col1"]
    df.at[0, "col1"]


def test_types_boolean_indexing() -> None:
    df = pd.DataFrame(data={"col1": [1, 2], "col2": [3, 4]})
    df[df > 1]
    df[~(df > 1.0)]


def test_types_head_tail() -> None:
    df = pd.DataFrame(data={"col1": [1, 2], "col2": [3, 4]})
    df.head(1)
    df.tail(1)


def test_types_assign() -> None:
    df = pd.DataFrame(data={"col1": [1, 2], "col2": [3, 4]})
    df.assign(col3=lambda frame: frame.sum(axis=1))
    df["col3"] = df.sum(axis=1)


def test_types_sample() -> None:
    df = pd.DataFrame(data={"col1": [1, 2], "col2": [3, 4]})
    df.sample(frac=0.5)
    df.sample(n=1)


def test_types_nlargest_nsmallest() -> None:
    df = pd.DataFrame(data={"col1": [1, 2], "col2": [3, 4]})
    df.nlargest(1, "col1")
    df.nsmallest(1, "col2")


def test_types_filter() -> None:
    df = pd.DataFrame(data={"col1": [1, 2], "col2": [3, 4]})
    df.filter(items=["col1"])
    df.filter(regex="co.*")
    df.filter(like="1")


def test_types_setting() -> None:
    df = pd.DataFrame(data={"col1": [1, 2], "col2": [3, 4]})
    df["col1"] = 1
    df[df == 1] = 7


def test_types_drop() -> None:
    df = pd.DataFrame(data={"col1": [1, 2], "col2": [3, 4]})
    res: pd.DataFrame = df.drop("col1", axis=1)
    res2: pd.DataFrame = df.drop(columns=["col1"])
    res3: pd.DataFrame = df.drop({0})
    res4: pd.DataFrame = df.drop(index={0})
    res5: pd.DataFrame = df.drop(columns={"col1"})
    res6: pd.DataFrame = df.drop(index=1)
    res7: pd.DataFrame = df.drop(labels=0)
    res8: None = df.drop([0, 0], inplace=True)


def test_types_dropna() -> None:
    df = pd.DataFrame(data={"col1": [np.nan, np.nan], "col2": [3, np.nan]})
    res: pd.DataFrame = df.dropna()
    res2: pd.DataFrame = df.dropna(axis=1, thresh=1)
    res3: None = df.dropna(axis=0, how="all", subset=["col1"], inplace=True)


def test_types_fillna() -> None:
    df = pd.DataFrame(data={"col1": [np.nan, np.nan], "col2": [3, np.nan]})
    res: pd.DataFrame = df.fillna(0)
    res2: None = df.fillna(method="pad", axis=1, inplace=True)


def test_types_sort_index() -> None:
    df = pd.DataFrame(data={"col1": [1, 2, 3, 4]}, index=[5, 1, 3, 2])
    df2 = pd.DataFrame(data={"col1": [1, 2, 3, 4]}, index=["a", "b", "c", "d"])
    res: pd.DataFrame = df.sort_index()
    level1 = (1, 2)
    res2: pd.DataFrame = df.sort_index(ascending=False, level=level1)
    level2: List[str] = ["a", "b", "c"]
    # error: Argument "level" to "sort_index" of "DataFrame" has incompatible type
    # "List[str]"; expected "Optional[Union[Hashable, int]]"
    res3: pd.DataFrame = df2.sort_index(level=level2)  # type: ignore[arg-type]
    res4: pd.DataFrame = df.sort_index(ascending=False, level=3)
    res5: None = df.sort_index(kind="mergesort", inplace=True)


# This was added in 1.1.0 https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
def test_types_sort_index_with_key() -> None:
    df = pd.DataFrame(data={"col1": [1, 2, 3, 4]}, index=["a", "b", "C", "d"])
    res: pd.DataFrame = df.sort_index(key=lambda k: k.str.lower())


def test_types_set_index() -> None:
    df = pd.DataFrame(
        data={"col1": [1, 2, 3, 4], "col2": ["a", "b", "c", "d"]}, index=[5, 1, 3, 2]
    )
    res: pd.DataFrame = df.set_index("col1")
    res2: pd.DataFrame = df.set_index("col1", drop=False)
    res3: pd.DataFrame = df.set_index("col1", append=True)
    res4: pd.DataFrame = df.set_index("col1", verify_integrity=True)
    res5: pd.DataFrame = df.set_index(["col1", "col2"])
    res6: None = df.set_index("col1", inplace=True)


def test_types_query() -> None:
    df = pd.DataFrame(data={"col1": [1, 2, 3, 4], "col2": [3, 0, 1, 7]})
    res: pd.DataFrame = df.query("col1 > col2")
    res2: None = df.query("col1 % col2 == 0", inplace=True)


def test_types_eval() -> None:
    df = pd.DataFrame(data={"col1": [1, 2, 3, 4], "col2": [3, 0, 1, 7]})
    df.eval("col1 > col2")
    res: None = df.eval("C = col1 % col2 == 0", inplace=True)


def test_types_sort_values() -> None:
    df = pd.DataFrame(data={"col1": [2, 1], "col2": [3, 4]})
    res: pd.DataFrame = df.sort_values("col1")
    res2: None = df.sort_values("col1", ascending=False, inplace=True)
    res3: pd.DataFrame = df.sort_values(by=["col1", "col2"], ascending=[True, False])


# This was added in 1.1.0 https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
def test_types_sort_values_with_key() -> None:
    df = pd.DataFrame(data={"col1": [2, 1], "col2": [3, 4]})
    res: pd.DataFrame = df.sort_values(by="col1", key=lambda k: -k)


def test_types_shift() -> None:
    df = pd.DataFrame(data={"col1": [1, 1], "col2": [3, 4]})
    df.shift()
    df.shift(1)
    df.shift(-1)


def test_types_rank() -> None:
    df = pd.DataFrame(data={"col1": [2, 1], "col2": [3, 4]})
    df.rank(axis=0, na_option="bottom")
    df.rank(method="min", pct=True)
    df.rank(method="dense", ascending=True)
    df.rank(method="first", numeric_only=True)


def test_types_mean() -> None:
    df = pd.DataFrame(data={"col1": [2, 1], "col2": [3, 4]})
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # float]", variable has type "Series")
    s1: pd.Series = df.mean()  # type: ignore[assignment]
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # float]", variable has type "Series")
    s2: pd.Series = df.mean(axis=0)  # type: ignore[assignment]
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # float]", variable has type "DataFrame")
    df2: pd.DataFrame = df.mean(level=0)  # type: ignore[assignment]
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # float]", variable has type "DataFrame")
    df3: pd.DataFrame = df.mean(axis=1, level=0)  # type: ignore[assignment]
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # float]", variable has type "DataFrame")
    df4: pd.DataFrame = df.mean(1, True, level=0)  # type: ignore[assignment]
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # float]", variable has type "Series")# error: Incompatible types in assignment
    # (expression has type "Union[Series, float]", variable has type "Series")
    s3: pd.Series = df.mean(  # type: ignore[assignment]
        axis=1, skipna=True, numeric_only=False
    )


def test_types_median() -> None:
    df = pd.DataFrame(data={"col1": [2, 1], "col2": [3, 4]})
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # float]", variable has type "Series")
    s1: pd.Series = df.median()  # type: ignore[assignment]
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # float]", variable has type "Series")
    s2: pd.Series = df.median(axis=0)  # type: ignore[assignment]
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # float]", variable has type "DataFrame")
    df2: pd.DataFrame = df.median(level=0)  # type: ignore[assignment]
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # float]", variable has type "DataFrame")
    df3: pd.DataFrame = df.median(axis=1, level=0)  # type: ignore[assignment]
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # float]", variable has type "DataFrame")
    df4: pd.DataFrame = df.median(1, True, level=0)  # type: ignore[assignment]
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # float]", variable has type "Series")
    s3: pd.Series = df.median(  # type: ignore[assignment]
        axis=1, skipna=True, numeric_only=False
    )


def test_types_itertuples() -> None:
    df = pd.DataFrame(data={"col1": [2, 1], "col2": [3, 4]})
    res1: Iterable[Tuple[Any, ...]] = df.itertuples()
    res2: Iterable[Tuple[Any, ...]] = df.itertuples(index=False, name="Foobar")
    res3: Iterable[Tuple[Any, ...]] = df.itertuples(index=False, name=None)


def test_types_sum() -> None:
    df = pd.DataFrame(data={"col1": [2, 1], "col2": [3, 4]})
    df.sum()
    df.sum(axis=1)


def test_types_cumsum() -> None:
    df = pd.DataFrame(data={"col1": [2, 1], "col2": [3, 4]})
    df.cumsum()
    df.sum(axis=0)


def test_types_min() -> None:
    df = pd.DataFrame(data={"col1": [2, 1], "col2": [3, 4]})
    df.min()
    df.min(axis=0)


def test_types_max() -> None:
    df = pd.DataFrame(data={"col1": [2, 1], "col2": [3, 4]})
    df.max()
    df.max(axis=0)


def test_types_quantile() -> None:
    df = pd.DataFrame(data={"col1": [2, 1], "col2": [3, 4]})
    df.quantile([0.25, 0.5])
    df.quantile(0.75)
    df.quantile()


def test_types_clip() -> None:
    df = pd.DataFrame(data={"col1": [20, 12], "col2": [3, 14]})
    df.clip(lower=5, upper=15)


def test_types_abs() -> None:
    df = pd.DataFrame(data={"col1": [-5, 1], "col2": [3, -14]})
    df.abs()


def test_types_var() -> None:
    df = pd.DataFrame(data={"col1": [2, 1], "col2": [1, 4]})
    df.var()
    df.var(axis=1, ddof=1)
    df.var(skipna=True, numeric_only=False)


def test_types_std() -> None:
    df = pd.DataFrame(data={"col1": [2, 1], "col2": [1, 4]})
    df.std()
    df.std(axis=1, ddof=1)
    df.std(skipna=True, numeric_only=False)


def test_types_idxmin() -> None:
    df = pd.DataFrame(data={"col1": [2, 1], "col2": [3, 4]})
    df.idxmin()
    df.idxmin(axis=0)


def test_types_idxmax() -> None:
    df = pd.DataFrame(data={"col1": [2, 1], "col2": [3, 4]})
    df.idxmax()
    df.idxmax(axis=0)


# This was added in 1.1.0 https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
def test_types_value_counts() -> None:
    df = pd.DataFrame(data={"col1": [1, 2], "col2": [1, 4]})
    df.value_counts()


def test_types_unique() -> None:
    # This is really more for of a Series test
    df = pd.DataFrame(data={"col1": [1, 2], "col2": [1, 4]})
    df["col1"].unique()


def test_types_apply() -> None:
    df = pd.DataFrame(data={"col1": [2, 1], "col2": [3, 4]})
    df.apply(lambda x: x ** 2)
    df.apply(np.exp)
    df.apply(str)


def test_types_applymap() -> None:
    df = pd.DataFrame(data={"col1": [2, 1], "col2": [3, 4]})
    df.applymap(lambda x: x ** 2)
    df.applymap(np.exp)
    df.applymap(str)
    # na_action parameter was added in 1.2.0
    # https://pandas.pydata.org/docs/whatsnew/v1.2.0.html
    df.applymap(np.exp, na_action="ignore")
    df.applymap(str, na_action=None)


def test_types_element_wise_arithmetic() -> None:
    df = pd.DataFrame(data={"col1": [2, 1], "col2": [3, 4]})
    df2 = pd.DataFrame(data={"col1": [10, 20], "col3": [3, 4]})

    df + df2
    df.add(df2, fill_value=0)

    df - df2
    df.sub(df2, fill_value=0)

    df * df2
    df.mul(df2, fill_value=0)

    df / df2
    df.div(df2, fill_value=0)

    df // df2
    df.floordiv(df2, fill_value=0)

    df % df2
    df.mod(df2, fill_value=0)

    # divmod operation was added in 1.2.0
    # https://pandas.pydata.org/docs/whatsnew/v1.2.0.html
    # noinspection PyTypeChecker
    divmod(df, df2)
    df.__divmod__(df2)
    df.__rdivmod__(df2)


def test_types_melt() -> None:
    df = pd.DataFrame(data={"col1": [1, 2], "col2": [3, 4]})
    df.melt()
    df.melt(id_vars=["col1"], value_vars=["col2"])
    df.melt(
        id_vars=["col1"],
        value_vars=["col2"],
        var_name="someVariable",
        value_name="someValue",
    )

    pd.melt(df)
    pd.melt(df, id_vars=["col1"], value_vars=["col2"])
    pd.melt(
        df,
        id_vars=["col1"],
        value_vars=["col2"],
        var_name="someVariable",
        value_name="someValue",
    )


def test_types_pivot() -> None:
    df = pd.DataFrame(
        data={
            "col1": ["first", "second", "third", "fourth"],
            "col2": [50, 70, 56, 111],
            "col3": ["A", "B", "B", "A"],
            "col4": [100, 102, 500, 600],
        }
    )
    df.pivot(index="col1", columns="col3", values="col2")
    df.pivot(index="col1", columns="col3")
    df.pivot(index="col1", columns="col3", values=["col2", "col4"])


def test_types_groupby() -> None:
    df = pd.DataFrame(data={"col1": [1, 1, 2], "col2": [3, 4, 5], "col3": [0, 1, 0]})
    df.index.name = "ind"
    df.groupby(by="col1")
    df.groupby(level="ind")
    df.groupby(by="col1", sort=False, as_index=True)
    df.groupby(by=["col1", "col2"])

    df1: pd.DataFrame = df.groupby(by="col1").agg("sum")
    df2: pd.DataFrame = df.groupby(level="ind").aggregate("sum")
    df3: pd.DataFrame = df.groupby(by="col1", sort=False, as_index=True).transform(
        lambda x: x.max()
    )
    # error: Incompatible types in assignment (expression has type "Union[Series,
    # DataFrame]", variable has type "DataFrame")
    df4: pd.DataFrame = df.groupby(  # type: ignore[assignment]
        by=["col1", "col2"]
    ).count()
    df5: pd.DataFrame = df.groupby(by=["col1", "col2"]).filter(lambda x: x["col1"] > 0)
    df6: pd.DataFrame = df.groupby(by=["col1", "col2"]).nunique()


# This was added in 1.1.0 https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
def test_types_group_by_with_dropna_keyword() -> None:
    df = pd.DataFrame(
        data={"col1": [1, 1, 2, 1], "col2": [2, None, 1, 2], "col3": [3, 4, 3, 2]}
    )
    df.groupby(by="col2", dropna=True).sum()
    df.groupby(by="col2", dropna=False).sum()
    df.groupby(by="col2").sum()


def test_types_merge() -> None:
    df = pd.DataFrame(data={"col1": [1, 1, 2], "col2": [3, 4, 5]})
    df2 = pd.DataFrame(data={"col1": [1, 1, 2], "col2": [0, 1, 0]})
    df.merge(df2)
    df.merge(df2, on="col1")
    df.merge(df2, on="col1", how="left")
    df.merge(df2, on=["col1", "col2"], how="left")
    df.merge(df2, on=("col1", "col2"), how="left")
    l: List[str] = ["col1", "col2"]
    df.merge(df2, on=l)


@td.skip_if_no("matplotlib")  # type: ignore[misc]
def test_types_plot() -> None:
    df = pd.DataFrame(data={"col1": [1, 1, 2], "col2": [3, 4, 5]})
    df.plot.hist()
    df.plot.scatter(x="col2", y="col1")


def test_types_window() -> None:
    df = pd.DataFrame(data={"col1": [1, 1, 2], "col2": [3, 4, 5]})
    df.expanding()
    df.expanding(axis=1, center=True)

    df.rolling(2)
    df.rolling(2, axis=1, center=True)


def test_types_cov() -> None:
    df = pd.DataFrame(data={"col1": [1, 1, 2], "col2": [3, 4, 5]})
    df.cov()
    df.cov(min_periods=1)
    # ddof param was added in 1.1.0
    # https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
    df.cov(ddof=2)


def test_types_to_numpy() -> None:
    df = pd.DataFrame(data={"col1": [1, 1, 2], "col2": [3, 4, 5]})
    df.to_numpy()
    df.to_numpy(dtype="str", copy=True)
    # na_value param was added in 1.1.0
    # https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
    df.to_numpy(na_value=0)


# error: Untyped decorator makes function "test_types_to_feather" untyped
@td.skip_if_no("tabulate")  # type: ignore[misc]
def test_to_markdown() -> None:
    df = pd.DataFrame(data={"col1": [1, 1, 2], "col2": [3, 4, 5]})
    df.to_markdown()
    df.to_markdown(buf=None, mode="wt")
    # index param was added in 1.1.0
    # https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
    df.to_markdown(index=False)


# error: Untyped decorator makes function "test_types_to_feather" untyped
@td.skip_if_no("pyarrow")  # type: ignore[misc]
def test_types_to_feather() -> None:
    df = pd.DataFrame(data={"col1": [1, 1, 2], "col2": [3, 4, 5]})
    df.to_feather("dummy_path")
    # kwargs for pyarrow.feather.write_feather added in 1.1.0
    # https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
    df.to_feather(
        "dummy_path",
        compression="zstd",
        compression_level=3,
        chunksize=2,
    )

    # to_feather has been able to accept a buffer since pandas 1.0.0
    # See https://pandas.pydata.org/docs/whatsnew/v1.0.0.html
    # Docstring and type were updated in 1.2.0.
    # https://github.com/pandas-dev/pandas/pull/35408
    with tempfile.TemporaryFile() as f:
        df.to_feather(f)


# compare() method added in 1.1.0
# https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
def test_types_compare() -> None:
    df1 = pd.DataFrame(
        data={"col1": [1, 1, 2, 1], "col2": [2, None, 1, 2], "col3": [3, 4, 3, 2]}
    )
    df2 = pd.DataFrame(
        data={"col1": [1, 2, 5, 6], "col2": [3, 4, 1, 1], "col3": [3, 4, 3, 2]}
    )
    df1.compare(df2)
    df2.compare(df1, align_axis=0, keep_shape=True, keep_equal=True)


def test_types_agg() -> None:
    df = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]], columns=["A", "B", "C"])
    df.agg("min")
    df.agg(x=("A", max), y=("B", "min"), z=("C", np.mean))
    df.agg("mean", axis=1)


def test_types_describe() -> None:
    df = pd.DataFrame(
        data={
            "col1": [1, 2, -4],
            "col2": [
                np.datetime64("2000-01-01"),
                np.datetime64("2010-01-01"),
                np.datetime64("2010-01-01"),
            ],
        }
    )
    df.describe()
    df.describe(percentiles=[0.5], include="all")
    df.describe(exclude=np.number)
    # datetime_is_numeric param added in 1.1.0
    # https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
    df.describe(datetime_is_numeric=True)


def test_types_to_string() -> None:
    df = pd.DataFrame(
        data={
            "col1": [1, None, -4],
            "col2": [
                np.datetime64("2000-01-01"),
                np.datetime64("2010-01-01"),
                np.datetime64("2010-01-01"),
            ],
        }
    )
    df.to_string(
        index=True,
        col_space=2,
        header=["a", "b"],
        na_rep="0",
        justify="left",
        max_rows=2,
        min_rows=0,
        max_cols=2,
        show_dimensions=True,
        line_width=3,
    )
    # col_space accepting list or dict added in 1.1.0
    # https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
    df.to_string(col_space=[1, 2])
    df.to_string(col_space={"col1": 1, "col2": 3})


def test_types_to_html() -> None:
    df = pd.DataFrame(
        data={
            "col1": [1, None, -4],
            "col2": [
                np.datetime64("2000-01-01"),
                np.datetime64("2010-01-01"),
                np.datetime64("2010-01-01"),
            ],
        }
    )
    df.to_html(
        index=True,
        col_space=2,
        header=["a", "b"],
        na_rep="0",
        justify="left",
        max_rows=2,
        max_cols=2,
        show_dimensions=True,
    )
    # col_space accepting list or dict added in 1.1.0
    # https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
    df.to_html(col_space=[1, 2])
    df.to_html(col_space={"col1": 1, "col2": 3})


def test_types_resample() -> None:
    df = pd.DataFrame({"values": [2, 11, 3, 13, 14, 18, 17, 19]})
    df["date"] = pd.date_range("01/01/2018", periods=8, freq="W")
    df.resample("M", on="date")
    # origin and offset params added in 1.1.0
    # https://pandas.pydata.org/docs/whatsnew/v1.1.0.html
    df.resample("20min", origin="epoch", offset=pd.Timedelta(2, "minutes"), on="date")


def test_types_from_dict() -> None:
    pd.DataFrame.from_dict({"col_1": [3, 2, 1, 0], "col_2": ["a", "b", "c", "d"]})
    pd.DataFrame.from_dict({1: [3, 2, 1, 0], 2: ["a", "b", "c", "d"]})
    pd.DataFrame.from_dict({"a": {1: 2}, "b": {3: 4, 1: 4}}, orient="index")
    pd.DataFrame.from_dict({"a": {"row1": 2}, "b": {"row2": 4, "row1": 4}})
    pd.DataFrame.from_dict({"a": (1, 2, 3), "b": (2, 4, 5)})
    pd.DataFrame.from_dict(
        data={"col_1": {"a": 1}, "col_2": {"a": 1, "b": 2}}, orient="columns"
    )


@td.skip_if_no("jinja")  # type: ignore[misc]
def test_pipe() -> None:
    def foo(df: pd.DataFrame) -> pd.DataFrame:
        return df

    df1: pd.DataFrame = pd.DataFrame({"a": [1]}).pipe(foo)

    df2: pd.DataFrame = (
        pd.DataFrame(
            {
                "price": [10, 11, 9, 13, 14, 18, 17, 19],
                "volume": [50, 60, 40, 100, 50, 100, 40, 50],
            }
        )
        .assign(week_starting=pd.date_range("01/01/2018", periods=8, freq="W"))
        .resample("M", on="week_starting")
        .pipe(foo)
    )

    df3: pd.DataFrame = pd.DataFrame({"a": [1], "b": [1]}).groupby("a").pipe(foo)

    df4: pd.DataFrame = pd.DataFrame({"a": [1], "b": [1]}).style.pipe(foo)


# set_flags() method added in 1.2.0
# https://pandas.pydata.org/docs/whatsnew/v1.2.0.html
def test_types_set_flags() -> None:
    pd.DataFrame([[1, 2], [8, 9]], columns=["A", "B"]).set_flags(
        allows_duplicate_labels=False
    )
    pd.DataFrame([[1, 2], [8, 9]], columns=["A", "A"]).set_flags(
        allows_duplicate_labels=True
    )
    pd.DataFrame([[1, 2], [8, 9]], columns=["A", "A"])


# error: Untyped decorator makes function "test_types_to_parquet" untyped
@td.skip_if_no("pyarrow")  # type: ignore[misc]
def test_types_to_parquet() -> None:
    df = pd.DataFrame([[1, 2], [8, 9]], columns=["A", "B"]).set_flags(
        allows_duplicate_labels=False
    )
    with tempfile.NamedTemporaryFile() as file:
        df.to_parquet(Path(file.name))
    # to_parquet() returns bytes when no path given since 1.2.0
    # https://pandas.pydata.org/docs/whatsnew/v1.2.0.html
    # error: Incompatible types in assignment (expression has type "Optional[bytes]",
    # variable has type "bytes")
    b: bytes = df.to_parquet()  # type: ignore[assignment]


def test_types_to_latex() -> None:
    df = pd.DataFrame([[1, 2], [8, 9]], columns=["A", "B"])
    df.to_latex(
        columns=["A"], label="some_label", caption="some_caption", multirow=True
    )
    df.to_latex(escape=False, decimal=",", column_format="r")
    # position param was added in 1.2.0
    # https://pandas.pydata.org/docs/whatsnew/v1.2.0.html
    df.to_latex(position="some")
    # caption param was extended to accept tuple in 1.2.0
    # https://pandas.pydata.org/docs/whatsnew/v1.2.0.html
    df.to_latex(caption=("cap1", "cap2"))


def test_types_explode() -> None:
    df = pd.DataFrame([[1, 2], [8, 9]], columns=["A", "B"])
    res1: pd.DataFrame = df.explode("A")
    res2: pd.DataFrame = df.explode("A", ignore_index=False)
    res3: pd.DataFrame = df.explode("A", ignore_index=True)


def test_types_rename() -> None:
    df = pd.DataFrame(columns=["a"])
    col_map = {"a": "b"}
    # error: Argument "columns" to "rename" of "DataFrame" has incompatible type
    # "Dict[str, str]"; expected "Optional[Union[Mapping[Hashable, Any],
    # Callable[[Hashable], Hashable]]]"
    df.rename(columns=col_map)  # type: ignore[arg-type]
    df.rename(columns={"a": "b"})
    df.rename(columns={1: "b"})
    # Apparently all of these calls are accepted by pandas
    df.rename(columns={None: "b"})
    df.rename(columns={"": "b"})
    df.rename(columns={(2, 1): "b"})


def test_types_eq() -> None:
    df1 = pd.DataFrame([[1, 2], [8, 9]], columns=["A", "B"])
    res1: pd.DataFrame = df1 == 1
    df2 = pd.DataFrame([[1, 2], [8, 9]], columns=["A", "B"])
    res2: pd.DataFrame = df1 == df2


def test_types_as_type() -> None:
    df1 = pd.DataFrame([[1, 2], [8, 9]], columns=["A", "B"])
    df2: pd.DataFrame = df1.astype({"A": "int32"})


def test_types_dot() -> None:
    df1 = pd.DataFrame([[0, 1, -2, -1], [1, 1, 1, 1]])
    df2 = pd.DataFrame([[0, 1], [1, 2], [-1, -1], [2, 0]])
    s1 = pd.Series([1, 1, 2, 1])
    np_array = np.array([[0, 1], [1, 2], [-1, -1], [2, 0]])
    # error: Incompatible types in assignment (expression has type "Union[DataFrame,
    # Series]", variable has type "DataFrame")
    df3: pd.DataFrame = df1 @ df2  # type: ignore[assignment]
    df4: pd.DataFrame = df1.dot(df2)
    # error: Incompatible types in assignment (expression has type "Union[DataFrame,
    # Series]", variable has type "DataFrame")
    df5: pd.DataFrame = df1 @ np_array  # type: ignore[assignment]
    df6: pd.DataFrame = df1.dot(np_array)
    df7: pd.Series = df1 @ s1
    df8: pd.Series = df1.dot(s1)
