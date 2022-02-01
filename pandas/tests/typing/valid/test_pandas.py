"""
Copyright (c) Virtus Lab sp. z o.o. (Ltd.)

Distributed under the terms of the MIT license.

The full license is in the STUBS_LICENSE file, distributed with this software.
"""
# flake8: noqa: F841
# TODO: many functions need return types annotations for pyright
# to run with reportGeneralTypeIssues = true
from typing import (
    Any,
    Dict,
    List,
    Union,
)

import pandas as pd
import pandas._testing as tm
from pandas.io.parsers import TextFileReader


def test_types_to_datetime() -> None:
    df = pd.DataFrame({"year": [2015, 2016], "month": [2, 3], "day": [4, 5]})
    # error: No overload variant of "to_datetime" matches argument type "DataFrame"
    pd.to_datetime(df)  # type: ignore[call-overload]
    # error: No overload variant of "to_datetime" matches argument types "DataFrame",
    # "str", "str", "bool"
    pd.to_datetime(  # type: ignore[call-overload]
        df, unit="s", origin="unix", infer_datetime_format=True
    )
    # error: No overload variant of "to_datetime" matches argument types "DataFrame",
    # "str", "bool", "None", "str", "bool"
    pd.to_datetime(  # type: ignore[call-overload]
        df, unit="ns", dayfirst=True, utc=None, format="%M:%D", exact=False
    )
    pd.to_datetime([1, 2], unit="D", origin=pd.Timestamp("01/01/2000"))
    pd.to_datetime([1, 2], unit="D", origin=3)


def test_types_concat() -> None:
    s = pd.Series([0, 1, -10])
    s2 = pd.Series([7, -5, 10])

    pd.concat([s, s2])
    pd.concat([s, s2], axis=1)
    pd.concat([s, s2], keys=["first", "second"], sort=True)
    pd.concat([s, s2], keys=["first", "second"], names=["source", "row"])

    # Depends on the axis
    # error: Argument 1 to "concat" has incompatible type "Dict[str, Series]"; expected
    # "Union[Iterable[DataFrame], Mapping[Hashable, DataFrame]]"
    rs1: Union[pd.Series, pd.DataFrame] = pd.concat(
        {"a": s, "b": s2}  # type:ignore[arg-type]
    )
    # error: Argument 1 to "concat" has incompatible type "Dict[str, Series]"; expected
    # "Union[Iterable[NDFrame], Mapping[Hashable, NDFrame]]"
    rs1a: Union[pd.Series, pd.DataFrame] = pd.concat(
        {"a": s, "b": s2}, axis=1  # type:ignore[arg-type]
    )
    # error: Argument 1 to "concat" has incompatible type "Dict[int, Series]"; expected
    # "Union[Iterable[DataFrame], Mapping[Hashable, DataFrame]]"
    rs2: Union[pd.Series, pd.DataFrame] = pd.concat(
        {1: s, 2: s2}  # type:ignore[arg-type]
    )
    # error: Argument 1 to "concat" has incompatible type "Dict[int, Series]"; expected
    # "Union[Iterable[NDFrame], Mapping[Hashable, NDFrame]]"
    rs2a: Union[pd.Series, pd.DataFrame] = pd.concat(
        {1: s, 2: s2}, axis=1  # type:ignore[arg-type]
    )
    # error: Argument 1 to "concat" has incompatible type "Dict[Optional[int], Series]";
    # expected "Union[Iterable[DataFrame], Mapping[Hashable, DataFrame]]"
    rs3: Union[pd.Series, pd.DataFrame] = pd.concat(
        {1: s, None: s2}  # type:ignore[arg-type]
    )
    # error: Argument 1 to "concat" has incompatible type "Dict[Optional[int], Series]";
    # expected "Union[Iterable[NDFrame], Mapping[Hashable, NDFrame]]"
    rs3a: Union[pd.Series, pd.DataFrame] = pd.concat(
        {1: s, None: s2}, axis=1  # type:ignore[arg-type]
    )

    df = pd.DataFrame(data={"col1": [1, 2], "col2": [3, 4]})
    df2 = pd.DataFrame(data={"col1": [10, 20], "col2": [30, 40]})

    pd.concat([df, df2])
    pd.concat([df, df2], axis=1)
    pd.concat([df, df2], keys=["first", "second"], sort=True)
    pd.concat([df, df2], keys=["first", "second"], names=["source", "row"])

    # error: Incompatible types in assignment (expression has type "Union[DataFrame,
    # Series]", variable has type "DataFrame")
    # error: Argument 1 to "concat" has incompatible type "Dict[str, DataFrame]";
    # expected "Union[Iterable[NDFrame], Mapping[Hashable, NDFrame]]"
    result: pd.DataFrame = pd.concat(  # type: ignore[assignment]
        {
            "a": pd.DataFrame([1, 2, 3]),
            "b": pd.DataFrame([4, 5, 6]),
        },  # type:ignore[arg-type]
        axis=1,
    )
    # error: Argument 1 to "concat" has incompatible type "Dict[str, Series]"; expected
    # "Union[Iterable[NDFrame], Mapping[Hashable, NDFrame]]"
    result2: Union[pd.DataFrame, pd.Series] = pd.concat(
        {
            "a": pd.Series([1, 2, 3]),
            "b": pd.Series([4, 5, 6]),
        },  # type:ignore[arg-type]
        axis=1,
    )

    # error: Argument 1 to "concat" has incompatible type "Dict[str, DataFrame]";
    # expected "Union[Iterable[DataFrame], Mapping[Hashable, DataFrame]]"
    rdf1: pd.DataFrame = pd.concat({"a": df, "b": df2})  # type:ignore[arg-type]
    # error: Argument 1 to "concat" has incompatible type "Dict[int, DataFrame]";
    # expected "Union[Iterable[DataFrame], Mapping[Hashable, DataFrame]]"
    rdf2: pd.DataFrame = pd.concat({1: df, 2: df2})  # type:ignore[arg-type]
    # error: Argument 1 to "concat" has incompatible type "Dict[Optional[int],
    # DataFrame]"; expected "Union[Iterable[DataFrame], Mapping[Hashable, DataFrame]]"
    rdf3: pd.DataFrame = pd.concat({1: df, None: df2})  # type:ignore[arg-type]

    rdf4: pd.DataFrame = pd.concat(list(map(lambda x: s2, ["some_value", 3])), axis=1)


def test_types_json_normalize() -> None:
    data1: List[Dict[str, Any]] = [
        {"id": 1, "name": {"first": "Coleen", "last": "Volk"}},
        {"name": {"given": "More", "family": "Regner"}},
        {"id": 2, "name": "Faye Raker"},
    ]
    df1: pd.DataFrame = pd.json_normalize(data=data1)
    df2: pd.DataFrame = pd.json_normalize(data=data1, max_level=0, sep=";")
    df3: pd.DataFrame = pd.json_normalize(
        data=data1, meta_prefix="id", record_prefix="name", errors="raise"
    )
    df4: pd.DataFrame = pd.json_normalize(data=data1, record_path=None, meta="id")
    data2: Dict[str, Any] = {"name": {"given": "More", "family": "Regner"}}
    df5: pd.DataFrame = pd.json_normalize(data=data2)


def test_types_read_csv() -> None:
    df = pd.DataFrame(data={"col1": [1, 2], "col2": [3, 4]})
    # error: Incompatible types in assignment (expression has type "Optional[str]",
    # variable has type "str")
    csv_df: str = df.to_csv()  # type: ignore[assignment]

    with tm.ensure_clean() as path:
        df.to_csv(path)
        df2: pd.DataFrame = pd.read_csv(path)
        df3: pd.DataFrame = pd.read_csv(path, sep="a", squeeze=False)
        df4: pd.DataFrame = pd.read_csv(
            path,
            header=None,
            prefix="b",
            mangle_dupe_cols=True,
            keep_default_na=False,
        )
        df5: pd.DataFrame = pd.read_csv(
            path, engine="python", true_values=[0, 1, 3], na_filter=False
        )
        df6: pd.DataFrame = pd.read_csv(
            path,
            skiprows=lambda x: x in [0, 2],
            skip_blank_lines=True,
            dayfirst=False,
        )
        df7: pd.DataFrame = pd.read_csv(path, nrows=2)
        tfr1: TextFileReader = pd.read_csv(path, nrows=2, iterator=True, chunksize=3)
        tfr2: TextFileReader = pd.read_csv(path, nrows=2, chunksize=1)
        tfr3: TextFileReader = pd.read_csv(path, nrows=2, iterator=False, chunksize=1)
        tfr4: TextFileReader = pd.read_csv(path, nrows=2, iterator=True)
