# Tests specifically aimed at detecting bad arguments.
# This file is organized by reason for exception.
#     1. always invalid argument values
#     2. missing column(s)
#     3. incompatible ops/dtype/args/kwargs
#     4. invalid result shape/type
# If your test does not fit into one of these categories, add to this list.

import re

import numpy as np
import pytest

from pandas import (
    DataFrame,
    Series,
    date_range,
    notna,
)
import pandas._testing as tm
from pandas.core.base import SpecificationError


@pytest.mark.parametrize("result_type", ["foo", 1])
def test_result_type_error(result_type, int_frame_const_col):
    # allowed result_type
    df = int_frame_const_col

    msg = (
        "invalid value for result_type, must be one of "
        "{None, 'reduce', 'broadcast', 'expand'}"
    )
    with pytest.raises(ValueError, match=msg):
        df.apply(lambda x: [1, 2, 3], axis=1, result_type=result_type)


def test_agg_raises():
    # GH 26513
    df = DataFrame({"A": [0, 1], "B": [1, 2]})
    msg = "Must provide"

    with pytest.raises(TypeError, match=msg):
        df.agg()


@pytest.mark.parametrize("box", [DataFrame, Series])
@pytest.mark.parametrize("method", ["apply", "agg", "transform"])
@pytest.mark.parametrize("func", [{"A": {"B": "sum"}}, {"A": {"B": ["sum"]}}])
def test_nested_renamer(box, method, func):
    # GH 35964
    obj = box({"A": [1]})
    match = "nested renamer is not supported"
    with pytest.raises(SpecificationError, match=match):
        getattr(obj, method)(func)


def test_transform_nested_renamer():
    # GH 35964
    match = "nested renamer is not supported"
    with pytest.raises(SpecificationError, match=match):
        Series([1]).transform({"A": {"B": ["sum"]}})


def test_agg_dict_nested_renaming_depr_agg():

    df = DataFrame({"A": range(5), "B": 5})

    # nested renaming
    msg = r"nested renamer is not supported"
    with pytest.raises(SpecificationError, match=msg):
        df.agg({"A": {"foo": "min"}, "B": {"bar": "max"}})


def test_agg_dict_nested_renaming_depr_transform():
    df = DataFrame({"A": range(5), "B": 5})

    # nested renaming
    msg = r"nested renamer is not supported"
    with pytest.raises(SpecificationError, match=msg):
        # mypy identifies the argument as an invalid type
        df.transform({"A": {"foo": "min"}, "B": {"bar": "max"}})


def test_apply_dict_depr():

    tsdf = DataFrame(
        np.random.randn(10, 3),
        columns=["A", "B", "C"],
        index=date_range("1/1/2000", periods=10),
    )
    msg = "nested renamer is not supported"
    with pytest.raises(SpecificationError, match=msg):
        tsdf.A.agg({"foo": ["sum", "mean"]})


@pytest.mark.parametrize("method", ["apply", "agg", "transform"])
@pytest.mark.parametrize("func", [{"B": "sum"}, {"B": ["sum"]}])
def test_missing_column(method, func):
    # GH 40004
    obj = DataFrame({"A": [1]})
    match = re.escape("Column(s) ['B'] do not exist")
    with pytest.raises(KeyError, match=match):
        getattr(obj, method)(func)


def test_transform_missing_columns(axis):
    # GH#35964
    df = DataFrame({"A": [1, 2], "B": [3, 4]})
    match = re.escape("Column(s) ['C'] do not exist")
    with pytest.raises(KeyError, match=match):
        df.transform({"C": "cumsum"})


def test_transform_mixed_column_name_dtypes():
    # GH39025
    df = DataFrame({"a": ["1"]})
    msg = r"Column\(s\) \[1, 'b'\] do not exist"
    with pytest.raises(KeyError, match=msg):
        df.transform({"a": int, 1: str, "b": int})


@pytest.mark.parametrize(
    "how, args", [("pct_change", ()), ("nsmallest", (1, ["a", "b"])), ("tail", 1)]
)
def test_apply_str_axis_1_raises(how, args):
    # GH 39211 - some ops don't support axis=1
    df = DataFrame({"a": [1, 2], "b": [3, 4]})
    msg = f"Operation {how} does not support axis=1"
    with pytest.raises(ValueError, match=msg):
        df.apply(how, axis=1, args=args)


def test_transform_axis_1_raises():
    # GH 35964
    msg = "No axis named 1 for object type Series"
    with pytest.raises(ValueError, match=msg):
        Series([1]).transform("sum", axis=1)


def test_apply_modify_traceback():
    data = DataFrame(
        {
            "A": [
                "foo",
                "foo",
                "foo",
                "foo",
                "bar",
                "bar",
                "bar",
                "bar",
                "foo",
                "foo",
                "foo",
            ],
            "B": [
                "one",
                "one",
                "one",
                "two",
                "one",
                "one",
                "one",
                "two",
                "two",
                "two",
                "one",
            ],
            "C": [
                "dull",
                "dull",
                "shiny",
                "dull",
                "dull",
                "shiny",
                "shiny",
                "dull",
                "shiny",
                "shiny",
                "shiny",
            ],
            "D": np.random.randn(11),
            "E": np.random.randn(11),
            "F": np.random.randn(11),
        }
    )

    data.loc[4, "C"] = np.nan

    def transform(row):
        if row["C"].startswith("shin") and row["A"] == "foo":
            row["D"] = 7
        return row

    def transform2(row):
        if notna(row["C"]) and row["C"].startswith("shin") and row["A"] == "foo":
            row["D"] = 7
        return row

    msg = "'float' object has no attribute 'startswith'"
    with pytest.raises(AttributeError, match=msg):
        data.apply(transform, axis=1)


@pytest.mark.parametrize(
    "df, func, expected",
    tm.get_cython_table_params(
        DataFrame([["a", "b"], ["b", "a"]]), [["cumprod", TypeError]]
    ),
)
def test_agg_cython_table_raises(df, func, expected, axis):
    # GH 21224
    msg = "can't multiply sequence by non-int of type 'str'"
    with pytest.raises(expected, match=msg):
        df.agg(func, axis=axis)


def test_transform_none_to_type():
    # GH#34377
    df = DataFrame({"a": [None]})
    msg = "Transform function failed"
    with pytest.raises(ValueError, match=msg):
        df.transform({"a": int})


def test_apply_broadcast_error(int_frame_const_col):
    df = int_frame_const_col

    # > 1 ndim
    msg = "too many dims to broadcast"
    with pytest.raises(ValueError, match=msg):
        df.apply(
            lambda x: np.array([1, 2]).reshape(-1, 2),
            axis=1,
            result_type="broadcast",
        )

    # cannot broadcast
    msg = "cannot broadcast result"
    with pytest.raises(ValueError, match=msg):
        df.apply(lambda x: [1, 2], axis=1, result_type="broadcast")

    with pytest.raises(ValueError, match=msg):
        df.apply(lambda x: Series([1, 2]), axis=1, result_type="broadcast")


def test_transform_and_agg_err_agg(axis, float_frame):
    # cannot both transform and agg
    msg = "cannot combine transform and aggregation operations"
    with pytest.raises(ValueError, match=msg):
        with np.errstate(all="ignore"):
            float_frame.agg(["max", "sqrt"], axis=axis)

    df = DataFrame({"A": range(5), "B": 5})

    def f():
        with np.errstate(all="ignore"):
            df.agg({"A": ["abs", "sum"], "B": ["mean", "max"]}, axis=axis)


def test_transform_and_agg_error_agg(string_series):
    # we are trying to transform with an aggregator
    msg = "cannot combine transform and aggregation"
    with pytest.raises(ValueError, match=msg):
        with np.errstate(all="ignore"):
            string_series.agg(["sqrt", "max"])

    msg = "cannot perform both aggregation and transformation"
    with pytest.raises(ValueError, match=msg):
        with np.errstate(all="ignore"):
            string_series.agg({"foo": np.sqrt, "bar": "sum"})


def test_transform_and_agg_err_transform(axis, float_frame):
    # GH 35964
    # cannot both transform and agg
    msg = "Function did not transform"
    with pytest.raises(ValueError, match=msg):
        float_frame.transform(["max", "min"], axis=axis)

    msg = "Function did not transform"
    with pytest.raises(ValueError, match=msg):
        float_frame.transform(["max", "sqrt"], axis=axis)


def test_transform_reducer_raises(all_reductions, frame_or_series):
    # GH 35964
    op = all_reductions

    obj = DataFrame({"A": [1, 2, 3]})
    if frame_or_series is not DataFrame:
        obj = obj["A"]

    msg = "Function did not transform"
    with pytest.raises(ValueError, match=msg):
        obj.transform(op)
    with pytest.raises(ValueError, match=msg):
        obj.transform([op])
    with pytest.raises(ValueError, match=msg):
        obj.transform({"A": op})
    with pytest.raises(ValueError, match=msg):
        obj.transform({"A": [op]})


def test_transform_wont_agg(string_series):
    # GH 35964
    # we are trying to transform with an aggregator
    msg = "Function did not transform"
    with pytest.raises(ValueError, match=msg):
        string_series.transform(["min", "max"])

    msg = "Function did not transform"
    with pytest.raises(ValueError, match=msg):
        with np.errstate(all="ignore"):
            string_series.transform(["sqrt", "max"])
