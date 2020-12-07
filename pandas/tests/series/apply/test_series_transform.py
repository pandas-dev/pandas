import numpy as np
import pytest

from pandas import DataFrame, Series, concat
import pandas._testing as tm
from pandas.core.base import SpecificationError
from pandas.core.groupby.base import transformation_kernels

# tshift only works on time index and is deprecated
# There is no Series.cumcount
series_kernels = [
    x for x in sorted(transformation_kernels) if x not in ["tshift", "cumcount"]
]


@pytest.mark.parametrize("op", series_kernels)
def test_transform_groupby_kernel(string_series, op):
    # GH 35964

    args = [0.0] if op == "fillna" else []
    ones = np.ones(string_series.shape[0])
    expected = string_series.groupby(ones).transform(op, *args)
    result = string_series.transform(op, 0, *args)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "ops, names",
    [
        ([np.sqrt], ["sqrt"]),
        ([np.abs, np.sqrt], ["absolute", "sqrt"]),
        (np.array([np.sqrt]), ["sqrt"]),
        (np.array([np.abs, np.sqrt]), ["absolute", "sqrt"]),
    ],
)
def test_transform_listlike(string_series, ops, names):
    # GH 35964
    with np.errstate(all="ignore"):
        expected = concat([op(string_series) for op in ops], axis=1)
        expected.columns = names
        result = string_series.transform(ops)
        tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("box", [dict, Series])
def test_transform_dictlike(string_series, box):
    # GH 35964
    with np.errstate(all="ignore"):
        expected = concat([np.sqrt(string_series), np.abs(string_series)], axis=1)
    expected.columns = ["foo", "bar"]
    result = string_series.transform(box({"foo": np.sqrt, "bar": np.abs}))
    tm.assert_frame_equal(result, expected)


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


def test_transform_none_to_type():
    # GH34377
    df = DataFrame({"a": [None]})
    msg = "Transform function failed"
    with pytest.raises(ValueError, match=msg):
        df.transform({"a": int})


def test_transform_axis_1_raises():
    # GH 35964
    msg = "No axis named 1 for object type Series"
    with pytest.raises(ValueError, match=msg):
        Series([1]).transform("sum", axis=1)


def test_transform_nested_renamer():
    # GH 35964
    match = "nested renamer is not supported"
    with pytest.raises(SpecificationError, match=match):
        Series([1]).transform({"A": {"B": ["sum"]}})
