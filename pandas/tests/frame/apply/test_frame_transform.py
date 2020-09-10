import operator
import re

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.tests.frame.common import zip_frames


def test_agg_transform(axis, float_frame):
    other_axis = 1 if axis in {0, "index"} else 0

    with np.errstate(all="ignore"):

from pandas import DataFrame, MultiIndex
import pandas._testing as tm
from pandas.core.base import SpecificationError
from pandas.core.groupby.base import transformation_kernels
from pandas.tests.frame.common import zip_frames


def test_transform(axis, float_frame):
    # GH 35964
    other_axis = 1 if axis in {0, "index"} else 0

    with np.errstate(all="ignore"):
        f_abs = np.abs(float_frame)
        f_sqrt = np.sqrt(float_frame)

        # ufunc
        result = float_frame.transform(np.sqrt, axis=axis)
        expected = f_sqrt.copy()
        tm.assert_frame_equal(result, expected)

        result = float_frame.transform(np.sqrt, axis=axis)
        tm.assert_frame_equal(result, expected)

        # list-like
        expected = f_sqrt.copy()
        if axis in {0, "index"}:
            expected.columns = pd.MultiIndex.from_product(
                [float_frame.columns, ["sqrt"]]
            )
        else:
            expected.index = pd.MultiIndex.from_product([float_frame.index, ["sqrt"]])
        result = float_frame.transform([np.sqrt], axis=axis)
        tm.assert_frame_equal(result, expected)

        # multiple items in list
        # these are in the order as if we are applying both
        # functions per series and then concatting
        expected = zip_frames([f_abs, f_sqrt], axis=other_axis)
        if axis in {0, "index"}:
            expected.columns = pd.MultiIndex.from_product(
                [float_frame.columns, ["absolute", "sqrt"]]
            )
        else:
            expected.index = pd.MultiIndex.from_product(
                [float_frame.index, ["absolute", "sqrt"]]
            )
        result = float_frame.transform([np.abs, "sqrt"], axis=axis)
        tm.assert_frame_equal(result, expected)


def test_transform_and_agg_err(axis, float_frame):
    # cannot both transform and agg
    msg = "transforms cannot produce aggregated results"
    with pytest.raises(ValueError, match=msg):
        float_frame.transform(["max", "min"], axis=axis)

    msg = "cannot combine transform and aggregation operations"
    with pytest.raises(ValueError, match=msg):
        with np.errstate(all="ignore"):
            float_frame.transform(["max", "sqrt"], axis=axis)


@pytest.mark.parametrize("method", ["abs", "shift", "pct_change", "cumsum", "rank"])
def test_transform_method_name(method):
    # GH 19760
    df = DataFrame({"A": [-1, 2]})
    result = df.transform(method)
    expected = operator.methodcaller(method)(df)
    tm.assert_frame_equal(result, expected)


def test_transform_and_agg_err(axis, float_frame):
    # GH 35964
    # cannot both transform and agg
    msg = "Function did not transform"
    with pytest.raises(ValueError, match=msg):
        float_frame.transform(["max", "min"], axis=axis)

    msg = "Function did not transform"
    with pytest.raises(ValueError, match=msg):
        float_frame.transform(["max", "sqrt"], axis=axis)


def test_agg_dict_nested_renaming_depr():
    df = DataFrame({"A": range(5), "B": 5})

    # nested renaming
    msg = r"nested renamer is not supported"
    with pytest.raises(SpecificationError, match=msg):
        df.transform({"A": {"foo": "min"}, "B": {"bar": "max"}})


def test_transform_reducer_raises(all_reductions):
    # GH 35964
    op = all_reductions
    df = DataFrame({"A": [1, 2, 3]})
    msg = "Function did not transform"
    with pytest.raises(ValueError, match=msg):
        df.transform(op)
    with pytest.raises(ValueError, match=msg):
        df.transform([op])
    with pytest.raises(ValueError, match=msg):
        df.transform({"A": op})
    with pytest.raises(ValueError, match=msg):
        df.transform({"A": [op]})


# mypy doesn't allow adding lists of different types
# https://github.com/python/mypy/issues/5492
@pytest.mark.parametrize("op", [*transformation_kernels, lambda x: x + 1])
def test_transform_bad_dtype(op):
    # GH 35964
    df = DataFrame({"A": 3 * [object]})  # DataFrame that will fail on most transforms
    if op in ("backfill", "shift", "pad", "bfill", "ffill"):
        pytest.xfail("Transform function works on any datatype")
    msg = "Transform function failed"
    with pytest.raises(ValueError, match=msg):
        df.transform(op)
    with pytest.raises(ValueError, match=msg):
        df.transform([op])
    with pytest.raises(ValueError, match=msg):
        df.transform({"A": op})
    with pytest.raises(ValueError, match=msg):
        df.transform({"A": [op]})


@pytest.mark.parametrize("op", transformation_kernels)
def test_transform_multi_dtypes(op):
    # GH 35964
    df = DataFrame({"A": ["a", "b", "c"], "B": [1, 2, 3]})

    # Determine which columns op will work on
    columns = []
    for column in df:
        try:
            df[column].transform(op)
            columns.append(column)
        except Exception:
            pass

    if len(columns) > 0:
        expected = df[columns].transform([op])
        result = df.transform([op])
        tm.assert_equal(result, expected)

        expected = df[columns].transform({column: op for column in columns})
        result = df.transform({column: op for column in columns})
        tm.assert_equal(result, expected)

        expected = df[columns].transform({column: [op] for column in columns})
        result = df.transform({column: [op] for column in columns})
        tm.assert_equal(result, expected)


@pytest.mark.parametrize("use_apply", [True, False])
def test_transform_passes_args(use_apply):
    # GH 35964
    # transform uses UDF either via apply or passing the entire DataFrame
    expected_args = [1, 2]
    expected_kwargs = {"c": 3}

    def f(x, a, b, c):
        # transform is using apply iff x is not a DataFrame
        if use_apply == isinstance(x, DataFrame):
            # Force transform to fallback
            raise ValueError
        assert [a, b] == expected_args
        assert c == expected_kwargs["c"]
        return x

    DataFrame([1]).transform(f, 0, *expected_args, **expected_kwargs)


def test_transform_missing_columns(axis):
    # GH 35964
    df = DataFrame({"A": [1, 2], "B": [3, 4]})
    match = re.escape("Column(s) ['C'] do not exist")
    with pytest.raises(SpecificationError, match=match):
        df.transform({"C": "cumsum"})
