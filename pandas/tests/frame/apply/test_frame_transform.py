import operator
import re

import numpy as np
import pytest

from pandas import DataFrame, MultiIndex, Series
import pandas._testing as tm
from pandas.core.base import SpecificationError
from pandas.core.groupby.base import transformation_kernels
from pandas.tests.frame.common import zip_frames


def test_transform_ufunc(axis, float_frame):
    # GH 35964
    with np.errstate(all="ignore"):
        f_sqrt = np.sqrt(float_frame)
    result = float_frame.transform(np.sqrt, axis=axis)
    expected = f_sqrt
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("op", transformation_kernels)
def test_transform_groupby_kernel(axis, float_frame, op):
    # GH 35964
    if op == "cumcount":
        pytest.xfail("DataFrame.cumcount does not exist")
    if op == "tshift":
        pytest.xfail("Only works on time index and is deprecated")
    if axis == 1 or axis == "columns":
        pytest.xfail("GH 36308: groupby.transform with axis=1 is broken")

    args = [0.0] if op == "fillna" else []
    if axis == 0 or axis == "index":
        ones = np.ones(float_frame.shape[0])
    else:
        ones = np.ones(float_frame.shape[1])
    expected = float_frame.groupby(ones, axis=axis).transform(op, *args)
    result = float_frame.transform(op, axis, *args)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "ops, names",
    [
        ([np.sqrt], ["sqrt"]),
        ([np.abs, np.sqrt], ["absolute", "sqrt"]),
        (np.array([np.sqrt]), ["sqrt"]),
        (np.array([np.abs, np.sqrt]), ["absolute", "sqrt"]),
    ],
)
def test_transform_listlike(axis, float_frame, ops, names):
    # GH 35964
    other_axis = 1 if axis in {0, "index"} else 0
    with np.errstate(all="ignore"):
        expected = zip_frames([op(float_frame) for op in ops], axis=other_axis)
    if axis in {0, "index"}:
        expected.columns = MultiIndex.from_product([float_frame.columns, names])
    else:
        expected.index = MultiIndex.from_product([float_frame.index, names])
    result = float_frame.transform(ops, axis=axis)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("ops", [[], np.array([])])
def test_transform_empty_listlike(float_frame, ops):
    with pytest.raises(ValueError, match="No transform functions were provided"):
        float_frame.transform(ops)


@pytest.mark.parametrize("box", [dict, Series])
def test_transform_dictlike(axis, float_frame, box):
    # GH 35964
    if axis == 0 or axis == "index":
        e = float_frame.columns[0]
        expected = float_frame[[e]].transform(np.abs)
    else:
        e = float_frame.index[0]
        expected = float_frame.iloc[[0]].transform(np.abs)
    result = float_frame.transform(box({e: np.abs}), axis=axis)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "ops",
    [
        {},
        {"A": []},
        {"A": [], "B": "cumsum"},
        {"A": "cumsum", "B": []},
        {"A": [], "B": ["cumsum"]},
        {"A": ["cumsum"], "B": []},
    ],
)
def test_transform_empty_dictlike(float_frame, ops):
    with pytest.raises(ValueError, match="No transform functions were provided"):
        float_frame.transform(ops)


@pytest.mark.parametrize("use_apply", [True, False])
def test_transform_udf(axis, float_frame, use_apply):
    # GH 35964
    # transform uses UDF either via apply or passing the entire DataFrame
    def func(x):
        # transform is using apply iff x is not a DataFrame
        if use_apply == isinstance(x, DataFrame):
            # Force transform to fallback
            raise ValueError
        return x + 1

    result = float_frame.transform(func, axis=axis)
    expected = float_frame + 1
    tm.assert_frame_equal(result, expected)


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
        # mypy identifies the argument as an invalid type
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

    # tshift is deprecated
    warn = None if op != "tshift" else FutureWarning
    with tm.assert_produces_warning(warn, check_stacklevel=False):
        with pytest.raises(ValueError, match=msg):
            df.transform(op)
        with pytest.raises(ValueError, match=msg):
            df.transform([op])
        with pytest.raises(ValueError, match=msg):
            df.transform({"A": op})
        with pytest.raises(ValueError, match=msg):
            df.transform({"A": [op]})


@pytest.mark.parametrize("op", transformation_kernels)
def test_transform_partial_failure(op):
    # GH 35964
    wont_fail = ["ffill", "bfill", "fillna", "pad", "backfill", "shift"]
    if op in wont_fail:
        pytest.xfail("Transform kernel is successful on all dtypes")
    if op == "cumcount":
        pytest.xfail("transform('cumcount') not implemented")
    if op == "tshift":
        pytest.xfail("Only works on time index; deprecated")

    # Using object makes most transform kernels fail
    df = DataFrame({"A": 3 * [object], "B": [1, 2, 3]})

    expected = df[["B"]].transform([op])
    result = df.transform([op])
    tm.assert_equal(result, expected)

    expected = df[["B"]].transform({"B": op})
    result = df.transform({"B": op})
    tm.assert_equal(result, expected)

    expected = df[["B"]].transform({"B": [op]})
    result = df.transform({"B": [op]})
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
