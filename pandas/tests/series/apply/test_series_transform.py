import numpy as np
import pytest

from pandas import DataFrame, Series, concat
import pandas._testing as tm
from pandas.core.base import SpecificationError
from pandas.core.groupby.base import transformation_kernels


def test_transform_ufunc(string_series):
    # GH 35964
    with np.errstate(all="ignore"):
        f_sqrt = np.sqrt(string_series)

    # ufunc
    result = string_series.transform(np.sqrt)
    expected = f_sqrt.copy()
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("op", transformation_kernels)
def test_transform_groupby_kernel(string_series, op):
    # GH 35964
    if op == "cumcount":
        pytest.xfail("Series.cumcount does not exist")
    if op == "tshift":
        pytest.xfail("Only works on time index and is deprecated")

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


@pytest.mark.parametrize("ops", [[], np.array([])])
def test_transform_empty_listlike(string_series, ops):
    with pytest.raises(ValueError, match="No transform functions were provided"):
        string_series.transform(ops)


@pytest.mark.parametrize("box", [dict, Series])
def test_transform_dictlike(string_series, box):
    # GH 35964
    with np.errstate(all="ignore"):
        expected = concat([np.sqrt(string_series), np.abs(string_series)], axis=1)
    expected.columns = ["foo", "bar"]
    result = string_series.transform(box({"foo": np.sqrt, "bar": np.abs}))
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "ops",
    [
        {},
        {"A": []},
        {"A": [], "B": ["cumsum"]},
        {"A": ["cumsum"], "B": []},
        {"A": [], "B": "cumsum"},
        {"A": "cumsum", "B": []},
    ],
)
def test_transform_empty_dictlike(string_series, ops):
    with pytest.raises(ValueError, match="No transform functions were provided"):
        string_series.transform(ops)


def test_transform_udf(axis, string_series):
    # GH 35964
    # via apply
    def func(x):
        if isinstance(x, Series):
            raise ValueError
        return x + 1

    result = string_series.transform(func)
    expected = string_series + 1
    tm.assert_series_equal(result, expected)

    # via map Series -> Series
    def func(x):
        if not isinstance(x, Series):
            raise ValueError
        return x + 1

    result = string_series.transform(func)
    expected = string_series + 1
    tm.assert_series_equal(result, expected)


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


def test_transform_reducer_raises(all_reductions):
    # GH 35964
    op = all_reductions
    s = Series([1, 2, 3])
    msg = "Function did not transform"
    with pytest.raises(ValueError, match=msg):
        s.transform(op)
    with pytest.raises(ValueError, match=msg):
        s.transform([op])
    with pytest.raises(ValueError, match=msg):
        s.transform({"A": op})
    with pytest.raises(ValueError, match=msg):
        s.transform({"A": [op]})


# mypy doesn't allow adding lists of different types
# https://github.com/python/mypy/issues/5492
@pytest.mark.parametrize("op", [*transformation_kernels, lambda x: x + 1])
def test_transform_bad_dtype(op):
    # GH 35964
    s = Series(3 * [object])  # Series that will fail on most transforms
    if op in ("backfill", "shift", "pad", "bfill", "ffill"):
        pytest.xfail("Transform function works on any datatype")

    msg = "Transform function failed"

    # tshift is deprecated
    warn = None if op != "tshift" else FutureWarning
    with tm.assert_produces_warning(warn, check_stacklevel=False):
        with pytest.raises(ValueError, match=msg):
            s.transform(op)
        with pytest.raises(ValueError, match=msg):
            s.transform([op])
        with pytest.raises(ValueError, match=msg):
            s.transform({"A": op})
        with pytest.raises(ValueError, match=msg):
            s.transform({"A": [op]})


@pytest.mark.parametrize("use_apply", [True, False])
def test_transform_passes_args(use_apply):
    # GH 35964
    # transform uses UDF either via apply or passing the entire Series
    expected_args = [1, 2]
    expected_kwargs = {"c": 3}

    def f(x, a, b, c):
        # transform is using apply iff x is not a Series
        if use_apply == isinstance(x, Series):
            # Force transform to fallback
            raise ValueError
        assert [a, b] == expected_args
        assert c == expected_kwargs["c"]
        return x

    Series([1]).transform(f, 0, *expected_args, **expected_kwargs)


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
