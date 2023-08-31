import numpy as np
import pytest

from pandas import (
    DataFrame,
    MultiIndex,
    Series,
    concat,
)
import pandas._testing as tm
from pandas.tests.apply.common import (
    series_transform_kernels,
    transform_obj,
)


@pytest.mark.parametrize(
    "args, kwargs, increment",
    [((), {}, 0), ((), {"a": 1}, 1), ((2, 3), {}, 32), ((1,), {"c": 2}, 201)],
)
def test_agg_args(args, kwargs, increment, series_ops_only):
    # GH 43357
    def f(x, a=0, b=0, c=0):
        return x + a + 10 * b + 100 * c

    s = Series([1, 2])
    result = transform_obj(s, f, *args, series_ops_only=series_ops_only, **kwargs)
    expected = s + increment
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
def test_transform_listlike(string_series, ops, names, series_ops_only):
    # GH 35964
    with np.errstate(all="ignore"):
        expected = concat([op(string_series) for op in ops], axis=1)
        expected.columns = names
        result = transform_obj(string_series, ops, series_ops_only=series_ops_only)
        tm.assert_frame_equal(result, expected)


def test_transform_listlike_func_with_args(series_ops_only):
    # GH 50624

    s = Series([1, 2, 3])

    def foo1(x, a=1, c=0):
        return x + a + c

    def foo2(x, b=2, c=0):
        return x + b + c

    msg = r"foo1\(\) got an unexpected keyword argument 'b'"
    with pytest.raises(TypeError, match=msg):
        transform_obj(s, [foo1, foo2], 3, series_ops_only=series_ops_only, b=3, c=4)

    result = transform_obj(s, [foo1, foo2], 3, series_ops_only=series_ops_only, c=4)
    expected = DataFrame({"foo1": [8, 9, 10], "foo2": [8, 9, 10]})
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("box", [dict, Series])
def test_transform_dictlike(string_series, box, series_ops_only):
    # GH 35964
    with np.errstate(all="ignore"):
        expected = concat([np.sqrt(string_series), np.abs(string_series)], axis=1)
    expected.columns = ["foo", "bar"]
    op = box({"foo": np.sqrt, "bar": np.abs})
    result = transform_obj(string_series, op, series_ops_only=series_ops_only)
    tm.assert_frame_equal(result, expected)


def test_transform_dictlike_mixed(series_ops_only):
    # GH 40018 - mix of lists and non-lists in values of a dictionary
    ser = Series([1, 4])
    op = {"b": ["sqrt", "abs"], "c": "sqrt"}
    result = transform_obj(ser, op, series_ops_only=series_ops_only)
    expected = DataFrame(
        [[1.0, 1, 1.0], [2.0, 4, 2.0]],
        columns=MultiIndex([("b", "c"), ("sqrt", "abs")], [(0, 0, 1), (0, 1, 0)]),
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("op", series_transform_kernels)
def test_transform_partial_failure(op, request, series_ops_only):
    # GH 35964
    if op in ("ffill", "bfill", "pad", "backfill", "shift"):
        request.node.add_marker(
            pytest.mark.xfail(reason=f"{op} is successful on any dtype")
        )

    # Using object makes most transform kernels fail
    ser = Series(3 * [object])

    if op == "ngroup" and series_ops_only:
        error = AttributeError
        msg = "'ngroup' is not a valid function for 'Series' object"
    elif op == "fillna" and series_ops_only:
        error = ValueError
        msg = "Must specify a fill 'value' or 'method'"
    elif op in ("fillna", "ngroup") and not series_ops_only:
        error = ValueError
        msg = "Transform function failed"
    else:
        error = TypeError
        msg = "|".join(
            [
                "not supported between instances of 'type' and 'type'",
                "unsupported operand type",
            ]
        )

    for func in [
        op,
        {"A": op, "B": "shift"},
        {"A": [op], "B": ["shift"]},
        {"A": [op, "shift"], "B": [op]},
    ]:
        with pytest.raises(error, match=msg):
            transform_obj(ser, func, series_ops_only=series_ops_only)


def test_transform_partial_failure_valueerror(series_ops_only):
    # GH 40211
    def noop(x):
        return x

    def raising_op(_):
        raise ValueError

    ser = Series(3 * [object])
    msg = "" if series_ops_only else "Transform function failed"

    for func in [
        [noop, raising_op],
        {"A": raising_op, "B": noop},
        {"A": [raising_op], "B": [noop]},
        {"A": [noop, raising_op], "B": [noop]},
    ]:
        with pytest.raises(ValueError, match=msg):
            transform_obj(ser, func, series_ops_only=series_ops_only)
