from __future__ import annotations

import operator

import numpy as np
import pytest

from pandas.errors import (
    Pandas4Warning,
)

from pandas.core.dtypes.dtypes import (
    ArrowDtype,
)

import pandas as pd
import pandas._testing as tm

pa = pytest.importorskip("pyarrow")


# GH#62423; matched instead of the leading clause, whose wording differs per warn site
depr_msg = "In a future version these will be treated as scalar-like"


class TestLogicalOps:
    """Various Series and DataFrame logical ops methods."""

    def test_kleene_or(self):
        a = pd.Series([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean[pyarrow]")
        b = pd.Series([True, False, None] * 3, dtype="boolean[pyarrow]")
        result = a | b
        expected = pd.Series(
            [True, True, True, True, False, None, True, None, None],
            dtype="boolean[pyarrow]",
        )
        tm.assert_series_equal(result, expected)

        result = b | a
        tm.assert_series_equal(result, expected)

        # ensure we haven't mutated anything inplace
        tm.assert_series_equal(
            a,
            pd.Series([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean[pyarrow]"),
        )
        tm.assert_series_equal(
            b, pd.Series([True, False, None] * 3, dtype="boolean[pyarrow]")
        )

    @pytest.mark.parametrize(
        "other, expected",
        [
            (None, [True, None, None]),
            (pd.NA, [True, None, None]),
            (True, [True, True, True]),
            (np.bool_(True), [True, True, True]),
            (False, [True, False, None]),
            (np.bool_(False), [True, False, None]),
        ],
    )
    def test_kleene_or_scalar(self, other, expected):
        a = pd.Series([True, False, None], dtype="boolean[pyarrow]")
        result = a | other
        expected = pd.Series(expected, dtype="boolean[pyarrow]")
        tm.assert_series_equal(result, expected)

        result = other | a
        tm.assert_series_equal(result, expected)

        # ensure we haven't mutated anything inplace
        tm.assert_series_equal(
            a, pd.Series([True, False, None], dtype="boolean[pyarrow]")
        )

    def test_kleene_and(self):
        a = pd.Series([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean[pyarrow]")
        b = pd.Series([True, False, None] * 3, dtype="boolean[pyarrow]")
        result = a & b
        expected = pd.Series(
            [True, False, None, False, False, False, None, False, None],
            dtype="boolean[pyarrow]",
        )
        tm.assert_series_equal(result, expected)

        result = b & a
        tm.assert_series_equal(result, expected)

        # ensure we haven't mutated anything inplace
        tm.assert_series_equal(
            a,
            pd.Series([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean[pyarrow]"),
        )
        tm.assert_series_equal(
            b, pd.Series([True, False, None] * 3, dtype="boolean[pyarrow]")
        )

    @pytest.mark.parametrize(
        "other, expected",
        [
            (None, [None, False, None]),
            (pd.NA, [None, False, None]),
            (True, [True, False, None]),
            (False, [False, False, False]),
            (np.bool_(True), [True, False, None]),
            (np.bool_(False), [False, False, False]),
        ],
    )
    def test_kleene_and_scalar(self, other, expected):
        a = pd.Series([True, False, None], dtype="boolean[pyarrow]")
        result = a & other
        expected = pd.Series(expected, dtype="boolean[pyarrow]")
        tm.assert_series_equal(result, expected)

        result = other & a
        tm.assert_series_equal(result, expected)

        # ensure we haven't mutated anything inplace
        tm.assert_series_equal(
            a, pd.Series([True, False, None], dtype="boolean[pyarrow]")
        )

    def test_kleene_xor(self):
        a = pd.Series([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean[pyarrow]")
        b = pd.Series([True, False, None] * 3, dtype="boolean[pyarrow]")
        result = a ^ b
        expected = pd.Series(
            [False, True, None, True, False, None, None, None, None],
            dtype="boolean[pyarrow]",
        )
        tm.assert_series_equal(result, expected)

        result = b ^ a
        tm.assert_series_equal(result, expected)

        # ensure we haven't mutated anything inplace
        tm.assert_series_equal(
            a,
            pd.Series([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean[pyarrow]"),
        )
        tm.assert_series_equal(
            b, pd.Series([True, False, None] * 3, dtype="boolean[pyarrow]")
        )

    @pytest.mark.parametrize(
        "other, expected",
        [
            (None, [None, None, None]),
            (pd.NA, [None, None, None]),
            (True, [False, True, None]),
            (np.bool_(True), [False, True, None]),
            (np.bool_(False), [True, False, None]),
        ],
    )
    def test_kleene_xor_scalar(self, other, expected):
        a = pd.Series([True, False, None], dtype="boolean[pyarrow]")
        result = a ^ other
        expected = pd.Series(expected, dtype="boolean[pyarrow]")
        tm.assert_series_equal(result, expected)

        result = other ^ a
        tm.assert_series_equal(result, expected)

        # ensure we haven't mutated anything inplace
        tm.assert_series_equal(
            a, pd.Series([True, False, None], dtype="boolean[pyarrow]")
        )

    @pytest.mark.parametrize(
        "op, exp",
        [
            ["__and__", True],
            ["__or__", True],
            ["__xor__", False],
        ],
    )
    def test_logical_masked_numpy(self, op, exp):
        # GH 52625
        data = [True, False, None]
        ser_masked = pd.Series(data, dtype="boolean")
        ser_pa = pd.Series(data, dtype="boolean[pyarrow]")
        result = getattr(ser_pa, op)(ser_masked)
        expected = pd.Series([exp, False, None], dtype=ArrowDtype(pa.bool_()))
        tm.assert_series_equal(result, expected)


def test_compare_range_len(data, comparison_op):
    # GH#63429 a range compares elementwise like the equivalent list.
    #  Note we can't go through _compare_other here: its pointwise
    #  expectation uses Series.combine, which treats the range as a scalar
    #  and so matches an all-False result too.
    ser = pd.Series(data)
    rng = range(len(ser))

    try:
        expected = comparison_op(ser, list(rng))
    except Exception as err:
        with tm.assert_produces_warning(Pandas4Warning, match=depr_msg):
            with pytest.raises(type(err)):
                comparison_op(ser, rng)
        return

    with tm.assert_produces_warning(Pandas4Warning, match=depr_msg):
        result = comparison_op(ser, rng)
    tm.assert_series_equal(result, expected)


def test_compare_range_mismatched_len(data, comparison_op):
    # GH#63429 the length check must not be bypassed for a range
    ser = pd.Series(data)
    rng = range(len(ser) + 1)

    with tm.assert_produces_warning(Pandas4Warning, match=depr_msg):
        with pytest.raises(ValueError, match="Lengths must match to compare"):
            comparison_op(ser, rng)


def test_compare_iterator_warns(data, comparison_op):
    # GH#31646 an iterator is a non-standard list-like (GH#62423)
    with tm.assert_produces_warning(Pandas4Warning, match=depr_msg):
        with pytest.raises(NotImplementedError, match="not implemented for"):
            comparison_op(pd.Series(data), iter(data))


def test_invalid_other_comp(data, comparison_op):
    # GH 48833
    with pytest.raises(
        NotImplementedError, match=".* not implemented for <class 'object'>"
    ):
        comparison_op(data, object())


@pytest.mark.parametrize("masked_dtype", ["boolean", "Int64", "Float64"])
def test_comp_masked_numpy(masked_dtype, comparison_op):
    # GH 52625
    data = [1, 0, None]
    ser_masked = pd.Series(data, dtype=masked_dtype)
    ser_pa = pd.Series(data, dtype=f"{masked_dtype.lower()}[pyarrow]")
    result = comparison_op(ser_pa, ser_masked)
    if comparison_op in [operator.lt, operator.gt, operator.ne]:
        exp = [False, False, None]
    else:
        exp = [True, True, None]
    expected = pd.Series(exp, dtype=ArrowDtype(pa.bool_()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("pa_type", tm.ALL_INT_PYARROW_DTYPES)
def test_bitwise(pa_type):
    # GH 54495
    dtype = ArrowDtype(pa_type)
    left = pd.Series([1, None, 3, 4], dtype=dtype)
    right = pd.Series([None, 3, 5, 4], dtype=dtype)

    result = left | right
    expected = pd.Series([None, None, 3 | 5, 4 | 4], dtype=dtype)
    tm.assert_series_equal(result, expected)

    result = left & right
    expected = pd.Series([None, None, 3 & 5, 4 & 4], dtype=dtype)
    tm.assert_series_equal(result, expected)

    result = left ^ right
    expected = pd.Series([None, None, 3 ^ 5, 4 ^ 4], dtype=dtype)
    tm.assert_series_equal(result, expected)

    result = ~left
    expected = ~(left.fillna(0).to_numpy())
    expected = pd.Series(expected, dtype=dtype).mask(left.isnull())
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "other", [np.array([1]), [1], pd.array([1], dtype="int64[pyarrow]")]
)
def test_cmp_length_mismatch_raises(other):
    # GH#62682 pyarrow otherwise raises "Array arguments must all be the
    #  same length"
    arr = pd.array([1, 2], dtype=ArrowDtype(pa.int64()))
    with pytest.raises(ValueError, match="Lengths must match to compare"):
        arr == other


@pytest.mark.parametrize("op", [operator.add, operator.eq, operator.and_])
@pytest.mark.parametrize(
    "pa_type, values", [(pa.int64(), [1, 2]), (pa.string(), ["a", "b"])]
)
def test_op_2d_ndarray_raises(op, pa_type, values):
    # GH#62682 match BaseMaskedArray instead of raising an opaque ArrowInvalid
    arr = pd.array(values, dtype=ArrowDtype(pa_type))
    other = np.array([[1, 2], [3, 4]])
    with pytest.raises(NotImplementedError, match="can only perform ops with 1-d"):
        op(arr, other)


@pytest.mark.parametrize("op", [operator.add, operator.eq, operator.and_])
@pytest.mark.parametrize(
    "other",
    [
        pd.arrays.IntegerArray(np.array([[1], [2]]), np.zeros((2, 1), dtype=bool)),
        pd.date_range("2020", periods=2)._data.reshape(2, 1),
        # ABCExtensionArray does not match NumpyExtensionArray
        pd.arrays.NumpyExtensionArray(np.array([[1, 2], [3, 4]])),
    ],
    ids=["masked", "datetimelike", "numpy_ea"],
)
def test_op_2d_extension_array_raises(other, op):
    # GH#62682 a 2-D EA operand leaked "Mask must be 1D array" for the
    #  masked case, and `==` silently compared all-False for the datetimelike one
    arr = pd.array([1, 2], dtype=ArrowDtype(pa.int64()))
    with pytest.raises(NotImplementedError, match="can only perform ops with 1-d"):
        op(arr, other)


@pytest.mark.parametrize("op", [operator.or_, operator.and_, operator.xor])
def test_logical_2d_ndarray_bool_raises(op):
    # GH#62682 the GH#60234 string-vs-bool arm returns before _evaluate_op_method,
    #  so it needs a 1-d check of its own; StringDtype is covered by test_logical_2d
    arr = pd.array(["a", "b"], dtype=ArrowDtype(pa.string()))
    other = np.array([[True, False], [True, False]])
    with pytest.raises(NotImplementedError, match="can only perform ops with 1-d"):
        op(other, arr)


def test_pow_missing_operand():
    # GH 55512
    k = pd.Series([2, None], dtype="int64[pyarrow]")
    result = k.pow(None, fill_value=3)
    expected = pd.Series([8, None], dtype="int64[pyarrow]")
    tm.assert_series_equal(result, expected)


def test_comparison_not_propagating_arrow_error():
    # GH#54944
    a = pd.Series([1 << 63], dtype="uint64[pyarrow]")
    b = pd.Series([None], dtype="int64[pyarrow]")
    with pytest.raises(pa.lib.ArrowInvalid, match="Integer value"):
        a < b


@pytest.mark.parametrize("box", [pd.Series, pd.Index])
def test_comparison_range_matches_numpy_backed(box, comparison_op):
    # GH#63429 an arrow-backed box compared to a range must give the same
    #  elementwise answer as the numpy-backed equivalent, not the all-False
    #  result that routing through ops.invalid_comparison produces.
    #  Series and Index gate both the GH#62423 warning and the length check
    #  separately before reaching the EA, so both need covering.
    values = [5, 1, 7]

    with tm.assert_produces_warning(Pandas4Warning, match=depr_msg):
        expected = comparison_op(box(values), range(3))
    with tm.assert_produces_warning(Pandas4Warning, match=depr_msg):
        result = comparison_op(box(values, dtype="int64[pyarrow]"), range(3))

    tm.assert_numpy_array_equal(np.asarray(result), np.asarray(expected))
    assert result.dtype == "bool[pyarrow]"

    with tm.assert_produces_warning(Pandas4Warning, match=depr_msg):
        with pytest.raises(ValueError, match="Lengths must match to compare"):
            comparison_op(box(values, dtype="int64[pyarrow]"), range(4))


def test_arrow_floordiv():
    # GH 55561
    a = pd.Series([-7], dtype="int64[pyarrow]")
    b = pd.Series([4], dtype="int64[pyarrow]")
    expected = pd.Series([-2], dtype="int64[pyarrow]")
    result = a // b
    tm.assert_series_equal(result, expected)


def test_arrow_floordiv_large_values():
    # GH 56645
    a = pd.Series([1425801600000000000], dtype="int64[pyarrow]")
    expected = pd.Series([1425801600000], dtype="int64[pyarrow]")
    result = a // 1_000_000
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("dtype", ["int64[pyarrow]", "uint64[pyarrow]"])
def test_arrow_floordiv_large_integral_result(dtype):
    # GH 56676
    a = pd.Series([18014398509481983], dtype=dtype)
    result = a // 1
    tm.assert_series_equal(result, a)


@pytest.mark.parametrize("pa_type", tm.SIGNED_INT_PYARROW_DTYPES)
def test_arrow_floordiv_larger_divisor(pa_type):
    # GH 56676
    dtype = ArrowDtype(pa_type)
    a = pd.Series([-23], dtype=dtype)
    result = a // 24
    expected = pd.Series([-1], dtype=dtype)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("pa_type", tm.SIGNED_INT_PYARROW_DTYPES)
def test_arrow_floordiv_integral_invalid(pa_type):
    # GH 56676
    min_value = np.iinfo(pa_type.to_pandas_dtype()).min
    a = pd.Series([min_value], dtype=ArrowDtype(pa_type))
    msg = "|".join(["overflow", "not in range"])
    with pytest.raises(pa.lib.ArrowInvalid, match=msg):
        a // -1
    with pytest.raises(pa.lib.ArrowInvalid, match="divide by zero"):
        a // 0


@pytest.mark.parametrize("dtype", tm.FLOAT_PYARROW_DTYPES_STR_REPR)
def test_arrow_floordiv_floating_0_divisor(dtype):
    # GH 56676
    a = pd.Series([2], dtype=dtype)
    result = a // 0
    expected = pd.Series([float("inf")], dtype=dtype)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("pa_type", tm.ALL_INT_PYARROW_DTYPES)
def test_arrow_integral_floordiv_large_values(pa_type):
    # GH 56676
    max_value = np.iinfo(pa_type.to_pandas_dtype()).max
    dtype = ArrowDtype(pa_type)
    a = pd.Series([max_value], dtype=dtype)
    b = pd.Series([1], dtype=dtype)
    result = a // b
    tm.assert_series_equal(result, a)


@pytest.mark.parametrize("dtype", ["int64[pyarrow]", "uint64[pyarrow]"])
def test_arrow_true_division_large_divisor(dtype):
    # GH 56706
    a = pd.Series([0], dtype=dtype)
    b = pd.Series([18014398509481983], dtype=dtype)
    expected = pd.Series([0], dtype="float64[pyarrow]")
    result = a / b
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("dtype", ["int64[pyarrow]", "uint64[pyarrow]"])
def test_arrow_floor_division_large_divisor(dtype):
    # GH 56706
    a = pd.Series([0], dtype=dtype)
    b = pd.Series([18014398509481983], dtype=dtype)
    expected = pd.Series([0], dtype=dtype)
    result = a // b
    tm.assert_series_equal(result, expected)


def test_ops_with_nan_is_na(using_nan_is_na):
    # GH#61732
    ser = pd.Series([-1, 0, 1], dtype="int64[pyarrow]")

    result = ser - np.nan
    if using_nan_is_na:
        assert result.isna().all()
    else:
        assert not result.isna().any()

    result = ser * np.nan
    if using_nan_is_na:
        assert result.isna().all()
    else:
        assert not result.isna().any()

    result = ser / 0
    if using_nan_is_na:
        assert result.isna()[1]
    else:
        assert not result.isna()[1]


def test_pow_with_all_na_float():
    # GH#62520

    s = pd.Series([None, None], dtype="float64[pyarrow]")
    result = s.pow(2)
    expected = pd.Series([pd.NA, pd.NA], dtype="float64[pyarrow]")
    tm.assert_series_equal(result, expected)


def test_mul_numpy_nullable_with_pyarrow_float():
    # GH#58602
    left = pd.Series(range(5), dtype="Float64")
    right = pd.Series(range(5), dtype="float64[pyarrow]")

    expected = pd.Series([0, 1, 4, 9, 16], dtype="float64[pyarrow]")

    result = left * right
    tm.assert_series_equal(result, expected)

    result2 = right * left
    tm.assert_series_equal(result2, expected)

    # while we're here, let's check __eq__
    result3 = left == right
    expected3 = pd.Series([True] * 5, dtype="bool[pyarrow]")
    tm.assert_series_equal(result3, expected3)

    result4 = right == left
    tm.assert_series_equal(result4, expected3)
