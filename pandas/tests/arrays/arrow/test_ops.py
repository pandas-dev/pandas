from __future__ import annotations

import operator
import re

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


@pytest.fixture(
    params=[
        pd.Categorical(["test"]),
        pd.offsets.Minute(3),
        pd.Interval(0, 1),
        object(),
    ],
    ids=["categorical", "dateoffset", "interval", "object"],
)
def unconvertible_object(request):
    """An object that pyarrow cannot convert to one of its own types."""
    return request.param


def object_array_of(value):
    # np.array([value, value]) would give a 2D array for array-like values
    #  such as Categorical
    result = np.empty(2, dtype=object)
    result.fill(value)
    return result


@pytest.mark.parametrize(
    "op", [operator.add, operator.sub, operator.mul, operator.and_, operator.or_]
)
def test_op_unconvertible_object_raises(unconvertible_object, op):
    # GH#62682 pyarrow cannot convert the operand, so we raise our own
    # TypeError instead of letting an ArrowInvalid/ArrowTypeError escape
    arr = pd.array([1, 2], dtype=ArrowDtype(pa.int64()))
    other = object_array_of(unconvertible_object)

    msg = f"operation '{op.__name__}' not supported for dtype 'int64[pyarrow]'"
    with pytest.raises(TypeError, match=re.escape(msg)):
        op(arr, other)


def test_cmp_unconvertible_object(unconvertible_object):
    # GH#62682 comparisons fall back to elementwise ops rather than raising
    # an ArrowInvalid/ArrowTypeError
    arr = pd.array([1, 2], dtype=ArrowDtype(pa.int64()))
    other = object_array_of(unconvertible_object)

    result = arr == other
    expected = pd.array([False, False], dtype=ArrowDtype(pa.bool_()))
    tm.assert_extension_array_equal(result, expected)

    msg = "|".join(
        ["not supported between", "Unordered Categoricals can only compare equality"]
    )
    with pytest.raises(TypeError, match=msg):
        arr < other


@pytest.mark.parametrize("other", [pd.offsets.Minute(3), pd.Interval(0, 1)])
def test_op_unconvertible_scalar(other):
    # GH#62682 the scalar path has its own _box_pa call; pyarrow raises
    #  ArrowTypeError for a DateOffset and ArrowInvalid for an Interval.
    #  Two non-NA entries are needed, or the reflected __ne__ below gets a
    #  length-1 array that its `not` accepts.
    arr = pd.array([1, 2, None], dtype=ArrowDtype(pa.int64()))

    msg = "operation 'add' not supported for dtype 'int64[pyarrow]'"
    with pytest.raises(TypeError, match=re.escape(msg)):
        arr + other

    result = arr == other
    expected = pd.array([False, False, None], dtype=ArrowDtype(pa.bool_()))
    tm.assert_extension_array_equal(result, expected)

    result = arr != other
    expected = pd.array([True, True, None], dtype=ArrowDtype(pa.bool_()))
    tm.assert_extension_array_equal(result, expected)

    with pytest.raises(TypeError, match="Invalid comparison"):
        arr < other


@pytest.mark.parametrize("op", [operator.eq, operator.ne, operator.lt, operator.ge])
def test_cmp_duration_offset_scalar(op):
    # GH#62682 pyarrow raises ArrowTypeError on a Tick; match the numpy-backed result
    # length 2 with differing values: `not` accepts a length-1 array, hiding the bug
    arr = pd.array([pd.Timedelta("1h"), pd.Timedelta("2h")])

    result = op(arr.astype("duration[ns][pyarrow]"), pd.offsets.Hour(1))
    tm.assert_numpy_array_equal(
        np.asarray(result, dtype=bool), np.asarray(op(arr, pd.offsets.Hour(1)))
    )


def test_cmp_partly_unconvertible_object(unconvertible_object):
    # GH#62682 the pointwise fallback still gives a convertible entry a real answer
    arr = pd.array([1, 2], dtype=ArrowDtype(pa.int64()))
    other = object_array_of(unconvertible_object)
    other[0] = 1

    result = arr == other
    expected = pd.array([True, False], dtype=ArrowDtype(pa.bool_()))
    tm.assert_extension_array_equal(result, expected)

    result = arr != other
    expected = pd.array([False, True], dtype=ArrowDtype(pa.bool_()))
    tm.assert_extension_array_equal(result, expected)


def test_cmp_array_valued_pointwise_result():
    # GH#62682 a scalar compared to a length-1 Categorical gives a length-1 bool
    #  array holding a real answer, which must not be discarded; a longer one is
    #  genuinely ambiguous and raises, both matching Int64
    arr = pd.array([1, 2], dtype=ArrowDtype(pa.int64()))

    result = arr == object_array_of(pd.Categorical([1]))
    expected = pd.array([True, False], dtype=ArrowDtype(pa.bool_()))
    tm.assert_extension_array_equal(result, expected)

    msg = "truth value of an array with more than one element is ambiguous"
    with pytest.raises(ValueError, match=msg):
        arr == object_array_of(pd.Categorical([1, 2]))


def test_cmp_list_of_array_likes():
    # GH#62682 isna on a LIST of array-likes returns a 2-D mask, which the
    #  object_array_of spelling above never produces
    other = [pd.Categorical(["a"]), pd.Categorical(["b"])]
    arr = pd.array([1, 2], dtype=ArrowDtype(pa.int64()))

    result = arr == other

    # Int64 gives the same answer
    expected = pd.array([False, False], dtype=ArrowDtype(pa.bool_()))
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize(
    "values, pa_type",
    [
        ([[1], [2]], pa.list_(pa.int64())),
        ([{"a": 1}, {"a": 2}], pa.struct([("a", pa.int64())])),
    ],
)
def test_cmp_nested_dtype_pointwise(values, pa_type):
    # GH#62682 nested dtypes still take the pointwise fallback
    arr = pd.array(values, dtype=ArrowDtype(pa_type))

    result = arr == [1, "a"]
    expected = pd.array([False, False], dtype=ArrowDtype(pa.bool_()))
    tm.assert_extension_array_equal(result, expected)


def test_cmp_list_dtype_numpy_scalar_elements():
    # GH#62682 numpy broadcasts [1] == np.int64(1) to array([True]); treat it
    #  like the other invalid list-likes
    arr = pd.array([[1], [2]], dtype=ArrowDtype(pa.list_(pa.int64())))
    other = pd.arrays.SparseArray([1, 2])

    result = arr == other
    expected = pd.array([False, False], dtype=ArrowDtype(pa.bool_()))
    tm.assert_extension_array_equal(result, expected)

    result = arr != other
    expected = pd.array([True, True], dtype=ArrowDtype(pa.bool_()))
    tm.assert_extension_array_equal(result, expected)

    with pytest.raises(TypeError, match="Invalid comparison"):
        arr < other


def test_cmp_ne_offset_array_with_na():
    # GH#62682 an NA pair must not evaluate `ne`: numpy defers to the offset,
    #  whose __ne__ is `not self == other`, and NA has no truth value
    arr = pd.array([pd.Timedelta("1h"), None], dtype="duration[ns][pyarrow]")
    other = object_array_of(pd.offsets.Hour(1))

    result = arr != other
    expected = pd.array([False, None], dtype=ArrowDtype(pa.bool_()))
    tm.assert_extension_array_equal(result, expected)


def test_cmp_all_na_unordered_still_raises(unconvertible_object):
    # GH#62682 an all-NA operand must not turn an unsupported ordered
    #  comparison into an all-NA answer, as it would if op were skipped
    arr = pd.array([None, None], dtype=ArrowDtype(pa.int64()))
    other = object_array_of(unconvertible_object)

    msg = "|".join(
        ["not supported between", "Unordered Categoricals can only compare equality"]
    )
    with pytest.raises(TypeError, match=msg):
        arr < other


def test_str_arith_object_fallback_not_wrapped():
    # GH#62682 the string path must keep falling back to object dtype rather
    #  than reporting the operand as unsupported
    class Radd:
        def __radd__(self, other):
            return other + "!"

    other = np.empty(2, dtype=object)
    other[:] = [Radd(), Radd()]
    arr = pd.array(["a", "b"], dtype=ArrowDtype(pa.string()))

    result = arr + other
    expected = pd.array(["a!", "b!"], dtype=ArrowDtype(pa.string()))
    tm.assert_extension_array_equal(result, expected)


def test_cmp_unconvertible_object_keeps_na():
    # GH#62682 an NA entry stays NA rather than becoming a concrete bool,
    #  on whichever side it appears
    arr = pd.array([1, None], dtype=ArrowDtype(pa.int64()))
    other = object_array_of(pd.Categorical(["test"]))

    result = arr == other
    expected = pd.array([False, None], dtype=ArrowDtype(pa.bool_()))
    tm.assert_extension_array_equal(result, expected)

    other[0] = np.nan
    result = pd.array([1, 2], dtype=ArrowDtype(pa.int64())) == other
    expected = pd.array([None, False], dtype=ArrowDtype(pa.bool_()))
    tm.assert_extension_array_equal(result, expected)


def test_cmp_mixed_object_keeps_na():
    # GH#60228 / GH#62682 the pointwise fallback must not turn an NA entry into
    #  a concrete bool
    arr = pd.array([1, None], dtype=ArrowDtype(pa.int64()))
    result = arr == np.array([1, "b"], dtype=object)
    expected = pd.array([True, None], dtype=ArrowDtype(pa.bool_()))
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize(
    "dtype, values", [(ArrowDtype(pa.int64()), [1, None]), (None, ["a", None])]
)
def test_arith_dataframe_of_unconvertible_objects(unconvertible_object, dtype, values):
    # GH#62682 the operand reaches ArrowExtensionArray as a DataFrame column, as
    #  in the issue; dtype=object is required, since inferred interval columns
    #  box as pyarrow structs and raise TypeError before reaching this path.
    arr = pd.array(values, dtype=dtype)
    df = pd.DataFrame([[unconvertible_object, unconvertible_object]], dtype=object)
    msg = "|".join(
        [
            "can only concatenate str",
            re.escape(f"operation 'add' not supported for dtype '{arr.dtype}'"),
        ]
    )
    with pytest.raises(TypeError, match=msg):
        arr + df


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
