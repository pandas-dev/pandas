"""
These the test the public routines exposed in types/common.py
related to inference and not otherwise tested in types/test_common.py

"""
import collections
from collections import namedtuple
from datetime import date, datetime, time, timedelta
from decimal import Decimal
from fractions import Fraction
from io import StringIO
from numbers import Number
import re

import numpy as np
import pytest
import pytz

from pandas._libs import iNaT, lib, missing as libmissing
import pandas.util._test_decorators as td

from pandas.core.dtypes import inference
from pandas.core.dtypes.common import (
    ensure_categorical,
    ensure_int32,
    is_bool,
    is_datetime64_any_dtype,
    is_datetime64_dtype,
    is_datetime64_ns_dtype,
    is_datetime64tz_dtype,
    is_float,
    is_integer,
    is_number,
    is_scalar,
    is_scipy_sparse,
    is_timedelta64_dtype,
    is_timedelta64_ns_dtype,
)

import pandas as pd
from pandas import (
    Categorical,
    DataFrame,
    DateOffset,
    DatetimeIndex,
    Index,
    Interval,
    Period,
    Series,
    Timedelta,
    TimedeltaIndex,
    Timestamp,
    isna,
)
import pandas._testing as tm
from pandas.core.arrays import IntegerArray


@pytest.fixture(params=[True, False], ids=str)
def coerce(request):
    return request.param


# collect all objects to be tested for list-like-ness; use tuples of objects,
# whether they are list-like or not (special casing for sets), and their ID
ll_params = [
    ([1], True, "list"),
    ([], True, "list-empty"),
    ((1,), True, "tuple"),
    (tuple(), True, "tuple-empty"),
    ({"a": 1}, True, "dict"),
    (dict(), True, "dict-empty"),
    ({"a", 1}, "set", "set"),
    (set(), "set", "set-empty"),
    (frozenset({"a", 1}), "set", "frozenset"),
    (frozenset(), "set", "frozenset-empty"),
    (iter([1, 2]), True, "iterator"),
    (iter([]), True, "iterator-empty"),
    ((x for x in [1, 2]), True, "generator"),
    ((_ for _ in []), True, "generator-empty"),
    (Series([1]), True, "Series"),
    (Series([], dtype=object), True, "Series-empty"),
    (Series(["a"]).str, True, "StringMethods"),
    (Series([], dtype="O").str, True, "StringMethods-empty"),
    (Index([1]), True, "Index"),
    (Index([]), True, "Index-empty"),
    (DataFrame([[1]]), True, "DataFrame"),
    (DataFrame(), True, "DataFrame-empty"),
    (np.ndarray((2,) * 1), True, "ndarray-1d"),
    (np.array([]), True, "ndarray-1d-empty"),
    (np.ndarray((2,) * 2), True, "ndarray-2d"),
    (np.array([[]]), True, "ndarray-2d-empty"),
    (np.ndarray((2,) * 3), True, "ndarray-3d"),
    (np.array([[[]]]), True, "ndarray-3d-empty"),
    (np.ndarray((2,) * 4), True, "ndarray-4d"),
    (np.array([[[[]]]]), True, "ndarray-4d-empty"),
    (np.array(2), False, "ndarray-0d"),
    (1, False, "int"),
    (b"123", False, "bytes"),
    (b"", False, "bytes-empty"),
    ("123", False, "string"),
    ("", False, "string-empty"),
    (str, False, "string-type"),
    (object(), False, "object"),
    (np.nan, False, "NaN"),
    (None, False, "None"),
]
objs, expected, ids = zip(*ll_params)


@pytest.fixture(params=zip(objs, expected), ids=ids)
def maybe_list_like(request):
    return request.param


def test_is_list_like(maybe_list_like):
    obj, expected = maybe_list_like
    expected = True if expected == "set" else expected
    assert inference.is_list_like(obj) == expected


def test_is_list_like_disallow_sets(maybe_list_like):
    obj, expected = maybe_list_like
    expected = False if expected == "set" else expected
    assert inference.is_list_like(obj, allow_sets=False) == expected


def test_is_sequence():
    is_seq = inference.is_sequence
    assert is_seq((1, 2))
    assert is_seq([1, 2])
    assert not is_seq("abcd")
    assert not is_seq(np.int64)

    class A:
        def __getitem__(self):
            return 1

    assert not is_seq(A())


def test_is_array_like():
    assert inference.is_array_like(Series([], dtype=object))
    assert inference.is_array_like(Series([1, 2]))
    assert inference.is_array_like(np.array(["a", "b"]))
    assert inference.is_array_like(Index(["2016-01-01"]))

    class DtypeList(list):
        dtype = "special"

    assert inference.is_array_like(DtypeList())

    assert not inference.is_array_like([1, 2, 3])
    assert not inference.is_array_like(tuple())
    assert not inference.is_array_like("foo")
    assert not inference.is_array_like(123)


@pytest.mark.parametrize(
    "inner",
    [
        [],
        [1],
        (1,),
        (1, 2),
        {"a": 1},
        {1, "a"},
        Series([1]),
        Series([], dtype=object),
        Series(["a"]).str,
        (x for x in range(5)),
    ],
)
@pytest.mark.parametrize("outer", [list, Series, np.array, tuple])
def test_is_nested_list_like_passes(inner, outer):
    result = outer([inner for _ in range(5)])
    assert inference.is_list_like(result)


@pytest.mark.parametrize(
    "obj",
    [
        "abc",
        [],
        [1],
        (1,),
        ["a"],
        "a",
        {"a"},
        [1, 2, 3],
        Series([1]),
        DataFrame({"A": [1]}),
        ([1, 2] for _ in range(5)),
    ],
)
def test_is_nested_list_like_fails(obj):
    assert not inference.is_nested_list_like(obj)


@pytest.mark.parametrize("ll", [{}, {"A": 1}, Series([1]), collections.defaultdict()])
def test_is_dict_like_passes(ll):
    assert inference.is_dict_like(ll)


@pytest.mark.parametrize(
    "ll",
    [
        "1",
        1,
        [1, 2],
        (1, 2),
        range(2),
        Index([1]),
        dict,
        collections.defaultdict,
        Series,
    ],
)
def test_is_dict_like_fails(ll):
    assert not inference.is_dict_like(ll)


@pytest.mark.parametrize("has_keys", [True, False])
@pytest.mark.parametrize("has_getitem", [True, False])
@pytest.mark.parametrize("has_contains", [True, False])
def test_is_dict_like_duck_type(has_keys, has_getitem, has_contains):
    class DictLike:
        def __init__(self, d):
            self.d = d

        if has_keys:

            def keys(self):
                return self.d.keys()

        if has_getitem:

            def __getitem__(self, key):
                return self.d.__getitem__(key)

        if has_contains:

            def __contains__(self, key) -> bool:
                return self.d.__contains__(key)

    d = DictLike({1: 2})
    result = inference.is_dict_like(d)
    expected = has_keys and has_getitem and has_contains

    assert result is expected


def test_is_file_like():
    class MockFile:
        pass

    is_file = inference.is_file_like

    data = StringIO("data")
    assert is_file(data)

    # No read / write attributes
    # No iterator attributes
    m = MockFile()
    assert not is_file(m)

    MockFile.write = lambda self: 0

    # Write attribute but not an iterator
    m = MockFile()
    assert not is_file(m)

    # gh-16530: Valid iterator just means we have the
    # __iter__ attribute for our purposes.
    MockFile.__iter__ = lambda self: self

    # Valid write-only file
    m = MockFile()
    assert is_file(m)

    del MockFile.write
    MockFile.read = lambda self: 0

    # Valid read-only file
    m = MockFile()
    assert is_file(m)

    # Iterator but no read / write attributes
    data = [1, 2, 3]
    assert not is_file(data)


test_tuple = collections.namedtuple("Test", ["a", "b", "c"])


@pytest.mark.parametrize("ll", [test_tuple(1, 2, 3)])
def test_is_names_tuple_passes(ll):
    assert inference.is_named_tuple(ll)


@pytest.mark.parametrize("ll", [(1, 2, 3), "a", Series({"pi": 3.14})])
def test_is_names_tuple_fails(ll):
    assert not inference.is_named_tuple(ll)


def test_is_hashable():

    # all new-style classes are hashable by default
    class HashableClass:
        pass

    class UnhashableClass1:
        __hash__ = None

    class UnhashableClass2:
        def __hash__(self):
            raise TypeError("Not hashable")

    hashable = (1, 3.14, np.float64(3.14), "a", tuple(), (1,), HashableClass())
    not_hashable = ([], UnhashableClass1())
    abc_hashable_not_really_hashable = (([],), UnhashableClass2())

    for i in hashable:
        assert inference.is_hashable(i)
    for i in not_hashable:
        assert not inference.is_hashable(i)
    for i in abc_hashable_not_really_hashable:
        assert not inference.is_hashable(i)

    # numpy.array is no longer collections.abc.Hashable as of
    # https://github.com/numpy/numpy/pull/5326, just test
    # is_hashable()
    assert not inference.is_hashable(np.array([]))


@pytest.mark.parametrize("ll", [re.compile("ad")])
def test_is_re_passes(ll):
    assert inference.is_re(ll)


@pytest.mark.parametrize("ll", ["x", 2, 3, object()])
def test_is_re_fails(ll):
    assert not inference.is_re(ll)


@pytest.mark.parametrize(
    "ll", [r"a", "x", r"asdf", re.compile("adsf"), r"\u2233\s*", re.compile(r"")]
)
def test_is_recompilable_passes(ll):
    assert inference.is_re_compilable(ll)


@pytest.mark.parametrize("ll", [1, [], object()])
def test_is_recompilable_fails(ll):
    assert not inference.is_re_compilable(ll)


class TestInference:
    def test_infer_dtype_bytes(self):
        compare = "bytes"

        # string array of bytes
        arr = np.array(list("abc"), dtype="S1")
        assert lib.infer_dtype(arr, skipna=True) == compare

        # object array of bytes
        arr = arr.astype(object)
        assert lib.infer_dtype(arr, skipna=True) == compare

        # object array of bytes with missing values
        assert lib.infer_dtype([b"a", np.nan, b"c"], skipna=True) == compare

    def test_isinf_scalar(self):
        # GH 11352
        assert libmissing.isposinf_scalar(float("inf"))
        assert libmissing.isposinf_scalar(np.inf)
        assert not libmissing.isposinf_scalar(-np.inf)
        assert not libmissing.isposinf_scalar(1)
        assert not libmissing.isposinf_scalar("a")

        assert libmissing.isneginf_scalar(float("-inf"))
        assert libmissing.isneginf_scalar(-np.inf)
        assert not libmissing.isneginf_scalar(np.inf)
        assert not libmissing.isneginf_scalar(1)
        assert not libmissing.isneginf_scalar("a")

    @pytest.mark.parametrize("maybe_int", [True, False])
    @pytest.mark.parametrize(
        "infinity", ["inf", "inF", "iNf", "Inf", "iNF", "InF", "INf", "INF"]
    )
    def test_maybe_convert_numeric_infinities(self, infinity, maybe_int):
        # see gh-13274
        na_values = {"", "NULL", "nan"}

        pos = np.array(["inf"], dtype=np.float64)
        neg = np.array(["-inf"], dtype=np.float64)

        msg = "Unable to parse string"

        out = lib.maybe_convert_numeric(
            np.array([infinity], dtype=object), na_values, maybe_int
        )
        tm.assert_numpy_array_equal(out, pos)

        out = lib.maybe_convert_numeric(
            np.array(["-" + infinity], dtype=object), na_values, maybe_int
        )
        tm.assert_numpy_array_equal(out, neg)

        out = lib.maybe_convert_numeric(
            np.array([infinity], dtype=object), na_values, maybe_int
        )
        tm.assert_numpy_array_equal(out, pos)

        out = lib.maybe_convert_numeric(
            np.array(["+" + infinity], dtype=object), na_values, maybe_int
        )
        tm.assert_numpy_array_equal(out, pos)

        # too many characters
        with pytest.raises(ValueError, match=msg):
            lib.maybe_convert_numeric(
                np.array(["foo_" + infinity], dtype=object), na_values, maybe_int
            )

    def test_maybe_convert_numeric_post_floatify_nan(self, coerce):
        # see gh-13314
        data = np.array(["1.200", "-999.000", "4.500"], dtype=object)
        expected = np.array([1.2, np.nan, 4.5], dtype=np.float64)
        nan_values = {-999, -999.0}

        out = lib.maybe_convert_numeric(data, nan_values, coerce)
        tm.assert_numpy_array_equal(out, expected)

    def test_convert_infs(self):
        arr = np.array(["inf", "inf", "inf"], dtype="O")
        result = lib.maybe_convert_numeric(arr, set(), False)
        assert result.dtype == np.float64

        arr = np.array(["-inf", "-inf", "-inf"], dtype="O")
        result = lib.maybe_convert_numeric(arr, set(), False)
        assert result.dtype == np.float64

    def test_scientific_no_exponent(self):
        # See PR 12215
        arr = np.array(["42E", "2E", "99e", "6e"], dtype="O")
        result = lib.maybe_convert_numeric(arr, set(), False, True)
        assert np.all(np.isnan(result))

    def test_convert_non_hashable(self):
        # GH13324
        # make sure that we are handing non-hashables
        arr = np.array([[10.0, 2], 1.0, "apple"], dtype=object)
        result = lib.maybe_convert_numeric(arr, set(), False, True)
        tm.assert_numpy_array_equal(result, np.array([np.nan, 1.0, np.nan]))

    def test_convert_numeric_uint64(self):
        arr = np.array([2 ** 63], dtype=object)
        exp = np.array([2 ** 63], dtype=np.uint64)
        tm.assert_numpy_array_equal(lib.maybe_convert_numeric(arr, set()), exp)

        arr = np.array([str(2 ** 63)], dtype=object)
        exp = np.array([2 ** 63], dtype=np.uint64)
        tm.assert_numpy_array_equal(lib.maybe_convert_numeric(arr, set()), exp)

        arr = np.array([np.uint64(2 ** 63)], dtype=object)
        exp = np.array([2 ** 63], dtype=np.uint64)
        tm.assert_numpy_array_equal(lib.maybe_convert_numeric(arr, set()), exp)

    @pytest.mark.parametrize(
        "arr",
        [
            np.array([2 ** 63, np.nan], dtype=object),
            np.array([str(2 ** 63), np.nan], dtype=object),
            np.array([np.nan, 2 ** 63], dtype=object),
            np.array([np.nan, str(2 ** 63)], dtype=object),
        ],
    )
    def test_convert_numeric_uint64_nan(self, coerce, arr):
        expected = arr.astype(float) if coerce else arr.copy()
        result = lib.maybe_convert_numeric(arr, set(), coerce_numeric=coerce)
        tm.assert_almost_equal(result, expected)

    def test_convert_numeric_uint64_nan_values(self, coerce):
        arr = np.array([2 ** 63, 2 ** 63 + 1], dtype=object)
        na_values = {2 ** 63}

        expected = (
            np.array([np.nan, 2 ** 63 + 1], dtype=float) if coerce else arr.copy()
        )
        result = lib.maybe_convert_numeric(arr, na_values, coerce_numeric=coerce)
        tm.assert_almost_equal(result, expected)

    @pytest.mark.parametrize(
        "case",
        [
            np.array([2 ** 63, -1], dtype=object),
            np.array([str(2 ** 63), -1], dtype=object),
            np.array([str(2 ** 63), str(-1)], dtype=object),
            np.array([-1, 2 ** 63], dtype=object),
            np.array([-1, str(2 ** 63)], dtype=object),
            np.array([str(-1), str(2 ** 63)], dtype=object),
        ],
    )
    def test_convert_numeric_int64_uint64(self, case, coerce):
        expected = case.astype(float) if coerce else case.copy()
        result = lib.maybe_convert_numeric(case, set(), coerce_numeric=coerce)
        tm.assert_almost_equal(result, expected)

    @pytest.mark.parametrize("value", [-(2 ** 63) - 1, 2 ** 64])
    def test_convert_int_overflow(self, value):
        # see gh-18584
        arr = np.array([value], dtype=object)
        result = lib.maybe_convert_objects(arr)
        tm.assert_numpy_array_equal(arr, result)

    def test_maybe_convert_objects_uint64(self):
        # see gh-4471
        arr = np.array([2 ** 63], dtype=object)
        exp = np.array([2 ** 63], dtype=np.uint64)
        tm.assert_numpy_array_equal(lib.maybe_convert_objects(arr), exp)

        # NumPy bug: can't compare uint64 to int64, as that
        # results in both casting to float64, so we should
        # make sure that this function is robust against it
        arr = np.array([np.uint64(2 ** 63)], dtype=object)
        exp = np.array([2 ** 63], dtype=np.uint64)
        tm.assert_numpy_array_equal(lib.maybe_convert_objects(arr), exp)

        arr = np.array([2, -1], dtype=object)
        exp = np.array([2, -1], dtype=np.int64)
        tm.assert_numpy_array_equal(lib.maybe_convert_objects(arr), exp)

        arr = np.array([2 ** 63, -1], dtype=object)
        exp = np.array([2 ** 63, -1], dtype=object)
        tm.assert_numpy_array_equal(lib.maybe_convert_objects(arr), exp)

    def test_maybe_convert_objects_datetime(self):
        # GH27438
        arr = np.array(
            [np.datetime64("2000-01-01"), np.timedelta64(1, "s")], dtype=object
        )
        exp = arr.copy()
        out = lib.maybe_convert_objects(arr, convert_datetime=1, convert_timedelta=1)
        tm.assert_numpy_array_equal(out, exp)

        arr = np.array([pd.NaT, np.timedelta64(1, "s")], dtype=object)
        exp = np.array([np.timedelta64("NaT"), np.timedelta64(1, "s")], dtype="m8[ns]")
        out = lib.maybe_convert_objects(arr, convert_datetime=1, convert_timedelta=1)
        tm.assert_numpy_array_equal(out, exp)

        arr = np.array([np.timedelta64(1, "s"), np.nan], dtype=object)
        exp = arr.copy()
        out = lib.maybe_convert_objects(arr, convert_datetime=1, convert_timedelta=1)
        tm.assert_numpy_array_equal(out, exp)

    @pytest.mark.parametrize(
        "exp",
        [
            IntegerArray(np.array([2, 0], dtype="i8"), np.array([False, True])),
            IntegerArray(np.array([2, 0], dtype="int64"), np.array([False, True])),
        ],
    )
    def test_maybe_convert_objects_nullable_integer(self, exp):
        # GH27335
        arr = np.array([2, np.NaN], dtype=object)
        result = lib.maybe_convert_objects(arr, convert_to_nullable_integer=1)

        tm.assert_extension_array_equal(result, exp)

    def test_mixed_dtypes_remain_object_array(self):
        # GH14956
        array = np.array([datetime(2015, 1, 1, tzinfo=pytz.utc), 1], dtype=object)
        result = lib.maybe_convert_objects(array, convert_datetime=1)
        tm.assert_numpy_array_equal(result, array)


class TestTypeInference:

    # Dummy class used for testing with Python objects
    class Dummy:
        pass

    def test_inferred_dtype_fixture(self, any_skipna_inferred_dtype):
        # see pandas/conftest.py
        inferred_dtype, values = any_skipna_inferred_dtype

        # make sure the inferred dtype of the fixture is as requested
        assert inferred_dtype == lib.infer_dtype(values, skipna=True)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_length_zero(self, skipna):
        result = lib.infer_dtype(np.array([], dtype="i4"), skipna=skipna)
        assert result == "integer"

        result = lib.infer_dtype([], skipna=skipna)
        assert result == "empty"

        # GH 18004
        arr = np.array([np.array([], dtype=object), np.array([], dtype=object)])
        result = lib.infer_dtype(arr, skipna=skipna)
        assert result == "empty"

    def test_integers(self):
        arr = np.array([1, 2, 3, np.int64(4), np.int32(5)], dtype="O")
        result = lib.infer_dtype(arr, skipna=True)
        assert result == "integer"

        arr = np.array([1, 2, 3, np.int64(4), np.int32(5), "foo"], dtype="O")
        result = lib.infer_dtype(arr, skipna=True)
        assert result == "mixed-integer"

        arr = np.array([1, 2, 3, 4, 5], dtype="i4")
        result = lib.infer_dtype(arr, skipna=True)
        assert result == "integer"

    @pytest.mark.parametrize(
        "arr, skipna",
        [
            (np.array([1, 2, np.nan, np.nan, 3], dtype="O"), False),
            (np.array([1, 2, np.nan, np.nan, 3], dtype="O"), True),
            (np.array([1, 2, 3, np.int64(4), np.int32(5), np.nan], dtype="O"), False),
            (np.array([1, 2, 3, np.int64(4), np.int32(5), np.nan], dtype="O"), True),
        ],
    )
    def test_integer_na(self, arr, skipna):
        # GH 27392
        result = lib.infer_dtype(arr, skipna=skipna)
        expected = "integer" if skipna else "integer-na"
        assert result == expected

    def test_infer_dtype_skipna_default(self):
        # infer_dtype `skipna` default deprecated in GH#24050,
        #  changed to True in GH#29876
        arr = np.array([1, 2, 3, np.nan], dtype=object)

        result = lib.infer_dtype(arr)
        assert result == "integer"

    def test_bools(self):
        arr = np.array([True, False, True, True, True], dtype="O")
        result = lib.infer_dtype(arr, skipna=True)
        assert result == "boolean"

        arr = np.array([np.bool_(True), np.bool_(False)], dtype="O")
        result = lib.infer_dtype(arr, skipna=True)
        assert result == "boolean"

        arr = np.array([True, False, True, "foo"], dtype="O")
        result = lib.infer_dtype(arr, skipna=True)
        assert result == "mixed"

        arr = np.array([True, False, True], dtype=bool)
        result = lib.infer_dtype(arr, skipna=True)
        assert result == "boolean"

        arr = np.array([True, np.nan, False], dtype="O")
        result = lib.infer_dtype(arr, skipna=True)
        assert result == "boolean"

        result = lib.infer_dtype(arr, skipna=False)
        assert result == "mixed"

    def test_floats(self):
        arr = np.array([1.0, 2.0, 3.0, np.float64(4), np.float32(5)], dtype="O")
        result = lib.infer_dtype(arr, skipna=True)
        assert result == "floating"

        arr = np.array([1, 2, 3, np.float64(4), np.float32(5), "foo"], dtype="O")
        result = lib.infer_dtype(arr, skipna=True)
        assert result == "mixed-integer"

        arr = np.array([1, 2, 3, 4, 5], dtype="f4")
        result = lib.infer_dtype(arr, skipna=True)
        assert result == "floating"

        arr = np.array([1, 2, 3, 4, 5], dtype="f8")
        result = lib.infer_dtype(arr, skipna=True)
        assert result == "floating"

    def test_decimals(self):
        # GH15690
        arr = np.array([Decimal(1), Decimal(2), Decimal(3)])
        result = lib.infer_dtype(arr, skipna=True)
        assert result == "decimal"

        arr = np.array([1.0, 2.0, Decimal(3)])
        result = lib.infer_dtype(arr, skipna=True)
        assert result == "mixed"

        arr = np.array([Decimal(1), Decimal("NaN"), Decimal(3)])
        result = lib.infer_dtype(arr, skipna=True)
        assert result == "decimal"

        arr = np.array([Decimal(1), np.nan, Decimal(3)], dtype="O")
        result = lib.infer_dtype(arr, skipna=True)
        assert result == "decimal"

    # complex is compatible with nan, so skipna has no effect
    @pytest.mark.parametrize("skipna", [True, False])
    def test_complex(self, skipna):
        # gets cast to complex on array construction
        arr = np.array([1.0, 2.0, 1 + 1j])
        result = lib.infer_dtype(arr, skipna=skipna)
        assert result == "complex"

        arr = np.array([1.0, 2.0, 1 + 1j], dtype="O")
        result = lib.infer_dtype(arr, skipna=skipna)
        assert result == "mixed"

        # gets cast to complex on array construction
        arr = np.array([1, np.nan, 1 + 1j])
        result = lib.infer_dtype(arr, skipna=skipna)
        assert result == "complex"

        arr = np.array([1.0, np.nan, 1 + 1j], dtype="O")
        result = lib.infer_dtype(arr, skipna=skipna)
        assert result == "mixed"

        # complex with nans stays complex
        arr = np.array([1 + 1j, np.nan, 3 + 3j], dtype="O")
        result = lib.infer_dtype(arr, skipna=skipna)
        assert result == "complex"

        # test smaller complex dtype; will pass through _try_infer_map fastpath
        arr = np.array([1 + 1j, np.nan, 3 + 3j], dtype=np.complex64)
        result = lib.infer_dtype(arr, skipna=skipna)
        assert result == "complex"

    def test_string(self):
        pass

    def test_unicode(self):
        arr = ["a", np.nan, "c"]
        result = lib.infer_dtype(arr, skipna=False)
        # This currently returns "mixed", but it's not clear that's optimal.
        # This could also return "string" or "mixed-string"
        assert result == "mixed"

        arr = ["a", np.nan, "c"]
        result = lib.infer_dtype(arr, skipna=True)
        assert result == "string"

        arr = ["a", "c"]
        result = lib.infer_dtype(arr, skipna=False)
        assert result == "string"

    @pytest.mark.parametrize(
        "dtype, missing, skipna, expected",
        [
            (float, np.nan, False, "floating"),
            (float, np.nan, True, "floating"),
            (object, np.nan, False, "floating"),
            (object, np.nan, True, "empty"),
            (object, None, False, "mixed"),
            (object, None, True, "empty"),
        ],
    )
    @pytest.mark.parametrize("box", [pd.Series, np.array])
    def test_object_empty(self, box, missing, dtype, skipna, expected):
        # GH 23421
        arr = box([missing, missing], dtype=dtype)

        result = lib.infer_dtype(arr, skipna=skipna)
        assert result == expected

    def test_datetime(self):

        dates = [datetime(2012, 1, x) for x in range(1, 20)]
        index = Index(dates)
        assert index.inferred_type == "datetime64"

    def test_infer_dtype_datetime(self):

        arr = np.array([Timestamp("2011-01-01"), Timestamp("2011-01-02")])
        assert lib.infer_dtype(arr, skipna=True) == "datetime"

        arr = np.array(
            [np.datetime64("2011-01-01"), np.datetime64("2011-01-01")], dtype=object
        )
        assert lib.infer_dtype(arr, skipna=True) == "datetime64"

        arr = np.array([datetime(2011, 1, 1), datetime(2012, 2, 1)])
        assert lib.infer_dtype(arr, skipna=True) == "datetime"

        # starts with nan
        for n in [pd.NaT, np.nan]:
            arr = np.array([n, pd.Timestamp("2011-01-02")])
            assert lib.infer_dtype(arr, skipna=True) == "datetime"

            arr = np.array([n, np.datetime64("2011-01-02")])
            assert lib.infer_dtype(arr, skipna=True) == "datetime64"

            arr = np.array([n, datetime(2011, 1, 1)])
            assert lib.infer_dtype(arr, skipna=True) == "datetime"

            arr = np.array([n, pd.Timestamp("2011-01-02"), n])
            assert lib.infer_dtype(arr, skipna=True) == "datetime"

            arr = np.array([n, np.datetime64("2011-01-02"), n])
            assert lib.infer_dtype(arr, skipna=True) == "datetime64"

            arr = np.array([n, datetime(2011, 1, 1), n])
            assert lib.infer_dtype(arr, skipna=True) == "datetime"

        # different type of nat
        arr = np.array(
            [np.timedelta64("nat"), np.datetime64("2011-01-02")], dtype=object
        )
        assert lib.infer_dtype(arr, skipna=False) == "mixed"

        arr = np.array(
            [np.datetime64("2011-01-02"), np.timedelta64("nat")], dtype=object
        )
        assert lib.infer_dtype(arr, skipna=False) == "mixed"

        # mixed datetime
        arr = np.array([datetime(2011, 1, 1), pd.Timestamp("2011-01-02")])
        assert lib.infer_dtype(arr, skipna=True) == "datetime"

        # should be datetime?
        arr = np.array([np.datetime64("2011-01-01"), pd.Timestamp("2011-01-02")])
        assert lib.infer_dtype(arr, skipna=True) == "mixed"

        arr = np.array([pd.Timestamp("2011-01-02"), np.datetime64("2011-01-01")])
        assert lib.infer_dtype(arr, skipna=True) == "mixed"

        arr = np.array([np.nan, pd.Timestamp("2011-01-02"), 1])
        assert lib.infer_dtype(arr, skipna=True) == "mixed-integer"

        arr = np.array([np.nan, pd.Timestamp("2011-01-02"), 1.1])
        assert lib.infer_dtype(arr, skipna=True) == "mixed"

        arr = np.array([np.nan, "2011-01-01", pd.Timestamp("2011-01-02")])
        assert lib.infer_dtype(arr, skipna=True) == "mixed"

    def test_infer_dtype_timedelta(self):

        arr = np.array([pd.Timedelta("1 days"), pd.Timedelta("2 days")])
        assert lib.infer_dtype(arr, skipna=True) == "timedelta"

        arr = np.array([np.timedelta64(1, "D"), np.timedelta64(2, "D")], dtype=object)
        assert lib.infer_dtype(arr, skipna=True) == "timedelta"

        arr = np.array([timedelta(1), timedelta(2)])
        assert lib.infer_dtype(arr, skipna=True) == "timedelta"

        # starts with nan
        for n in [pd.NaT, np.nan]:
            arr = np.array([n, Timedelta("1 days")])
            assert lib.infer_dtype(arr, skipna=True) == "timedelta"

            arr = np.array([n, np.timedelta64(1, "D")])
            assert lib.infer_dtype(arr, skipna=True) == "timedelta"

            arr = np.array([n, timedelta(1)])
            assert lib.infer_dtype(arr, skipna=True) == "timedelta"

            arr = np.array([n, pd.Timedelta("1 days"), n])
            assert lib.infer_dtype(arr, skipna=True) == "timedelta"

            arr = np.array([n, np.timedelta64(1, "D"), n])
            assert lib.infer_dtype(arr, skipna=True) == "timedelta"

            arr = np.array([n, timedelta(1), n])
            assert lib.infer_dtype(arr, skipna=True) == "timedelta"

        # different type of nat
        arr = np.array([np.datetime64("nat"), np.timedelta64(1, "D")], dtype=object)
        assert lib.infer_dtype(arr, skipna=False) == "mixed"

        arr = np.array([np.timedelta64(1, "D"), np.datetime64("nat")], dtype=object)
        assert lib.infer_dtype(arr, skipna=False) == "mixed"

    def test_infer_dtype_period(self):
        # GH 13664
        arr = np.array([pd.Period("2011-01", freq="D"), pd.Period("2011-02", freq="D")])
        assert lib.infer_dtype(arr, skipna=True) == "period"

        arr = np.array([pd.Period("2011-01", freq="D"), pd.Period("2011-02", freq="M")])
        assert lib.infer_dtype(arr, skipna=True) == "period"

        # starts with nan
        for n in [pd.NaT, np.nan]:
            arr = np.array([n, pd.Period("2011-01", freq="D")])
            assert lib.infer_dtype(arr, skipna=True) == "period"

            arr = np.array([n, pd.Period("2011-01", freq="D"), n])
            assert lib.infer_dtype(arr, skipna=True) == "period"

        # different type of nat
        arr = np.array(
            [np.datetime64("nat"), pd.Period("2011-01", freq="M")], dtype=object
        )
        assert lib.infer_dtype(arr, skipna=False) == "mixed"

        arr = np.array(
            [pd.Period("2011-01", freq="M"), np.datetime64("nat")], dtype=object
        )
        assert lib.infer_dtype(arr, skipna=False) == "mixed"

    @pytest.mark.parametrize(
        "data",
        [
            [datetime(2017, 6, 12, 19, 30), datetime(2017, 3, 11, 1, 15)],
            [Timestamp("20170612"), Timestamp("20170311")],
            [
                Timestamp("20170612", tz="US/Eastern"),
                Timestamp("20170311", tz="US/Eastern"),
            ],
            [date(2017, 6, 12), Timestamp("20170311", tz="US/Eastern")],
            [np.datetime64("2017-06-12"), np.datetime64("2017-03-11")],
            [np.datetime64("2017-06-12"), datetime(2017, 3, 11, 1, 15)],
        ],
    )
    def test_infer_datetimelike_array_datetime(self, data):
        assert lib.infer_datetimelike_array(data) == "datetime"

    @pytest.mark.parametrize(
        "data",
        [
            [timedelta(2017, 6, 12), timedelta(2017, 3, 11)],
            [timedelta(2017, 6, 12), date(2017, 3, 11)],
            [np.timedelta64(2017, "D"), np.timedelta64(6, "s")],
            [np.timedelta64(2017, "D"), timedelta(2017, 3, 11)],
        ],
    )
    def test_infer_datetimelike_array_timedelta(self, data):
        assert lib.infer_datetimelike_array(data) == "timedelta"

    def test_infer_datetimelike_array_date(self):
        arr = [date(2017, 6, 12), date(2017, 3, 11)]
        assert lib.infer_datetimelike_array(arr) == "date"

    @pytest.mark.parametrize(
        "data",
        [
            ["2017-06-12", "2017-03-11"],
            [20170612, 20170311],
            [20170612.5, 20170311.8],
            [Dummy(), Dummy()],
            [Timestamp("20170612"), Timestamp("20170311", tz="US/Eastern")],
            [Timestamp("20170612"), 20170311],
            [timedelta(2017, 6, 12), Timestamp("20170311", tz="US/Eastern")],
        ],
    )
    def test_infer_datetimelike_array_mixed(self, data):
        assert lib.infer_datetimelike_array(data) == "mixed"

    @pytest.mark.parametrize(
        "first, expected",
        [
            [[None], "mixed"],
            [[np.nan], "mixed"],
            [[pd.NaT], "nat"],
            [[datetime(2017, 6, 12, 19, 30), pd.NaT], "datetime"],
            [[np.datetime64("2017-06-12"), pd.NaT], "datetime"],
            [[date(2017, 6, 12), pd.NaT], "date"],
            [[timedelta(2017, 6, 12), pd.NaT], "timedelta"],
            [[np.timedelta64(2017, "D"), pd.NaT], "timedelta"],
        ],
    )
    @pytest.mark.parametrize("second", [None, np.nan])
    def test_infer_datetimelike_array_nan_nat_like(self, first, second, expected):
        first.append(second)
        assert lib.infer_datetimelike_array(first) == expected

    def test_infer_dtype_all_nan_nat_like(self):
        arr = np.array([np.nan, np.nan])
        assert lib.infer_dtype(arr, skipna=True) == "floating"

        # nan and None mix are result in mixed
        arr = np.array([np.nan, np.nan, None])
        assert lib.infer_dtype(arr, skipna=True) == "empty"
        assert lib.infer_dtype(arr, skipna=False) == "mixed"

        arr = np.array([None, np.nan, np.nan])
        assert lib.infer_dtype(arr, skipna=True) == "empty"
        assert lib.infer_dtype(arr, skipna=False) == "mixed"

        # pd.NaT
        arr = np.array([pd.NaT])
        assert lib.infer_dtype(arr, skipna=False) == "datetime"

        arr = np.array([pd.NaT, np.nan])
        assert lib.infer_dtype(arr, skipna=False) == "datetime"

        arr = np.array([np.nan, pd.NaT])
        assert lib.infer_dtype(arr, skipna=False) == "datetime"

        arr = np.array([np.nan, pd.NaT, np.nan])
        assert lib.infer_dtype(arr, skipna=False) == "datetime"

        arr = np.array([None, pd.NaT, None])
        assert lib.infer_dtype(arr, skipna=False) == "datetime"

        # np.datetime64(nat)
        arr = np.array([np.datetime64("nat")])
        assert lib.infer_dtype(arr, skipna=False) == "datetime64"

        for n in [np.nan, pd.NaT, None]:
            arr = np.array([n, np.datetime64("nat"), n])
            assert lib.infer_dtype(arr, skipna=False) == "datetime64"

            arr = np.array([pd.NaT, n, np.datetime64("nat"), n])
            assert lib.infer_dtype(arr, skipna=False) == "datetime64"

        arr = np.array([np.timedelta64("nat")], dtype=object)
        assert lib.infer_dtype(arr, skipna=False) == "timedelta"

        for n in [np.nan, pd.NaT, None]:
            arr = np.array([n, np.timedelta64("nat"), n])
            assert lib.infer_dtype(arr, skipna=False) == "timedelta"

            arr = np.array([pd.NaT, n, np.timedelta64("nat"), n])
            assert lib.infer_dtype(arr, skipna=False) == "timedelta"

        # datetime / timedelta mixed
        arr = np.array([pd.NaT, np.datetime64("nat"), np.timedelta64("nat"), np.nan])
        assert lib.infer_dtype(arr, skipna=False) == "mixed"

        arr = np.array([np.timedelta64("nat"), np.datetime64("nat")], dtype=object)
        assert lib.infer_dtype(arr, skipna=False) == "mixed"

    def test_is_datetimelike_array_all_nan_nat_like(self):
        arr = np.array([np.nan, pd.NaT, np.datetime64("nat")])
        assert lib.is_datetime_array(arr)
        assert lib.is_datetime64_array(arr)
        assert not lib.is_timedelta_or_timedelta64_array(arr)

        arr = np.array([np.nan, pd.NaT, np.timedelta64("nat")])
        assert not lib.is_datetime_array(arr)
        assert not lib.is_datetime64_array(arr)
        assert lib.is_timedelta_or_timedelta64_array(arr)

        arr = np.array([np.nan, pd.NaT, np.datetime64("nat"), np.timedelta64("nat")])
        assert not lib.is_datetime_array(arr)
        assert not lib.is_datetime64_array(arr)
        assert not lib.is_timedelta_or_timedelta64_array(arr)

        arr = np.array([np.nan, pd.NaT])
        assert lib.is_datetime_array(arr)
        assert lib.is_datetime64_array(arr)
        assert lib.is_timedelta_or_timedelta64_array(arr)

        arr = np.array([np.nan, np.nan], dtype=object)
        assert not lib.is_datetime_array(arr)
        assert not lib.is_datetime64_array(arr)
        assert not lib.is_timedelta_or_timedelta64_array(arr)

        assert lib.is_datetime_with_singletz_array(
            np.array(
                [
                    pd.Timestamp("20130101", tz="US/Eastern"),
                    pd.Timestamp("20130102", tz="US/Eastern"),
                ],
                dtype=object,
            )
        )
        assert not lib.is_datetime_with_singletz_array(
            np.array(
                [
                    pd.Timestamp("20130101", tz="US/Eastern"),
                    pd.Timestamp("20130102", tz="CET"),
                ],
                dtype=object,
            )
        )

    @pytest.mark.parametrize(
        "func",
        [
            "is_datetime_array",
            "is_datetime64_array",
            "is_bool_array",
            "is_timedelta_or_timedelta64_array",
            "is_date_array",
            "is_time_array",
            "is_interval_array",
            "is_period_array",
        ],
    )
    def test_other_dtypes_for_array(self, func):
        func = getattr(lib, func)
        arr = np.array(["foo", "bar"])
        assert not func(arr)

        arr = np.array([1, 2])
        assert not func(arr)

    def test_date(self):

        dates = [date(2012, 1, day) for day in range(1, 20)]
        index = Index(dates)
        assert index.inferred_type == "date"

        dates = [date(2012, 1, day) for day in range(1, 20)] + [np.nan]
        result = lib.infer_dtype(dates, skipna=False)
        assert result == "mixed"

        result = lib.infer_dtype(dates, skipna=True)
        assert result == "date"

    def test_is_numeric_array(self):

        assert lib.is_float_array(np.array([1, 2.0]))
        assert lib.is_float_array(np.array([1, 2.0, np.nan]))
        assert not lib.is_float_array(np.array([1, 2]))

        assert lib.is_integer_array(np.array([1, 2]))
        assert not lib.is_integer_array(np.array([1, 2.0]))

    def test_is_string_array(self):

        assert lib.is_string_array(np.array(["foo", "bar"]))
        assert not lib.is_string_array(
            np.array(["foo", "bar", pd.NA], dtype=object), skipna=False
        )
        assert lib.is_string_array(
            np.array(["foo", "bar", pd.NA], dtype=object), skipna=True
        )
        # NaN is not valid for string array, just NA
        assert not lib.is_string_array(
            np.array(["foo", "bar", np.nan], dtype=object), skipna=True
        )

        assert not lib.is_string_array(np.array([1, 2]))

    def test_to_object_array_tuples(self):
        r = (5, 6)
        values = [r]
        lib.to_object_array_tuples(values)

        # make sure record array works
        record = namedtuple("record", "x y")
        r = record(5, 6)
        values = [r]
        lib.to_object_array_tuples(values)

    def test_object(self):

        # GH 7431
        # cannot infer more than this as only a single element
        arr = np.array([None], dtype="O")
        result = lib.infer_dtype(arr, skipna=False)
        assert result == "mixed"
        result = lib.infer_dtype(arr, skipna=True)
        assert result == "empty"

    def test_to_object_array_width(self):
        # see gh-13320
        rows = [[1, 2, 3], [4, 5, 6]]

        expected = np.array(rows, dtype=object)
        out = lib.to_object_array(rows)
        tm.assert_numpy_array_equal(out, expected)

        expected = np.array(rows, dtype=object)
        out = lib.to_object_array(rows, min_width=1)
        tm.assert_numpy_array_equal(out, expected)

        expected = np.array(
            [[1, 2, 3, None, None], [4, 5, 6, None, None]], dtype=object
        )
        out = lib.to_object_array(rows, min_width=5)
        tm.assert_numpy_array_equal(out, expected)

    def test_is_period(self):
        assert lib.is_period(pd.Period("2011-01", freq="M"))
        assert not lib.is_period(pd.PeriodIndex(["2011-01"], freq="M"))
        assert not lib.is_period(pd.Timestamp("2011-01"))
        assert not lib.is_period(1)
        assert not lib.is_period(np.nan)

    def test_categorical(self):

        # GH 8974
        arr = Categorical(list("abc"))
        result = lib.infer_dtype(arr, skipna=True)
        assert result == "categorical"

        result = lib.infer_dtype(Series(arr), skipna=True)
        assert result == "categorical"

        arr = Categorical(list("abc"), categories=["cegfab"], ordered=True)
        result = lib.infer_dtype(arr, skipna=True)
        assert result == "categorical"

        result = lib.infer_dtype(Series(arr), skipna=True)
        assert result == "categorical"

    def test_interval(self):
        idx = pd.IntervalIndex.from_breaks(range(5), closed="both")
        inferred = lib.infer_dtype(idx, skipna=False)
        assert inferred == "interval"

        inferred = lib.infer_dtype(idx._data, skipna=False)
        assert inferred == "interval"

        inferred = lib.infer_dtype(pd.Series(idx), skipna=False)
        assert inferred == "interval"

    @pytest.mark.parametrize("klass", [pd.array, pd.Series])
    @pytest.mark.parametrize("skipna", [True, False])
    @pytest.mark.parametrize("data", [["a", "b", "c"], ["a", "b", pd.NA]])
    def test_string_dtype(self, data, skipna, klass):
        # StringArray
        val = klass(data, dtype="string")
        inferred = lib.infer_dtype(val, skipna=skipna)
        assert inferred == "string"

    @pytest.mark.parametrize("klass", [pd.array, pd.Series])
    @pytest.mark.parametrize("skipna", [True, False])
    @pytest.mark.parametrize("data", [[True, False, True], [True, False, pd.NA]])
    def test_boolean_dtype(self, data, skipna, klass):
        # BooleanArray
        val = klass(data, dtype="boolean")
        inferred = lib.infer_dtype(val, skipna=skipna)
        assert inferred == "boolean"


class TestNumberScalar:
    def test_is_number(self):

        assert is_number(True)
        assert is_number(1)
        assert is_number(1.1)
        assert is_number(1 + 3j)
        assert is_number(np.bool(False))
        assert is_number(np.int64(1))
        assert is_number(np.float64(1.1))
        assert is_number(np.complex128(1 + 3j))
        assert is_number(np.nan)

        assert not is_number(None)
        assert not is_number("x")
        assert not is_number(datetime(2011, 1, 1))
        assert not is_number(np.datetime64("2011-01-01"))
        assert not is_number(Timestamp("2011-01-01"))
        assert not is_number(Timestamp("2011-01-01", tz="US/Eastern"))
        assert not is_number(timedelta(1000))
        assert not is_number(Timedelta("1 days"))

        # questionable
        assert not is_number(np.bool_(False))
        assert is_number(np.timedelta64(1, "D"))

    def test_is_bool(self):
        assert is_bool(True)
        assert is_bool(np.bool(False))
        assert is_bool(np.bool_(False))

        assert not is_bool(1)
        assert not is_bool(1.1)
        assert not is_bool(1 + 3j)
        assert not is_bool(np.int64(1))
        assert not is_bool(np.float64(1.1))
        assert not is_bool(np.complex128(1 + 3j))
        assert not is_bool(np.nan)
        assert not is_bool(None)
        assert not is_bool("x")
        assert not is_bool(datetime(2011, 1, 1))
        assert not is_bool(np.datetime64("2011-01-01"))
        assert not is_bool(Timestamp("2011-01-01"))
        assert not is_bool(Timestamp("2011-01-01", tz="US/Eastern"))
        assert not is_bool(timedelta(1000))
        assert not is_bool(np.timedelta64(1, "D"))
        assert not is_bool(Timedelta("1 days"))

    def test_is_integer(self):
        assert is_integer(1)
        assert is_integer(np.int64(1))

        assert not is_integer(True)
        assert not is_integer(1.1)
        assert not is_integer(1 + 3j)
        assert not is_integer(np.bool(False))
        assert not is_integer(np.bool_(False))
        assert not is_integer(np.float64(1.1))
        assert not is_integer(np.complex128(1 + 3j))
        assert not is_integer(np.nan)
        assert not is_integer(None)
        assert not is_integer("x")
        assert not is_integer(datetime(2011, 1, 1))
        assert not is_integer(np.datetime64("2011-01-01"))
        assert not is_integer(Timestamp("2011-01-01"))
        assert not is_integer(Timestamp("2011-01-01", tz="US/Eastern"))
        assert not is_integer(timedelta(1000))
        assert not is_integer(Timedelta("1 days"))
        assert not is_integer(np.timedelta64(1, "D"))

    def test_is_float(self):
        assert is_float(1.1)
        assert is_float(np.float64(1.1))
        assert is_float(np.nan)

        assert not is_float(True)
        assert not is_float(1)
        assert not is_float(1 + 3j)
        assert not is_float(np.bool(False))
        assert not is_float(np.bool_(False))
        assert not is_float(np.int64(1))
        assert not is_float(np.complex128(1 + 3j))
        assert not is_float(None)
        assert not is_float("x")
        assert not is_float(datetime(2011, 1, 1))
        assert not is_float(np.datetime64("2011-01-01"))
        assert not is_float(Timestamp("2011-01-01"))
        assert not is_float(Timestamp("2011-01-01", tz="US/Eastern"))
        assert not is_float(timedelta(1000))
        assert not is_float(np.timedelta64(1, "D"))
        assert not is_float(Timedelta("1 days"))

    def test_is_datetime_dtypes(self):

        ts = pd.date_range("20130101", periods=3)
        tsa = pd.date_range("20130101", periods=3, tz="US/Eastern")

        assert is_datetime64_dtype("datetime64")
        assert is_datetime64_dtype("datetime64[ns]")
        assert is_datetime64_dtype(ts)
        assert not is_datetime64_dtype(tsa)

        assert not is_datetime64_ns_dtype("datetime64")
        assert is_datetime64_ns_dtype("datetime64[ns]")
        assert is_datetime64_ns_dtype(ts)
        assert is_datetime64_ns_dtype(tsa)

        assert is_datetime64_any_dtype("datetime64")
        assert is_datetime64_any_dtype("datetime64[ns]")
        assert is_datetime64_any_dtype(ts)
        assert is_datetime64_any_dtype(tsa)

        assert not is_datetime64tz_dtype("datetime64")
        assert not is_datetime64tz_dtype("datetime64[ns]")
        assert not is_datetime64tz_dtype(ts)
        assert is_datetime64tz_dtype(tsa)

        for tz in ["US/Eastern", "UTC"]:
            dtype = f"datetime64[ns, {tz}]"
            assert not is_datetime64_dtype(dtype)
            assert is_datetime64tz_dtype(dtype)
            assert is_datetime64_ns_dtype(dtype)
            assert is_datetime64_any_dtype(dtype)

    def test_is_timedelta(self):
        assert is_timedelta64_dtype("timedelta64")
        assert is_timedelta64_dtype("timedelta64[ns]")
        assert not is_timedelta64_ns_dtype("timedelta64")
        assert is_timedelta64_ns_dtype("timedelta64[ns]")

        tdi = TimedeltaIndex([1e14, 2e14], dtype="timedelta64[ns]")
        assert is_timedelta64_dtype(tdi)
        assert is_timedelta64_ns_dtype(tdi)
        assert is_timedelta64_ns_dtype(tdi.astype("timedelta64[ns]"))

        # Conversion to Int64Index:
        assert not is_timedelta64_ns_dtype(tdi.astype("timedelta64"))
        assert not is_timedelta64_ns_dtype(tdi.astype("timedelta64[h]"))


class TestIsScalar:
    def test_is_scalar_builtin_scalars(self):
        assert is_scalar(None)
        assert is_scalar(True)
        assert is_scalar(False)
        assert is_scalar(Fraction())
        assert is_scalar(0.0)
        assert is_scalar(1)
        assert is_scalar(complex(2))
        assert is_scalar(float("NaN"))
        assert is_scalar(np.nan)
        assert is_scalar("foobar")
        assert is_scalar(b"foobar")
        assert is_scalar(datetime(2014, 1, 1))
        assert is_scalar(date(2014, 1, 1))
        assert is_scalar(time(12, 0))
        assert is_scalar(timedelta(hours=1))
        assert is_scalar(pd.NaT)
        assert is_scalar(pd.NA)

    def test_is_scalar_builtin_nonscalars(self):
        assert not is_scalar({})
        assert not is_scalar([])
        assert not is_scalar([1])
        assert not is_scalar(())
        assert not is_scalar((1,))
        assert not is_scalar(slice(None))
        assert not is_scalar(Ellipsis)

    def test_is_scalar_numpy_array_scalars(self):
        assert is_scalar(np.int64(1))
        assert is_scalar(np.float64(1.0))
        assert is_scalar(np.int32(1))
        assert is_scalar(np.complex64(2))
        assert is_scalar(np.object_("foobar"))
        assert is_scalar(np.str_("foobar"))
        assert is_scalar(np.unicode_("foobar"))
        assert is_scalar(np.bytes_(b"foobar"))
        assert is_scalar(np.datetime64("2014-01-01"))
        assert is_scalar(np.timedelta64(1, "h"))

    def test_is_scalar_numpy_zerodim_arrays(self):
        for zerodim in [
            np.array(1),
            np.array("foobar"),
            np.array(np.datetime64("2014-01-01")),
            np.array(np.timedelta64(1, "h")),
            np.array(np.datetime64("NaT")),
        ]:
            assert not is_scalar(zerodim)
            assert is_scalar(lib.item_from_zerodim(zerodim))

    @pytest.mark.filterwarnings("ignore::PendingDeprecationWarning")
    def test_is_scalar_numpy_arrays(self):
        assert not is_scalar(np.array([]))
        assert not is_scalar(np.array([[]]))
        assert not is_scalar(np.matrix("1; 2"))

    def test_is_scalar_pandas_scalars(self):
        assert is_scalar(Timestamp("2014-01-01"))
        assert is_scalar(Timedelta(hours=1))
        assert is_scalar(Period("2014-01-01"))
        assert is_scalar(Interval(left=0, right=1))
        assert is_scalar(DateOffset(days=1))

    def test_is_scalar_pandas_containers(self):
        assert not is_scalar(Series(dtype=object))
        assert not is_scalar(Series([1]))
        assert not is_scalar(DataFrame())
        assert not is_scalar(DataFrame([[1]]))
        assert not is_scalar(Index([]))
        assert not is_scalar(Index([1]))

    def test_is_scalar_number(self):
        # Number() is not recognied by PyNumber_Check, so by extension
        #  is not recognized by is_scalar, but instances of non-abstract
        #  subclasses are.

        class Numeric(Number):
            def __init__(self, value):
                self.value = value

            def __int__(self):
                return self.value

        num = Numeric(1)
        assert is_scalar(num)


def test_datetimeindex_from_empty_datetime64_array():
    for unit in ["ms", "us", "ns"]:
        idx = DatetimeIndex(np.array([], dtype=f"datetime64[{unit}]"))
        assert len(idx) == 0


def test_nan_to_nat_conversions():

    df = DataFrame(
        dict({"A": np.asarray(range(10), dtype="float64"), "B": Timestamp("20010101")})
    )
    df.iloc[3:6, :] = np.nan
    result = df.loc[4, "B"].value
    assert result == iNaT

    s = df["B"].copy()
    s._data = s._data.setitem(indexer=tuple([slice(8, 9)]), value=np.nan)
    assert isna(s[8])

    assert s[8].value == np.datetime64("NaT").astype(np.int64)


@td.skip_if_no_scipy
@pytest.mark.filterwarnings("ignore::PendingDeprecationWarning")
def test_is_scipy_sparse(spmatrix):  # noqa: F811
    assert is_scipy_sparse(spmatrix([[0, 1]]))
    assert not is_scipy_sparse(np.array([1]))


def test_ensure_int32():
    values = np.arange(10, dtype=np.int32)
    result = ensure_int32(values)
    assert result.dtype == np.int32

    values = np.arange(10, dtype=np.int64)
    result = ensure_int32(values)
    assert result.dtype == np.int32


def test_ensure_categorical():
    values = np.arange(10, dtype=np.int32)
    result = ensure_categorical(values)
    assert result.dtype == "category"

    values = Categorical(values)
    result = ensure_categorical(values)
    tm.assert_categorical_equal(result, values)
