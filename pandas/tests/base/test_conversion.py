import numpy as np
import pytest

from pandas.core.dtypes.common import is_datetime64_dtype, is_timedelta64_dtype
from pandas.core.dtypes.dtypes import DatetimeTZDtype

import pandas as pd
from pandas import CategoricalIndex, Series, Timedelta, Timestamp
from pandas.core.arrays import DatetimeArray, PandasArray, TimedeltaArray
import pandas.util.testing as tm


class TestToIterable:
    # test that we convert an iterable to python types

    dtypes = [
        ("int8", int),
        ("int16", int),
        ("int32", int),
        ("int64", int),
        ("uint8", int),
        ("uint16", int),
        ("uint32", int),
        ("uint64", int),
        ("float16", float),
        ("float32", float),
        ("float64", float),
        ("datetime64[ns]", Timestamp),
        ("datetime64[ns, US/Eastern]", Timestamp),
        ("timedelta64[ns]", Timedelta),
    ]

    @pytest.mark.parametrize("dtype, rdtype", dtypes)
    @pytest.mark.parametrize(
        "method",
        [
            lambda x: x.tolist(),
            lambda x: x.to_list(),
            lambda x: list(x),
            lambda x: list(x.__iter__()),
        ],
        ids=["tolist", "to_list", "list", "iter"],
    )
    @pytest.mark.filterwarnings("ignore:\\n    Passing:FutureWarning")
    # TODO(GH-24559): Remove the filterwarnings
    def test_iterable(self, index_or_series, method, dtype, rdtype):
        # gh-10904
        # gh-13258
        # coerce iteration to underlying python / pandas types
        typ = index_or_series
        s = typ([1], dtype=dtype)
        result = method(s)[0]
        assert isinstance(result, rdtype)

    @pytest.mark.parametrize(
        "dtype, rdtype, obj",
        [
            ("object", object, "a"),
            ("object", int, 1),
            ("category", object, "a"),
            ("category", int, 1),
        ],
    )
    @pytest.mark.parametrize(
        "method",
        [
            lambda x: x.tolist(),
            lambda x: x.to_list(),
            lambda x: list(x),
            lambda x: list(x.__iter__()),
        ],
        ids=["tolist", "to_list", "list", "iter"],
    )
    def test_iterable_object_and_category(
        self, index_or_series, method, dtype, rdtype, obj
    ):
        # gh-10904
        # gh-13258
        # coerce iteration to underlying python / pandas types
        typ = index_or_series
        s = typ([obj], dtype=dtype)
        result = method(s)[0]
        assert isinstance(result, rdtype)

    @pytest.mark.parametrize("dtype, rdtype", dtypes)
    def test_iterable_items(self, dtype, rdtype):
        # gh-13258
        # test if items yields the correct boxed scalars
        # this only applies to series
        s = Series([1], dtype=dtype)
        _, result = list(s.items())[0]
        assert isinstance(result, rdtype)

        _, result = list(s.items())[0]
        assert isinstance(result, rdtype)

    @pytest.mark.parametrize(
        "dtype, rdtype", dtypes + [("object", int), ("category", int)]
    )
    @pytest.mark.filterwarnings("ignore:\\n    Passing:FutureWarning")
    # TODO(GH-24559): Remove the filterwarnings
    def test_iterable_map(self, index_or_series, dtype, rdtype):
        # gh-13236
        # coerce iteration to underlying python / pandas types
        typ = index_or_series
        s = typ([1], dtype=dtype)
        result = s.map(type)[0]
        if not isinstance(rdtype, tuple):
            rdtype = tuple([rdtype])
        assert result in rdtype

    @pytest.mark.parametrize(
        "method",
        [
            lambda x: x.tolist(),
            lambda x: x.to_list(),
            lambda x: list(x),
            lambda x: list(x.__iter__()),
        ],
        ids=["tolist", "to_list", "list", "iter"],
    )
    def test_categorial_datetimelike(self, method):
        i = CategoricalIndex([Timestamp("1999-12-31"), Timestamp("2000-12-31")])

        result = method(i)[0]
        assert isinstance(result, Timestamp)

    def test_iter_box(self):
        vals = [Timestamp("2011-01-01"), Timestamp("2011-01-02")]
        s = Series(vals)
        assert s.dtype == "datetime64[ns]"
        for res, exp in zip(s, vals):
            assert isinstance(res, Timestamp)
            assert res.tz is None
            assert res == exp

        vals = [
            Timestamp("2011-01-01", tz="US/Eastern"),
            Timestamp("2011-01-02", tz="US/Eastern"),
        ]
        s = Series(vals)

        assert s.dtype == "datetime64[ns, US/Eastern]"
        for res, exp in zip(s, vals):
            assert isinstance(res, Timestamp)
            assert res.tz == exp.tz
            assert res == exp

        # timedelta
        vals = [Timedelta("1 days"), Timedelta("2 days")]
        s = Series(vals)
        assert s.dtype == "timedelta64[ns]"
        for res, exp in zip(s, vals):
            assert isinstance(res, Timedelta)
            assert res == exp

        # period
        vals = [pd.Period("2011-01-01", freq="M"), pd.Period("2011-01-02", freq="M")]
        s = Series(vals)
        assert s.dtype == "Period[M]"
        for res, exp in zip(s, vals):
            assert isinstance(res, pd.Period)
            assert res.freq == "M"
            assert res == exp


@pytest.mark.parametrize(
    "array, expected_type, dtype",
    [
        (np.array([0, 1], dtype=np.int64), np.ndarray, "int64"),
        (np.array(["a", "b"]), np.ndarray, "object"),
        (pd.Categorical(["a", "b"]), pd.Categorical, "category"),
        (
            pd.DatetimeIndex(["2017", "2018"], tz="US/Central"),
            DatetimeArray,
            "datetime64[ns, US/Central]",
        ),
        (
            pd.PeriodIndex([2018, 2019], freq="A"),
            pd.core.arrays.PeriodArray,
            pd.core.dtypes.dtypes.PeriodDtype("A-DEC"),
        ),
        (
            pd.IntervalIndex.from_breaks([0, 1, 2]),
            pd.core.arrays.IntervalArray,
            "interval",
        ),
        # This test is currently failing for datetime64[ns] and timedelta64[ns].
        # The NumPy type system is sufficient for representing these types, so
        # we just use NumPy for Series / DataFrame columns of these types (so
        # we get consolidation and so on).
        # However, DatetimeIndex and TimedeltaIndex use the DateLikeArray
        # abstraction to for code reuse.
        # At the moment, we've judged that allowing this test to fail is more
        # practical that overriding Series._values to special case
        # Series[M8[ns]] and Series[m8[ns]] to return a DateLikeArray.
        pytest.param(
            pd.DatetimeIndex(["2017", "2018"]),
            np.ndarray,
            "datetime64[ns]",
            marks=[pytest.mark.xfail(reason="datetime _values", strict=True)],
        ),
        pytest.param(
            pd.TimedeltaIndex([10 ** 10]),
            np.ndarray,
            "m8[ns]",
            marks=[pytest.mark.xfail(reason="timedelta _values", strict=True)],
        ),
    ],
)
def test_values_consistent(array, expected_type, dtype):
    l_values = pd.Series(array)._values
    r_values = pd.Index(array)._values
    assert type(l_values) is expected_type
    assert type(l_values) is type(r_values)

    tm.assert_equal(l_values, r_values)


@pytest.mark.parametrize(
    "array, expected",
    [
        (np.array([0, 1], dtype=np.int64), np.array([0, 1], dtype=np.int64)),
        (np.array(["0", "1"]), np.array(["0", "1"], dtype=object)),
        (pd.Categorical(["a", "a"]), np.array([0, 0], dtype="int8")),
        (
            pd.DatetimeIndex(["2017-01-01T00:00:00"]),
            np.array(["2017-01-01T00:00:00"], dtype="M8[ns]"),
        ),
        (
            pd.DatetimeIndex(["2017-01-01T00:00:00"], tz="US/Eastern"),
            np.array(["2017-01-01T05:00:00"], dtype="M8[ns]"),
        ),
        (pd.TimedeltaIndex([10 ** 10]), np.array([10 ** 10], dtype="m8[ns]")),
        (
            pd.PeriodIndex(["2017", "2018"], freq="D"),
            np.array([17167, 17532], dtype=np.int64),
        ),
    ],
)
def test_ndarray_values(array, expected):
    l_values = pd.Series(array)._ndarray_values
    r_values = pd.Index(array)._ndarray_values
    tm.assert_numpy_array_equal(l_values, r_values)
    tm.assert_numpy_array_equal(l_values, expected)


@pytest.mark.parametrize("arr", [np.array([1, 2, 3])])
def test_numpy_array(arr):
    ser = pd.Series(arr)
    result = ser.array
    expected = PandasArray(arr)
    tm.assert_extension_array_equal(result, expected)


def test_numpy_array_all_dtypes(any_numpy_dtype):
    ser = pd.Series(dtype=any_numpy_dtype)
    result = ser.array
    if is_datetime64_dtype(any_numpy_dtype):
        assert isinstance(result, DatetimeArray)
    elif is_timedelta64_dtype(any_numpy_dtype):
        assert isinstance(result, TimedeltaArray)
    else:
        assert isinstance(result, PandasArray)


@pytest.mark.parametrize(
    "array, attr",
    [
        (pd.Categorical(["a", "b"]), "_codes"),
        (pd.core.arrays.period_array(["2000", "2001"], freq="D"), "_data"),
        (pd.core.arrays.integer_array([0, np.nan]), "_data"),
        (pd.core.arrays.IntervalArray.from_breaks([0, 1]), "_left"),
        (pd.SparseArray([0, 1]), "_sparse_values"),
        (DatetimeArray(np.array([1, 2], dtype="datetime64[ns]")), "_data"),
        # tz-aware Datetime
        (
            DatetimeArray(
                np.array(
                    ["2000-01-01T12:00:00", "2000-01-02T12:00:00"], dtype="M8[ns]"
                ),
                dtype=DatetimeTZDtype(tz="US/Central"),
            ),
            "_data",
        ),
    ],
)
def test_array(array, attr, index_or_series):
    box = index_or_series
    if array.dtype.name in ("Int64", "Sparse[int64, 0]") and box is pd.Index:
        pytest.skip("No index type for {}".format(array.dtype))
    result = box(array, copy=False).array

    if attr:
        array = getattr(array, attr)
        result = getattr(result, attr)

    assert result is array


def test_array_multiindex_raises():
    idx = pd.MultiIndex.from_product([["A"], ["a", "b"]])
    with pytest.raises(ValueError, match="MultiIndex"):
        idx.array


@pytest.mark.parametrize(
    "array, expected",
    [
        (np.array([1, 2], dtype=np.int64), np.array([1, 2], dtype=np.int64)),
        (pd.Categorical(["a", "b"]), np.array(["a", "b"], dtype=object)),
        (
            pd.core.arrays.period_array(["2000", "2001"], freq="D"),
            np.array([pd.Period("2000", freq="D"), pd.Period("2001", freq="D")]),
        ),
        (
            pd.core.arrays.integer_array([0, np.nan]),
            np.array([0, pd.NA], dtype=object),
        ),
        (
            pd.core.arrays.IntervalArray.from_breaks([0, 1, 2]),
            np.array([pd.Interval(0, 1), pd.Interval(1, 2)], dtype=object),
        ),
        (pd.SparseArray([0, 1]), np.array([0, 1], dtype=np.int64)),
        # tz-naive datetime
        (
            DatetimeArray(np.array(["2000", "2001"], dtype="M8[ns]")),
            np.array(["2000", "2001"], dtype="M8[ns]"),
        ),
        # tz-aware stays tz`-aware
        (
            DatetimeArray(
                np.array(
                    ["2000-01-01T06:00:00", "2000-01-02T06:00:00"], dtype="M8[ns]"
                ),
                dtype=DatetimeTZDtype(tz="US/Central"),
            ),
            np.array(
                [
                    pd.Timestamp("2000-01-01", tz="US/Central"),
                    pd.Timestamp("2000-01-02", tz="US/Central"),
                ]
            ),
        ),
        # Timedelta
        (
            TimedeltaArray(np.array([0, 3600000000000], dtype="i8"), freq="H"),
            np.array([0, 3600000000000], dtype="m8[ns]"),
        ),
    ],
)
def test_to_numpy(array, expected, index_or_series):
    box = index_or_series
    thing = box(array)

    if array.dtype.name in ("Int64", "Sparse[int64, 0]") and box is pd.Index:
        pytest.skip("No index type for {}".format(array.dtype))

    result = thing.to_numpy()
    tm.assert_numpy_array_equal(result, expected)


@pytest.mark.parametrize("as_series", [True, False])
@pytest.mark.parametrize(
    "arr", [np.array([1, 2, 3], dtype="int64"), np.array(["a", "b", "c"], dtype=object)]
)
def test_to_numpy_copy(arr, as_series):
    obj = pd.Index(arr, copy=False)
    if as_series:
        obj = pd.Series(obj.values, copy=False)

    # no copy by default
    result = obj.to_numpy()
    assert np.shares_memory(arr, result) is True

    result = obj.to_numpy(copy=False)
    assert np.shares_memory(arr, result) is True

    # copy=True
    result = obj.to_numpy(copy=True)
    assert np.shares_memory(arr, result) is False


@pytest.mark.parametrize("as_series", [True, False])
def test_to_numpy_dtype(as_series):
    tz = "US/Eastern"
    obj = pd.DatetimeIndex(["2000", "2001"], tz=tz)
    if as_series:
        obj = pd.Series(obj)

    # preserve tz by default
    result = obj.to_numpy()
    expected = np.array(
        [pd.Timestamp("2000", tz=tz), pd.Timestamp("2001", tz=tz)], dtype=object
    )
    tm.assert_numpy_array_equal(result, expected)

    result = obj.to_numpy(dtype="object")
    tm.assert_numpy_array_equal(result, expected)

    result = obj.to_numpy(dtype="M8[ns]")
    expected = np.array(["2000-01-01T05", "2001-01-01T05"], dtype="M8[ns]")
    tm.assert_numpy_array_equal(result, expected)
