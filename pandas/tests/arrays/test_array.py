import datetime
import decimal

import numpy as np
import pytest
import pytz

from pandas.core.dtypes.dtypes import registry

import pandas as pd
import pandas._testing as tm
from pandas.api.extensions import register_extension_dtype
from pandas.api.types import is_scalar
from pandas.core.arrays import PandasArray, integer_array, period_array
from pandas.tests.extension.decimal import DecimalArray, DecimalDtype, to_decimal


@pytest.mark.parametrize(
    "data, dtype, expected",
    [
        # Basic NumPy defaults.
        ([1, 2], None, pd.arrays.IntegerArray._from_sequence([1, 2])),
        ([1, 2], object, PandasArray(np.array([1, 2], dtype=object))),
        (
            [1, 2],
            np.dtype("float32"),
            PandasArray(np.array([1.0, 2.0], dtype=np.dtype("float32"))),
        ),
        (
            np.array([1, 2], dtype="int64"),
            None,
            pd.arrays.IntegerArray._from_sequence([1, 2]),
        ),
        # String alias passes through to NumPy
        ([1, 2], "float32", PandasArray(np.array([1, 2], dtype="float32"))),
        # Period alias
        (
            [pd.Period("2000", "D"), pd.Period("2001", "D")],
            "Period[D]",
            period_array(["2000", "2001"], freq="D"),
        ),
        # Period dtype
        (
            [pd.Period("2000", "D")],
            pd.PeriodDtype("D"),
            period_array(["2000"], freq="D"),
        ),
        # Datetime (naive)
        (
            [1, 2],
            np.dtype("datetime64[ns]"),
            pd.arrays.DatetimeArray._from_sequence(
                np.array([1, 2], dtype="datetime64[ns]")
            ),
        ),
        (
            np.array([1, 2], dtype="datetime64[ns]"),
            None,
            pd.arrays.DatetimeArray._from_sequence(
                np.array([1, 2], dtype="datetime64[ns]")
            ),
        ),
        (
            pd.DatetimeIndex(["2000", "2001"]),
            np.dtype("datetime64[ns]"),
            pd.arrays.DatetimeArray._from_sequence(["2000", "2001"]),
        ),
        (
            pd.DatetimeIndex(["2000", "2001"]),
            None,
            pd.arrays.DatetimeArray._from_sequence(["2000", "2001"]),
        ),
        (
            ["2000", "2001"],
            np.dtype("datetime64[ns]"),
            pd.arrays.DatetimeArray._from_sequence(["2000", "2001"]),
        ),
        # Datetime (tz-aware)
        (
            ["2000", "2001"],
            pd.DatetimeTZDtype(tz="CET"),
            pd.arrays.DatetimeArray._from_sequence(
                ["2000", "2001"], dtype=pd.DatetimeTZDtype(tz="CET")
            ),
        ),
        # Timedelta
        (
            ["1H", "2H"],
            np.dtype("timedelta64[ns]"),
            pd.arrays.TimedeltaArray._from_sequence(["1H", "2H"]),
        ),
        (
            pd.TimedeltaIndex(["1H", "2H"]),
            np.dtype("timedelta64[ns]"),
            pd.arrays.TimedeltaArray._from_sequence(["1H", "2H"]),
        ),
        (
            pd.TimedeltaIndex(["1H", "2H"]),
            None,
            pd.arrays.TimedeltaArray._from_sequence(["1H", "2H"]),
        ),
        # Category
        (["a", "b"], "category", pd.Categorical(["a", "b"])),
        (
            ["a", "b"],
            pd.CategoricalDtype(None, ordered=True),
            pd.Categorical(["a", "b"], ordered=True),
        ),
        # Interval
        (
            [pd.Interval(1, 2), pd.Interval(3, 4)],
            "interval",
            pd.arrays.IntervalArray.from_tuples([(1, 2), (3, 4)]),
        ),
        # Sparse
        ([0, 1], "Sparse[int64]", pd.SparseArray([0, 1], dtype="int64")),
        # IntegerNA
        ([1, None], "Int16", integer_array([1, None], dtype="Int16")),
        (pd.Series([1, 2]), None, PandasArray(np.array([1, 2], dtype=np.int64))),
        # String
        (["a", None], "string", pd.arrays.StringArray._from_sequence(["a", None])),
        (
            ["a", None],
            pd.StringDtype(),
            pd.arrays.StringArray._from_sequence(["a", None]),
        ),
        # Boolean
        ([True, None], "boolean", pd.arrays.BooleanArray._from_sequence([True, None])),
        (
            [True, None],
            pd.BooleanDtype(),
            pd.arrays.BooleanArray._from_sequence([True, None]),
        ),
        # Index
        (pd.Index([1, 2]), None, PandasArray(np.array([1, 2], dtype=np.int64))),
        # Series[EA] returns the EA
        (
            pd.Series(pd.Categorical(["a", "b"], categories=["a", "b", "c"])),
            None,
            pd.Categorical(["a", "b"], categories=["a", "b", "c"]),
        ),
        # "3rd party" EAs work
        ([decimal.Decimal(0), decimal.Decimal(1)], "decimal", to_decimal([0, 1])),
        # pass an ExtensionArray, but a different dtype
        (
            period_array(["2000", "2001"], freq="D"),
            "category",
            pd.Categorical([pd.Period("2000", "D"), pd.Period("2001", "D")]),
        ),
    ],
)
def test_array(data, dtype, expected):
    result = pd.array(data, dtype=dtype)
    tm.assert_equal(result, expected)


def test_array_copy():
    a = np.array([1, 2])
    # default is to copy
    b = pd.array(a, dtype=a.dtype)
    assert np.shares_memory(a, b._ndarray) is False

    # copy=True
    b = pd.array(a, dtype=a.dtype, copy=True)
    assert np.shares_memory(a, b._ndarray) is False

    # copy=False
    b = pd.array(a, dtype=a.dtype, copy=False)
    assert np.shares_memory(a, b._ndarray) is True


cet = pytz.timezone("CET")


@pytest.mark.parametrize(
    "data, expected",
    [
        # period
        (
            [pd.Period("2000", "D"), pd.Period("2001", "D")],
            period_array(["2000", "2001"], freq="D"),
        ),
        # interval
        (
            [pd.Interval(0, 1), pd.Interval(1, 2)],
            pd.arrays.IntervalArray.from_breaks([0, 1, 2]),
        ),
        # datetime
        (
            [pd.Timestamp("2000"), pd.Timestamp("2001")],
            pd.arrays.DatetimeArray._from_sequence(["2000", "2001"]),
        ),
        (
            [datetime.datetime(2000, 1, 1), datetime.datetime(2001, 1, 1)],
            pd.arrays.DatetimeArray._from_sequence(["2000", "2001"]),
        ),
        (
            np.array([1, 2], dtype="M8[ns]"),
            pd.arrays.DatetimeArray(np.array([1, 2], dtype="M8[ns]")),
        ),
        (
            np.array([1, 2], dtype="M8[us]"),
            pd.arrays.DatetimeArray(np.array([1000, 2000], dtype="M8[ns]")),
        ),
        # datetimetz
        (
            [pd.Timestamp("2000", tz="CET"), pd.Timestamp("2001", tz="CET")],
            pd.arrays.DatetimeArray._from_sequence(
                ["2000", "2001"], dtype=pd.DatetimeTZDtype(tz="CET")
            ),
        ),
        (
            [
                datetime.datetime(2000, 1, 1, tzinfo=cet),
                datetime.datetime(2001, 1, 1, tzinfo=cet),
            ],
            pd.arrays.DatetimeArray._from_sequence(["2000", "2001"], tz=cet),
        ),
        # timedelta
        (
            [pd.Timedelta("1H"), pd.Timedelta("2H")],
            pd.arrays.TimedeltaArray._from_sequence(["1H", "2H"]),
        ),
        (
            np.array([1, 2], dtype="m8[ns]"),
            pd.arrays.TimedeltaArray(np.array([1, 2], dtype="m8[ns]")),
        ),
        (
            np.array([1, 2], dtype="m8[us]"),
            pd.arrays.TimedeltaArray(np.array([1000, 2000], dtype="m8[ns]")),
        ),
        # integer
        ([1, 2], pd.arrays.IntegerArray._from_sequence([1, 2])),
        ([1, None], pd.arrays.IntegerArray._from_sequence([1, None])),
        # string
        (["a", "b"], pd.arrays.StringArray._from_sequence(["a", "b"])),
        (["a", None], pd.arrays.StringArray._from_sequence(["a", None])),
        # Boolean
        ([True, False], pd.arrays.BooleanArray._from_sequence([True, False])),
        ([True, None], pd.arrays.BooleanArray._from_sequence([True, None])),
    ],
)
def test_array_inference(data, expected):
    result = pd.array(data)
    tm.assert_equal(result, expected)


@pytest.mark.parametrize(
    "data",
    [
        # mix of frequencies
        [pd.Period("2000", "D"), pd.Period("2001", "A")],
        # mix of closed
        [pd.Interval(0, 1, closed="left"), pd.Interval(1, 2, closed="right")],
        # Mix of timezones
        [pd.Timestamp("2000", tz="CET"), pd.Timestamp("2000", tz="UTC")],
        # Mix of tz-aware and tz-naive
        [pd.Timestamp("2000", tz="CET"), pd.Timestamp("2000")],
        np.array([pd.Timestamp("2000"), pd.Timestamp("2000", tz="CET")]),
    ],
)
def test_array_inference_fails(data):
    result = pd.array(data)
    expected = PandasArray(np.array(data, dtype=object))
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize("data", [np.array([[1, 2], [3, 4]]), [[1, 2], [3, 4]]])
def test_nd_raises(data):
    with pytest.raises(ValueError, match="PandasArray must be 1-dimensional"):
        pd.array(data, dtype="int64")


def test_scalar_raises():
    with pytest.raises(ValueError, match="Cannot pass scalar '1'"):
        pd.array(1)


# ---------------------------------------------------------------------------
# A couple dummy classes to ensure that Series and Indexes are unboxed before
# getting to the EA classes.


@register_extension_dtype
class DecimalDtype2(DecimalDtype):
    name = "decimal2"

    @classmethod
    def construct_array_type(cls):
        """
        Return the array type associated with this dtype.

        Returns
        -------
        type
        """
        return DecimalArray2


class DecimalArray2(DecimalArray):
    @classmethod
    def _from_sequence(cls, scalars, dtype=None, copy=False):
        if isinstance(scalars, (pd.Series, pd.Index)):
            raise TypeError

        return super()._from_sequence(scalars, dtype=dtype, copy=copy)


def test_array_unboxes(index_or_series):
    box = index_or_series

    data = box([decimal.Decimal("1"), decimal.Decimal("2")])
    # make sure it works
    with pytest.raises(TypeError):
        DecimalArray2._from_sequence(data)

    result = pd.array(data, dtype="decimal2")
    expected = DecimalArray2._from_sequence(data.values)
    tm.assert_equal(result, expected)


@pytest.fixture
def registry_without_decimal():
    idx = registry.dtypes.index(DecimalDtype)
    registry.dtypes.pop(idx)
    yield
    registry.dtypes.append(DecimalDtype)


def test_array_not_registered(registry_without_decimal):
    # check we aren't on it
    assert registry.find("decimal") is None
    data = [decimal.Decimal("1"), decimal.Decimal("2")]

    result = pd.array(data, dtype=DecimalDtype)
    expected = DecimalArray._from_sequence(data)
    tm.assert_equal(result, expected)


class TestArrayAnalytics:
    def test_searchsorted(self, string_dtype):
        arr = pd.array(["a", "b", "c"], dtype=string_dtype)

        result = arr.searchsorted("a", side="left")
        assert is_scalar(result)
        assert result == 0

        result = arr.searchsorted("a", side="right")
        assert is_scalar(result)
        assert result == 1

    def test_searchsorted_numeric_dtypes_scalar(self, any_real_dtype):
        arr = pd.array([1, 3, 90], dtype=any_real_dtype)
        result = arr.searchsorted(30)
        assert is_scalar(result)
        assert result == 2

        result = arr.searchsorted([30])
        expected = np.array([2], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_searchsorted_numeric_dtypes_vector(self, any_real_dtype):
        arr = pd.array([1, 3, 90], dtype=any_real_dtype)
        result = arr.searchsorted([2, 30])
        expected = np.array([1, 2], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize(
        "arr, val",
        [
            [
                pd.date_range("20120101", periods=10, freq="2D"),
                pd.Timestamp("20120102"),
            ],
            [
                pd.date_range("20120101", periods=10, freq="2D", tz="Asia/Hong_Kong"),
                pd.Timestamp("20120102", tz="Asia/Hong_Kong"),
            ],
            [
                pd.timedelta_range(start="1 day", end="10 days", periods=10),
                pd.Timedelta("2 days"),
            ],
        ],
    )
    def test_search_sorted_datetime64_scalar(self, arr, val):
        arr = pd.array(arr)
        result = arr.searchsorted(val)
        assert is_scalar(result)
        assert result == 1

    def test_searchsorted_sorter(self, any_real_dtype):
        arr = pd.array([3, 1, 2], dtype=any_real_dtype)
        result = arr.searchsorted([0, 3], sorter=np.argsort(arr))
        expected = np.array([0, 2], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)
