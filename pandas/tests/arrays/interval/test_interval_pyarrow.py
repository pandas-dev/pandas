import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays import IntervalArray


def test_arrow_extension_type():
    pa = pytest.importorskip("pyarrow")

    from pandas.core.arrays.arrow.extension_types import ArrowIntervalType

    p1 = ArrowIntervalType(pa.int64(), "left")
    p2 = ArrowIntervalType(pa.int64(), "left")
    p3 = ArrowIntervalType(pa.int64(), "right")

    assert p1.closed == "left"
    assert p1 == p2
    assert p1 != p3
    assert hash(p1) == hash(p2)
    assert hash(p1) != hash(p3)


def test_arrow_array():
    pa = pytest.importorskip("pyarrow")

    from pandas.core.arrays.arrow.extension_types import ArrowIntervalType

    intervals = pd.interval_range(1, 5, freq=1).array

    result = pa.array(intervals)
    assert isinstance(result.type, ArrowIntervalType)
    assert result.type.closed == intervals.closed
    assert result.type.subtype == pa.int64()
    assert result.storage.field("left").equals(pa.array([1, 2, 3, 4], type="int64"))
    assert result.storage.field("right").equals(pa.array([2, 3, 4, 5], type="int64"))

    expected = pa.array([{"left": i, "right": i + 1} for i in range(1, 5)])
    assert result.storage.equals(expected)

    # convert to its storage type
    result = pa.array(intervals, type=expected.type)
    assert result.equals(expected)

    # unsupported conversions
    with pytest.raises(TypeError, match="Not supported to convert IntervalArray"):
        pa.array(intervals, type="float64")

    with pytest.raises(TypeError, match="Not supported to convert IntervalArray"):
        pa.array(intervals, type=ArrowIntervalType(pa.float64(), "left"))


def test_arrow_array_missing(using_nan_is_na):
    pa = pytest.importorskip("pyarrow")

    from pandas.core.arrays.arrow.extension_types import ArrowIntervalType

    arr = IntervalArray.from_breaks([0.0, 1.0, 2.0, 3.0])
    arr[1] = None

    result = pa.array(arr)
    assert isinstance(result.type, ArrowIntervalType)
    assert result.type.closed == arr.closed
    assert result.type.subtype == pa.float64()

    # fields have missing values (not NaN)
    left = pa.array([0.0, None, 2.0], type="float64")
    right = pa.array([1.0, None, 3.0], type="float64")
    assert result.storage.field("left").equals(left)
    assert result.storage.field("right").equals(right)

    # structarray itself also has missing values on the array level
    vals = [
        {"left": 0.0, "right": 1.0},
        {"left": None, "right": None},
        {"left": 2.0, "right": 3.0},
    ]
    expected = pa.StructArray.from_pandas(vals, mask=np.array([False, True, False]))
    assert result.storage.equals(expected)


@pytest.mark.filterwarnings(
    "ignore:Passing a BlockManager to DataFrame:DeprecationWarning"
)
@pytest.mark.parametrize(
    "breaks",
    [[0.0, 1.0, 2.0, 3.0], pd.date_range("2017", periods=4, freq="D")],
    ids=["float", "datetime64[ns]"],
)
def test_arrow_table_roundtrip(breaks):
    pa = pytest.importorskip("pyarrow")

    from pandas.core.arrays.arrow.extension_types import ArrowIntervalType

    arr = IntervalArray.from_breaks(breaks)
    arr[1] = None
    df = pd.DataFrame({"a": arr})

    table = pa.table(df)
    assert isinstance(table.field("a").type, ArrowIntervalType)
    result = table.to_pandas()
    assert isinstance(result["a"].dtype, pd.IntervalDtype)
    tm.assert_frame_equal(result, df)

    table2 = pa.concat_tables([table, table])
    result = table2.to_pandas()
    expected = pd.concat([df, df], ignore_index=True)
    tm.assert_frame_equal(result, expected)

    # GH#41040
    table = pa.table(
        [pa.chunked_array([], type=table.column(0).type)], schema=table.schema
    )
    result = table.to_pandas()
    tm.assert_frame_equal(result, expected[0:0])


@pytest.mark.filterwarnings(
    "ignore:Passing a BlockManager to DataFrame:DeprecationWarning"
)
@pytest.mark.parametrize(
    "breaks",
    [[0.0, 1.0, 2.0, 3.0], pd.date_range("2017", periods=4, freq="D")],
    ids=["float", "datetime64[ns]"],
)
def test_arrow_table_roundtrip_without_metadata(breaks):
    pa = pytest.importorskip("pyarrow")

    arr = IntervalArray.from_breaks(breaks)
    arr[1] = None
    df = pd.DataFrame({"a": arr})

    table = pa.table(df)
    # remove the metadata
    table = table.replace_schema_metadata()
    assert table.schema.metadata is None

    result = table.to_pandas()
    assert isinstance(result["a"].dtype, pd.IntervalDtype)
    tm.assert_frame_equal(result, df)


def test_from_arrow_from_raw_struct_array():
    # in case pyarrow lost the Interval extension type (eg on parquet roundtrip
    # with datetime64[ns] subtype, see GH-45881), still allow conversion
    # from arrow to IntervalArray
    pa = pytest.importorskip("pyarrow")

    arr = pa.array([{"left": 0, "right": 1}, {"left": 1, "right": 2}])
    dtype = pd.IntervalDtype(np.dtype("int64"), closed="neither")

    result = dtype.__from_arrow__(arr)
    expected = IntervalArray.from_breaks(
        np.array([0, 1, 2], dtype="int64"), closed="neither"
    )
    tm.assert_extension_array_equal(result, expected)

    result = dtype.__from_arrow__(pa.chunked_array([arr]))
    tm.assert_extension_array_equal(result, expected)


def test_from_arrow_datetime_tz():
    # GH#64297 - IntervalDtype with DatetimeTZDtype subtype must convert
    # from Arrow, not raise TypeError on np.asarray(..., dtype=DatetimeTZDtype)
    pa = pytest.importorskip("pyarrow")

    from pandas.core.arrays.arrow.extension_types import ArrowIntervalType

    tz = "Europe/Brussels"
    left = pd.Timestamp("2012-01-01", tz=tz)
    right = pd.Timestamp("2013-01-01", tz=tz)
    arrow_type = ArrowIntervalType(pa.timestamp("us", tz=tz), "right")
    arr = pa.array([{"left": left, "right": right}], type=arrow_type)

    dtype = pd.IntervalDtype(pd.DatetimeTZDtype("us", tz), closed="right")
    result = dtype.__from_arrow__(arr)

    expected = IntervalArray.from_arrays(
        pd.DatetimeIndex([left]),
        pd.DatetimeIndex([right]),
        closed="right",
    )
    tm.assert_extension_array_equal(result, expected)

    # chunked and empty
    result = dtype.__from_arrow__(pa.chunked_array([arr]))
    tm.assert_extension_array_equal(result, expected)

    empty = dtype.__from_arrow__(pa.chunked_array([], type=arrow_type))
    assert len(empty) == 0
    assert empty.dtype == dtype


def test_arrow_array_datetime_tz_roundtrip():
    # GH#64297 - full Arrow round-trip for tz-aware interval bounds
    pa = pytest.importorskip("pyarrow")

    from pandas.core.arrays.arrow.extension_types import ArrowIntervalType

    breaks = pd.date_range("2017-01-01", periods=4, freq="D", tz="US/Eastern")
    arr = IntervalArray.from_breaks(breaks)
    arr[1] = None

    result = pa.array(arr)
    assert isinstance(result.type, ArrowIntervalType)
    assert result.type.closed == arr.closed
    assert result.type.subtype == pa.timestamp("us", tz="US/Eastern")

    back = arr.dtype.__from_arrow__(result)
    tm.assert_extension_array_equal(back, arr)


def test_from_arrow_nullable_integer_subtype_falls_back_to_numpy():
    # GH#64297 - IntervalArray does not generally retain nullable integer
    # subtypes, but conversion should no longer raise TypeError.
    pa = pytest.importorskip("pyarrow")

    from pandas.core.arrays.arrow.extension_types import ArrowIntervalType

    arr = pa.array(
        [{"left": 1, "right": 2}], type=ArrowIntervalType(pa.int64(), "right")
    )
    dtype = pd.IntervalDtype(pd.Int64Dtype(), closed="right")
    result = dtype.__from_arrow__(arr)
    expected = IntervalArray.from_arrays([1], [2], closed="right")
    tm.assert_extension_array_equal(result, expected)
