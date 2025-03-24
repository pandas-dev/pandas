# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.

import datetime
import decimal
import pytest
import weakref

try:
    import numpy as np
except ImportError:
    np = None

import pyarrow as pa
import pyarrow.compute as pc


@pytest.mark.parametrize(['value', 'ty', 'klass'], [
    (False, None, pa.BooleanScalar),
    (True, None, pa.BooleanScalar),
    (1, None, pa.Int64Scalar),
    (-1, None, pa.Int64Scalar),
    (1, pa.int8(), pa.Int8Scalar),
    (1, pa.uint8(), pa.UInt8Scalar),
    (1, pa.int16(), pa.Int16Scalar),
    (1, pa.uint16(), pa.UInt16Scalar),
    (1, pa.int32(), pa.Int32Scalar),
    (1, pa.uint32(), pa.UInt32Scalar),
    (1, pa.int64(), pa.Int64Scalar),
    (1, pa.uint64(), pa.UInt64Scalar),
    (1.0, None, pa.DoubleScalar),
    (1.0, pa.float32(), pa.FloatScalar),
    (decimal.Decimal("1.123"), None, pa.Decimal128Scalar),
    (decimal.Decimal("1.1234567890123456789012345678901234567890"),
     None, pa.Decimal256Scalar),
    ("string", None, pa.StringScalar),
    (b"bytes", None, pa.BinaryScalar),
    ("largestring", pa.large_string(), pa.LargeStringScalar),
    (b"largebytes", pa.large_binary(), pa.LargeBinaryScalar),
    ("string_view", pa.string_view(), pa.StringViewScalar),
    (b"bytes_view", pa.binary_view(), pa.BinaryViewScalar),
    (b"abc", pa.binary(3), pa.FixedSizeBinaryScalar),
    ([1, 2, 3], None, pa.ListScalar),
    ([1, 2, 3, 4], pa.large_list(pa.int8()), pa.LargeListScalar),
    ([1, 2, 3, 4, 5], pa.list_(pa.int8(), 5), pa.FixedSizeListScalar),
    ([1, 2, 3], pa.list_view(pa.int8()), pa.ListViewScalar),
    ([1, 2, 3, 4], pa.large_list_view(pa.int8()), pa.LargeListViewScalar),
    (datetime.date.today(), None, pa.Date32Scalar),
    (datetime.date.today(), pa.date64(), pa.Date64Scalar),
    (datetime.datetime.now(), None, pa.TimestampScalar),
    (datetime.datetime.now().time().replace(microsecond=0), pa.time32('s'),
     pa.Time32Scalar),
    (datetime.datetime.now().time(), None, pa.Time64Scalar),
    (datetime.timedelta(days=1), None, pa.DurationScalar),
    (pa.MonthDayNano([1, -1, -10100]), None,
     pa.MonthDayNanoIntervalScalar),
    ({'a': 1, 'b': [1, 2]}, None, pa.StructScalar),
    ([('a', 1), ('b', 2)], pa.map_(pa.string(), pa.int8()), pa.MapScalar),
])
def test_basics(value, ty, klass, pickle_module):
    s = pa.scalar(value, type=ty)
    s.validate()
    s.validate(full=True)
    assert isinstance(s, klass)
    assert s.as_py() == value
    assert s == pa.scalar(value, type=ty)
    assert s != value
    assert s != "else"
    assert hash(s) == hash(s)
    assert s.is_valid is True
    assert s != None  # noqa: E711

    s = pa.scalar(None, type=s.type)
    assert s.is_valid is False
    assert s.as_py() is None
    assert s != pa.scalar(value, type=ty)

    # test pickle roundtrip
    restored = pickle_module.loads(pickle_module.dumps(s))
    assert s.equals(restored)

    # test that scalars are weak-referenceable
    wr = weakref.ref(s)
    assert wr() is not None
    del s
    assert wr() is None


# This test is a copy of test_basics but only for float16 (HalfFloatScalar)
# which currently requires a numpy scalar to create it. The test collection
# fails if numpy is used on the parametrization when not present.
@pytest.mark.numpy
def test_basics_np_required(pickle_module):
    value, ty, klass = np.float16(1.0), pa.float16(), pa.HalfFloatScalar
    s = pa.scalar(value, type=ty)
    s.validate()
    s.validate(full=True)
    assert isinstance(s, klass)
    assert s.as_py() == value
    assert s == pa.scalar(value, type=ty)
    assert s != value
    assert s != "else"
    assert hash(s) == hash(s)
    assert s.is_valid is True
    assert s != None  # noqa: E711

    s = pa.scalar(None, type=s.type)
    assert s.is_valid is False
    assert s.as_py() is None
    assert s != pa.scalar(value, type=ty)

    # test pickle roundtrip
    restored = pickle_module.loads(pickle_module.dumps(s))
    assert s.equals(restored)

    # test that scalars are weak-referenceable
    wr = weakref.ref(s)
    assert wr() is not None
    del s
    assert wr() is None


def test_invalid_scalar():
    s = pc.cast(pa.scalar(b"\xff"), pa.string(), safe=False)
    s.validate()
    with pytest.raises(ValueError,
                       match="string scalar contains invalid UTF8 data"):
        s.validate(full=True)


def test_null_singleton():
    with pytest.raises(RuntimeError):
        pa.NullScalar()


def test_nulls(pickle_module):
    null = pa.scalar(None)
    assert null is pa.NA
    assert null.as_py() is None
    assert null != "something"
    assert (null == pa.scalar(None)) is True
    assert (null == 0) is False
    assert pa.NA == pa.NA
    assert pa.NA not in [5]

    arr = pa.array([None, None])
    for v in arr:
        assert v is pa.NA
        assert v.as_py() is None

    # test pickle roundtrip
    restored = pickle_module.loads(pickle_module.dumps(null))
    assert restored.equals(null)

    # test that scalars are weak-referenceable
    wr = weakref.ref(null)
    assert wr() is not None
    del null
    assert wr() is not None  # singleton


def test_hashing():
    # ARROW-640
    values = list(range(500))
    arr = pa.array(values + values)
    set_from_array = set(arr)
    assert isinstance(set_from_array, set)
    assert len(set_from_array) == 500


def test_hashing_struct_scalar():
    # GH-35360
    a = pa.array([[{'a': 5}, {'a': 6}], [{'a': 7}, None]])
    b = pa.array([[{'a': 7}, None]])
    hash1 = hash(a[1])
    hash2 = hash(b[0])
    assert hash1 == hash2


@pytest.mark.timezone_data
def test_timestamp_scalar():
    a = repr(pa.scalar("0000-01-01").cast(pa.timestamp("s")))
    assert a == "<pyarrow.TimestampScalar: '0000-01-01T00:00:00'>"
    b = repr(pa.scalar(datetime.datetime(2015, 1, 1), type=pa.timestamp('s', tz='UTC')))
    assert b == "<pyarrow.TimestampScalar: '2015-01-01T00:00:00+0000'>"
    c = repr(pa.scalar(datetime.datetime(2015, 1, 1), type=pa.timestamp('us')))
    assert c == "<pyarrow.TimestampScalar: '2015-01-01T00:00:00.000000'>"
    d = repr(pc.assume_timezone(
        pa.scalar("2000-01-01").cast(pa.timestamp("s")), "America/New_York"))
    assert d == "<pyarrow.TimestampScalar: '2000-01-01T00:00:00-0500'>"


def test_bool():
    false = pa.scalar(False)
    true = pa.scalar(True)

    assert isinstance(false, pa.BooleanScalar)
    assert isinstance(true, pa.BooleanScalar)

    assert repr(true) == "<pyarrow.BooleanScalar: True>"
    assert str(true) == "True"
    assert repr(false) == "<pyarrow.BooleanScalar: False>"
    assert str(false) == "False"

    assert true.as_py() is True
    assert false.as_py() is False


def test_numerics():
    # int64
    s = pa.scalar(1)
    assert isinstance(s, pa.Int64Scalar)
    assert repr(s) == "<pyarrow.Int64Scalar: 1>"
    assert str(s) == "1"
    assert s.as_py() == 1

    with pytest.raises(OverflowError):
        pa.scalar(-1, type='uint8')

    # float64
    s = pa.scalar(1.5)
    assert isinstance(s, pa.DoubleScalar)
    assert repr(s) == "<pyarrow.DoubleScalar: 1.5>"
    assert str(s) == "1.5"
    assert s.as_py() == 1.5

    if np is not None:
        # float16
        s = pa.scalar(np.float16(0.5), type='float16')
        assert isinstance(s, pa.HalfFloatScalar)
        # on numpy2 repr(np.float16(0.5)) == "np.float16(0.5)"
        # on numpy1 repr(np.float16(0.5)) == "0.5"
        assert repr(s) == f"<pyarrow.HalfFloatScalar: {np.float16(0.5)!r}>"
        assert str(s) == "0.5"
        assert s.as_py() == 0.5


def test_decimal128():
    v = decimal.Decimal("1.123")
    s = pa.scalar(v)
    assert isinstance(s, pa.Decimal128Scalar)
    assert s.as_py() == v
    assert s.type == pa.decimal128(4, 3)

    v = decimal.Decimal("1.1234")
    with pytest.raises(pa.ArrowInvalid):
        pa.scalar(v, type=pa.decimal128(4, scale=3))
    with pytest.raises(pa.ArrowInvalid):
        pa.scalar(v, type=pa.decimal128(5, scale=3))

    s = pa.scalar(v, type=pa.decimal128(5, scale=4))
    assert isinstance(s, pa.Decimal128Scalar)
    assert s.as_py() == v


def test_decimal256():
    v = decimal.Decimal("1234567890123456789012345678901234567890.123")
    s = pa.scalar(v)
    assert isinstance(s, pa.Decimal256Scalar)
    assert s.as_py() == v
    assert s.type == pa.decimal256(43, 3)

    v = decimal.Decimal("1.1234")
    with pytest.raises(pa.ArrowInvalid):
        pa.scalar(v, type=pa.decimal256(4, scale=3))
    with pytest.raises(pa.ArrowInvalid):
        pa.scalar(v, type=pa.decimal256(5, scale=3))

    s = pa.scalar(v, type=pa.decimal256(5, scale=4))
    assert isinstance(s, pa.Decimal256Scalar)
    assert s.as_py() == v


def test_date():
    # ARROW-5125
    d1 = datetime.date(3200, 1, 1)
    d2 = datetime.date(1960, 1, 1)

    for ty in [pa.date32(), pa.date64()]:
        for d in [d1, d2]:
            s = pa.scalar(d, type=ty)
            assert s.as_py() == d


def test_date_cast():
    # ARROW-10472 - casting fo scalars doesn't segfault
    scalar = pa.scalar(datetime.datetime(2012, 1, 1), type=pa.timestamp("us"))
    expected = datetime.date(2012, 1, 1)
    for ty in [pa.date32(), pa.date64()]:
        result = scalar.cast(ty)
        assert result.as_py() == expected


def test_time_from_datetime_time():
    t1 = datetime.time(18, 0)
    t2 = datetime.time(21, 0)

    types = [pa.time32('s'), pa.time32('ms'), pa.time64('us'), pa.time64('ns')]
    for ty in types:
        for t in [t1, t2]:
            s = pa.scalar(t, type=ty)
            assert s.as_py() == t


@pytest.mark.parametrize(['value', 'time_type'], [
    (1, pa.time32("s")),
    (2**30, pa.time32("s")),
    (None, pa.time32("s")),
    (1, pa.time32("ms")),
    (2**30, pa.time32("ms")),
    (None, pa.time32("ms")),
    (1, pa.time64("us")),
    (2**62, pa.time64("us")),
    (None, pa.time64("us")),
    (1, pa.time64("ns")),
    (2**62, pa.time64("ns")),
    (None, pa.time64("ns")),
    (1, pa.date32()),
    (2**30, pa.date32()),
    (None, pa.date32()),
    (1, pa.date64()),
    (2**62, pa.date64()),
    (None, pa.date64()),
    (1, pa.timestamp("ns")),
    (2**62, pa.timestamp("ns")),
    (None, pa.timestamp("ns")),
    (1, pa.duration("ns")),
    (2**62, pa.duration("ns")),
    (None, pa.duration("ns")),
    ((1, 2, -3), pa.month_day_nano_interval()),
    (None, pa.month_day_nano_interval()),
])
def test_temporal_values(value, time_type: pa.DataType):
    time_scalar = pa.scalar(value, type=time_type)
    time_scalar.validate(full=True)
    assert time_scalar.value == value


def test_cast():
    val = pa.scalar(5, type='int8')
    assert val.cast('int64') == pa.scalar(5, type='int64')
    assert val.cast('uint32') == pa.scalar(5, type='uint32')
    assert val.cast('string') == pa.scalar('5', type='string')
    with pytest.raises(ValueError):
        pa.scalar('foo').cast('int32')


@pytest.mark.timezone_data
def test_cast_timestamp_to_string():
    # GH-35370
    pytest.importorskip("pytz")
    import pytz
    dt = datetime.datetime(2000, 1, 1, 0, 0, 0, tzinfo=pytz.utc)
    ts = pa.scalar(dt, type=pa.timestamp("ns", tz="UTC"))
    assert ts.cast(pa.string()) == pa.scalar('2000-01-01 00:00:00.000000000Z')


def test_cast_float_to_int():
    # GH-35040
    float_scalar = pa.scalar(1.5, type=pa.float64())
    unsafe_cast = float_scalar.cast(pa.int64(), safe=False)
    expected_unsafe_cast = pa.scalar(1, type=pa.int64())
    assert unsafe_cast == expected_unsafe_cast
    with pytest.raises(pa.ArrowInvalid):
        float_scalar.cast(pa.int64())  # verify default is safe cast


def test_cast_int_to_float():
    # GH-34901
    int_scalar = pa.scalar(18014398509481983, type=pa.int64())
    unsafe_cast = int_scalar.cast(pa.float64(), safe=False)
    expected_unsafe_cast = pa.scalar(18014398509481983.0, type=pa.float64())
    assert unsafe_cast == expected_unsafe_cast
    with pytest.raises(pa.ArrowInvalid):
        int_scalar.cast(pa.float64())  # verify default is safe cast


@pytest.mark.parametrize("typ", [pa.date32(), pa.date64()])
def test_cast_string_to_date(typ):
    scalar = pa.scalar('2021-01-01')
    result = scalar.cast(typ)
    assert result == pa.scalar(datetime.date(2021, 1, 1), type=typ)


@pytest.mark.pandas
def test_timestamp():
    import pandas as pd
    arr = pd.date_range('2000-01-01 12:34:56', periods=10).values

    units = ['ns', 'us', 'ms', 's']

    for i, unit in enumerate(units):
        dtype = 'datetime64[{}]'.format(unit)
        arrow_arr = pa.Array.from_pandas(arr.astype(dtype))
        expected = pd.Timestamp('2000-01-01 12:34:56')

        assert arrow_arr[0].as_py() == expected
        assert arrow_arr[0].value * 1000**i == expected.value

        tz = 'America/New_York'
        arrow_type = pa.timestamp(unit, tz=tz)

        dtype = 'datetime64[{}]'.format(unit)
        arrow_arr = pa.Array.from_pandas(arr.astype(dtype), type=arrow_type)
        expected = (pd.Timestamp('2000-01-01 12:34:56')
                    .tz_localize('utc')
                    .tz_convert(tz))

        assert arrow_arr[0].as_py() == expected
        assert arrow_arr[0].value * 1000**i == expected.value


@pytest.mark.nopandas
def test_timestamp_nanos_nopandas():
    # ARROW-5450
    pytest.importorskip("pytz")
    import pytz
    tz = 'America/New_York'
    ty = pa.timestamp('ns', tz=tz)

    # 2000-01-01 00:00:00 + 1 microsecond
    s = pa.scalar(946684800000000000 + 1000, type=ty)

    tzinfo = pytz.timezone(tz)
    expected = datetime.datetime(2000, 1, 1, microsecond=1, tzinfo=tzinfo)
    expected = tzinfo.fromutc(expected)
    result = s.as_py()
    assert result == expected
    assert result.year == 1999
    assert result.hour == 19

    # Non-zero nanos yields ValueError
    s = pa.scalar(946684800000000001, type=ty)
    with pytest.raises(ValueError):
        s.as_py()


def test_timestamp_no_overflow():
    # ARROW-5450
    pytest.importorskip("pytz")
    import pytz

    timestamps = [
        datetime.datetime(1, 1, 1, 0, 0, 0, tzinfo=pytz.utc),
        datetime.datetime(9999, 12, 31, 23, 59, 59, 999999, tzinfo=pytz.utc),
        datetime.datetime(1970, 1, 1, 0, 0, 0, tzinfo=pytz.utc),
    ]
    for ts in timestamps:
        s = pa.scalar(ts, type=pa.timestamp("us", tz="UTC"))
        assert s.as_py() == ts


def test_timestamp_fixed_offset_print():
    # ARROW-13896
    pytest.importorskip("pytz")
    arr = pa.array([0], pa.timestamp('s', tz='+02:00'))
    assert str(arr[0]) == "1970-01-01 02:00:00+02:00"


@pytest.mark.numpy
def test_duration():
    arr = np.array([0, 3600000000000], dtype='timedelta64[ns]')

    units = ['us', 'ms', 's']

    for i, unit in enumerate(units):
        dtype = 'timedelta64[{}]'.format(unit)
        arrow_arr = pa.array(arr.astype(dtype))
        expected = datetime.timedelta(seconds=60*60)
        assert isinstance(arrow_arr[1].as_py(), datetime.timedelta)
        assert arrow_arr[1].as_py() == expected
        assert (arrow_arr[1].value * 1000**(i+1) ==
                expected.total_seconds() * 1e9)


@pytest.mark.pandas
def test_duration_nanos_pandas():
    import pandas as pd
    arr = pa.array([0, 3600000000000], type=pa.duration('ns'))
    expected = pd.Timedelta('1 hour')
    assert isinstance(arr[1].as_py(), pd.Timedelta)
    assert arr[1].as_py() == expected
    assert arr[1].value == expected.value

    # Non-zero nanos work fine
    arr = pa.array([946684800000000001], type=pa.duration('ns'))
    assert arr[0].as_py() == pd.Timedelta(946684800000000001, unit='ns')


@pytest.mark.nopandas
def test_duration_nanos_nopandas():
    arr = pa.array([0, 3600000000000], pa.duration('ns'))
    expected = datetime.timedelta(seconds=60*60)
    assert isinstance(arr[1].as_py(), datetime.timedelta)
    assert arr[1].as_py() == expected
    assert arr[1].value == expected.total_seconds() * 1e9

    # Non-zero nanos yields ValueError
    arr = pa.array([946684800000000001], type=pa.duration('ns'))
    with pytest.raises(ValueError):
        arr[0].as_py()


def test_month_day_nano_interval():
    triple = pa.MonthDayNano([-3600, 1800, -50])
    arr = pa.array([triple])
    assert isinstance(arr[0].as_py(), pa.MonthDayNano)
    assert arr[0].as_py() == triple
    assert arr[0].value == triple


@pytest.mark.parametrize('value', ['foo', 'mañana'])
@pytest.mark.parametrize(('ty', 'scalar_typ'), [
    (pa.string(), pa.StringScalar),
    (pa.large_string(), pa.LargeStringScalar),
    (pa.string_view(), pa.StringViewScalar),
])
def test_string(value, ty, scalar_typ):
    s = pa.scalar(value, type=ty)
    assert isinstance(s, scalar_typ)
    assert s.as_py() == value
    assert s.as_py() != 'something'
    assert repr(value) in repr(s)
    assert str(s) == str(value)

    buf = s.as_buffer()
    assert isinstance(buf, pa.Buffer)
    assert buf.to_pybytes() == value.encode()


@pytest.mark.parametrize('value', [b'foo', b'bar'])
@pytest.mark.parametrize(('ty', 'scalar_typ'), [
    (pa.binary(), pa.BinaryScalar),
    (pa.large_binary(), pa.LargeBinaryScalar),
    (pa.binary_view(), pa.BinaryViewScalar),
])
def test_binary(value, ty, scalar_typ):
    s = pa.scalar(value, type=ty)
    assert isinstance(s, scalar_typ)
    assert s.as_py() == value
    assert str(s) == str(value)
    assert repr(value) in repr(s)
    assert s.as_py() == value
    assert s != b'xxxxx'

    buf = s.as_buffer()
    assert isinstance(buf, pa.Buffer)
    assert buf.to_pybytes() == value


def test_fixed_size_binary():
    s = pa.scalar(b'foof', type=pa.binary(4))
    assert isinstance(s, pa.FixedSizeBinaryScalar)
    assert s.as_py() == b'foof'

    with pytest.raises(pa.ArrowInvalid):
        pa.scalar(b'foof5', type=pa.binary(4))


@pytest.mark.parametrize(('ty', 'klass'), [
    (pa.list_(pa.string()), pa.ListScalar),
    (pa.large_list(pa.string()), pa.LargeListScalar),
    (pa.list_view(pa.string()), pa.ListViewScalar),
    (pa.large_list_view(pa.string()), pa.LargeListViewScalar)
])
def test_list(ty, klass):
    v = ['foo', None]
    s = pa.scalar(v, type=ty)
    assert s.type == ty
    assert len(s) == 2
    assert isinstance(s.values, pa.Array)
    assert s.values.to_pylist() == v
    assert isinstance(s, klass)
    assert repr(v) in repr(s)
    assert s.as_py() == v
    assert s[0].as_py() == 'foo'
    assert s[1].as_py() is None
    assert s[-1] == s[1]
    assert s[-2] == s[0]
    with pytest.raises(IndexError):
        s[-3]
    with pytest.raises(IndexError):
        s[2]


@pytest.mark.numpy
@pytest.mark.parametrize('ty', [
    pa.list_(pa.int64()),
    pa.large_list(pa.int64()),
    pa.list_view(pa.int64()),
    pa.large_list_view(pa.int64()),
    None
])
def test_list_from_numpy(ty):
    s = pa.scalar(np.array([1, 2, 3], dtype=np.int64()), type=ty)
    if ty is None:
        ty = pa.list_(pa.int64())  # expected inferred type
    assert s.type == ty
    assert s.as_py() == [1, 2, 3]


@pytest.mark.pandas
@pytest.mark.parametrize('factory', [
    pa.list_,
    pa.large_list,
    pa.list_view,
    pa.large_list_view
])
def test_list_from_pandas(factory):
    import pandas as pd

    s = pa.scalar(pd.Series([1, 2, 3]))
    assert s.as_py() == [1, 2, 3]

    cases = [
        (np.nan, 'null'),
        (['string', np.nan], factory(pa.binary())),
        (['string', np.nan], factory(pa.utf8())),
        ([b'string', np.nan], factory(pa.binary(6))),
        ([True, np.nan], factory(pa.bool_())),
        ([decimal.Decimal('0'), np.nan], factory(pa.decimal128(12, 2))),
    ]
    for case, ty in cases:
        # Both types of exceptions are raised. May want to clean that up
        with pytest.raises((ValueError, TypeError)):
            pa.scalar(case, type=ty)

        # from_pandas option suppresses failure
        s = pa.scalar(case, type=ty, from_pandas=True)


def test_fixed_size_list():
    s = pa.scalar([1, None, 3], type=pa.list_(pa.int64(), 3))

    assert len(s) == 3
    assert isinstance(s, pa.FixedSizeListScalar)
    assert repr(s) == "<pyarrow.FixedSizeListScalar: [1, None, 3]>"
    assert s.as_py() == [1, None, 3]
    assert s[0].as_py() == 1
    assert s[1].as_py() is None
    assert s[-1] == s[2]
    with pytest.raises(IndexError):
        s[-4]
    with pytest.raises(IndexError):
        s[3]


def test_struct():
    ty = pa.struct([
        pa.field('x', pa.int16()),
        pa.field('y', pa.float32())
    ])

    v = {'x': 2, 'y': 3.5}
    s = pa.scalar(v, type=ty)
    assert list(s) == list(s.keys()) == ['x', 'y']
    assert list(s.values()) == [
        pa.scalar(2, type=pa.int16()),
        pa.scalar(3.5, type=pa.float32())
    ]
    assert list(s.items()) == [
        ('x', pa.scalar(2, type=pa.int16())),
        ('y', pa.scalar(3.5, type=pa.float32()))
    ]
    assert 'x' in s
    assert 'y' in s
    assert 'z' not in s
    assert 0 not in s

    assert s.as_py() == v
    assert repr(s) != repr(v)
    assert repr(s.as_py()) == repr(v)
    assert len(s) == 2
    assert isinstance(s['x'], pa.Int16Scalar)
    assert isinstance(s['y'], pa.FloatScalar)
    assert s['x'].as_py() == 2
    assert s['y'].as_py() == 3.5

    with pytest.raises(KeyError):
        s['nonexistent']

    s = pa.scalar(None, type=ty)
    assert list(s) == list(s.keys()) == ['x', 'y']
    assert s.as_py() is None
    assert 'x' in s
    assert 'y' in s
    assert isinstance(s['x'], pa.Int16Scalar)
    assert isinstance(s['y'], pa.FloatScalar)
    assert s['x'].is_valid is False
    assert s['y'].is_valid is False
    assert s['x'].as_py() is None
    assert s['y'].as_py() is None


def test_struct_duplicate_fields():
    ty = pa.struct([
        pa.field('x', pa.int16()),
        pa.field('y', pa.float32()),
        pa.field('x', pa.int64()),
    ])
    s = pa.scalar([('x', 1), ('y', 2.0), ('x', 3)], type=ty)

    assert list(s) == list(s.keys()) == ['x', 'y', 'x']
    assert len(s) == 3
    assert s == s
    assert list(s.items()) == [
        ('x', pa.scalar(1, pa.int16())),
        ('y', pa.scalar(2.0, pa.float32())),
        ('x', pa.scalar(3, pa.int64()))
    ]

    assert 'x' in s
    assert 'y' in s
    assert 'z' not in s
    assert 0 not in s

    # getitem with field names fails for duplicate fields, works for others
    with pytest.raises(KeyError):
        s['x']

    assert isinstance(s['y'], pa.FloatScalar)
    assert s['y'].as_py() == 2.0

    # getitem with integer index works for all fields
    assert isinstance(s[0], pa.Int16Scalar)
    assert s[0].as_py() == 1
    assert isinstance(s[1], pa.FloatScalar)
    assert s[1].as_py() == 2.0
    assert isinstance(s[2], pa.Int64Scalar)
    assert s[2].as_py() == 3

    assert "pyarrow.StructScalar" in repr(s)

    with pytest.raises(ValueError, match="duplicate field names"):
        s.as_py()


def test_map(pickle_module):
    ty = pa.map_(pa.string(), pa.int8())
    v = [('a', 1), ('b', 2)]
    s = pa.scalar(v, type=ty)

    assert len(s) == 2
    assert isinstance(s, pa.MapScalar)
    assert isinstance(s.values, pa.Array)
    assert repr(s) == "<pyarrow.MapScalar: [('a', 1), ('b', 2)]>"
    assert s.values.to_pylist() == [
        {'key': 'a', 'value': 1},
        {'key': 'b', 'value': 2}
    ]

    # test iteration
    for i, j in zip(s, v):
        assert i == j

    # test iteration with missing values
    for _ in pa.scalar(None, type=ty):
        pass

    assert s.as_py() == v
    assert s[1] == (
        pa.scalar('b', type=pa.string()),
        pa.scalar(2, type=pa.int8())
    )
    assert s[-1] == s[1]
    assert s[-2] == s[0]
    with pytest.raises(IndexError):
        s[-3]
    with pytest.raises(IndexError):
        s[2]

    restored = pickle_module.loads(pickle_module.dumps(s))
    assert restored.equals(s)


def test_dictionary(pickle_module):
    indices = pa.array([2, None, 1, 2, 0, None])
    dictionary = pa.array(['foo', 'bar', 'baz'])

    arr = pa.DictionaryArray.from_arrays(indices, dictionary)
    expected = ['baz', None, 'bar', 'baz', 'foo', None]
    assert arr.to_pylist() == expected

    for j, (i, v) in enumerate(zip(indices, expected)):
        s = arr[j]

        assert s.as_py() == v
        assert s.value.as_py() == v
        assert s.index.equals(i)
        assert s.dictionary.equals(dictionary)

        restored = pickle_module.loads(pickle_module.dumps(s))
        assert restored.equals(s)


def test_run_end_encoded():
    run_ends = [3, 5, 10, 12, 19]
    values = [1, 2, 1, None, 3]
    arr = pa.RunEndEncodedArray.from_arrays(run_ends, values)

    scalar = arr[0]
    assert isinstance(scalar, pa.RunEndEncodedScalar)
    assert isinstance(scalar.value, pa.Int64Scalar)
    assert scalar.value == pa.array(values)[0]
    assert scalar.as_py() == 1

    # null -> .value is still a scalar, as_py returns None
    scalar = arr[10]
    assert isinstance(scalar.value, pa.Int64Scalar)
    assert scalar.as_py() is None

    # constructing a scalar directly doesn't work yet
    with pytest.raises(NotImplementedError):
        pa.scalar(1, pa.run_end_encoded(pa.int64(), pa.int64()))


def test_union(pickle_module):
    # sparse
    arr = pa.UnionArray.from_sparse(
        pa.array([0, 0, 1, 1], type=pa.int8()),
        [
            pa.array(["a", "b", "c", "d"]),
            pa.array([1, 2, 3, 4])
        ]
    )
    for s in arr:
        s.validate(full=True)
        assert isinstance(s, pa.UnionScalar)
        assert s.type.equals(arr.type)
        assert s.is_valid is True
        with pytest.raises(pa.ArrowNotImplementedError):
            pickle_module.loads(pickle_module.dumps(s))

    assert arr[0].type_code == 0
    assert arr[0].as_py() == "a"
    assert arr[1].type_code == 0
    assert arr[1].as_py() == "b"
    assert arr[2].type_code == 1
    assert arr[2].as_py() == 3
    assert arr[3].type_code == 1
    assert arr[3].as_py() == 4

    # dense
    arr = pa.UnionArray.from_dense(
        types=pa.array([0, 1, 0, 0, 1, 1, 0], type='int8'),
        value_offsets=pa.array([0, 0, 2, 1, 1, 2, 3], type='int32'),
        children=[
            pa.array([b'a', b'b', b'c', b'd'], type='binary'),
            pa.array([1, 2, 3], type='int64')
        ]
    )
    for s in arr:
        s.validate(full=True)
        assert isinstance(s, pa.UnionScalar)
        assert s.type.equals(arr.type)
        assert s.is_valid is True
        with pytest.raises(pa.ArrowNotImplementedError):
            pickle_module.loads(pickle_module.dumps(s))

    assert arr[0].type_code == 0
    assert arr[0].as_py() == b'a'
    assert arr[5].type_code == 1
    assert arr[5].as_py() == 3


def test_map_scalar_as_py_with_custom_field_name():
    """
    Check we can call `MapScalar.as_py` with custom field names

    See https://github.com/apache/arrow/issues/36809
    """
    assert pa.scalar(
        [("foo", "bar")],
        pa.map_(
            pa.string(),
            pa.string()
        ),
    ).as_py() == [("foo", "bar")]

    assert pa.scalar(
        [("foo", "bar")],
        pa.map_(
            pa.field("custom_key", pa.string(), nullable=False),
            pa.field("custom_value", pa.string()),
        ),
    ).as_py() == [("foo", "bar")]
