from datetime import (
    datetime,
    timedelta,
)
import operator
import zoneinfo

import numpy as np
import pytest

from pandas._libs.tslibs import iNaT
from pandas.errors import Pandas4Warning

import pandas as pd
import pandas._testing as tm
from pandas.core import roperator
from pandas.core.arrays import (
    DatetimeArray,
    PeriodArray,
    TimedeltaArray,
)


class TestNaTDeprecations:
    @pytest.mark.parametrize(
        "old_attr, new_attr",
        [
            ("dayofweek", "day_of_week"),
            ("dayofyear", "day_of_year"),
            ("daysinmonth", "days_in_month"),
        ],
    )
    def test_deprecated_day_attrs(self, old_attr, new_attr):
        # GH#46768
        msg = f"NaTType.{old_attr} is deprecated"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = getattr(pd.NaT, old_attr)
        assert np.isnan(result)


class TestNaTFormatting:
    def test_repr(self):
        assert repr(pd.NaT) == "NaT"

    def test_str(self):
        assert str(pd.NaT) == "NaT"

    def test_isoformat(self):
        assert pd.NaT.isoformat() == "NaT"


@pytest.mark.parametrize(
    "nat,idx",
    [
        (pd.Timestamp("NaT"), DatetimeArray),
        (pd.Timedelta("NaT"), TimedeltaArray),
        (pd.Period("NaT", freq="M"), PeriodArray),
    ],
)
def test_nat_fields(nat, idx):
    for field in idx._field_ops:
        # weekday is a property of DTI, but a method
        # on NaT/Timestamp for compat with datetime
        if field == "weekday":
            continue

        result = getattr(pd.NaT, field)
        assert np.isnan(result)

        result = getattr(nat, field)
        assert np.isnan(result)

    for field in idx._bool_ops:
        result = getattr(pd.NaT, field)
        assert result is False

        result = getattr(nat, field)
        assert result is False


def test_nat_vector_field_access():
    idx = pd.DatetimeIndex(["1/1/2000", None, None, "1/4/2000"])

    for field in DatetimeArray._field_ops:
        # weekday is a property of DTI, but a method
        # on NaT/Timestamp for compat with datetime
        if field == "weekday":
            continue

        result = getattr(idx, field)
        expected = pd.Index([getattr(x, field) for x in idx])
        tm.assert_index_equal(result, expected)

    ser = pd.Series(idx)

    for field in DatetimeArray._field_ops:
        # weekday is a property of DTI, but a method
        # on NaT/Timestamp for compat with datetime
        if field == "weekday":
            continue

        result = getattr(ser.dt, field)
        expected = [getattr(x, field) for x in idx]
        tm.assert_series_equal(result, pd.Series(expected))

    for field in DatetimeArray._bool_ops:
        result = getattr(ser.dt, field)
        expected = [getattr(x, field) for x in idx]
        tm.assert_series_equal(result, pd.Series(expected))


@pytest.mark.parametrize("klass", [pd.Timestamp, pd.Timedelta, pd.Period])
@pytest.mark.parametrize(
    "value", [None, np.nan, iNaT, float("nan"), pd.NaT, "NaT", "nat", "", "NAT"]
)
def test_identity(klass, value):
    assert klass(value) is pd.NaT


@pytest.mark.parametrize("klass", [pd.Timestamp, pd.Timedelta])
@pytest.mark.parametrize("method", ["round", "floor", "ceil"])
@pytest.mark.parametrize("freq", ["s", "5s", "min", "5min", "h", "5h"])
def test_round_nat(klass, method, freq):
    # see gh-14940
    ts = klass("nat")

    round_method = getattr(ts, method)
    assert round_method(freq) is ts


@pytest.mark.parametrize(
    "method",
    [
        "astimezone",
        "combine",
        "ctime",
        "dst",
        "fromordinal",
        "fromtimestamp",
        "fromisocalendar",
        "isocalendar",
        "strftime",
        "strptime",
        "time",
        "timestamp",
        "timetuple",
        "timetz",
        "toordinal",
        "tzname",
        "utcfromtimestamp",
        "utcnow",
        "utcoffset",
        "utctimetuple",
    ],
)
def test_nat_methods_raise(method):
    # see gh-9513, gh-17329
    msg = f"NaTType does not support {method}"

    with pytest.raises(ValueError, match=msg):
        getattr(pd.NaT, method)()


@pytest.mark.parametrize("method", ["weekday", "isoweekday"])
def test_nat_methods_nan(method):
    # see gh-9513, gh-17329
    assert np.isnan(getattr(pd.NaT, method)())


@pytest.mark.parametrize(
    "method", ["date", "now", "replace", "today", "tz_convert", "tz_localize"]
)
def test_nat_methods_nat(method):
    # see gh-8254, gh-9513, gh-17329
    assert getattr(pd.NaT, method)() is pd.NaT


@pytest.mark.parametrize(
    "get_nat", [lambda x: pd.NaT, lambda x: pd.Timedelta(x), lambda x: pd.Timestamp(x)]
)
def test_nat_iso_format(get_nat):
    # see gh-12300
    assert get_nat("NaT").isoformat() == "NaT"
    assert get_nat("NaT").isoformat(timespec="nanoseconds") == "NaT"


@pytest.mark.parametrize(
    "klass,expected",
    [
        (pd.Timestamp, ["normalize", "to_julian_date", "to_period", "unit"]),
        (
            pd.Timedelta,
            [
                "components",
                "resolution_string",
                "to_pytimedelta",
                "to_timedelta64",
                "unit",
                "view",
            ],
        ),
    ],
)
def test_missing_public_nat_methods(klass, expected):
    # see gh-17327
    #
    # NaT should have *most* of the Timestamp and Timedelta methods.
    # Here, we check which public methods NaT does not have. We
    # ignore any missing private methods.
    nat_names = dir(pd.NaT)
    klass_names = dir(klass)

    missing = [x for x in klass_names if x not in nat_names and not x.startswith("_")]
    missing.sort()

    assert missing == expected


def _get_overlap_public_nat_methods(klass, as_tuple=False):
    """
    Get overlapping public methods between NaT and another class.

    Parameters
    ----------
    klass : type
        The class to compare with NaT
    as_tuple : bool, default False
        Whether to return a list of tuples of the form (klass, method).

    Returns
    -------
    overlap : list
    """
    nat_names = dir(pd.NaT)
    klass_names = dir(klass)

    overlap = [
        x
        for x in nat_names
        if x in klass_names and not x.startswith("_") and callable(getattr(klass, x))
    ]

    # Timestamp takes precedence over Timedelta in terms of overlap.
    if klass is pd.Timedelta:
        ts_names = dir(pd.Timestamp)
        overlap = [x for x in overlap if x not in ts_names]

    if as_tuple:
        overlap = [(klass, method) for method in overlap]

    overlap.sort()
    return overlap


@pytest.mark.parametrize(
    "klass,expected",
    [
        (
            pd.Timestamp,
            [
                "as_unit",
                "astimezone",
                "ceil",
                "combine",
                "ctime",
                "date",
                "day_name",
                "dst",
                "floor",
                "fromisocalendar",
                "fromisoformat",
                "fromordinal",
                "fromtimestamp",
                "isocalendar",
                "isoformat",
                "isoweekday",
                "month_name",
                "now",
                "replace",
                "round",
                "strftime",
                "strptime",
                "time",
                "timestamp",
                "timetuple",
                "timetz",
                "to_datetime64",
                "to_numpy",
                "to_pydatetime",
                "today",
                "toordinal",
                "tz_convert",
                "tz_localize",
                "tzname",
                "utcfromtimestamp",
                "utcnow",
                "utcoffset",
                "utctimetuple",
                "weekday",
            ],
        ),
        (pd.Timedelta, ["total_seconds"]),
    ],
)
def test_overlap_public_nat_methods(klass, expected):
    # see gh-17327
    #
    # NaT should have *most* of the Timestamp and Timedelta methods.
    # In case when Timestamp, Timedelta, and NaT are overlap, the overlap
    # is considered to be with Timestamp and NaT, not Timedelta.
    assert _get_overlap_public_nat_methods(klass) == expected


@pytest.mark.parametrize(
    "compare",
    (
        _get_overlap_public_nat_methods(pd.Timestamp, True)
        + _get_overlap_public_nat_methods(pd.Timedelta, True)
    ),
    ids=lambda x: f"{x[0].__name__}.{x[1]}",
)
def test_nat_doc_strings(compare):
    # see gh-17327
    #
    # The docstrings for overlapping methods should match.
    klass, method = compare
    klass_doc = getattr(klass, method).__doc__

    if klass == pd.Timestamp and method == "isoformat":
        pytest.skip(
            "Ignore differences with Timestamp.isoformat() as they're intentional"
        )

    if method == "to_numpy":
        # GH#44460 can return either dt64 or td64 depending on dtype,
        #  different docstring is intentional
        pytest.skip(f"different docstring for {method} is intentional")

    nat_doc = getattr(pd.NaT, method).__doc__
    assert klass_doc == nat_doc


_ops = {
    "left_plus_right": lambda a, b: a + b,
    "right_plus_left": lambda a, b: b + a,
    "left_minus_right": lambda a, b: a - b,
    "right_minus_left": lambda a, b: b - a,
    "left_times_right": lambda a, b: a * b,
    "right_times_left": lambda a, b: b * a,
    "left_div_right": lambda a, b: a / b,
    "right_div_left": lambda a, b: b / a,
}


@pytest.mark.parametrize("op_name", list(_ops.keys()))
@pytest.mark.parametrize(
    "value,val_type",
    [
        (2, "scalar"),
        (1.5, "floating"),
        (np.nan, "floating"),
        ("foo", "str"),
        (timedelta(3600), "timedelta"),
        (pd.Timedelta("5s"), "timedelta"),
        (datetime(2014, 1, 1), "timestamp"),
        (pd.Timestamp("2014-01-01"), "timestamp"),
        (pd.Timestamp("2014-01-01", tz="UTC"), "timestamp"),
        (pd.Timestamp("2014-01-01", tz="US/Eastern"), "timestamp"),
        (datetime(2014, 1, 1).astimezone(zoneinfo.ZoneInfo("Asia/Tokyo")), "timestamp"),
    ],
)
def test_nat_arithmetic_scalar(op_name, value, val_type):
    # see gh-6873
    invalid_ops = {
        "scalar": {"right_div_left"},
        "floating": {
            "right_div_left",
            "left_minus_right",
            "right_minus_left",
            "left_plus_right",
            "right_plus_left",
        },
        "str": set(_ops.keys()),
        "timedelta": {"left_times_right", "right_times_left"},
        "timestamp": {
            "left_times_right",
            "right_times_left",
            "left_div_right",
            "right_div_left",
        },
    }

    op = _ops[op_name]

    if op_name in invalid_ops.get(val_type, set()):
        if (
            val_type == "timedelta"
            and "times" in op_name
            and isinstance(value, pd.Timedelta)
        ):
            typs = "(Timedelta|NaTType)"
            msg = rf"unsupported operand type\(s\) for \*: '{typs}' and '{typs}'"
        elif val_type == "str":
            # un-specific check here because the message comes from str
            #  and varies by method
            msg = "|".join(
                [
                    "can only concatenate str",
                    "unsupported operand type",
                    "can't multiply sequence",
                    "Can't convert 'NaTType'",
                    "must be str, not NaTType",
                ]
            )
        else:
            msg = "unsupported operand type"

        with pytest.raises(TypeError, match=msg):
            op(pd.NaT, value)
    else:
        if val_type == "timedelta" and "div" in op_name:
            expected = np.nan
        else:
            expected = pd.NaT

        assert op(pd.NaT, value) is expected


@pytest.mark.parametrize(
    "val,expected",
    [(np.nan, pd.NaT), (pd.NaT, np.nan), (np.timedelta64("NaT", "ns"), np.nan)],
)
def test_nat_rfloordiv_timedelta(val, expected):
    # see gh-#18846
    #
    # See also test_timedelta.TestTimedeltaArithmetic.test_floordiv
    td = pd.Timedelta(hours=3, minutes=4)
    assert td // val is expected


@pytest.mark.parametrize(
    "op_name",
    ["left_plus_right", "right_plus_left", "left_minus_right", "right_minus_left"],
)
@pytest.mark.parametrize(
    "value",
    [
        pd.DatetimeIndex(["2011-01-01", "2011-01-02"], dtype="M8[ns]", name="x"),
        pd.DatetimeIndex(
            ["2011-01-01", "2011-01-02"], dtype="M8[ns, US/Eastern]", name="x"
        ),
        DatetimeArray._from_sequence(["2011-01-01", "2011-01-02"], dtype="M8[ns]"),
        DatetimeArray._from_sequence(
            ["2011-01-01", "2011-01-02"], dtype=pd.DatetimeTZDtype(tz="US/Pacific")
        ),
        pd.TimedeltaIndex(["1 day", "2 day"], name="x"),
    ],
)
def test_nat_arithmetic_index(op_name, value):
    # see gh-11718
    exp_name = "x"
    exp_data = [pd.NaT] * 2

    if value.dtype.kind == "M" and "plus" in op_name:
        expected = pd.DatetimeIndex(exp_data, tz=value.tz, name=exp_name)
    else:
        expected = pd.TimedeltaIndex(exp_data, name=exp_name)
    expected = expected.as_unit(value.unit)

    if not isinstance(value, pd.Index):
        expected = expected.array

    op = _ops[op_name]
    result = op(pd.NaT, value)
    tm.assert_equal(result, expected)


@pytest.mark.parametrize(
    "op_name",
    ["left_plus_right", "right_plus_left", "left_minus_right", "right_minus_left"],
)
@pytest.mark.parametrize(
    "box", [pd.TimedeltaIndex, pd.Series, TimedeltaArray._from_sequence]
)
def test_nat_arithmetic_td64_vector(op_name, box):
    # see gh-19124
    vec = box(["1 day", "2 day"], dtype="timedelta64[ns]")
    box_nat = box([pd.NaT, pd.NaT], dtype="timedelta64[ns]")
    tm.assert_equal(_ops[op_name](vec, pd.NaT), box_nat)


@pytest.mark.parametrize(
    "dtype,op,out_dtype",
    [
        ("datetime64[ns]", operator.add, "datetime64[ns]"),
        ("datetime64[ns]", roperator.radd, "datetime64[ns]"),
        ("datetime64[ns]", operator.sub, "timedelta64[ns]"),
        ("datetime64[ns]", roperator.rsub, "timedelta64[ns]"),
        ("timedelta64[ns]", operator.add, "datetime64[ns]"),
        ("timedelta64[ns]", roperator.radd, "datetime64[ns]"),
        ("timedelta64[ns]", operator.sub, "datetime64[ns]"),
        ("timedelta64[ns]", roperator.rsub, "timedelta64[ns]"),
    ],
)
def test_nat_arithmetic_ndarray(dtype, op, out_dtype):
    other = np.arange(10).astype(dtype)
    result = op(pd.NaT, other)

    expected = np.empty(other.shape, dtype=out_dtype)
    expected.fill("NaT")
    tm.assert_numpy_array_equal(result, expected)


def test_nat_pinned_docstrings():
    # see gh-17327
    assert pd.NaT.ctime.__doc__ == pd.Timestamp.ctime.__doc__


def test_to_numpy_alias():
    # GH 24653: alias .to_numpy() for scalars
    expected = pd.NaT.to_datetime64()
    result = pd.NaT.to_numpy()

    assert pd.isna(expected) and pd.isna(result)

    # GH#44460
    result = pd.NaT.to_numpy("M8[s]")
    assert isinstance(result, np.datetime64)
    assert result.dtype == "M8[s]"

    result = pd.NaT.to_numpy("m8[ns]")
    assert isinstance(result, np.timedelta64)
    assert result.dtype == "m8[ns]"

    result = pd.NaT.to_numpy("m8[s]")
    assert isinstance(result, np.timedelta64)
    assert result.dtype == "m8[s]"

    with pytest.raises(ValueError, match="NaT.to_numpy dtype must be a "):
        pd.NaT.to_numpy(np.int64)


@pytest.mark.parametrize(
    "other",
    [
        pd.Timedelta(0),
        pd.Timedelta(0).to_pytimedelta(),
        pd.Timedelta(0).to_timedelta64(),
        pd.Timestamp(0),
        pd.Timestamp(0).to_pydatetime(),
        pd.Timestamp(0).to_datetime64(),
        pd.Timestamp(0).tz_localize("UTC"),
        pd.NaT,
    ],
)
def test_nat_comparisons(compare_operators_no_eq_ne, other):
    # GH 26039
    opname = compare_operators_no_eq_ne

    assert getattr(pd.NaT, opname)(other) is False

    op = getattr(operator, opname.strip("_"))
    assert op(pd.NaT, other) is False
    assert op(other, pd.NaT) is False


@pytest.mark.parametrize("other_and_type", [("foo", "str"), (2, "int"), (2.0, "float")])
@pytest.mark.parametrize(
    "symbol_and_op",
    [("<=", operator.le), ("<", operator.lt), (">=", operator.ge), (">", operator.gt)],
)
def test_nat_comparisons_invalid(other_and_type, symbol_and_op):
    # GH#35585
    other, other_type = other_and_type
    symbol, op = symbol_and_op

    assert not pd.NaT == other
    assert not other == pd.NaT

    assert pd.NaT != other
    assert other != pd.NaT

    msg = f"'{symbol}' not supported between instances of 'NaTType' and '{other_type}'"
    with pytest.raises(TypeError, match=msg):
        op(pd.NaT, other)

    msg = f"'{symbol}' not supported between instances of '{other_type}' and 'NaTType'"
    with pytest.raises(TypeError, match=msg):
        op(other, pd.NaT)


@pytest.mark.parametrize(
    "other",
    [
        np.array(["foo"] * 2, dtype=object),
        np.array([2, 3], dtype="int64"),
        np.array([2.0, 3.5], dtype="float64"),
    ],
    ids=["str", "int", "float"],
)
def test_nat_comparisons_invalid_ndarray(other):
    # GH#40722
    expected = np.array([False, False])
    result = pd.NaT == other
    tm.assert_numpy_array_equal(result, expected)
    result = other == pd.NaT
    tm.assert_numpy_array_equal(result, expected)

    expected = np.array([True, True])
    result = pd.NaT != other
    tm.assert_numpy_array_equal(result, expected)
    result = other != pd.NaT
    tm.assert_numpy_array_equal(result, expected)

    for symbol, op in [
        ("<=", operator.le),
        ("<", operator.lt),
        (">=", operator.ge),
        (">", operator.gt),
    ]:
        msg = f"'{symbol}' not supported between"

        with pytest.raises(TypeError, match=msg):
            op(pd.NaT, other)

        if other.dtype == np.dtype("object"):
            # uses the reverse operator, so symbol changes
            msg = None
        with pytest.raises(TypeError, match=msg):
            op(other, pd.NaT)


def test_compare_date(fixed_now_ts):
    # GH#39151 comparing NaT with date object is deprecated
    # See also: tests.scalar.timestamps.test_comparisons::test_compare_date

    dt = fixed_now_ts.to_pydatetime().date()

    msg = "Cannot compare NaT with datetime.date object"
    for left, right in [(pd.NaT, dt), (dt, pd.NaT)]:
        assert not left == right
        assert left != right

        with pytest.raises(TypeError, match=msg):
            left < right
        with pytest.raises(TypeError, match=msg):
            left <= right
        with pytest.raises(TypeError, match=msg):
            left > right
        with pytest.raises(TypeError, match=msg):
            left >= right


@pytest.mark.parametrize(
    "obj",
    [
        pd.offsets.YearEnd(2),
        pd.offsets.YearBegin(2),
        pd.offsets.MonthBegin(1),
        pd.offsets.MonthEnd(2),
        pd.offsets.MonthEnd(12),
        pd.offsets.Day(2),
        pd.offsets.Day(5),
        pd.offsets.Hour(24),
        pd.offsets.Hour(3),
        pd.offsets.Minute(),
        np.timedelta64(3, "h"),
        np.timedelta64(4, "h"),
        np.timedelta64(3200, "s"),
        np.timedelta64(3600, "s"),
        np.timedelta64(3600 * 24, "s"),
        np.timedelta64(2, "D"),
        np.timedelta64(365, "D"),
        timedelta(-2),
        timedelta(365),
        timedelta(minutes=120),
        timedelta(days=4, minutes=180),
        timedelta(hours=23),
        timedelta(hours=23, minutes=30),
        timedelta(hours=48),
    ],
)
def test_nat_addsub_tdlike_scalar(obj):
    assert pd.NaT + obj is pd.NaT
    assert obj + pd.NaT is pd.NaT
    assert pd.NaT - obj is pd.NaT


def test_pickle(temp_file):
    # GH#4606
    p = tm.round_trip_pickle(pd.NaT, temp_file)
    assert p is pd.NaT
