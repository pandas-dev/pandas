import numpy as np
import pytest

from pandas._libs.tslibs import (
    Period,
    Resolution,
    to_offset,
)
from pandas._libs.tslibs.dtypes import _attrname_to_abbrevs

import pandas._testing as tm


@pytest.mark.parametrize(
    "freqstr,exp_freqstr",
    [("D", "D"), ("W", "D"), ("ME", "D"), ("s", "s"), ("min", "s"), ("h", "s")],
)
def test_get_to_timestamp_base(freqstr, exp_freqstr):
    off = to_offset(freqstr)
    per = Period._from_ordinal(1, off)
    exp_code = to_offset(exp_freqstr)._period_dtype_code

    result_code = per._dtype._get_to_timestamp_base()
    assert result_code == exp_code


@pytest.mark.parametrize(
    "freqstr,expected",
    [
        ("Y", "year"),
        ("Q", "quarter"),
        ("M", "month"),
        ("D", "day"),
        ("h", "hour"),
        ("min", "minute"),
        ("s", "second"),
        ("ms", "millisecond"),
        ("us", "microsecond"),
        ("ns", "nanosecond"),
    ],
)
def test_get_attrname_from_abbrev(freqstr, expected):
    assert Resolution.get_reso_from_freqstr(freqstr).attrname == expected


@pytest.mark.parametrize("freq", ["D", "h", "min", "s", "ms", "us", "ns"])
def test_get_freq_roundtrip2(freq):
    obj = Resolution.get_reso_from_freqstr(freq)
    result = _attrname_to_abbrevs[obj.attrname]
    assert freq == result


@pytest.mark.parametrize(
    "args,expected",
    [
        ((1.5, "min"), (90, "s")),
        ((62.4, "min"), (3744, "s")),
        ((1.04, "h"), (3744, "s")),
        ((1, "D"), (1, "D")),
        ((0.342931, "h"), (1234551600, "us")),
        ((1.2345, "D"), (106660800, "ms")),
    ],
)
def test_resolution_bumping(args, expected):
    # see gh-14378
    off = to_offset(str(args[0]) + args[1])
    assert off.n == expected[0]
    assert off._prefix == expected[1]


@pytest.mark.parametrize(
    "args",
    [
        (0.5, "ns"),
        # Too much precision in the input can prevent.
        (0.3429324798798269273987982, "h"),
    ],
)
def test_cat(args):
    msg = "Invalid frequency"

    with pytest.raises(ValueError, match=msg):
        to_offset(str(args[0]) + args[1])


@pytest.mark.parametrize(
    "freqstr,expected",
    [
        ("1h", "2021-01-01T09:00:00"),
        ("1D", "2021-01-02T08:00:00"),
        ("1W", "2021-01-03T08:00:00"),
        ("1ME", "2021-01-31T08:00:00"),
        ("1Y", "2021-12-31T08:00:00"),
    ],
)
def test_compatibility(freqstr, expected):
    ts_np = np.datetime64("2021-01-01T08:00:00.00")
    do = to_offset(freqstr)
    assert ts_np + do == np.datetime64(expected)


@pytest.mark.parametrize("freq", ["A", "H", "T", "S", "L", "U", "N"])
def test_units_A_H_T_S_L_U_N_deprecated_from_attrname_to_abbrevs(freq):
    # GH#52536
    msg = f"'{freq}' is deprecated and will be removed in a future version."

    with tm.assert_produces_warning(FutureWarning, match=msg):
        Resolution.get_reso_from_freqstr(freq)
