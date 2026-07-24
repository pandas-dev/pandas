# coding: utf-8

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import pytest

from contextlib import nullcontext
import datetime
import itertools

import pytz
from traitlets import TraitError

from ..widget_datetime import DatetimePicker


dt_1442 = datetime.datetime(1442, 1, 1, tzinfo=pytz.utc)
dt_1664 = datetime.datetime(1664, 1, 1, tzinfo=pytz.utc)
dt_1994 = datetime.datetime(1994, 1, 1, tzinfo=pytz.utc)
dt_2002 = datetime.datetime(2002, 2, 20, 13, 37, 42, 7, tzinfo=pytz.utc)
dt_2056 = datetime.datetime(2056, 1, 1, tzinfo=pytz.utc)

def test_time_creation_blank():
    w = DatetimePicker()
    assert w.value is None


def test_time_creation_value():
    dt = datetime.datetime.now(pytz.utc)
    w = DatetimePicker(value=dt)
    assert w.value is dt


def test_datetime_validate_value_none():
    dt = dt_2002
    dt_min = dt_1442
    dt_max = dt_2056
    w = DatetimePicker(value=dt, min=dt_min, max=dt_max)
    w.value = None
    assert w.value is None


def _permuted_dts():
    ret = []
    combos = list(itertools.product([None, dt_1442, dt_2002, dt_2056], repeat=3))
    for vals in combos:
        expected = vals[0]
        if vals[1] and vals[2] and vals[1] > vals[2]:
            expected = TraitError
        elif vals[0] is None:
            pass
        elif vals[1] and vals[1] > vals[0]:
            expected = vals[1]
        elif vals[2] and vals[2] < vals[0]:
            expected = vals[2]
        ret.append(vals + (expected,))
    return ret


@pytest.mark.parametrize(
    "input_value,input_min,input_max,expected",
    _permuted_dts()
)
def test_datetime_cross_validate_value_min_max(
    input_value,
    input_min,
    input_max,
    expected,
):
    w = DatetimePicker(value=dt_2002, min=dt_2002, max=dt_2002)
    should_raise = expected is TraitError
    with pytest.raises(expected) if should_raise else nullcontext():
        with w.hold_trait_notifications():
            w.value = input_value
            w.min = input_min
            w.max = input_max
    if not should_raise:
        assert w.value is expected


def test_datetime_validate_value_vs_min():
    dt = dt_2002
    dt_min = datetime.datetime(2019, 1, 1, tzinfo=pytz.utc)
    dt_max = dt_2056
    w = DatetimePicker(min=dt_min, max=dt_max)
    w.value = dt
    assert w.value.year == 2019


def test_datetime_validate_value_vs_max():
    dt = dt_2002
    dt_min = dt_1664
    dt_max = dt_1994
    w = DatetimePicker(min=dt_min, max=dt_max)
    w.value = dt
    assert w.value.year == 1994


def test_datetime_validate_min_vs_value():
    dt = dt_2002
    dt_min = datetime.datetime(2019, 1, 1, tzinfo=pytz.utc)
    dt_max = dt_2056
    w = DatetimePicker(value=dt, max=dt_max)
    w.min = dt_min
    assert w.value.year == 2019


def test_datetime_validate_min_vs_max():
    dt = dt_2002
    dt_min = datetime.datetime(2112, 1, 1, tzinfo=pytz.utc)
    dt_max = dt_2056
    w = DatetimePicker(value=dt, max=dt_max)
    with pytest.raises(TraitError):
        w.min = dt_min


def test_datetime_validate_max_vs_value():
    dt = dt_2002
    dt_min = dt_1664
    dt_max = dt_1994
    w = DatetimePicker(value=dt, min=dt_min)
    w.max = dt_max
    assert w.value.year == 1994


def test_datetime_validate_max_vs_min():
    dt = dt_2002
    dt_min = dt_1664
    dt_max = datetime.datetime(1337, 1, 1, tzinfo=pytz.utc)
    w = DatetimePicker(value=dt, min=dt_min)
    with pytest.raises(TraitError):
        w.max = dt_max


def test_datetime_validate_naive():
    dt = dt_2002
    dt_min = dt_1442
    dt_max = dt_2056

    w = DatetimePicker(value=dt, min=dt_min, max=dt_max)
    with pytest.raises(TraitError):
        w.max = dt_max.replace(tzinfo=None)
    with pytest.raises(TraitError):
        w.min = dt_min.replace(tzinfo=None)
    with pytest.raises(TraitError):
        w.value = dt.replace(tzinfo=None)


def test_datetime_tzinfo():
    tz = pytz.timezone('Australia/Sydney')
    dt = datetime.datetime(2002, 2, 20, 13, 37, 42, 7, tzinfo=tz)
    w = DatetimePicker(value=dt)
    assert w.value == dt
    # tzinfo only changes upon input from user
    assert w.value.tzinfo == tz
