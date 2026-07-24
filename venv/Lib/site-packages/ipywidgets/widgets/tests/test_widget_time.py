# coding: utf-8

# Copyright (c) Vidar Tonaas Fauske.
# Distributed under the terms of the Modified BSD License.

import pytest

import datetime

from traitlets import TraitError

from ..widget_time import TimePicker


def test_time_creation_blank():
    w = TimePicker()
    assert w.value is None


def test_time_creation_value():
    t = datetime.time()
    w = TimePicker(value=t)
    assert w.value is t


def test_time_cross_validate_value_min_max():
    w = TimePicker(value=datetime.time(2), min=datetime.time(2), max=datetime.time(2))
    with w.hold_trait_notifications():
        w.value = None
        w.min = datetime.time(4)
        w.max = datetime.time(6)
    assert w.value is None
    with w.hold_trait_notifications():
        w.value = datetime.time(4)
        w.min = None
        w.max = None
    assert w.value == datetime.time(4)


def test_time_validate_value_none():
    t = datetime.time(13, 37, 42, 7)
    t_min = datetime.time(2)
    t_max = datetime.time(22)
    w = TimePicker(value=t, min=t_min, max=t_max)
    w.value = None
    assert w.value is None


def test_time_validate_value_vs_min():
    t = datetime.time(13, 37, 42, 7)
    t_min = datetime.time(14)
    t_max = datetime.time(22)
    w = TimePicker(min=t_min, max=t_max)
    w.value = t
    assert w.value.hour == 14


def test_time_validate_value_vs_max():
    t = datetime.time(13, 37, 42, 7)
    t_min = datetime.time(2)
    t_max = datetime.time(12)
    w = TimePicker(min=t_min, max=t_max)
    w.value = t
    assert w.value.hour == 12


def test_time_validate_min_vs_value():
    t = datetime.time(13, 37, 42, 7)
    t_min = datetime.time(14)
    t_max = datetime.time(22)
    w = TimePicker(value=t, max=t_max)
    w.min = t_min
    assert w.value.hour == 14


def test_time_validate_min_vs_max():
    t = datetime.time(13, 37, 42, 7)
    t_min = datetime.time(14)
    t_max = datetime.time(12)
    w = TimePicker(value=t, max=t_max)
    with pytest.raises(TraitError):
        w.min = t_min


def test_time_validate_max_vs_value():
    t = datetime.time(13, 37, 42, 7)
    t_min = datetime.time(2)
    t_max = datetime.time(12)
    w = TimePicker(value=t, min=t_min)
    w.max = t_max
    assert w.value.hour == 12


def test_time_validate_max_vs_min():
    t = datetime.time(13, 37, 42, 7)
    t_min = datetime.time(2)
    t_max = datetime.time(1)
    w = TimePicker(value=t, min=t_min)
    with pytest.raises(TraitError):
        w.max = t_max
