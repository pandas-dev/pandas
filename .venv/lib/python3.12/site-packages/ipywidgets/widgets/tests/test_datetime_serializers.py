# coding: utf-8

# Copyright (c) Vidar Tonaas Fauske.
# Distributed under the terms of the Modified BSD License.

import pytest

import datetime
import pytz

from traitlets import TraitError

from ..trait_types import (
    time_to_json,
    time_from_json,
    datetime_to_json,
    datetime_from_json,
)


def test_time_serialize_none():
    assert time_to_json(None, None) == None


def test_time_serialize_value():
    t = datetime.time(13, 37, 42, 7000)
    assert time_to_json(t, None) == dict(
        hours=13, minutes=37, seconds=42, milliseconds=7
    )


def test_time_deserialize_none():
    assert time_from_json(None, None) == None


def test_time_deserialize_value():
    v = dict(hours=13, minutes=37, seconds=42, milliseconds=7)
    assert time_from_json(v, None) == datetime.time(13, 37, 42, 7000)


def test_datetime_serialize_none():
    assert datetime_to_json(None, None) == None


def test_datetime_serialize_value():
    t = datetime.datetime(2002, 2, 20, 13, 37, 42, 7000, pytz.utc)
    assert datetime_to_json(t, None) == dict(
        year=2002,
        month=1,  # Months are 0-based indices in JS
        date=20,
        hours=13,
        minutes=37,
        seconds=42,
        milliseconds=7,
    )


def test_datetime_serialize_non_utz():
    # Non-existent timezone, so it will never be the local one:
    tz = pytz.FixedOffset(42)
    t = datetime.datetime(2002, 2, 20, 13, 37, 42, 7000, tz)
    assert datetime_to_json(t, None) == dict(
        year=2002,
        month=1,  # Months are 0-based indices in JS
        date=20,
        hours=12,
        minutes=55,
        seconds=42,
        milliseconds=7,
    )


def test_datetime_deserialize_none():
    assert datetime_from_json(None, None) == None


def test_datetime_deserialize_value():
    tz = pytz.FixedOffset(42)
    v = dict(
        year=2002,
        month=1,  # Months are 0-based indices in JS
        date=20,
        hours=13,
        minutes=37,
        seconds=42,
        milliseconds=7,
    )
    assert datetime_from_json(v, None) == datetime.datetime(
        2002, 2, 20, 14, 19, 42, 7000, tz
    )
