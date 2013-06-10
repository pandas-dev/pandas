#!/usr/bin/env python
# coding: utf-8

from pytest import raises
from pandas.msgpack import packb, unpackb

import datetime


class DummyException(Exception):
    pass


def test_raise_on_find_unsupported_value():
    with raises(TypeError):
        packb(datetime.datetime.now())


def test_raise_from_object_hook():
    def hook(obj):
        raise DummyException
    raises(DummyException, unpackb, packb({}), object_hook=hook)
    raises(DummyException, unpackb, packb({'fizz': 'buzz'}), object_hook=hook)
    raises(DummyException, unpackb, packb({'fizz': 'buzz'}), object_pairs_hook=hook)
    raises(DummyException, unpackb, packb({'fizz': {'buzz': 'spam'}}), object_hook=hook)
    raises(DummyException, unpackb, packb({'fizz': {'buzz': 'spam'}}), object_pairs_hook=hook)


def test_invalidvalue():
    with raises(ValueError):
        unpackb(b'\xd9\x97#DL_')
