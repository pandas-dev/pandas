#!/usr/bin/env python
# coding: utf-8

import unittest
import nose

import datetime
from pandas.msgpack import packb, unpackb

class DummyException(Exception):
    pass

class TestExceptions(unittest.TestCase):

    def test_raise_on_find_unsupported_value(self):
        import datetime
        self.assertRaises(TypeError, packb, datetime.datetime.now())

    def test_raise_from_object_hook(self):
        def hook(obj):
            raise DummyException
        self.assertRaises(DummyException, unpackb, packb({}), object_hook=hook)
        self.assertRaises(DummyException, unpackb, packb({'fizz': 'buzz'}), object_hook=hook)
        self.assertRaises(DummyException, unpackb, packb({'fizz': 'buzz'}), object_pairs_hook=hook)
        self.assertRaises(DummyException, unpackb, packb({'fizz': {'buzz': 'spam'}}), object_hook=hook)
        self.assertRaises(DummyException, unpackb, packb({'fizz': {'buzz': 'spam'}}), object_pairs_hook=hook)

    def test_invalidvalue(self):
        self.assertRaises(ValueError, unpackb, b'\xd9\x97#DL_')
