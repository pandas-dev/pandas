#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

An example/recipe/test for implementing custom accessors.

"""
import unittest
import pandas.util.testing as tm

import pandas as pd

from pandas.core.accessors import (wrap_delegate_names,
                                   PandasDelegate, AccessorProperty)

# Example 1:
# An accessor for attributes of custom class in a Series with object dtype.


class State(object):
    """
    A dummy class for which only two states have the attributes implemented.
    """
    def __repr__(self):
        return repr(self.name)

    def __init__(self, name):
        self.name = name
        self._abbrev_dict = {'California': 'CA', 'Alabama': 'AL'}

    @property
    def abbrev(self):
        return self._abbrev_dict[self.name]

    @abbrev.setter
    def abbrev(self, value):
        self._abbrev_dict[self.name] = value

    def fips(self):
        return {'California': 6, 'Alabama': 1}[self.name]


@wrap_delegate_names(delegate=State,
                     accessors=["fips"],
                     typ="method")
@wrap_delegate_names(delegate=State,
                     accessors=["abbrev"],
                     typ="property")
class StateDelegate(PandasDelegate):

    def __init__(self, values):
        self.values = values

    @classmethod
    def _make_accessor(cls, data):
        """
        When implementing custom accessors, `_make_accessor` is the place
        to do validation that the attributes be accessed will actually be
        present in the underlying data.
        """
        if not isinstance(data, pd.Series):
            raise ValueError('Input must be a Series of States')
        elif not data.apply(lambda x: isinstance(x, State)).all():
            raise ValueError('All entries must be State objects')
        return StateDelegate(data)

    def _delegate_method(self, name, *args, **kwargs):
        state_method = lambda x: getattr(x, name)(*args, **kwargs)
        return self.values.apply(state_method)

    def _delegate_property_get(self, name):
        state_property = lambda x: getattr(x, name)
        return self.values.apply(state_property)

    def _delegate_property_set(self, name, new_values):
        """
        Setting properties via accessors is permitted but discouraged.
        """
        for (obj, val) in zip(self.values, new_values):
            setattr(obj, name, val)


class TestVectorizedAccessor(unittest.TestCase):

    @classmethod
    def setup_class(cls):
        pd.Series.state = AccessorProperty(StateDelegate)

        cls.ser = pd.Series([State('Alabama'), State('California')])

    @classmethod
    def teardown_class(cls):
        del pd.Series.state
        # TODO: is there a nicer way to do this with `mock`?

    def test_method(self):
        ser = self.ser
        fips = pd.Series([1, 6])
        tm.assert_series_equal(ser.state.fips(), fips)

    def test_property_get(self):
        ser = self.ser
        abbrev = pd.Series(['AL', 'CA'])
        tm.assert_series_equal(ser.state.abbrev, abbrev)

    def test_property_set(self):
        ser = self.ser.copy()

        ser.state.abbrev = ['Foo', 'Bar']
        new_abbrev = pd.Series(['Foo', 'Bar'])
        tm.assert_series_equal(ser.state.abbrev, new_abbrev)


@wrap_delegate_names(delegate=pd.Series,
                     accessors=["real", "imag"],
                     typ="property")
@wrap_delegate_names(delegate=pd.Series,
                     accessors=["abs"],
                     typ="method")
class ForgotToOverride(PandasDelegate):
    # A case where the relevant methods were not overridden.  Everything
    # should raise NotImplementedError or TypeError
    @classmethod
    def _make_accessor(cls, data):
        return cls(data)


class TestUnDelegated(unittest.TestCase):

    @classmethod
    def setup_class(cls):
        pd.Series.forgot = AccessorProperty(ForgotToOverride)

        cls.ser = pd.Series(range(-2, 2))

    @classmethod
    def teardown_class(cls):
        del pd.Series.forgot

    def test_get_fails(self):
        forgot = self.ser.forgot
        with self.assertRaises(TypeError):
            forgot.real

        with self.assertRaises(TypeError):
            forgot.imag

    def test_set_fails(self):
        forgot = self.ser.forgot
        with self.assertRaises(TypeError):
            forgot.real = range(5)

        # Check that the underlying hasn't been affected
        tm.assert_series_equal(self.ser, pd.Series(range(-2, 2)))

    def test_method_fails(self):
        forgot = self.ser.forgot
        with self.assertRaises(TypeError):
            forgot.abs()
