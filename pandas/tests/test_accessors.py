#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

An example/recipe/test for implementing custom accessors.

"""

import pandas as pd

from pandas.core.accessors import (wrap_delegate_names,
								   PandasDelegate, AccessorProperty)

class State(object):
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
        #self._freeze()

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





def test_geo_state_accessor():
	import pandas.util.testing as tm

	pd.Series.state = AccessorProperty(StateDelegate)

	ser = pd.Series([State('Alabama'), State('California')])

	abbrev = pd.Series(['AL', 'CA'])
	tm.assert_series_equal(ser.state.abbrev, abbrev)

	fips = pd.Series([1, 6])
	tm.assert_series_equal(ser.state.fips(), fips)



	ser.state.abbrev = ['Foo', 'Bar']

	new_abbrev = pd.Series(['Foo', 'Bar'])
	tm.assert_series_equal(ser.state.abbrev, new_abbrev)



