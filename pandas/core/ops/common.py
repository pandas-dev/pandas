"""
Boilerplate ops code shared by many methods.
"""
from functools import wraps

from pandas._libs import lib

from pandas.core.dtypes.generic import ABCDataFrame


def unpack_and_defer(method):

	@wraps(method)
	def new_method(self, other, *args, **kwargs):
		# args/kwargs included for compat with flex methods

		if isinstance(other, ABCDataFrame) and not isinstance(self, ABCDataFrame):
			return NotImplemented

		other = lib.item_from_zerodim(other)

		return method(self, other, *args, **kwargs)

	return new_method
