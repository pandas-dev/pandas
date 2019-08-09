"""
Boilerplate functions used in defining binary operations.
"""
from functools import wraps

import numpy as np

from pandas._libs.lib import item_from_zerodim

from pandas.core.dtypes.generic import ABCDataFrame, ABCSeries, ABCIndexClass
from pandas.core.dtypes.common import is_list_like


def unpack_arraylike(obj, match_len: int):
	obj = item_from_zerodim(obj)

	if is_list_like(obj) and len(obj) != match_len:
		# NB: we must have already checked for if we need to
		#  return NotImplemented for DataFrame
		raise ValueError("Lengths must match")

	if isinstance(obj, list):
		# TODO: what about tuple?
		obj = np.asarray(obj)
	return obj


def unpack_and_defer(name):
	def wrapper(method):
		return _unpack_and_defer(method, name)

	return wrapper


def _unpack_and_defer(method, name):

	is_cmp = name.strip('__') in {'eq', 'ne', 'lt', 'le', 'gt', 'ge'}

	@wraps(method)
	def new_method(self, other):

		for cls in [ABCDataFrame, ABCSeries, ABCIndexClass]:
			if isinstance(self, cls):
				break
			if isinstance(other, cls):
				if is_cmp and cls is ABCSeries and isinstance(self, ABCIndexClass):
					# For comparison ops, Index does *not* defer to Series
					break
				else:
					return NotImplemented

		return method(self, other)

	return new_method
