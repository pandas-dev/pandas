"""
Implementations of high-level numpy functions that are ExtensionArray-compatible. 
"""
import numpy as np

from pandas._typing import ArrayLike


def tile(arr: ArrayLike, shape) -> ArrayLike:
	raise NotImplementedError


def broadcast_to(arr: ArrayLike, shape) -> ArrayLike:
	if isinstance(arr, np.ndarray):
		return np.broadcast_to(arr, shape)

	values = arr._values_for_factorize()[0]

	btvalues = np.broadcast_to(values, shape)
	result = type(arr)._from_factorized(btvalues, arr)
	return result
