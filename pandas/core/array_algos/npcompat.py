"""
Implementations of high-level numpy functions that are ExtensionArray-compatible. 
"""
import numpy as np

from pandas._typing import ArrayLike


def tile(arr: ArrayLike, shape) -> ArrayLike:
	raise NotImplementedError


def broadcast_to(arr: ArrayLike, shape) -> ArrayLike:
	raise NotImplementedError
