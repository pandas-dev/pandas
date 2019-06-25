"""
ExtensionArray subclasses with compatibility for 2-dimensional arrays
"""
from typing import Any, Tuple, Union

import numpy as np

from pandas._libs.lib import is_integer
from pandas.errors import AbstractMethodError

from pandas.core.arrays.base import ExtensionArray
from pandas.core.dtypes.generic import ABCPandasArray


class ReshapeableArray(ExtensionArray):
    """
    ReshapeableArray holds a non-reshape-able ExtensionArray and supports
    reshaping methods.
    """
    _allows_2d = True

    def __init__(self, values: ExtensionArray, shape: Tuple[int, ...]):
        assert (isinstance(values, ExtensionArray)
                and not values._allows_2d), type(values)
        assert not isinstance(values, ABCPandasArray)
        self._1dvalues = values

        assert np.prod(shape) == values.size, (np.prod(shape), values.size)
        self._shape = shape

    def __len__(self):
        return self.shape[0]

    @property
    def shape(self) -> Tuple[int, ...]:
        return self._shape

    # --------------------------------------------------
    # Direct pass-through attributes

    @property
    def dtype(self):
        return self._1dvalues.dtype

    @property
    def size(self) -> int:
        return self._1dvalues.size

    @property
    def nbytes(self) -> int:
        return self._1dvalues.nbytes

    def copy(self, deep: bool = False):
        result = self._1dvalues.copy(deep=deep)
        return type(self)(result, shape=self.shape)

    def _formatting_values(self):
        # TODO: should this be reshaped?
        return self._1dvalues._formatting_values()

    # NB: Not a classmethod since we need access to self._1dvalues
    def _from_factorized(self, values, original):
        result = self._1dvalues._from_factorized(values, original)
        shape = (result.size,)
        return type(self)(result, shape=shape)

    # NB: Not a classmethod since we need access to self._1dvalues
    def _from_sequence(self, scalars, dtype=None, copy=False):
        result = self._1dvalues._from_sequence(scalars, dtype=dtype, copy=copy)
        shape = (result.size,)
        return type(self)(result, shape=shape)

    # NB: Not a classmethod since we need access to self._1dvalues
    def _concat_same_type(self, to_concat):
        result = self._1dvalues._concat_same_type(to_concat)
        shape = (result.size,)
        return type(self)(result, shape=shape)

    def shift(self, periods: int = 1, fill_value: object = None):
        # FIXME: technically wrong to allow if we dont have ndim == 1

        result = self._1dvalues.shift(periods, fill_value=fill_value)
        return type(self)(result, shape=self.shape)

    # --------------------------------------------------
    # Lightly Modified pass-through methods

    def __repr__(self):
        head = ('<{cls}> shape={shape} Wrapping:\n'
                .format(cls=type(self).__name__, shape=self.shape))
        result = head + repr(self._1dvalues)
        return result

    def __iter__(self):
        if self.ndim == 1:
            for item in self._1dvalues:
                yield item
        else:
            for n in range(len(self)):
                yield self[n]

    def isna(self):
        result = self._1dvalues.isna()
        if isinstance(result, np.ndarray):
            result = result.reshape(self.shape)
        else:
            result = type(self)(result, shape=self.shape)
        return result

    def astype(self, dtype, copy=True):
        result = self._1dvalues.astype(dtype=dtype, copy=copy)
        if isinstance(result, np.ndarray):
            result = result.reshape(self.shape)
        else:
            result = type(self)(result, shape=self.shape)
        return result

    def fillna(self, value=None, method=None, limit=None):
        result = self._1dvalues.fillna(value=value, method=method, limit=limit)
        return type(self)(result, shape=self.shape)

    def __sub__(self, other):
        assert isinstance(other, type(self))
        assert other.shape == self.shape
        result = self._1dvalues - other._1dvalues
        return type(self)(result, shape=self.shape)

    def __array__(self, dtype=None):
        result = np.array(self._1dvalues, dtype=dtype)
        return result.reshape(self.shape)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        # implementing for sparse tests
        invals = list(inputs)
        invals = [x if x is not self else self._1dvalues for x in invals]
        invals = tuple(invals)
        result = getattr(ufunc, method)(*invals, **kwargs)
        if (isinstance(result, type(self._1dvalues))
                and result.size == self.size):
            return type(self)(result, shape=self.shape)
        return result

    # TODO: implement this for other comparisons; this one is needed
    #  for Categorical.replace to work in a pytables test.
    def __eq__(self, other):
        if np.ndim(other) == 0:
            # scalars, dont need to worry about alignment
            pass
        elif other.shape == self.shape:
            pass
        elif self.ndim > 1:
            # TODO: should we allow for the NotImplemented before this?
            raise NotImplementedError(self.shape, other.shape)

        result = self._1dvalues.__eq__(other)
        if result is NotImplemented:
            return result
        assert (isinstance(result, np.ndarray)
                and result.dtype == np.bool_), result
        return result.reshape(self.shape)

    def __ne__(self, other):
        eq = self.__eq__(other)
        if eq is NotImplemented:
            return NotImplemented
        return ~eq

    # --------------------------------------------------
    # Heavily-Modified pass-through methods

    def __getitem__(self, key):
        if self.ndim == 1:
            result = self._1dvalues[key]
            if np.ndim(result) == 0:
                # i.e. scalar
                return result
            shape = (result.size,)
            return type(self)(result, shape=shape)

        assert self.ndim == 2

        if isinstance(key, slice) and key == slice(None):
            # Note: we make a shallow copy
            return type(self)(self._1dvalues, shape=self.shape)

        if is_integer(key) and key == 0 and self.shape[0] == 1:
            # squeeze
            shape = (self.size,)
            return type(self)(self._1dvalues, shape=shape)

        if (isinstance(key, np.ndarray) and key.dtype == np.bool_
                and key.shape == (len(self),) and key.all()):
            return type(self)(self._1dvalues, shape=self.shape)

        if self.shape[0] != 1:
            raise NotImplementedError(key, self.shape)

        if not isinstance(key, tuple) or len(key) != 2:
            raise NotImplementedError(key, self.shape)

        if key[0] is Ellipsis:
            key = (slice(None), key[1])

        if key[0] == 0:
            result = self._1dvalues[key[1]]
            if np.ndim(result) == 0:
                return result
            if not isinstance(result, type(self._1dvalues)):
                # e.g. for object dtype
                # pandas/tests/sparse/test_indexing.py::test_frame_indexing_single
                return result
            shape = (result.size,)
            return type(self)(result, shape=shape)

        if key[0] == slice(None) and isinstance(key[1], slice):
            result = self._1dvalues[key[1]]
            shape = (1, result.size,)
            return type(self)(result, shape=shape)

        if key[0] == slice(None):
            # FIXME: in some places using tuple fails
            #  (e.g. DateTimearray, in others we get numpy warnings)
            result = self._1dvalues[[key[1]]]
            if np.ndim(result) == 0:
                return result
            if not isinstance(result, type(self._1dvalues)):
                # e.g. for object dtype
                #  pandas/tests/sparse/test_indexing.py::test_frame_indexing_single
                return result
            shape = (1, result.size)
            return type(self)(result, shape=shape)

        raise NotImplementedError(key, self.shape)

    def __setitem__(self, key: Union[int, np.ndarray], value: Any) -> None:
        if self.ndim == 1:
            # TODO: do we need to unpack value if it is wrapped in type(self)?
            self._1dvalues[key] = value
            return

        assert self.ndim == 2

        if (isinstance(key, tuple) and len(key) == 2
                and key[0] == 0 and self.shape[0] == 1):
            # TODO: Do we need to squeeze value?
            self._1dvalues[key[1]] = value
            return

        if (isinstance(key, np.ndarray) and key.dtype == np.bool_
                and key.shape == self.shape):
            if self.shape[0] == 1:
                key1 = key[0, :]
                if isinstance(value, np.ndarray) and value.shape == key.shape:
                    value = value[0, :]
                self._1dvalues[key1] = value
                return

        if isinstance(key, slice) and key == slice(None):
            if (isinstance(value, np.ndarray) and value.shape == self.shape
                    and self.shape[0] == 1):
                value = value[0, :]
            self._1dvalues[key] = value
            return

        raise NotImplementedError(key, self.shape)

    def take(self, indices, allow_fill=False, fill_value=None, axis=0):
        if self.ndim == 1 and axis == 0:
            result = self._1dvalues.take(indices, allow_fill=allow_fill,
                                         fill_value=fill_value)
            shape = (result.size,)
            return type(self)(result, shape=shape)

        assert self.ndim == 2
        if axis == 1 and self.shape[0] == 1:
            result = self._1dvalues.take(indices, allow_fill=allow_fill,
                                         fill_value=fill_value)
            shape = (1, result.size)
            return type(self)(result, shape)

        if axis == 0 and self.shape[1] == 1:
            result = self.T.take(indices, allow_fill=allow_fill,
                                 fill_value=fill_value, axis=1)
            return result.T

        raise NotImplementedError(indices, self.shape, axis)

    # --------------------------------------------------
    # Magic

    def __dir__(self):
        own = object.__dir__(self)
        inherited = dir(self._1dvalues)
        result = set(own).union(inherited)
        return list(result)

    def __getattr__(self, key):
        if key in object.__dir__(self):
            # TODO: why cant we do object.__hasattr__?
            # TODO: avoid getting method from base class
            return object.__getattribute__(self, key)

        values = object.__getattribute__(self, "_1dvalues")
        result = getattr(values, key)

        if isinstance(result, ExtensionArray):
            raise NotImplementedError(key)
        if isinstance(result, np.ndarray) and result.size == self.size:
            # FIXME: you need to wrap callables...
            return result.reshape(self.shape)
        return result

    # --------------------------------------------------
    # Reshape Methods

    def _copy_with_shape(self, shape):
        # NB: copy is _never_ deep
        shape = _tuplify_shape(self.size, shape)
        return type(self)(self._1dvalues, shape=shape)

    def reshape(self, *shape):
        # numpy accepts either a single tuple or an expanded tuple
        return self._copy_with_shape(shape)

    def transpose(self, axes):
        raise NotImplementedError(axes)

    @property
    def T(self):
        if self.ndim == 1:
            return self.copy(deep=False)
        if self.ndim == 2:
            shape = self.shape[::-1]
            return type(self)(self._1dvalues, shape=shape)
        raise NotImplementedError

    def ravel(self, order=None):
        if order is not None:
            raise NotImplementedError
        shape = (self.size,)
        return self._copy_with_shape(shape)

    def swapaxes(self, *axes):
        if axes == (0, 1) and self.ndim == 2:
            return self.T

        if axes == (1, 2) and self.shape[2] == 1 and self.ndim == 3:
            # pandas/core/reshape/reshape.py::get_new_values
            # TODO: uh check we're doing this right
            shape = (self.shape[0], 1, self.shape[1])
            return type(self)(self._1dvalues, shape=shape)
        raise NotImplementedError(axes, self.shape)


class ReshapeMixin:
    """
    Mixin for ExtensionArray subclasses that define `reshape` and related
    methods.

    Subclass must implement _wrap_data property.

    Notes
    -----
    - We assume that the constructor will accept:
        type(self)(self._wrap_data.reshape(shape), dtype=self.dtype)
      If not, then the methods below will need to be overriden.
    - We assume that the only 2D shapes taken will be (N, 1) and (1, N).
      This ensures that we can reshape, transpose, and ravel without worrying
      about column-order/row-order.
    """
    _allows_2d = True

    @property
    def _wrap_data(self) -> np.ndarray:
        """
        The underlying reshape-able array that we are wrapping.
        """
        raise AbstractMethodError(self)

    # --------------------------------------------------
    # Shape Attributes

    @property
    def shape(self) -> Tuple[int, ...]:
        """
        Return a tuple of the array dimensions.
        """
        return self._wrap_data.shape

    def __len__(self) -> int:
        return self.shape[0]

    # --------------------------------------------------
    # Reshape Methods

    def reshape(self, *shape):
        # numpy accepts either a single tuple or an expanded tuple
        data = self._wrap_data.reshape(*shape)
        return type(self)(data, dtype=self.dtype)

    def transpose(self, axes):
        data = self._wrap_data.transpose(axes)
        return type(self)(data, dtype=self.dtype)

    @property
    def T(self):
        data = self._wrap_data.T
        return type(self)(data, dtype=self.dtype)

    def ravel(self, order=None):
        data = self._wrap_data.ravel(order=order)
        return type(self)(data, dtype=self.dtype)

    def swapaxes(self, *axes):
        data = self._wrap_data.swapaxes(*axes)
        return type(self)(data, dtype=self.dtype)


def _tuplify_shape(size: int, shape) -> Tuple[int, ...]:
    """
    Convert a passed shape into a valid tuple.

    Following ndarray.reshape, we accept either `reshape(a, b)` or
    `reshape((a, b))`, the latter being canonical.

    Parameters
    ----------
    size : int
    shape : tuple

    Returns
    -------
    tuple[int, ...]
    """
    if len(shape) == 0:
        raise ValueError("shape must be a non-empty tuple of integers",
                         shape)

    if len(shape) == 1:
        if is_integer(shape[0]):
            pass
        else:
            shape = shape[0]
            if not isinstance(shape, tuple):
                raise ValueError("shape must be a non-empty tuple of integers",
                                 shape)

    if not all(is_integer(x) for x in shape):
        raise ValueError("shape must be a non-empty tuple of integers", shape)

    if any(x < -1 for x in shape):
        raise ValueError("Invalid shape {shape}".format(shape=shape))

    if -1 in shape:
        if shape.count(-1) != 1:
            raise ValueError("Invalid shape {shape}".format(shape=shape))
        idx = shape.index(-1)
        others = [n for n in shape if n != -1]
        prod = np.prod(others)
        dim = size // prod
        shape = shape[:idx] + (dim,) + shape[idx + 1:]

    if np.prod(shape) != size:
        raise ValueError("Product of shape ({shape}) must match "
                         "size ({size})".format(shape=shape,
                                                size=size))
    return shape


def unwrap_reshapeable(values):
    if isinstance(values, ReshapeableArray):
        # TODO: require we are only working with 1D?
        return values._1dvalues
    return values
