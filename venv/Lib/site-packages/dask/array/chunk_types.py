from __future__ import annotations

import numpy as np

# Start list of valid chunk types, to be added to with guarded imports
_HANDLED_CHUNK_TYPES = [np.ndarray, np.ma.MaskedArray]


def register_chunk_type(type):
    """Register the given type as a valid chunk and downcast array type

    Parameters
    ----------
    type : type
        Duck array type to be registered as a type Dask can safely wrap as a chunk and
        to which Dask does not defer in arithmetic operations and NumPy
        functions/ufuncs.

    Notes
    -----
    A :py:class:`dask.array.Array` can contain any sufficiently "NumPy-like" array in
    its chunks. These are also referred to as "duck arrays" since they match the most
    important parts of NumPy's array API, and so, behave the same way when relying on
    duck typing.

    However, for multiple duck array types to interoperate properly, they need to
    properly defer to each other in arithmetic operations and NumPy functions/ufuncs
    according to a well-defined type casting hierarchy (
    `see NEP 13 <https://numpy.org/neps/nep-0013-ufunc-overrides.html#type-casting-hierarchy>`__
    ). In an effort to maintain this hierarchy, Dask defers to all other duck array
    types except those in its internal registry. By default, this registry contains

    * :py:class:`numpy.ndarray`
    * :py:class:`numpy.ma.MaskedArray`
    * :py:class:`cupy.ndarray`
    * :py:class:`sparse.SparseArray`
    * :py:class:`scipy.sparse.spmatrix`

    This function exists to append any other types to this registry. If a type is not
    in this registry, and yet is a downcast type (it comes below
    :py:class:`dask.array.Array` in the type casting hierarchy), a ``TypeError`` will
    be raised due to all operand types returning ``NotImplemented``.

    Examples
    --------
    Using a mock ``FlaggedArray`` class as an example chunk type unknown to Dask with
    minimal duck array API:

    >>> import numpy.lib.mixins
    >>> class FlaggedArray(numpy.lib.mixins.NDArrayOperatorsMixin):
    ...     def __init__(self, a, flag=False):
    ...         self.a = a
    ...         self.flag = flag
    ...     def __repr__(self):
    ...         return f"Flag: {self.flag}, Array: " + repr(self.a)
    ...     def __array__(self):
    ...         return np.asarray(self.a)
    ...     def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
    ...         if method == '__call__':
    ...             downcast_inputs = []
    ...             flag = False
    ...             for input in inputs:
    ...                 if isinstance(input, self.__class__):
    ...                     flag = flag or input.flag
    ...                     downcast_inputs.append(input.a)
    ...                 elif isinstance(input, np.ndarray):
    ...                     downcast_inputs.append(input)
    ...                 else:
    ...                     return NotImplemented
    ...             return self.__class__(ufunc(*downcast_inputs, **kwargs), flag)
    ...         else:
    ...             return NotImplemented
    ...     @property
    ...     def shape(self):
    ...         return self.a.shape
    ...     @property
    ...     def ndim(self):
    ...         return self.a.ndim
    ...     @property
    ...     def dtype(self):
    ...         return self.a.dtype
    ...     def __getitem__(self, key):
    ...         return type(self)(self.a[key], self.flag)
    ...     def __setitem__(self, key, value):
    ...         self.a[key] = value

    Before registering ``FlaggedArray``, both types will attempt to defer to the
    other:

    >>> import dask.array as da
    >>> da.ones(5) - FlaggedArray(np.ones(5), True)
    Traceback (most recent call last):
    ...
    TypeError: operand type(s) all returned NotImplemented ...

    However, once registered, Dask will be able to handle operations with this new
    type:

    >>> da.register_chunk_type(FlaggedArray)
    >>> x = da.ones(5) - FlaggedArray(np.ones(5), True)
    >>> x
    dask.array<sub, shape=(5,), dtype=float64, chunksize=(5,), chunktype=dask.FlaggedArray>
    >>> x.compute()
    Flag: True, Array: array([0., 0., 0., 0., 0.])
    """
    _HANDLED_CHUNK_TYPES.append(type)


try:
    import cupy

    register_chunk_type(cupy.ndarray)
except ImportError:
    pass

try:
    from cupyx.scipy.sparse import spmatrix

    register_chunk_type(spmatrix)
except ImportError:
    pass

try:
    import sparse

    register_chunk_type(sparse.SparseArray)
except ImportError:
    pass

try:
    import scipy.sparse

    register_chunk_type(scipy.sparse.spmatrix)
except ImportError:
    pass

try:
    from scipy.sparse import sparray

    register_chunk_type(sparray)
except ImportError:
    pass


def is_valid_chunk_type(type):
    """Check if given type is a valid chunk and downcast array type"""
    try:
        return type in _HANDLED_CHUNK_TYPES or issubclass(
            type, tuple(_HANDLED_CHUNK_TYPES)
        )
    except TypeError:
        return False


def is_valid_array_chunk(array):
    """Check if given array is of a valid type to operate with"""
    return array is None or isinstance(array, tuple(_HANDLED_CHUNK_TYPES))
