#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################

from __future__ import annotations

import builtins
import inspect
import math
import tempfile
from abc import abstractmethod
from collections import OrderedDict, namedtuple
from functools import reduce
from itertools import product
from typing import TYPE_CHECKING, Any, NamedTuple, Protocol, runtime_checkable

from numpy.exceptions import ComplexWarning

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence

from dataclasses import asdict

import ndindex
import numpy as np

import blosc2
from blosc2 import SpecialValue, blosc2_ext, compute_chunks_blocks
from blosc2.info import InfoReporter
from blosc2.schunk import SChunk

from .linalg import matmul
from .utils import (
    _get_local_slice,
    _get_selection,
    get_chunks_idx,
    npbinvert,
    nplshift,
    nprshift,
    process_key,
    slice_to_chunktuple,
)

# These functions in ufunc_map in ufunc_map_1param are implemented in numexpr and so we call
# those instead (since numexpr uses multithreading it is faster)
ufunc_map = {
    np.add: "+",
    np.subtract: "-",
    np.multiply: "*",
    np.divide: "/",
    np.true_divide: "/",
    np.floor_divide: "//",
    np.power: "**",
    np.less: "<",
    np.less_equal: "<=",
    np.greater: ">",
    np.greater_equal: ">=",
    np.equal: "==",
    np.not_equal: "!=",
    np.bitwise_and: "&",
    np.bitwise_or: "|",
    np.bitwise_xor: "^",
    np.arctan2: "arctan2",
    nplshift: "<<",  # nplshift selected above according to numpy version
    nprshift: ">>",  # nprshift selected above according to numpy version
    np.remainder: "%",
    np.nextafter: "nextafter",
    np.copysign: "copysign",
    np.hypot: "hypot",
    np.maximum: "maximum",
    np.minimum: "minimum",
}

# implemented in numexpr
ufunc_map_1param = {
    np.sqrt: "sqrt",
    np.sin: "sin",
    np.cos: "cos",
    np.tan: "tan",
    np.arcsin: "arcsin",
    np.arccos: "arccos",
    np.arctan: "arctan",
    np.sinh: "sinh",
    np.cosh: "cosh",
    np.tanh: "tanh",
    np.arcsinh: "arcsinh",
    np.arccosh: "arccosh",
    np.arctanh: "arctanh",
    np.exp: "exp",
    np.expm1: "expm1",
    np.log: "log",
    np.log10: "log10",
    np.log1p: "log1p",
    np.log2: "log2",
    np.abs: "abs",
    np.conj: "conj",
    np.real: "real",
    np.imag: "imag",
    npbinvert: "~",  # npbinvert selected above according to numpy version
    np.isnan: "isnan",
    np.isfinite: "isfinite",
    np.isinf: "isinf",
    np.floor: "floor",
    np.ceil: "ceil",
    np.trunc: "trunc",
    np.signbit: "signbit",
    np.round: "round",
}


@runtime_checkable
class Array(Protocol):
    """
    A typing protocol for array-like objects with basic array interface.

    This protocol describes the basic interface required by blosc2 arrays.
    It is implemented by blosc2 classes (:ref:`NDArray`, :ref:`NDField`,
    :ref:`LazyArray`, :ref:`C2Array`, :ref:`ProxyNDSource`...)
    and is compatible with NumPy arrays and other array-like containers
    (e.g., PyTorch, TensorFlow, Dask, Zarr, ...).
    """

    @property
    def dtype(self) -> Any:
        """The data type of the array."""
        ...

    @property
    def shape(self) -> tuple[int, ...]:
        """The shape of the array."""
        ...

    def __len__(self) -> int:
        """The length of the array."""
        ...

    def __getitem__(self, key: Any) -> Any:
        """Get items from the array."""
        ...


def is_documented_by(original):
    def wrapper(target):
        target.__doc__ = original.__doc__
        return target

    return wrapper


def is_inside_new_expr() -> bool:
    """
    Whether the current code is being executed during the creation of new expression.
    """
    # Get the current call stack
    stack = inspect.stack()
    return builtins.any(frame_info.function in {"_new_expr", "_open_lazyarray"} for frame_info in stack)


def make_key_hashable(key):
    if isinstance(key, slice):
        return (key.start, key.stop, key.step)
    elif isinstance(key, tuple | list):
        return tuple(make_key_hashable(k) for k in key)
    elif isinstance(key, np.ndarray):
        return tuple(key.tolist())
    else:
        return key


def get_ndarray_start_stop(ndim, key, shape):
    # key should be Nones and slices
    none_mask, start, stop, step = [], [], [], []
    for i, s in enumerate(key):
        none_mask.append(s is None)
        if s is not None:
            start.append(s.start if s.start is not None else 0)
            stop.append(s.stop if s.stop is not None else shape[i - np.sum(none_mask)])
            step.append(s.step if s.step is not None else 1)
    # Check that start and stop values do not exceed the shape
    for i in range(ndim):
        if start[i] < 0:
            start[i] = shape[i] + start[i]
        if start[i] > shape[i]:
            start[i] = shape[i]
        if stop[i] < 0:
            stop[i] = shape[i] + stop[i]
        if stop[i] > shape[i]:
            stop[i] = shape[i]

    return start, stop, tuple(step), none_mask


def are_partitions_aligned(shape, chunks, blocks):
    """
    Check if the partitions defined by chunks and blocks are aligned with the shape.

    This function verifies that the shape is aligned with the chunks and the chunks are aligned
    with the blocks.

    Returns
    -------
    bool
        True if the partitions are aligned, False otherwise.
    """
    # Check alignment
    alignment_shape_chunks = builtins.all(s % c == 0 for s, c in zip(shape, chunks, strict=True))
    if not alignment_shape_chunks:
        return False
    return builtins.all(c % b == 0 for c, b in zip(chunks, blocks, strict=True))


def are_partitions_behaved(shape, chunks, blocks):
    """
    Check if the partitions defined by chunks and blocks are well-behaved with respect to the shape.

    This function verifies that partitions are C-contiguous with respect the outer container.

    Returns
    -------
    bool
        True if the partitions are well-behaved, False otherwise.
    """

    # Check C-contiguity among partitions
    def check_contiguity(shape, part):
        ndims = len(shape)
        inner_dim = ndims - 1
        for i, size, unit in zip(reversed(range(ndims)), reversed(shape), reversed(part), strict=True):
            if size > unit:
                if i < inner_dim:
                    if size % unit != 0:
                        return False
                else:
                    if size != unit:
                        return False
                inner_dim = i
        return True

    # Check C-contiguity for blocks inside chunks
    if not check_contiguity(chunks, blocks):
        return False

    # Check C-contiguity for chunks inside shape
    return check_contiguity(shape, chunks)


def get_flat_slices_orig(shape: tuple[int], s: tuple[slice, ...]) -> list[slice]:
    """
    From array with `shape`, get the flattened list of slices corresponding to `s`.

    Parameters
    ----------
    shape: tuple[int]
        The shape of the array.
    s: tuple[slice]
        The slice we want to flatten.

    Returns
    -------
    list[slice]
        A list of slices that correspond to the slice `s`.
    """
    # Note: this has been rewritten to use cython, see get_flat_slices
    # It is kept here for reference
    #
    # Process the slice s to get start and stop indices
    key = np.index_exp[s]
    start = [k.start if k.start is not None else 0 for k in key]
    # For stop, cap the values to the shape (shape may not be an exact multiple of the chunks)
    stop = [builtins.min(k.stop if k.stop is not None else shape[i], shape[i]) for i, k in enumerate(key)]

    # Calculate the strides for each dimension
    strides = np.cumprod((1,) + shape[::-1][:-1])[::-1]

    # Generate the 1-dimensional slices
    slices = []
    current_slice_start = None
    current_slice_end = None
    for idx in np.ndindex(*[stop[i] - start[i] for i in range(len(shape))]):
        flat_idx = builtins.sum((start[i] + idx[i]) * strides[i] for i in range(len(shape)))
        if current_slice_start is None:
            current_slice_start = flat_idx
            current_slice_end = flat_idx
        elif flat_idx == current_slice_end + 1:
            current_slice_end = flat_idx
        else:
            slices.append(slice(current_slice_start, current_slice_end + 1))
            current_slice_start = flat_idx
            current_slice_end = flat_idx

    if current_slice_start is not None:
        slices.append(slice(current_slice_start, current_slice_end + 1))

    return slices


def get_flat_slices(
    shape: tuple[int],
    s: tuple[slice, ...],
    c_order: bool = True,
) -> list[slice]:
    """
    From array with `shape`, get the flattened list of slices corresponding to `s`.

    Parameters
    ----------
    shape: tuple
        The shape of the array.
    s: tuple
        The slice we want to flatten.
    c_order: bool
        Whether to flatten the slices in C order (row-major) or just plain order.
        Default is C order.

    Returns
    -------
    list
        A list of slices that correspond to the slice `s`.
    """
    ndim = len(shape)
    if ndim == 0:
        # this will likely cause failure since expected output is tuple of slices
        # however, the list conversion in the last line causes the process to be killed for some reason if shape = ()
        return ()
    start = [s[i].start if s[i].start is not None else 0 for i in range(ndim)]
    stop = [builtins.min(s[i].stop if s[i].stop is not None else shape[i], shape[i]) for i in range(ndim)]
    # Steps are not used in the computation, so raise an error if they are not None or 1
    if builtins.any(s[i].step not in (None, 1) for i in range(ndim)):
        raise ValueError("steps are not supported in slices")

    # Calculate the strides for each dimension
    # Both methods are equivalent
    # strides = np.cumprod((1,) + shape[::-1][:-1])[::-1]
    strides = [reduce(lambda x, y: x * y, shape[i + 1 :], 1) for i in range(ndim)]

    # Convert lists to numpy arrays
    start = np.array(start, dtype=np.int64)
    stop = np.array(stop, dtype=np.int64)
    strides = np.array(strides, dtype=np.int64)

    if not c_order:
        # Generate just a single 1-dimensional slice
        flat_start = np.sum(start * strides)
        # Compute the size of the slice
        flat_size = math.prod(stop - start)
        return [slice(flat_start, flat_start + flat_size)]

    # Generate and return the 1-dimensional slices in C order
    return list(blosc2_ext.slice_flatter(start, stop, strides))


def reshape(
    src: blosc2.Array,
    shape: tuple | list,
    c_order: bool = True,
    **kwargs: Any,
) -> NDArray:
    """Returns an array containing the same data with a new shape.

    This only works when src.shape is 1-dimensional. Multidim case for src is
    interesting, but not supported yet.

    Parameters
    ----------
    src: :ref:`NDArray` or :ref:`NDField` or :ref:`LazyArray` or :ref:`C2Array`
        The input array.
    shape : tuple or list
        The new shape of the array. It should have the same number of elements
        as the current shape.
    c_order: bool
        Whether to reshape the array in C order (row-major) or insertion order.
        Insertion order means that values will be stored in the array
        following the order of chunks in the source array.
        Default is C order.
    kwargs : dict, optional
        Additional keyword arguments supported by the :func:`empty` constructor.

    Returns
    -------
    out: :ref:`NDArray`
        A new array with the requested shape.

    Examples
    --------
    >>> import blosc2
    >>> import numpy as np
    >>> shape = [23 * 11]
    >>> a = np.arange(np.prod(shape))
    >>> # Create an array
    >>> b = blosc2.asarray(a)
    >>> # Reshape the array
    >>> c = blosc2.reshape(b, (11, 23))
    >>> print(c.shape)
    (11, 23)
    """

    if src.ndim != 1:
        raise ValueError("reshape only works when src.shape is 1-dimensional")
    # Check if the new shape is valid
    if math.prod(shape) != math.prod(src.shape):
        raise ValueError("total size of new array must be unchanged")

    # Create the new array
    dst = empty(shape, dtype=src.dtype, **kwargs)

    if is_inside_new_expr() or 0 in shape:
        # We already have the dtype and shape, so return immediately
        return dst

    if shape == ():  # get_flat_slices fails for this case so just return directly
        dst[()] = src[()] if src.shape == () else src[0]
        return dst

    # Copy the data chunk by chunk
    for dst_chunk in dst.iterchunks_info():
        dst_slice = tuple(
            slice(c * s, (c + 1) * s) for c, s in zip(dst_chunk.coords, dst.chunks, strict=False)
        )
        # Cap the stop indices in dst_slices to the dst.shape, and create a new list of slices
        dst_slice = tuple(
            slice(s.start, builtins.min(s.stop, sh)) for s, sh in zip(dst_slice, dst.shape, strict=False)
        )
        size_dst_slice = math.prod([s.stop - s.start for s in dst_slice])
        # Find the series of slices in source array that correspond to the destination chunk
        # (assuming the source array is 1-dimensional here)
        # t0 = time()
        # src_slices = get_flat_slices_orig(dst.shape, dst_slice)
        # Use the get_flat_slices which uses a much faster iterator in cython
        src_slices = get_flat_slices(dst.shape, dst_slice, c_order)
        # print(f"Time to get slices: {time() - t0:.3f} s")
        # Compute the size for slices in the source array
        size_src_slices = builtins.sum(s.stop - s.start for s in src_slices)
        if size_src_slices != size_dst_slice:
            raise ValueError("source slice size is not equal to the destination chunk size")
        # Now, assemble the slices for assignment in the destination array
        dst_buf = np.empty(size_dst_slice, dtype=src.dtype)
        dst_buf_len = 0
        for src_slice in src_slices:
            slice_size = src_slice.stop - src_slice.start
            dst_buf_slice = slice(dst_buf_len, dst_buf_len + slice_size)
            dst_buf_len += slice_size
            if hasattr(src, "res_getitem"):
                # Fast path for lazy UDFs (important for e.g. arange or linspace)
                # This essentially avoids the need to create a new,
                # potentially large NumPy array in memory.
                # This is not critical for Linux, but it is for Windows/Mac.
                dst_buf[dst_buf_slice] = src.res_getitem[src_slice]
            else:
                dst_buf[dst_buf_slice] = src[src_slice]
        # Compute the shape of dst_slice
        dst_slice_shape = tuple(s.stop - s.start for s in dst_slice)
        # ... and assign the buffer to the destination array
        dst[dst_slice] = dst_buf.reshape(dst_slice_shape)

    return dst


def _check_allowed_dtypes(
    value: bool | int | float | str | blosc2.Array,
):
    def _is_array_like(v: Any) -> bool:
        try:
            # Try Protocol runtime check first (works when possible)
            if isinstance(v, blosc2.Array):
                return True
        except Exception:
            # Some runtime contexts may raise (or return False) — fall back to duck typing
            pass
        # Structural fallback: common minimal array interface
        return hasattr(v, "shape") and hasattr(v, "dtype") and callable(getattr(v, "__getitem__", None))

    if not (_is_array_like(value) or np.isscalar(value)):
        raise RuntimeError(
            f"Expected blosc2.Array or scalar instances and you provided a '{type(value)}' instance"
        )


def sum(
    ndarr: blosc2.Array,
    axis: int | tuple[int] | None = None,
    dtype: np.dtype | str = None,
    keepdims: bool = False,
    **kwargs: Any,
) -> blosc2.Array | int | float | complex | bool:
    """
    Return the sum of array elements over a given axis.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array or expression.
    axis: int or tuple of ints, optional
        Axis or axes along which a sum is performed. By default, axis=None,
        sums all the elements of the input array. If axis is negative,
        it counts from the last to the first axis.
    dtype: np.dtype or list str, optional
        The type of the returned array and of the accumulator in which the
        elements are summed. The dtype of :paramref:`ndarr` is used by default unless it has
        an integer dtype of less precision than the default platform integer.
    keepdims: bool, optional
        If set to True, the reduced axes are left in the result
        as dimensions with size one. With this option, the result will broadcast
        correctly against the input array.
    fp_accuracy: :ref:`blosc2.FPAccuracy`, optional
        Specifies the floating-point accuracy for reductions on :ref:`LazyExpr`.
        Passed to :func:`LazyExpr.compute` when :paramref:`ndarr` is a LazyExpr.
    kwargs: dict, optional
        Additional keyword arguments supported by the :func:`empty` constructor.

    Returns
    -------
    sum_along_axis: np.ndarray or :ref:`NDArray` or scalar
        The sum of the elements along the axis.

    References
    ----------
    `np.sum <https://numpy.org/doc/stable/reference/generated/numpy.sum.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> # Example array
    >>> array = np.array([[1, 2, 3], [4, 5, 6]])
    >>> nd_array = blosc2.asarray(array)
    >>> # Sum all elements in the array (axis=None)
    >>> total_sum = blosc2.sum(nd_array)
    >>> print("Sum of all elements:", total_sum)
    21
    >>> # Sum along axis 0 (columns)
    >>> sum_axis_0 = blosc2.sum(nd_array, axis=0)
    >>> print("Sum along axis 0 (columns):", sum_axis_0)
    Sum along axis 0 (columns): [5 7 9]
    """
    return ndarr.sum(axis=axis, dtype=dtype, keepdims=keepdims, **kwargs)


def mean(
    ndarr: blosc2.Array,
    axis: int | tuple[int] | None = None,
    dtype: np.dtype | str = None,
    keepdims: bool = False,
    **kwargs: Any,
) -> blosc2.Array | int | float | complex | bool:
    """
    Return the arithmetic mean along the specified axis.

    The parameters are documented in the :func:`sum <blosc2.sum>`.

    Returns
    -------
    mean_along_axis: np.ndarray or :ref:`NDArray` or scalar
        The mean of the elements along the axis.

    References
    ----------
    `np.mean <https://numpy.org/doc/stable/reference/generated/numpy.mean.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> # Example array
    >>> array = np.array([[1, 2, 3], [4, 5, 6]]
    >>> nd_array = blosc2.asarray(array)
    >>> # Compute the mean of all elements in the array (axis=None)
    >>> overall_mean = blosc2.mean(nd_array)
    >>> print("Mean of all elements:", overall_mean)
    Mean of all elements: 3.5
    """
    return ndarr.mean(axis=axis, dtype=dtype, keepdims=keepdims, **kwargs)


def std(
    ndarr: blosc2.Array,
    axis: int | tuple[int] | None = None,
    dtype: np.dtype | str = None,
    ddof: int = 0,
    keepdims: bool = False,
    **kwargs: Any,
) -> blosc2.Array | int | float | bool:
    """
    Return the standard deviation along the specified axis.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array or expression.
    axis: int or tuple of ints, optional
        Axis or axes along which the standard deviation is computed. By default, `axis=None`
        computes the standard deviation of the flattened array.
    dtype: np.dtype or list str, optional
        Type to use in computing the standard deviation. For integer inputs, the
        default is float32; for floating point inputs, it is the same as the input dtype.
    ddof: int, optional
        Means Delta Degrees of Freedom. The divisor used in calculations is N - ddof,
        where N represents the number of elements. By default, ddof is zero.
    keepdims: bool, optional
        If set to True, the reduced axes are left in the result as
        dimensions with size one. This ensures that the result will broadcast correctly
        against the input array.
    fp_accuracy: :ref:`blosc2.FPAccuracy`, optional
        Specifies the floating-point accuracy for reductions on :ref:`LazyExpr`.
        Passed to :func:`LazyExpr.compute` when :paramref:`ndarr` is a LazyExpr.
    kwargs: dict, optional
        Additional keyword arguments that are supported by the :func:`empty` constructor.

    Returns
    -------
    std_along_axis: np.ndarray or :ref:`NDArray` or scalar
        The standard deviation of the elements along the axis.

    References
    ----------
    `np.std <https://numpy.org/doc/stable/reference/generated/numpy.std.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> # Create an instance of NDArray with some data
    >>> array = np.array([[1, 2, 3], [4, 5, 6]])
    >>> nd_array = blosc2.asarray(array)
    >>> # Compute the standard deviation of the entire array
    >>> std_all = blosc2.std(nd_array)
    >>> print("Standard deviation of the entire array:", std_all)
    Standard deviation of the entire array: 1.707825127659933
    >>> # Compute the standard deviation along axis 0 (columns)
    >>> std_axis0 = blosc2.std(nd_array, axis=0)
    >>> print("Standard deviation along axis 0:", std_axis0)
    Standard deviation along axis 0: [1.5 1.5 1.5]
    """
    return ndarr.std(axis=axis, dtype=dtype, ddof=ddof, keepdims=keepdims, **kwargs)


def var(
    ndarr: blosc2.Array,
    axis: int | tuple[int] | None = None,
    dtype: np.dtype | str = None,
    ddof: int = 0,
    keepdims: bool = False,
    **kwargs: Any,
) -> blosc2.Array | int | float | bool:
    """
    Return the variance along the specified axis.

    The parameters are documented in the :func:`std <blosc2.std>`.

    Returns
    -------
    var_along_axis: np.ndarray or :ref:`NDArray` or scalar
        The variance of the elements along the axis.

    References
    ----------
    `np.var <https://numpy.org/doc/stable/reference/generated/numpy.var.html>`_


    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> # Create an instance of NDArray with some data
    >>> array = np.array([[1, 2, 3], [4, 5, 6]])
    >>> nd_array = blosc2.asarray(array)
    >>> # Compute the variance of the entire array
    >>> var_all = blosc2.var(nd_array)
    >>> print("Variance of the entire array:", var_all)
    Variance of the entire array: 2.9166666666666665
    >>> # Compute the variance along axis 0 (columns)
    >>> var_axis0 = blosc2.var(nd_array, axis=0)
    >>> print("Variance along axis 0:", var_axis0)
    Variance along axis 0: [2.25 2.25 2.25]
    """
    return ndarr.var(axis=axis, dtype=dtype, ddof=ddof, keepdims=keepdims, **kwargs)


def prod(
    ndarr: blosc2.Array,
    axis: int | tuple[int] | None = None,
    dtype: np.dtype | str = None,
    keepdims: bool = False,
    **kwargs: Any,
) -> blosc2.Array | int | float | complex | bool:
    """
    Return the product of array elements over a given axis.

    The parameters are documented in the :func:`sum <blosc2.sum>`.

    Returns
    -------
    product_along_axis: np.ndarray or :ref:`NDArray` or scalar
        The product of the elements along the axis.

    References
    ----------
    `np.prod <https://numpy.org/doc/stable/reference/generated/numpy.prod.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> # Create an instance of NDArray with some data
    >>> array = np.array([[11, 22, 33], [4, 15, 36]])
    >>> nd_array = blosc2.asarray(array)
    >>> # Compute the product of all elements in the array
    >>> prod_all = blosc2.prod(nd_array)
    >>> print("Product of all elements in the array:", prod_all)
    Product of all elements in the array: 17249760
    >>> # Compute the product along axis 1 (rows)
    >>> prod_axis1 = blosc2.prod(nd_array, axis=1)
    >>> print("Product along axis 1:", prod_axis1)
    Product along axis 1: [7986 2160]
    """
    return ndarr.prod(axis=axis, dtype=dtype, keepdims=keepdims, **kwargs)


def min(
    ndarr: blosc2.Array,
    axis: int | tuple[int] | None = None,
    keepdims: bool = False,
    **kwargs: Any,
) -> blosc2.Array | int | float | complex | bool:
    """
    Return the minimum along a given axis.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array or expression.
    axis: int or tuple of ints, optional
        Axis or axes along which to operate. By default, flattened input is used.
    keepdims: bool, optional
        If set to True, the axes which are reduced are left in the result as
        dimensions with size one. With this option, the result will broadcast correctly
        against the input array.
    fp_accuracy: :ref:`blosc2.FPAccuracy`, optional
        Specifies the floating-point accuracy for reductions on :ref:`LazyExpr`.
        Passed to :func:`LazyExpr.compute` when :paramref:`ndarr` is a LazyExpr.
    kwargs: dict, optional
        Keyword arguments that are supported by the :func:`empty` constructor.

    Returns
    -------
    min_along_axis: np.ndarray or :ref:`NDArray` or scalar
        The minimum of the elements along the axis.

    References
    ----------
    `np.min <https://numpy.org/doc/stable/reference/generated/numpy.min.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> array = np.array([1, 3, 7, 8, 9, 31])
    >>> nd_array = blosc2.asarray(array)
    >>> min_all = blosc2.min(nd_array)
    >>> print("Minimum of all elements in the array:", min_all)
    Minimum of all elements in the array: 1
    >>> # Compute the minimum along axis 0 with keepdims=True
    >>> min_keepdims = blosc2.min(nd_array, axis=0, keepdims=True)
    >>> print("Minimum along axis 0 with keepdims=True:", min_keepdims)
    Minimum along axis 0 with keepdims=True:  [1]
    """
    return ndarr.min(axis=axis, keepdims=keepdims, **kwargs)


def max(
    ndarr: blosc2.Array,
    axis: int | tuple[int] | None = None,
    keepdims: bool = False,
    **kwargs: Any,
) -> blosc2.Array | int | float | complex | bool:
    """
    Return the maximum along a given axis.

    The parameters are documented in the :func:`min <blosc2.min>`.

    Returns
    -------
    max_along_axis: np.ndarray or :ref:`NDArray` or scalar
        The maximum of the elements along the axis.

    References
    ----------
    `np.max <https://numpy.org/doc/stable/reference/generated/numpy.max.html>`_

    Examples
    --------
    >>> import blosc2
    >>> import numpy as np
    >>> data = np.array([[11, 2, 36, 24, 5, 69], [73, 81, 49, 6, 73, 0]])
    >>> ndarray = blosc2.asarray(data)
    >>> print("NDArray data:", ndarray[:])
    NDArray data:  [[11  2 36 24  5 69]
                    [73 81 49  6 73  0]]
    >>> # Compute the maximum along axis 0 and 1
    >>> max_along_axis_0 = blosc2.max(ndarray, axis=0)
    >>> print("Maximum along axis 0:", max_along_axis_0)
    Maximum along axis 0: [73 81 49 24 73 69]
    >>> max_along_axis_1 = blosc2.max(ndarray, axis=1)
    >>> print("Maximum along axis 1:", max_along_axis_1)
    Maximum along axis 1: [69 81]
    >>> max_flattened = blosc2.max(ndarray)
    >>> print("Maximum of the flattened array:", max_flattened)
    Maximum of the flattened array: 81
    """
    return ndarr.max(axis=axis, keepdims=keepdims, **kwargs)


def any(
    ndarr: blosc2.Array,
    axis: int | tuple[int] | None = None,
    keepdims: bool = False,
    **kwargs: Any,
) -> blosc2.Array | bool:
    """
    Test whether any array element along a given axis evaluates to True.

    The parameters are documented in the :func:`min <blosc2.min>`.

    Returns
    -------
    any_along_axis: np.ndarray or :ref:`NDArray` or scalar
        The result of the evaluation along the axis.

    References
    ----------
    `np.any <https://numpy.org/doc/stable/reference/generated/numpy.any.html>`_

    Examples
    --------
    >>> import blosc2
    >>> import numpy as np
    >>> data = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0]])
    >>> # Convert the NumPy array to a Blosc2 NDArray
    >>> ndarray = blosc2.asarray(data)
    >>> print("NDArray data:", ndarray[:])
    NDArray data: [[1 0 0]
                    [0 1 0]
                    [0 0 0]]
    >>> any_along_axis_0 = blosc2.any(ndarray, axis=0)
    >>> print("Any along axis 0:", any_along_axis_0)
    Any along axis 0: [True True False]
    >>> any_flattened = blosc2.any(ndarray)
    >>> print("Any in the flattened array:", any_flattened)
    Any in the flattened array: True
    """
    return ndarr.any(axis=axis, keepdims=keepdims, **kwargs)


def argmin(
    ndarr: blosc2.Array, axis: int | None = None, keepdims: bool = False, **kwargs
) -> blosc2.Array | int:
    """
    Returns the indices of the minimum values along a specified axis.

    When the minimum value occurs multiple times, only the indices corresponding to the first occurrence are returned.

    Parameters
    ----------
    x: blosc2.Array
        Input array. Should have a real-valued data type.

    axis: int | None
        Axis along which to search. If None, return index of the minimum value of flattened array. Default: None.

    keepdims: bool
        If True, reduced axis included in the result as singleton dimension. Otherwise, axis not included in the result. Default: False.
    fp_accuracy: :ref:`blosc2.FPAccuracy`, optional
        Specifies the floating-point accuracy for reductions on :ref:`LazyExpr`.
        Passed to :func:`LazyExpr.compute` when :paramref:`ndarr` is a LazyExpr.

    Returns
    -------
    out: blosc2.Array
        If axis is None, a zero-dimensional array containing the index of the first occurrence of the minimum value; otherwise, a non-zero-dimensional array containing the indices of the minimum values.
    """
    return ndarr.argmin(axis=axis, keepdims=keepdims, **kwargs)


def argmax(
    ndarr: blosc2.Array, axis: int | None = None, keepdims: bool = False, **kwargs
) -> blosc2.Array | int:
    """
    Returns the indices of the maximum values along a specified axis.

    When the maximum value occurs multiple times, only the indices corresponding to the first occurrence are returned.

    Parameters
    ----------
    x: blosc2.Array
        Input array. Should have a real-valued data type.

    axis: int | None
        Axis along which to search. If None, return index of the maximum value of flattened array. Default: None.

    keepdims: bool
        If True, reduced axis included in the result as singleton dimension. Otherwise, axis not included in the result. Default: False.
    fp_accuracy: :ref:`blosc2.FPAccuracy`, optional
        Specifies the floating-point accuracy for reductions on :ref:`LazyExpr`.
        Passed to :func:`LazyExpr.compute` when :paramref:`ndarr` is a LazyExpr.

    Returns
    -------
    out: blosc2.Array
        If axis is None, a zero-dimensional array containing the index of the first occurrence of the maximum value; otherwise, a non-zero-dimensional array containing the indices of the maximum values.
    """
    return ndarr.argmax(axis=axis, keepdims=keepdims, **kwargs)


def all(
    ndarr: blosc2.Array,
    axis: int | tuple[int] | None = None,
    keepdims: bool = False,
    **kwargs: Any,
) -> blosc2.Array | bool:
    """
    Test whether all array elements along a given axis evaluate to True.

    The parameters are documented in the :func:`min <blosc2.min>`.

    Returns
    -------
    all_along_axis: np.ndarray or :ref:`NDArray` or scalar
        The result of the evaluation along the axis.

    References
    ----------
    `np.all <https://numpy.org/doc/stable/reference/generated/numpy.all.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> data = np.array([True, True, False, True, True, True])
    >>> ndarray = blosc2.asarray(data)
    >>> # Test if all elements are True along the default axis (flattened array)
    >>> result_flat = blosc2.all(ndarray)
    >>> print("All elements are True (flattened):", result_flat)
    All elements are True (flattened): False
    """
    return ndarr.all(axis=axis, keepdims=keepdims, **kwargs)


def sin(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Compute the trigonometric sine, element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array containing angles in radians.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing the sine of the input angles. The result can be evaluated.

    References
    ----------
    `np.sin <https://numpy.org/doc/stable/reference/generated/numpy.sin.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> angles = np.array([0, np.pi/6, np.pi/4, np.pi/2, np.pi])
    >>> nd_array = blosc2.asarray(angles)
    >>> result_ = blosc2.sin(nd_array)
    >>> result = result_[:]
    >>> print("Angles in radians:", angles)
    Angles in radians: [0.         0.52359878 0.78539816 1.57079633 3.14159265]
    >>> print("Sine of the angles:", result)
    Sine of the angles: [0.00000000e+00 5.00000000e-01 7.07106781e-01 1.00000000e+00
    1.22464680e-16]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "sin", None))


def cos(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Trigonometric cosine, element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array containing angles in radians.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing the cosine of the input angles. The result can be evaluated.

    References
    ----------
    `np.cos <https://numpy.org/doc/stable/reference/generated/numpy.cos.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> angles = np.array([0, np.pi/6, np.pi/4, np.pi/2, np.pi])
    >>> nd_array = blosc2.asarray(angles)
    >>> result_ = blosc2.cos(nd_array)
    >>> result = result_[:]
    >>> print("Angles in radians:", angles)
    Angles in radians: [0.         0.52359878 0.78539816 1.57079633 3.14159265]
    >>> print("Cosine of the angles:", result)
    Cosine of the angles: [ 1.00000000e+00  8.66025404e-01  7.07106781e-01  6.12323400e-17
    -1.00000000e+00]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "cos", None))


def tan(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Compute the trigonometric tangent, element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array containing angles in radians.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing the tangent of the input angles.
        The result can be evaluated.

    References
    ----------
    `np.tan <https://numpy.org/doc/stable/reference/generated/numpy.tan.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> angles = np.array([0, np.pi/6, np.pi/4, np.pi/2, np.pi])
    >>> nd_array = blosc2.asarray(angles)
    >>> result_ = blosc2.tan(nd_array)
    >>> result = result_[:]
    >>> print("Angles in radians:", angles)
    Angles in radians: [0.         0.52359878 0.78539816 1.57079633 3.14159265]
    >>> print("Tangent of the angles:", result)
    Tangent of the angles: [ 0.00000000e+00  5.77350269e-01  1.00000000e+00  1.63312394e+16
    -1.22464680e-16]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "tan", None))


def sqrt(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Return the non-negative square-root of an array, element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing the square root of the input array.
        The result can be evaluated.

    References
    ----------
    `np.sqrt <https://numpy.org/doc/stable/reference/generated/numpy.sqrt.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> data = np.array([0, np.pi/6, np.pi/4, np.pi/2, np.pi])
    >>> nd_array = blosc2.asarray(data)
    >>> result_ = blosc2.sqrt(nd_array)
    >>> result = result_[:]
    >>> print("Original numbers:", data)
    Original numbers: [ 0  1  4  9 16 25]
    >>> print("Square roots:", result)
    Square roots: [0. 1. 2. 3. 4. 5.]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "sqrt", None))


def sinh(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Hyperbolic sine, element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing the hyperbolic sine of the input array.
        The result can be evaluated.

    References
    ----------
    `np.sinh <https://numpy.org/doc/stable/reference/generated/numpy.sinh.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> numbers = np.array([-2, -1, 0, 1, 2])
    >>> ndarray = blosc2.asarray(numbers)
    >>> result_lazy = blosc2.sinh(ndarray)
    >>> result = result_lazy[:]
    >>> print("Original numbers:", numbers)
    Original numbers: [-2 -1  0  1  2]
    >>> print("Hyperbolic sine:", result)
    Hyperbolic sine: [-3.62686041 -1.17520119  0.          1.17520119  3.62686041]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "sinh", None))


def cosh(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Compute the hyperbolic cosine, element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing the hyperbolic cosine of the input array.
        The result can be evaluated.

    References
    ----------
    `np.cosh <https://numpy.org/doc/stable/reference/generated/numpy.cosh.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> numbers = np.array([-2, -1, 0, 1, 2])
    >>> ndarray = blosc2.asarray(numbers)
    >>> result_lazy = blosc2.cosh(ndarray)
    >>> result = result_lazy[:]
    >>> print("Original numbers:", numbers)
    Original numbers: [-2 -1  0  1  2]
    >>> print("Hyperbolic cosine:", result)
    Hyperbolic cosine: [3.76219569 1.54308063 1.         1.54308063 3.76219569]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "cosh", None))


def tanh(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Compute the hyperbolic tangent, element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing the hyperbolic tangent of the input array.
        The result can be evaluated.

    References
    ----------
    `np.tanh <https://numpy.org/doc/stable/reference/generated/numpy.tanh.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> numbers = np.array([-2, -1, 0, 1, 2])
    >>> ndarray = blosc2.asarray(numbers)
    >>> result_lazy = blosc2.tanh(ndarray)
    >>> result = result_lazy[:]
    >>> print("Original numbers:", numbers)
    Original numbers: [-2 -1  0  1  2]
    >>> print("Hyperbolic tangent:", result)
    Hyperbolic tangent: [-0.96402758 -0.76159416  0.          0.76159416  0.96402758]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "tanh", None))


def arcsin(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Compute the inverse sine, element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing the inverse sine of the input array.
        The result can be evaluated.

    References
    ----------
    `np.arcsin <https://numpy.org/doc/stable/reference/generated/numpy.arcsin.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> numbers = np.array([-1, -0.5, 0, 0.5, 1])
    >>> ndarray = blosc2.asarray(numbers)
    >>> result_lazy = blosc2.arcsin(ndarray)
    >>> result = result_lazy[:]
    >>> print("Original numbers:", numbers)
    Original numbers: [-1.  -0.5  0.   0.5  1. ]
    >>> print("Arcsin:", result)
    Arcsin: [-1.57079633 -0.52359878  0.          0.52359878  1.57079633]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "arcsin", None))


asin = arcsin  # alias


def arccos(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Compute the inverse cosine, element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing the inverse cosine of the input array.
        The result can be evaluated.

    References
    ----------
    `np.arccos <https://numpy.org/doc/stable/reference/generated/numpy.arccos.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> numbers = np.array([-1, -0.5, 0, 0.5, 1])
    >>> ndarray = blosc2.asarray(numbers)
    >>> result_lazy = blosc2.arccos(ndarray)
    >>> result = result_lazy[:]
    >>> print("Original numbers:", numbers)
    Original numbers: [-1.  -0.5  0.   0.5  1. ]
    >>> print("Arccos:", result)
    Arccos: [3.14159265 2.0943951  1.57079633 1.04719755 0.        ]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "arccos", None))


acos = arccos  # alias


def arctan(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Compute the inverse tangent, element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing the inverse tangent of the input array.
        The result can be evaluated.

    References
    ----------
    `np.arctan <https://numpy.org/doc/stable/reference/generated/numpy.arctan.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> numbers = np.array([-1, -0.5, 0, 0.5, 1])
    >>> ndarray = blosc2.asarray(numbers)
    >>> result_lazy = blosc2.arctan(ndarray)
    >>> result = result_lazy[:]
    >>> print("Original numbers:", numbers)
    Original numbers: [-1.  -0.5  0.   0.5  1. ]
    >>> print("Arctan:", result)
    Arctan: [-0.78539816 -0.46364761  0.          0.46364761  0.78539816]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "arctan", None))


atan = arctan  # alias


def arctan2(ndarr1: blosc2.Array, ndarr2: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Compute the element-wise arc tangent of ``ndarr1 / ndarr2`` choosing the quadrant correctly.

    Parameters
    ----------
    ndarr1: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array`
        The first input array.
    ndarr2: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array`
        The second input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing the element-wise arc tangent of ``ndarr1 / ndarr2``.
        The result can be evaluated.

    References
    ----------
    `np.arctan2 <https://numpy.org/doc/stable/reference/generated/numpy.arctan2.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> y = np.array([0, 1, 0, -1, 1])
    >>> x = np.array([1, 1, -1, -1, 0])
    >>> ndarray_y = blosc2.asarray(y)
    >>> ndarray_x = blosc2.asarray(x)
    >>> result_lazy = blosc2.arctan2(ndarray_y, ndarray_x)
    >>> result = result_lazy[:]
    >>> print("y:", y)
    y: [ 0  1  0 -1  1]
    >>> print("x:", x)
    x: [ 1  1 -1 -1  0]
    >>> print("Arctan2(y, x):", result)
    Arctan2(y, x): [ 0.          0.78539816  3.14159265 -2.35619449  1.57079633]
    """
    return blosc2.LazyExpr(new_op=(ndarr1, "arctan2", ndarr2))


atan2 = arctan2  # alias


def arcsinh(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Compute the inverse hyperbolic sine, element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing the inverse hyperbolic sine of the input array.
        The result can be evaluated.

    References
    ----------
    `np.arcsinh <https://numpy.org/doc/stable/reference/generated/numpy.arcsinh.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> values = np.array([-2, -1, 0, 1, 2])
    >>> ndarray = blosc2.asarray(values)
    >>> result_lazy = blosc2.arcsinh(ndarray)
    >>> result = result_lazy[:]
    >>> print("Original values:", values)
    Original values: [-2 -1  0  1  2]
    >>> print("Arcsinh:", result)
    Arcsinh: [-1.44363548 -0.88137359  0.          0.88137359  1.44363548]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "arcsinh", None))


asinh = arcsinh  # alias


def arccosh(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Compute the inverse hyperbolic cosine, element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing the inverse hyperbolic cosine of the input array.
        The result can be evaluated.

    References
    ----------
    `np.arccosh <https://numpy.org/doc/stable/reference/generated/numpy.arccosh.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> values = np.array([1, 2, 3, 4, 5])
    >>> ndarray = blosc2.asarray(values)
    >>> result_lazy = blosc2.arccosh(ndarray)
    >>> result = result_lazy[:]
    >>> print("Original values:", values)
    Original values: [1 2 3 4 5]
    >>> print("Arccosh:", result)
    Arccosh: [0.         1.3169579  1.76274717 2.06343707 2.29243167]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "arccosh", None))


acosh = arccosh  # alias


def arctanh(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Compute the inverse hyperbolic tangent, element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing the inverse hyperbolic tangent of the input array.
        The result can be evaluated.

    References
    ----------
    `np.arctanh <https://numpy.org/doc/stable/reference/generated/numpy.arctanh.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> values = np.array([-0.9, -0.5, 0, 0.5, 0.9])
    >>> ndarray = blosc2.asarray(values)
    >>> result_lazy = blosc2.arctanh(ndarray)
    >>> result = result_lazy[:]
    >>> print("Original values:", values)
    Original values: [-0.9 -0.5  0.   0.5  0.9]
    >>> print("Arctanh:", result)
    Arctanh: [-1.47221949 -0.54930614  0.          0.54930614  1.47221949]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "arctanh", None))


atanh = arctanh  # alias


def exp(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Calculate the exponential of all elements in the input array.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing the exponential of the input array.
        The result can be evaluated.

    References
    ----------
    `np.exp <https://numpy.org/doc/stable/reference/generated/numpy.exp.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> values = np.array([0, 1, 2, 3, 4])
    >>> ndarray = blosc2.asarray(values)
    >>> result_lazy = blosc2.exp(ndarray)
    >>> result = result_lazy[:]
    >>> print("Original values:", values)
    Original values: [0 1 2 3 4]
    >>> print("Exponential:", result)
    Exponential: [ 1.          2.71828183  7.3890561  20.08553692 54.59815003]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "exp", None))


def expm1(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Calculate ``exp(ndarr) - 1`` for all elements in the array.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing ``exp(ndarr) - 1`` of the input array.
        The result can be evaluated.

    References
    ----------
    `np.expm1 <https://numpy.org/doc/stable/reference/generated/numpy.expm1.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> values = np.array([-1, -0.5, 0, 0.5, 1])
    >>> ndarray = blosc2.asarray(values)
    >>> result_lazy = blosc2.expm1(ndarray)
    >>> result = result_lazy[:]
    >>> print("Original values:", values)
    Original values: [-1.  -0.5  0.   0.5  1. ]
    >>> print("Expm1:", result)
    Expm1: [-0.63212056 -0.39346934  0.          0.64872127  1.71828183]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "expm1", None))


def log(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Compute the natural logarithm, element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing the natural logarithm of the input array

    References
    ----------
    `np.log <https://numpy.org/doc/stable/reference/generated/numpy.log.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> values = np.array([1, 2, 3, 4, 5])
    >>> ndarray = blosc2.asarray(values)
    >>> result_lazy = blosc2.log(ndarray)
    >>> result = result_lazy[:]
    >>> print("Original values:", values)
    Original values: [1 2 3 4 5]
    >>> print("Logarithm (base e):", result)
    Logarithm (base e): [0.         0.69314718 1.09861229 1.38629436 1.60943791]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "log", None))


def log10(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Return the base 10 logarithm of the input array, element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing the base 10 logarithm of the input array.

    References
    ----------
    `np.log10 <https://numpy.org/doc/stable/reference/generated/numpy.log10.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> values = np.array([1, 10, 100, 1000, 10000])
    >>> ndarray = blosc2.asarray(values)
    >>> result_lazy = blosc2.log10(ndarray)
    >>> result = result_lazy[:]
    >>> print("Original values:", values)
    Original values: [    1    10   100  1000 10000]
    >>> print("Logarithm (base 10):", result)
    Logarithm (base 10): [0. 1. 2. 3. 4.]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "log10", None))


def log1p(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Return the natural logarithm of one plus the input array, element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing the natural logarithm of one plus the input array.

    References
    ----------
    `np.log1p <https://numpy.org/doc/stable/reference/generated/numpy.log1p.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> values = np.array([-0.9, -0.5, 0, 0.5, 0.9])
    >>> ndarray = blosc2.asarray(values)
    >>> result_lazy = blosc2.log1p(ndarray)
    >>> result = result_lazy[:]
    >>> print("Original values:", values)
    Original values: [-0.9 -0.5  0.   0.5  0.9]
    >>> print("Log1p (log(1 + x)):", result)
    Log1p (log(1 + x)): [-2.30258509 -0.69314718  0.          0.40546511  0.64185389]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "log1p", None))


def log2(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Return the base 2 logarithm of the input array, element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing the base 2 logarithm of the input array.

    References
    ----------
    `np.log2 <https://numpy.org/doc/stable/reference/generated/numpy.log2.html>`_

    """
    return blosc2.LazyExpr(new_op=(ndarr, "log2", None))


def conj(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Return the complex conjugate, element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing the complex conjugate of the input array.

    References
    ----------
    `np.conj <https://numpy.org/doc/stable/reference/generated/numpy.conj.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> values = np.array([1+2j, 3-4j, -5+6j, 7-8j])
    >>> ndarray = blosc2.asarray(values)
    >>> result_ = blosc2.conj(ndarray)
    >>> result = result_[:]
    >>> print("Original values:", values)
    Original values: [ 1.+2.j  3.-4.j -5.+6.j  7.-8.j]
    >>> print("Complex conjugates:", result)
    Complex conjugates: [ 1.-2.j  3.+4.j -5.-6.j  7.+8.j]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "conj", None))


def real(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Return the real part of the complex array, element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing the real part of the input array.

    References
    ----------
    `np.real <https://numpy.org/doc/stable/reference/generated/numpy.real.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> complex_values = np.array([1+2j, 3-4j, -5+6j, 7-8j])
    >>> ndarray = blosc2.asarray(complex_values)
    >>> result_ = blosc2.real(ndarray)
    >>> result = result_[:]
    >>> print("Original complex values:", complex_values)
    Original values: [ 1.+2.j  3.-4.j -5.+6.j  7.-8.j]
    >>> print("Real parts:", result)
    Real parts: [ 1.  3. -5.  7.]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "real", None))


def imag(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Return the imaginary part of the complex array, element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression representing the imaginary part of the input array.

    References
    ----------
    `np.imag <https://numpy.org/doc/stable/reference/generated/numpy.imag.html>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> complex_values = np.array([2+3j, -1+4j, 0-2j, 5+6j])
    >>> ndarray = blosc2.asarray(complex_values)
    >>> result_ = blosc2.imag(ndarray)
    >>> result = result_[:]
    >>> print("Original complex values:", complex_values)
    Original complex values: [ 2.+3.j -1.+4.j  0.-2.j  5.+6.j]
    >>> print("Imaginary parts:", result)
    Imaginary parts: [ 3.  4. -2.  6.]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "imag", None))


def contains(ndarr: blosc2.Array, value: str | bytes | blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Check if the array contains a specified value.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array`
        The input array.
    value: str or bytes or :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array`
        The value to be checked.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression that can be evaluated to check if the value
        is contained in the array.

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> values = np.array([b"apple", b"xxbananaxxx", b"cherry", b"date"])
    >>> text_values = blosc2.asarray(values)
    >>> value_to_check = b"banana"
    >>> expr = blosc2.contains(text_values, value_to_check)
    >>> result = expr.compute()
    >>> print("Contains 'banana':", result[:])
    Contains 'banana': [False  True False False]
    """
    if not isinstance(value, str | bytes | NDArray):
        raise TypeError("value should be a string, bytes or a NDArray!")
    return blosc2.LazyExpr(new_op=(ndarr, "contains", value))


def abs(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Calculate the absolute value element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression that can be evaluated to get the absolute values.

    References
    ----------
    `np.abs <https://numpy.org/doc/stable/reference/generated/numpy.absolute.html#numpy.absolute>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> values = np.array([-5, -3, 0, 2, 4])
    >>> ndarray = blosc2.asarray(values)
    >>> result_ = blosc2.abs(ndarray)
    >>> result = result_[:]
    >>> print("Original values:", values)
    Original values: [-5 -3  0  2  4]
    >>> print("Absolute values:", result)
    Absolute values: [5. 3. 0. 2. 4.]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "abs", None))


def isnan(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Return True/False for not-a-number values element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression that can be evaluated to get the True/False array of results.

    References
    ----------
    `np.isnan <https://numpy.org/doc/stable/reference/generated/numpy.isnan.html#numpy.isnan>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> values = np.array([-5, -3, np.nan, 2, 4])
    >>> ndarray = blosc2.asarray(values)
    >>> result_ = blosc2.isnan(ndarray)
    >>> result = result_[:]
    >>> print("isnan:", result)
    isnan: [False, False, True, False, False]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "isnan", None))


def isfinite(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Return True/False for finite values element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression that can be evaluated to get the True/False array of results.

    References
    ----------
    `np.isfinite <https://numpy.org/doc/stable/reference/generated/numpy.isfinite.html#numpy.isfinite>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> values = np.array([-5, -3, np.inf, 2, 4])
    >>> ndarray = blosc2.asarray(values)
    >>> result_ = blosc2.isfinite(ndarray)
    >>> result = result_[:]
    >>> print("isfinite:", result)
    isfinite: [True, True, False, True, True]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "isfinite", None))


def isinf(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
    """
    Return True/False for infinite values element-wise.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    Returns
    -------
    out: :ref:`LazyExpr`
        A lazy expression that can be evaluated to get the True/False array of results.

    References
    ----------
    `np.isinf <https://numpy.org/doc/stable/reference/generated/numpy.isinf.html#numpy.isinf>`_

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> values = np.array([-5, -3, np.inf, 2, 4])
    >>> ndarray = blosc2.asarray(values)
    >>> result_ = blosc2.isinf(ndarray)
    >>> result = result_[:]
    >>> print("isinf:", result)
    isinf: [False, False, True, False, False]
    """
    return blosc2.LazyExpr(new_op=(ndarr, "isinf", None))


# def nonzero(ndarr: blosc2.Array, /) -> blosc2.LazyExpr:
#     """
#     Return indices of nonzero values.

#     Parameters
#     ----------
#     ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
#         The input array.

#     Returns
#     -------
#     out: :ref:`LazyExpr`
#         A lazy expression that can be evaluated to get the array of results.

#     References
#     ----------
#     `np.nonzero <https://numpy.org/doc/stable/reference/generated/numpy.nonzero.html#numpy.nonzero>`_
#     """
#     # FIXME: This is not correct
#     return ndarr.__ne__(0)


def count_nonzero(ndarr: blosc2.Array, axis: int | Sequence[int] | None = None) -> int:
    """
    Return number of nonzero values along axes.

    Parameters
    ----------
    ndarr: :ref:`NDArray` or :ref:`NDField` or :ref:`C2Array` or :ref:`LazyExpr`
        The input array.

    axis: int | Sequence[int] | None
        Axes along which to count nonzero entries. If None, sum over whole array. Default: None.

    Returns
    -------
    out: int
        Number of nonzero elements.

    References
    ----------
    `np.count_nonzero <https://numpy.org/doc/stable/reference/generated/numpy.count_nonzero.html#numpy.count_nonzero>`_
    """
    # TODO: Optimise this
    return sum(ndarr.__ne__(0), axis=axis)


def equal(
    x1: blosc2.Array,
    x2: blosc2.Array,
) -> blosc2.LazyExpr:
    """
    Computes the truth value of x1_i == x2_i for each element x1_i of the input array x1
    with the respective element x2_i of the input array x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array. May have any data type.

    x2:blosc2.Array
        Second input array. Must be compatible with x1. May have any data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.equal <https://numpy.org/doc/stable/reference/generated/numpy.equal.html#numpy.equal>`_
    """
    return x1.__eq__(x2)


def not_equal(
    x1: blosc2.Array,
    x2: blosc2.Array,
) -> blosc2.LazyExpr:
    """
    Computes the truth value of x1_i != x2_i for each element x1_i of the input array x1
    with the respective element x2_i of the input array x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array. May have any data type.

    x2:blosc2.Array
        Second input array. Must be compatible with x1. May have any data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.not_equal <https://numpy.org/doc/stable/reference/generated/numpy.not_equal.html#numpy.not_equal>`_
    """
    return x1.__ne__(x2)


def less_equal(
    x1: blosc2.Array,
    x2: blosc2.Array,
) -> blosc2.LazyExpr:
    """
    Computes the truth value of x1_i <= x2_i for each element x1_i of the input array x1
    with the respective element x2_i of the input array x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array. May have any data type.

    x2:blosc2.Array
        Second input array. Must be compatible with x1. May have any data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.less_equal <https://numpy.org/doc/stable/reference/generated/numpy.less_equal.html#numpy.less_equal>`_
    """
    return x1.__le__(x2)


def less(
    x1: blosc2.Array,
    x2: blosc2.Array,
) -> blosc2.LazyExpr:
    """
    Computes the truth value of x1_i < x2_i for each element x1_i of the input array x1
    with the respective element x2_i of the input array x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array. May have any data type.

    x2:blosc2.Array
        Second input array. Must be compatible with x1. May have any data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.less <https://numpy.org/doc/stable/reference/generated/numpy.less.html#numpy.less>`_
    """
    return x1.__lt__(x2)


def greater_equal(
    x1: blosc2.Array,
    x2: blosc2.Array,
) -> blosc2.LazyExpr:
    """
    Computes the truth value of x1_i >= x2_i for each element x1_i of the input array x1
    with the respective element x2_i of the input array x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array. May have any data type.

    x2:blosc2.Array
        Second input array. Must be compatible with x1. May have any data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.greater_equal <https://numpy.org/doc/stable/reference/generated/numpy.greater_equal.html#numpy.greater_equal>`_
    """
    return x1.__ge__(x2)


def greater(
    x1: blosc2.Array,
    x2: blosc2.Array,
) -> blosc2.LazyExpr:
    """
    Computes the truth value of x1_i > x2_i for each element x1_i of the input array x1
    with the respective element x2_i of the input array x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array. May have any data type.

    x2:blosc2.Array
        Second input array. Must be compatible with x1. May have any data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.greater <https://numpy.org/doc/stable/reference/generated/numpy.greater.html#numpy.greater>`_
    """
    return x1.__gt__(x2)


def multiply(
    x1: blosc2.Array,
    x2: blosc2.Array,
) -> blosc2.LazyExpr:
    """
    Computes the value of x1_i * x2_i for each element x1_i of the input array x1
    with the respective element x2_i of the input array x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array. May have any data type.

    x2:blosc2.Array
        Second input array. Must be compatible with x1. May have any data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.multiply <https://numpy.org/doc/stable/reference/generated/numpy.multiply.html#numpy.multiply>`_
    """
    return x1 * x2


def divide(
    x1: blosc2.Array,
    x2: blosc2.Array,
) -> blosc2.LazyExpr:
    """
    Computes the value of x1_i / x2_i for each element x1_i of the input array x1
    with the respective element x2_i of the input array x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array. May have any data type.

    x2:blosc2.Array
        Second input array. Must be compatible with x1. May have any data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.divide <https://numpy.org/doc/stable/reference/generated/numpy.divide.html#numpy.divide>`_
    """
    return x1 / x2


def nextafter(
    x1: blosc2.Array,
    x2: blosc2.Array,
) -> blosc2.LazyExpr:
    """
    Returns the next representable floating-point value for each element x1_i of the input
    array x1 in the direction of the respective element x2_i of the input array x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array. Real-valued floating point dtype.

    x2:blosc2.Array
        Second input array. Must be compatible with x1 and have same data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.nextafter <https://numpy.org/doc/stable/reference/generated/numpy.nextafter.html#numpy.nextafter>`_
    """
    return blosc2.LazyExpr(new_op=(x1, "nextafter", x2))


def hypot(
    x1: blosc2.Array,
    x2: blosc2.Array,
) -> blosc2.LazyExpr:
    """
    Computes the square root of the sum of squares for each element x1_i of the input array
    x1 with the respective element x2_i of the input array x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array. Real-valued floating point dtype.

    x2:blosc2.Array
        Second input array. Must be compatible with x1. Real-valued floating point dtype.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.hypot <https://numpy.org/doc/stable/reference/generated/numpy.hypot.html#numpy.hypot>`_
    """
    return blosc2.LazyExpr(new_op=(x1, "hypot", x2))


def copysign(
    x1: blosc2.Array,
    x2: blosc2.Array,
) -> blosc2.LazyExpr:
    """
    Composes a floating-point value with the magnitude of x1_i and the sign of x2_i
    for each element of the input array x1.

    Parameters
    ----------
    x1: blosc2.Array
        First input array. Real-valued floating point dtype.

    x2:blosc2.Array
        Second input array. Must be compatible with x1. Real-valued floating point dtype.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.copysign <https://numpy.org/doc/stable/reference/generated/numpy.copysign.html#numpy.copysign>`_
    """
    return blosc2.LazyExpr(new_op=(x1, "copysign", x2))


def maximum(
    x1: blosc2.Array,
    x2: blosc2.Array,
) -> blosc2.LazyExpr:
    """
    Computes the maximum value for each element x1_i of the input array x1 relative to the
    respective element x2_i of the input array x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array. Real-valued dtype.

    x2:blosc2.Array
        Second input array. Must be compatible with x1. Real-valued dtype.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.maximum <https://numpy.org/doc/stable/reference/generated/numpy.maximum.html#numpy.maximum>`_
    """
    return blosc2.LazyExpr(new_op=(x1, "maximum", x2))


def minimum(
    x1: blosc2.Array,
    x2: blosc2.Array,
) -> blosc2.LazyExpr:
    """
    Computes the minimum value for each element x1_i of the input array x1 relative to the
    respective element x2_i of the input array x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array. Real-valued dtype.

    x2:blosc2.Array
        Second input array. Must be compatible with x1. Real-valued dtype.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.minimum <https://numpy.org/doc/stable/reference/generated/numpy.minimum.html#numpy.minimum>`_
    """
    return blosc2.LazyExpr(new_op=(x1, "minimum", x2))


def reciprocal(x: blosc2.Array) -> blosc2.LazyExpr:
    """
    Computes the value of 1/x1_i for each element x1_i of the input array x1.

    Parameters
    ----------
    x: blosc2.Array
        First input array, floating-point data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.reciprocal <https://numpy.org/doc/stable/reference/generated/numpy.reciprocal.html#numpy.reciprocal>`_
    """
    return 1.0 / x


def floor(x: blosc2.Array) -> blosc2.LazyExpr:
    """
    Rounds each element x_i of the input array x to the greatest (i.e., closest to +infinity)
    integer-valued number that is not greater than x_i.

    Parameters
    ----------
    x: blosc2.Array
        First input array. May have any real-valued data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.floor <https://numpy.org/doc/stable/reference/generated/numpy.floor.html#numpy.floor>`_
    """
    return blosc2.LazyExpr(new_op=(x, "floor", None))


def ceil(x: blosc2.Array) -> blosc2.LazyExpr:
    """
    Rounds each element x_i of the input array x to the smallest (i.e., closest to -infinity)
    integer-valued number that is not smaller than x_i.

    Parameters
    ----------
    x: blosc2.Array
        First input array. May have any real-valued data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.ceil <https://numpy.org/doc/stable/reference/generated/numpy.ceil.html#numpy.ceil>`_
    """
    return blosc2.LazyExpr(new_op=(x, "ceil", None))


def trunc(x: blosc2.Array) -> blosc2.LazyExpr:
    """
    Rounds each element x_i of the input array x to the closest to 0
    integer-valued number.

    Parameters
    ----------
    x: blosc2.Array
        First input array. May have any real-valued data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.trunc <https://numpy.org/doc/stable/reference/generated/numpy.trunc.html#numpy.trunc>`_
    """
    return blosc2.LazyExpr(new_op=(x, "trunc", None))


def signbit(x: blosc2.Array) -> blosc2.LazyExpr:
    """
    Determines whether the sign bit is set for each element x_i of the input array x.

    The sign bit of a real-valued floating-point number x_i is set whenever x_i is either -0,
    less than zero, or a signed NaN (i.e., a NaN value whose sign bit is 1).

    Parameters
    ----------
    x: blosc2.Array
        First input array. May have any real-valued floating-point data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.signbit <https://numpy.org/doc/stable/reference/generated/numpy.signbit.html#numpy.signbit>`_
    """
    return blosc2.LazyExpr(new_op=(x, "signbit", None))


def sign(x: blosc2.Array) -> blosc2.LazyExpr:
    """
    Returns an indication of the sign of a number for each element x_i of the input array x.

    Parameters
    ----------
    x: blosc2.Array
        First input array. May have any numeric data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results (-1, 0 or 1).

    References
    ----------
    `np.sign <https://numpy.org/doc/stable/reference/generated/numpy.sign.html#numpy.sign>`_
    """
    return blosc2.LazyExpr(new_op=(x, "sign", None))


def round(x: blosc2.Array) -> blosc2.LazyExpr:
    """
    Rounds each element x_i of the input array x to the nearest integer-valued number.

    Parameters
    ----------
    x: blosc2.Array
        First input array. May have any numeric data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results (-1, 0 or 1).

    References
    ----------
    `np.round <https://numpy.org/doc/stable/reference/generated/numpy.round.html#numpy.round>`_
    """
    return blosc2.LazyExpr(new_op=(x, "round", None))


def floor_divide(
    x1: blosc2.Array,
    x2: blosc2.Array,
) -> blosc2.LazyExpr:
    """
    Computes the value of x1_i // x2_i for each element x1_i of the input array x1
    with the respective element x2_i of the input array x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array. May have any real-valued data type.

    x2:blosc2.Array
        Second input array. Must be compatible with x1. May have any real-valued data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.floor_divide <https://numpy.org/doc/stable/reference/generated/numpy.floor_divide.html#numpy.floor_divide>`_
    """
    return x1 // x2


def add(
    x1: blosc2.Array,
    x2: blosc2.Array,
) -> blosc2.LazyExpr:
    """
    Computes the value of x1_i + x2_i for each element x1_i of the input array x1
    with the respective element x2_i of the input array x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array. May have any data type.

    x2:blosc2.Array
        Second input array. Must be compatible with x1. May have any data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.add <https://numpy.org/doc/stable/reference/generated/numpy.add.html#numpy.add>`_
    """
    return x1 + x2


def subtract(
    x1: blosc2.Array,
    x2: blosc2.Array,
) -> blosc2.LazyExpr:
    """
    Computes the value of x1_i - x2_i for each element x1_i of the input array x1
    with the respective element x2_i of the input array x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array. May have any data type.

    x2:blosc2.Array
        Second input array. Must be compatible with x1. May have any data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.subtract <https://numpy.org/doc/stable/reference/generated/numpy.subtract.html#numpy.subtract>`_
    """
    return x1 - x2


def square(x1: blosc2.Array) -> blosc2.LazyExpr:
    """
    Computes the value of x1_i**2 for each element x1_i of the input array x1.

    Parameters
    ----------
    x1: blosc2.Array
        First input array. May have any data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.square <https://numpy.org/doc/stable/reference/generated/numpy.square.html#numpy.square>`_
    """
    return x1 * x1


def pow(
    x1: blosc2.Array | int | float | complex,
    x2: blosc2.Array | int | float | complex,
) -> blosc2.LazyExpr:
    """
    Computes the value of x1_i**x2_i for each element x1_i of the input array x1 and x2_i
    of x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array. May have any data type.

    x2:blosc2.Array
        Second input array. Must be compatible with x1. May have any data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.pow <https://numpy.org/doc/stable/reference/generated/numpy.pow.html#numpy.pow>`_
    """
    return x1**x2


def logical_xor(
    x1: blosc2.Array | int | float | complex,
    x2: blosc2.Array | int | float | complex,
) -> blosc2.LazyExpr:
    """
    Computes the value of x1_i ^ x2_i for each element x1_i of the input array x1 and x2_i
    of x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array, boolean.

    x2:blosc2.Array
        Second input array. Must be compatible with x1, boolean.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.logical_xor <https://numpy.org/doc/stable/reference/generated/numpy.logical_xor.html#numpy.logical_xor>`_
    """
    if blosc2.result_type(x1, x2) != blosc2.bool_:
        raise TypeError("Both operands must be boolean types for logical ops.")
    return x1 ^ x2


def logical_and(
    x1: blosc2.Array | int | float | complex,
    x2: blosc2.Array | int | float | complex,
) -> blosc2.LazyExpr:
    """
    Computes the value of x1_i & x2_i for each element x1_i of the input array x1 and x2_i
    of x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array, boolean.

    x2:blosc2.Array
        Second input array. Must be compatible with x1. Boolean.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.logical_and <https://numpy.org/doc/stable/reference/generated/numpy.logical_and.html#numpy.logical_and>`_
    """
    if blosc2.result_type(x1, x2) != blosc2.bool_:
        raise TypeError("Both operands must be boolean types for logical ops.")
    return x1 & x2


def logical_or(
    x1: blosc2.Array | int | float | complex,
    x2: blosc2.Array | int | float | complex,
) -> blosc2.LazyExpr:
    """
    Computes the value of x1_i | x2_i for each element x1_i of the input array x1 and x2_i
    of x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array, boolean.

    x2: blosc2.Array
        Second input array. Must be compatible with x1, boolean.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.logical_or <https://numpy.org/doc/stable/reference/generated/numpy.logical_or.html#numpy.logical_or>`_
    """
    if blosc2.result_type(x1, x2) != blosc2.bool_:
        raise TypeError("Both operands must be boolean types for logical ops.")
    return x1 | x2


def logical_not(
    x1: blosc2.Array | int | float | complex,
) -> blosc2.LazyExpr:
    """
    Computes the value of ~x1_i for each element x1_i of the input array x1.

    Parameters
    ----------
    x1: blosc2.Array
        Input array, boolean.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.logical_not <https://numpy.org/doc/stable/reference/generated/numpy.logical_not.html#numpy.logical_not>`_
    """
    if blosc2.result_type(x1) != blosc2.bool_:
        raise TypeError("Operand must be boolean type for logical ops.")
    return ~x1


def bitwise_xor(
    x1: blosc2.Array | int | float | complex,
    x2: blosc2.Array | int | float | complex,
) -> blosc2.LazyExpr:
    """
    Computes the value of x1_i ^ x2_i for each element x1_i of the input array x1 and x2_i
    of x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array, integer or boolean.

    x2:blosc2.Array
        Second input array. Must be compatible with x1, integer or boolean.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.bitwise_xor <https://numpy.org/doc/stable/reference/generated/numpy.bitwise_xor.html#numpy.bitwise_xor>`_
    """
    return x1 ^ x2


def bitwise_and(
    x1: blosc2.Array | int | float | complex,
    x2: blosc2.Array | int | float | complex,
) -> blosc2.LazyExpr:
    """
    Computes the value of x1_i & x2_i for each element x1_i of the input array x1 and x2_i
    of x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array, integer or boolean.

    x2:blosc2.Array
        Second input array. Must be compatible with x1. Integer or boolean.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.bitwise_and <https://numpy.org/doc/stable/reference/generated/numpy.bitwise_and.html#numpy.bitwise_and>`_
    """
    return x1 & x2


def bitwise_or(
    x1: blosc2.Array | int | float | complex,
    x2: blosc2.Array | int | float | complex,
) -> blosc2.LazyExpr:
    """
    Computes the value of x1_i | x2_i for each element x1_i of the input array x1 and x2_i
    of x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array, integer or boolean.

    x2: blosc2.Array
        Second input array. Must be compatible with x1, integer or boolean.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.bitwise_or <https://numpy.org/doc/stable/reference/generated/numpy.bitwise_or.html#numpy.bitwise_or>`_
    """
    return x1 | x2


def bitwise_invert(
    x1: blosc2.Array | int | float | complex,
) -> blosc2.LazyExpr:
    """
    Computes the value of ~x1_i for each element x1_i of the input array x1.

    Parameters
    ----------
    x1: blosc2.Array
        Input array, integer or boolean.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.bitwise_invert <https://numpy.org/doc/stable/reference/generated/numpy.bitwise_invert.html#numpy.bitwise_invert>`_
    """
    return ~x1


def bitwise_right_shift(
    x1: blosc2.Array | int | float | complex,
    x2: blosc2.Array | int | float | complex,
) -> blosc2.LazyExpr:
    """
    Shifts the bits of each element x1_i of the input array x1 to the right according to
    the respective element x2_i of the input array x2.

    Note: This operation is an arithmetic shift (i.e., sign-propagating) and thus equivalent to
    floor division by a power of two.

    Parameters
    ----------
    x1: blosc2.Array
        First input array, integer.

    x2: blosc2.Array
        Second input array. Must be compatible with x1, integer.
    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.bitwise_right_shift <https://numpy.org/doc/stable/reference/generated/numpy.bitwise_right_shift.html#numpy.bitwise_right_shift>`_
    """
    return x1.__rshift__(x2)


def bitwise_left_shift(
    x1: blosc2.Array | int | float | complex,
    x2: blosc2.Array | int | float | complex,
) -> blosc2.LazyExpr:
    """
    Shifts the bits of each element x1_i of the input array x1 to the left by appending x2_i
    (i.e., the respective element in the input array x2) zeros to the right of x1_i.

    Note: this operation is equivalent to multiplying x1 by 2**x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array, integer.

    x2: blosc2.Array
        Second input array. Must be compatible with x1, integer.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.bitwise_left_shift <https://numpy.org/doc/stable/reference/generated/numpy.bitwise_left_shift.html#numpy.bitwise_left_shift>`_
    """
    return x1.__lshift__(x2)


def positive(
    x1: blosc2.Array | int | float | complex,
) -> blosc2.LazyExpr:
    """
    Computes the numerical positive of each element x_i (i.e., out_i = +x_i) of the input array x.

    Parameters
    ----------
    x1: blosc2.Array
        First input array. May have any data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.positive <https://numpy.org/doc/stable/reference/generated/numpy.positive.html#numpy.positive>`_
    """
    return blosc2.LazyExpr(new_op=(0, "+", x1))


def negative(
    x1: blosc2.Array | int | float | complex,
) -> blosc2.LazyExpr:
    """
    Computes the numerical negative of each element x_i (i.e., out_i = -x_i) of the input array x.

    Parameters
    ----------
    x1: blosc2.Array
        First input array. May have any data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.negative <https://numpy.org/doc/stable/reference/generated/numpy.negative.html#numpy.negative>`_
    """
    return blosc2.LazyExpr(new_op=(0, "-", x1))


def remainder(
    x1: blosc2.Array | int | float | complex,
    x2: blosc2.Array | int | float | complex,
) -> blosc2.LazyExpr:
    """
    Returns the remainder of division for each element x1_i of the input array x1 and the
    respective element x2_i of the input array x2.

    Note: This function is equivalent to the Python modulus operator x1_i % x2_i.

    Parameters
    ----------
    x1: blosc2.Array
        First input array. May have any data type.

    x2: blosc2.Array
        Second input array. Must be compatible with x1. May have any data type.

    Returns
    -------
    out: LazyExpr
        A LazyArray containing the element-wise results.

    References
    ----------
    `np.remainder <https://numpy.org/doc/stable/reference/generated/numpy.remainder.html#numpy.remainder>`_
    """
    return blosc2.LazyExpr(new_op=(x1, "%", x2))


def clip(
    x: blosc2.Array,
    min: int | float | blosc2.Array | None = None,
    max: int | float | blosc2.Array | None = None,
    **kwargs: Any,
) -> NDArray:
    """
    Clamps each element x_i of the input array x to the range [min, max].

    Parameters
    ----------
    x: blosc2.Array
        Input array. Should have a real-valued data type.

    min: int | float | blosc2.Array | None
        Lower-bound of the range to which to clamp. If None, no lower bound must be applied.
        Default: None.

    max: int | float | blosc2.Array | None
        Upper-bound of the range to which to clamp. If None, no upper bound must be applied.
        Default: None.

    kwargs: Any
        kwargs accepted by the :func:`empty` constructor

    Returns
    -------
    out: NDArray
        An array containing element-wise results.

    """

    def chunkwise_clip(inputs, output, offset):
        x, min, max = inputs
        output[:] = np.clip(x, min, max)

    dtype = blosc2.result_type(x)
    return blosc2.lazyudf(chunkwise_clip, (x, min, max), dtype=dtype, shape=x.shape, **kwargs)


def logaddexp(x1: int | float | blosc2.Array, x2: int | float | blosc2.Array, **kwargs: Any) -> NDArray:
    """
    Calculates the logarithm of the sum of exponentiations log(exp(x1) + exp(x2)) for
    each element x1_i of the input array x1 with the respective element x2_i of the
    input array x2.

    Parameters
    ----------
    x1: blosc2.Array
        First input array. May have any real-valued floating-point data type.

    x2: blosc2.Array
        Second input array. Must be compatible with x1. May have any
        real-valued floating-point data type.

    kwargs: Any
        kwargs accepted by the :func:`empty` constructor

    Returns
    -------
    out: NDArray
        An array containing element-wise results.

    """

    def chunkwise_logaddexp(inputs, output, offset):
        x1, x2 = inputs
        output[:] = np.logaddexp(x1, x2)

    dtype = blosc2.result_type(x1, x2)
    if dtype == blosc2.bool_:
        raise TypeError("logaddexp doesn't accept boolean arguments.")

    if np.issubdtype(dtype, np.integer):
        dtype = blosc2.float32
    return blosc2.lazyudf(chunkwise_logaddexp, (x1, x2), dtype=dtype, shape=x1.shape, **kwargs)


# implemented in python-blosc2
local_ufunc_map = {
    np.logaddexp: logaddexp,
    np.logical_not: logical_not,
    np.logical_and: logical_and,
    np.logical_or: logical_or,
    np.logical_xor: logical_xor,
    np.matmul: matmul,
}


class Operand:
    """Base class for all operands in expressions."""

    _device = "cpu"

    def __array_namespace__(self, api_version: str | None = None) -> Any:
        """Return an object with all the functions and attributes of the module."""
        return blosc2

    # Provide minimal __array_interface__ to allow NumPy to work with this object
    @property
    def __array_interface__(self):
        return {
            "shape": self.shape,
            "typestr": self.dtype.str,
            "data": self[()],
            "version": 3,
        }

    @property
    @abstractmethod
    def dtype(self) -> np.dtype:
        """
        Get the data type of the :ref:`Operand`.

        Returns
        -------
        out: np.dtype
            The data type of the :ref:`Operand`.
        """
        pass

    @property
    @abstractmethod
    def shape(self) -> tuple[int]:
        """
        Get the shape of the :ref:`Operand`.

        Returns
        -------
        out: tuple
                The shape of the :ref:`Operand`.
        """
        pass

    @property
    @abstractmethod
    def ndim(self) -> int:
        """
        Get the number of dimensions of the :ref:`Operand`.

        Returns
        -------
        out: int
            The number of dimensions of the :ref:`Operand`.
        """
        pass

    @property
    @abstractmethod
    def info(self) -> InfoReporter:
        """
        Get information about the :ref:`Operand`.

        Returns
        -------
        out: InfoReporter
            A printable class with information about the :ref:`Operand`.
        """
        pass

    @property
    def device(self):
        "Hardware device the array data resides on. Always equal to 'cpu'."
        return self._device

    def to_device(self: NDArray, device: str):
        """
        Copy the array from the device on which it currently resides to the specified device.

        Parameters
        ----------
        self: NDArray
            Array instance.

        device: str
            Device to move array object to. Returns error except when device=='cpu'.

        Returns
        -------
        out: NDArray
            If device='cpu', the same array; else raises an Error.
        """
        if device != "cpu":
            raise ValueError(f"Unsupported device: {device}. Only 'cpu' is accepted.")
        return self

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        # Handle operations at the array level
        if method != "__call__":
            return NotImplemented

        if ufunc in local_ufunc_map:
            return local_ufunc_map[ufunc](*inputs)

        if ufunc in ufunc_map:
            value = inputs[0] if inputs[1] is self else inputs[1]
            _check_allowed_dtypes(value)
            return blosc2.LazyExpr(new_op=(inputs[0], ufunc_map[ufunc], inputs[1]))

        if ufunc in ufunc_map_1param:
            value = inputs[0]
            _check_allowed_dtypes(value)
            return blosc2.LazyExpr(new_op=(value, ufunc_map_1param[ufunc], None))

        return NotImplemented  # if not implemented in numexpr will default to NumPy

    def __add__(self, value: int | float | blosc2.Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "+", value))

    def __iadd__(self, value: int | float | blosc2.Array, /) -> blosc2.LazyExpr:
        return self.__add__(value)

    @is_documented_by(negative)
    def __neg__(self) -> blosc2.LazyExpr:
        return negative(self)

    @is_documented_by(positive)
    def __pos__(self) -> blosc2.LazyExpr:
        return positive(self)

    @is_documented_by(remainder)
    def __mod__(self, other) -> blosc2.LazyExpr:
        return remainder(self, other)

    def __radd__(self, value: int | float | blosc2.Array, /) -> blosc2.LazyExpr:
        return self.__add__(value)

    def __sub__(self, value: int | float | blosc2.Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "-", value))

    def __isub__(self, value: int | float | blosc2.Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "-", value))

    def __rsub__(self, value: int | float | blosc2.Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(value, "-", self))

    @is_documented_by(multiply)
    def __mul__(self, value: int | float | blosc2.Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "*", value))

    def __imul__(self, value: int | float | blosc2.Array, /) -> blosc2.LazyExpr:
        return self.__mul__(value)

    def __rmul__(self, value: int | float | blosc2.Array, /) -> blosc2.LazyExpr:
        return self.__mul__(value)

    def __truediv__(self, value: int | float | blosc2.Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "/", value))

    def __itruediv__(self, value: int | float | blosc2.Array, /) -> blosc2.LazyExpr:
        return self.__truediv__(value)

    def __rtruediv__(self, value: int | float | blosc2.Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(value, "/", self))

    @is_documented_by(floor_divide)
    def __floordiv__(self, value: int | float | blosc2.Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "//", value))

    def __lt__(self, value: int | float | blosc2.Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "<", value))

    def __le__(self, value: int | float | blosc2.Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "<=", value))

    def __gt__(self, value: int | float | blosc2.Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, ">", value))

    def __ge__(self, value: int | float | blosc2.Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, ">=", value))

    def __eq__(self, value: int | float | blosc2.Array, /):
        _check_allowed_dtypes(value)
        if blosc2._disable_overloaded_equal:
            return self is value
        return blosc2.LazyExpr(new_op=(self, "==", value))

    def __ne__(self, value: int | float | blosc2.Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "!=", value))

    def __pow__(self, value: int | float | blosc2.Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "**", value))

    def __ipow__(self, value: int | float | blosc2.Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "**", value))

    def __rpow__(self, value: int | float | blosc2.Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(value, "**", self))

    @is_documented_by(abs)
    def __abs__(self) -> blosc2.LazyExpr:
        return abs(self)

    @is_documented_by(bitwise_and)
    def __and__(self, value: int | float | blosc2.Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "&", value))

    @is_documented_by(bitwise_xor)
    def __xor__(self, other) -> blosc2.LazyExpr:
        return blosc2.LazyExpr(new_op=(self, "^", other))

    @is_documented_by(bitwise_or)
    def __or__(self, other) -> blosc2.LazyExpr:
        return blosc2.LazyExpr(new_op=(self, "|", other))

    @is_documented_by(bitwise_invert)
    def __invert__(self) -> blosc2.LazyExpr:
        return blosc2.LazyExpr(new_op=(self, "~", None))

    @is_documented_by(bitwise_right_shift)
    def __rshift__(self, other) -> blosc2.LazyExpr:
        return blosc2.LazyExpr(new_op=(self, ">>", other))

    @is_documented_by(bitwise_left_shift)
    def __lshift__(self, other) -> blosc2.LazyExpr:
        return blosc2.LazyExpr(new_op=(self, "<<", other))

    def __bool__(self) -> bool:
        if math.prod(self.shape) != 1:
            raise ValueError(f"The truth value of an array of shape {self.shape} is ambiguous.")
        return bool(self[()])

    def __float__(self) -> float:
        if math.prod(self.shape) != 1:
            raise ValueError(f"Cannot convert array of shape {self.shape} to float.")
        return float(self[()])

    def __int__(self) -> bool:
        if math.prod(self.shape) != 1:
            raise ValueError(f"Cannot convert array of shape {self.shape} to int.")
        return int(self[()])

    def __index__(self) -> bool:
        if not np.issubdtype(self.dtype, np.integer):
            raise ValueError(
                f"Cannot convert array of dtype {self.dtype} to index array (must have dtype int)."
            )
        return self.__int__()

    def __complex__(self) -> complex:
        if math.prod(self.shape) != 1:
            raise ValueError(f"Cannot convert array of shape {self.shape} to complex float.")
        return complex(self[()])

    def item(self) -> float | bool | complex | int:
        """
        Copy an element of an array to a standard Python scalar and return it.
        """
        return self[()].item()

    def where(self, value1=None, value2=None):
        """
        Select ``value1`` or ``value2`` values based on ``True``/``False`` for ``self``.

        Parameters
        ----------
        value1: array_like, optional
            The value to select when element of ``self`` is True.
        value2: array_like, optional
            The value to select when element of ``self`` is False.

        Returns
        -------
        out: LazyExpr
            A new expression with the where condition applied.
        """
        expr = blosc2.LazyExpr._new_expr("o0", {"o0": self}, guess=False)
        return expr.where(value1, value2)

    @is_documented_by(sum)
    def sum(self, axis=None, dtype=None, keepdims=False, **kwargs):
        expr = blosc2.LazyExpr(new_op=(self, None, None))
        return expr.sum(axis=axis, dtype=dtype, keepdims=keepdims, **kwargs)

    @is_documented_by(mean)
    def mean(self, axis=None, dtype=None, keepdims=False, **kwargs):
        expr = blosc2.LazyExpr(new_op=(self, None, None))
        return expr.mean(axis=axis, dtype=dtype, keepdims=keepdims, **kwargs)

    @is_documented_by(std)
    def std(self, axis=None, dtype=None, ddof=0, keepdims=False, **kwargs):
        expr = blosc2.LazyExpr(new_op=(self, None, None))
        return expr.std(axis=axis, dtype=dtype, ddof=ddof, keepdims=keepdims, **kwargs)

    @is_documented_by(var)
    def var(self, axis=None, dtype=None, ddof=0, keepdims=False, **kwargs):
        expr = blosc2.LazyExpr(new_op=(self, None, None))
        return expr.var(axis=axis, dtype=dtype, ddof=ddof, keepdims=keepdims, **kwargs)

    @is_documented_by(prod)
    def prod(self, axis=None, dtype=None, keepdims=False, **kwargs):
        expr = blosc2.LazyExpr(new_op=(self, None, None))
        return expr.prod(axis=axis, dtype=dtype, keepdims=keepdims, **kwargs)

    @is_documented_by(min)
    def min(self, axis=None, keepdims=False, **kwargs):
        expr = blosc2.LazyExpr(new_op=(self, None, None))
        return expr.min(axis=axis, keepdims=keepdims, **kwargs)

    @is_documented_by(max)
    def max(self, axis=None, keepdims=False, **kwargs):
        expr = blosc2.LazyExpr(new_op=(self, None, None))
        return expr.max(axis=axis, keepdims=keepdims, **kwargs)

    @is_documented_by(argmax)
    def argmax(self, axis=None, keepdims=False, **kwargs):
        expr = blosc2.LazyExpr(new_op=(self, None, None))
        return expr.argmax(axis=axis, keepdims=keepdims, **kwargs)

    @is_documented_by(argmin)
    def argmin(self, axis=None, keepdims=False, **kwargs):
        expr = blosc2.LazyExpr(new_op=(self, None, None))
        return expr.argmin(axis=axis, keepdims=keepdims, **kwargs)

    @is_documented_by(any)
    def any(self, axis=None, keepdims=False, **kwargs):
        expr = blosc2.LazyExpr(new_op=(self, None, None))
        return expr.any(axis=axis, keepdims=keepdims, **kwargs)

    @is_documented_by(all)
    def all(self, axis=None, keepdims=False, **kwargs):
        expr = blosc2.LazyExpr(new_op=(self, None, None))
        return expr.all(axis=axis, keepdims=keepdims, **kwargs)


class LimitedSizeDict(OrderedDict):
    def __init__(self, max_entries, *args, **kwargs):
        self.max_entries = max_entries
        super().__init__(*args, **kwargs)

    def __setitem__(self, key, value):
        if len(self) >= self.max_entries:
            self.popitem(last=False)
        super().__setitem__(key, value)


def detect_aligned_chunks(
    key: Sequence[slice], shape: Sequence[int], chunks: Sequence[int], consecutive: bool = False
) -> list[int]:
    """
    Detect whether a multidimensional slice is aligned with chunk boundaries.

    Parameters
    ----------
    key : Sequence of slice
        The multidimensional slice to check.
    shape : Sequence of int
        Shape of the NDArray.
    chunks : Sequence of int
        Chunk shape of the NDArray.
    consecutive : bool, default=False
        If True, check if the chunks are consecutive in storage order.
        If False, only check for chunk boundary alignment.

    Returns
    -------
    list[int]
        List of chunk indices (in C-order) that the slice overlaps with.
        If the slice isn't aligned with chunk boundaries, returns an empty list.
        If consecutive=True and chunks aren't consecutive, returns an empty list.
    """
    if len(key) != len(shape):
        return []

    # Check that slice boundaries are exact multiple of chunk boundaries
    for i, s in enumerate(key):
        if s.start is not None and s.start % chunks[i] != 0:
            return []
        if s.stop is not None and s.stop % chunks[i] != 0:
            return []

    # Parse the slice boundaries
    start_indices = []
    end_indices = []
    n_chunks = []

    for i, s in enumerate(key):
        start = s.start if s.start is not None else 0
        stop = s.stop if s.stop is not None else shape[i]
        chunk_size = chunks[i]
        start_idx = start // chunk_size
        end_idx = stop // chunk_size
        start_indices.append(start_idx)
        end_indices.append(end_idx)
        n_chunks.append(shape[i] // chunk_size)

    # Get all chunk combinations in the slice
    indices = [range(start, end) for start, end in zip(start_indices, end_indices, strict=False)]
    result = []

    for combination in product(*indices):
        flat_index = 0
        multiplier = 1
        for idx, n in zip(reversed(range(len(n_chunks))), reversed(n_chunks), strict=False):
            flat_index += combination[idx] * multiplier
            multiplier *= n
        result.append(flat_index)

    # Check if chunks are consecutive if requested
    if consecutive and result:
        sorted_result = sorted(result)
        if sorted_result[-1] - sorted_result[0] + 1 != len(sorted_result):
            return []

        # The array of indices must be consecutive
        for i in range(len(sorted_result) - 1):
            if sorted_result[i + 1] - sorted_result[i] != 1:
                return []

    return sorted(result)


class NDOuterIterator:
    def __init__(self, ndarray: NDArray | NDField, cache_size=1):
        self.ndarray = ndarray
        self.outer_dim_size = ndarray.shape[0]
        self.inner_shape = ndarray.shape[1:]
        self.current_index = 0
        # Cache for 1D arrays; for higher dimensions, the implementation should be more involved
        self.chunk_size = ndarray.chunks[0] if len(ndarray.shape) == 1 else None
        self.cache = {} if len(ndarray.shape) == 1 else None
        self.cache_size = cache_size

    def __iter__(self):
        return self

    def __next__(self):
        if self.current_index >= self.outer_dim_size:
            raise StopIteration

        outer_index = self.current_index
        self.current_index += 1

        if self.cache is not None:
            chunk_index = outer_index // self.chunk_size
            local_index = outer_index % self.chunk_size

            if chunk_index not in self.cache:
                if len(self.cache) >= self.cache_size:
                    self.cache.pop(next(iter(self.cache)))
                self.cache[chunk_index] = self.ndarray[
                    chunk_index * self.chunk_size : (chunk_index + 1) * self.chunk_size
                ]

            return self.cache[chunk_index][local_index]
        else:
            return self.ndarray[outer_index]


class NDArray(blosc2_ext.NDArray, Operand):
    def __init__(self, **kwargs):
        self._schunk = SChunk(_schunk=kwargs["_schunk"], _is_view=True)  # SChunk Python instance
        self._keep_last_read = False
        # Where to store the last read data
        self._last_read = {}
        base = kwargs.pop("_base", None)
        super().__init__(kwargs["_array"], base=base)
        # Accessor to fields
        self._fields = {}
        if self.dtype.fields:
            for field in self.dtype.fields:
                self._fields[field] = NDField(self, field)

    @property
    def cparams(self) -> blosc2.CParams:
        """The compression parameters used by the array."""
        return self.schunk.cparams

    @property
    def dparams(self) -> blosc2.DParams:
        """The decompression parameters used by the array."""
        return self.schunk.dparams

    @property
    def nbytes(self) -> int:
        """The number of bytes used by the array."""
        return self.schunk.nbytes

    @property
    def cbytes(self) -> int:
        """The number of compressed bytes used by the array."""
        return self.schunk.cbytes

    @property
    def cratio(self) -> float:
        """The compression ratio of the array."""
        return self.schunk.cratio

    # TODO: Uncomment when blosc2.Storage is available
    # @property
    # def storage(self) -> blosc2.Storage:
    #     """The storage of the array."""
    #     return self.schunk.storage

    @property
    def urlpath(self) -> str:
        """The URL path of the array."""
        return self.schunk.urlpath

    @property
    def meta(self) -> dict:
        """The metadata of the array."""
        return self.schunk.meta

    @property
    def vlmeta(self) -> dict:
        """The variable-length metadata of the array."""
        return self.schunk.vlmeta

    @property
    def fields(self) -> dict:
        """
        Dictionary with the fields of the structured array.

        Returns
        -------
        fields: dict
            A dictionary with the fields of the structured array.

        See Also
        --------
        :ref:`NDField`

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> shape = (10,)
        >>> dtype = np.dtype([('a', np.int32), ('b', np.float64)])
        >>> # Create a structured array
        >>> sa = blosc2.zeros(shape, dtype=dtype)
        >>> # Check that fields are equal
        >>> assert sa.fields['a'] == sa.fields['b']
        """
        return self._fields

    @property
    def keep_last_read(self) -> bool:
        """Indicates whether the last read data should be kept in memory."""
        return self._keep_last_read

    @keep_last_read.setter
    def keep_last_read(self, value: bool) -> None:
        """Set whether the last read data should be kept in memory.

        This always clears the last read data (if any).
        """
        if not isinstance(value, bool):
            raise TypeError("keep_last_read should be a boolean")
        # Reset last read data
        self._last_read.clear()
        self._keep_last_read = value

    @property
    def info(self) -> InfoReporter:
        """
        Print information about this array.

        Examples
        --------
        >>> import numpy as np
        >>> import blosc2
        >>> my_array = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        >>> array = blosc2.asarray(my_array)
        >>> print(array.info)
        type    : NDArray
        shape   : (10,)
        chunks  : (10,)
        blocks  : (10,)
        dtype   : int64
        cratio  : 0.73
        cparams : {'blocksize': 80,
        'clevel': 1,
        'codec': <Codec.ZSTD: 5>,
        'codec_meta': 0,
        'filters': [<Filter.NOFILTER: 0>,
                <Filter.NOFILTER: 0>,
                <Filter.NOFILTER: 0>,
                <Filter.NOFILTER: 0>,
                <Filter.NOFILTER: 0>,
                <Filter.SHUFFLE: 1>],
        'filters_meta': [0, 0, 0, 0, 0, 0],
        'nthreads': 4,
        'splitmode': <SplitMode.ALWAYS_SPLIT: 1>,
        'typesize': 8,
        'use_dict': 0}
        dparams : {'nthreads': 4}
        """
        return InfoReporter(self)

    @property
    def info_items(self) -> list:
        """A list of tuples with the information about this array.
        Each tuple contains the name of the attribute and its value.
        """
        items = []
        items += [("type", f"{self.__class__.__name__}")]
        items += [("shape", self.shape)]
        items += [("chunks", self.chunks)]
        items += [("blocks", self.blocks)]
        items += [("dtype", self.dtype)]
        items += [("nbytes", self.nbytes)]
        items += [("cbytes", self.cbytes)]
        items += [("cratio", f"{self.cratio:.2f}")]
        items += [("cparams", self.cparams)]
        items += [("dparams", self.dparams)]
        return items

    @property
    def schunk(self) -> blosc2.SChunk:
        """
        The :ref:`SChunk <SChunk>` reference of the :ref:`NDArray`.
        All the attributes from the :ref:`SChunk <SChunk>` can be accessed through
        this instance as `self.schunk`.

        See Also
        --------
        :ref:`SChunk Attributes <SChunkAttributes>`
        """
        return self._schunk

    @property
    def shape(self) -> tuple[int]:
        """Returns the data shape of this container.

        If the shape is a multiple of each dimension of :attr:`chunks`,
        it will be the same as :attr:`ext_shape`.

        See Also
        --------
        :attr:`ext_shape`
        """
        return super().shape

    @property
    def ext_shape(self) -> tuple[int]:
        """The padded data shape.

        The padded data is filled with zeros to make the real data fit into blocks and chunks, but it
        will never be retrieved as actual data (so the user can ignore this).
        In case :attr:`shape` is multiple in each dimension of :attr:`chunks` it will be the same
        as :attr:`shape`.

        See Also
        --------
        :attr:`shape`
        :attr:`chunks`
        """
        return super().ext_shape

    @property
    def chunks(self) -> tuple[int]:
        """Returns the data chunk shape of this container.

        If the chunk shape is a multiple of each dimension of :attr:`blocks`,
        it will be the same as :attr:`ext_chunks`.

        See Also
        --------
        :attr:`ext_chunks`
        """
        return super().chunks

    @property
    def ext_chunks(self) -> tuple[int]:
        """
        Returns the padded chunk shape which defines the chunksize in the associated schunk.

        This will be the chunk shape used to store each chunk, filling the extra positions
        with zeros (padding). If the :attr:`chunks` is a multiple of
        each dimension of :attr:`blocks` it will be the same as :attr:`chunks`.

        See Also
        --------
        :attr:`chunks`
        """
        return super().ext_chunks

    @property
    def blocks(self) -> tuple[int]:
        """The block shape of this container."""
        return super().blocks

    @property
    def ndim(self) -> int:
        """The number of dimensions of this container."""
        return super().ndim

    @property
    def size(self) -> int:
        """The size (in elements) for this container."""
        return super().size

    @property
    def chunksize(self) -> int:
        """Returns the data chunk size (in bytes) for this container.

        This will not be the same as
        :attr:`SChunk.chunksize <blosc2.schunk.SChunk.chunksize>`
        in case :attr:`chunks` is not multiple in
        each dimension of :attr:`blocks` (or equivalently, if :attr:`chunks` is
        not the same as :attr:`ext_chunks`).

        See Also
        --------
        :attr:`chunks`
        :attr:`ext_chunks`
        """
        return super().chunksize

    @property
    def dtype(self) -> np.dtype:
        """
        Data-type of the array's elements.
        """
        return super().dtype

    @property
    def blocksize(self) -> int:
        """The block size (in bytes) for this container.

        This is a shortcut to
        :attr:`SChunk.blocksize <blosc2.schunk.SChunk.blocksize>` and can be accessed
        through the :attr:`schunk` attribute as well.

        See Also
        --------
        :attr:`schunk`

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> array = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        >>> ndarray = blosc2.asarray(array)
        >>> print("Block size:", ndarray.blocksize)
        Block size: 80
        """
        return self._schunk.blocksize

    @property
    def oindex(self) -> OIndex:
        """Shortcut for orthogonal (outer) indexing, see :func:`get_oselection_numpy`"""
        return OIndex(self)

    # @property
    # def vindex(self) -> VIndex:
    #     """Shortcut for vectorised indexing. Not yet supported."""
    #     return VIndex(self)

    @property
    def T(self):
        """Return the transpose of a 2-dimensional array."""
        if self.ndim != 2:
            raise ValueError("This property only works for 2-dimensional arrays.")
        return blosc2.linalg.permute_dims(self)

    @property
    def mT(self):
        """Transpose of a matrix (or a stack of matrices)."""
        if self.ndim < 2:
            raise ValueError("This property only works for N-dimensional arrays with N>=2.")
        axes = np.arange(self.ndim)
        axes[-1] = self.ndim - 2
        axes[-2] = self.ndim - 1
        return blosc2.linalg.permute_dims(self, axes=axes)

    def get_fselection_numpy(self, key: list | np.ndarray) -> np.ndarray:
        """
        Select a slice from the array using a fancy index.
        Closely matches NumPy fancy indexing behaviour, except in
        some edge cases which are not supported by ndindex.
        Array indices separated by slice object - e.g. arr[0, :10, [0,1]] - are NOT supported.
        See https://www.blosc.org/posts/blosc2-fancy-indexing for more details.

        Parameters
        ----------
        key: list or np.ndarray

        Returns
        -------
        out: np.ndarray

        """
        # TODO: Make this faster and avoid running out of memory - avoid broadcasting keys

        ## Can't do this because ndindex doesn't support all the same indexing cases as Numpy
        # if math.prod(self.shape) * self.dtype.itemsize < blosc2.MAX_FAST_PATH_SIZE:
        #     return self[:][key]  # load into memory for smallish arrays
        shape = self.shape
        chunks = self.chunks

        # TODO: try to optimise and avoid this expand which seems to copy - maybe np.broadcast
        _slice = ndindex.ndindex(key).expand(shape)  # handles negative indices -> positive internally
        out_shape = _slice.newshape(shape)
        _slice = _slice.raw
        # now all indices are slices or arrays of integers (or booleans)
        # # moreover, all arrays are consecutive (otherwise an error is raised)

        if np.all([isinstance(s, (slice, np.ndarray)) for s in _slice]) and np.all(
            [s.dtype is not bool for s in _slice if isinstance(s, np.ndarray)]
        ):
            chunks = np.array(chunks)
            #       |------|
            # ------| arrs |------
            arridxs = [i for i, s in enumerate(_slice) if isinstance(s, np.ndarray)]
            begin, end = arridxs[0], arridxs[-1] + 1

            start, stop, step, _ = get_ndarray_start_stop(begin, _slice[:begin], self.shape[:begin])
            prior_tuple = tuple(
                slice(s, st, stp) for s, st, stp in zip(start, stop, step, strict=True)
            )  # convert to start and stop +ve
            start, stop, step, _ = get_ndarray_start_stop(
                len(self.shape[end:]), _slice[end:], self.shape[end:]
            )
            post_tuple = tuple(
                slice(s, st, stp) for s, st, stp in zip(start, stop, step, strict=True)
            )  # convert to start and stop +ve

            flat_shape = tuple(
                (i.stop - i.start - i.step // builtins.abs(i.step)) // i.step + 1 for i in prior_tuple
            )
            idx_dim = np.prod(_slice[begin].shape, dtype=np.int32)

            # TODO: find a nicer way to do the copy maybe
            arr = np.empty((idx_dim, end - begin), dtype=_slice[begin].dtype)
            for i, s in enumerate(_slice[begin:end]):
                arr[:, i] = s.reshape(-1)  # have to do a copy

            flat_shape += (idx_dim,)
            flat_shape += tuple(
                (i.stop - i.start - i.step // builtins.abs(i.step)) // i.step + 1 for i in post_tuple
            )
            # out_shape could have new dims if indexing arrays are not all 1D
            # (we have just flattened them so need to handle accordingly)
            divider = chunks[begin:end]
            chunked_arr = arr // divider
            if arr.shape[-1] == 1:  # 1D chunks, can avoid loading whole chunks
                idx_order = np.argsort(arr.squeeze(axis=1), axis=-1)  # sort by real index
                chunk_nitems = np.bincount(chunked_arr.reshape(-1), minlength=self.schunk.nchunks)
                unique_chunks = np.nonzero(chunk_nitems)[0][:, None]  # add dummy axis
                chunk_nitems = chunk_nitems[unique_chunks]
            else:
                chunked_arr = np.ascontiguousarray(
                    chunked_arr
                )  # ensure C-order memory to allow structured dtype view
                # TODO: check that avoids sort and copy (alternative: maybe do a bincount with structured data types?)
                _, row_ids, idx_inv, chunk_nitems = np.unique(
                    chunked_arr.view([("", chunked_arr.dtype)] * chunked_arr.shape[1]),
                    return_counts=True,
                    return_index=True,
                    return_inverse=True,
                )
                # In some versions of Numpy, output of np.unique has dummy dimension
                idx_inv = idx_inv if len(idx_inv.shape) == 1 else idx_inv.squeeze(-1)
                unique_chunks = chunked_arr[row_ids]
                # sort by chunks (can't sort by index since larger index could belong to lower chunk)
                # e.g. chunks of (100, 10) means (50, 15) has chunk idx (0,1) but (60,5) has (0, 0)
                idx_order = np.argsort(idx_inv)
            sorted_idxs = arr[idx_order]
            out = np.empty(flat_shape, dtype=self.dtype)
            shape = np.array(shape)

            chunk_nitems_cumsum = np.cumsum(chunk_nitems)
            cprior_slices = [
                slice_to_chunktuple(s, c) for s, c in zip(prior_tuple, chunks[:begin], strict=True)
            ]
            cpost_slices = [slice_to_chunktuple(s, c) for s, c in zip(post_tuple, chunks[end:], strict=True)]
            # TODO: rewrite to allow interleaved slices/array indexes
            for chunk_i, chunk_idx in enumerate(unique_chunks):
                start = 0 if chunk_i == 0 else chunk_nitems_cumsum[chunk_i - 1]
                stop = chunk_nitems_cumsum[chunk_i]
                selection = sorted_idxs[start:stop]
                out_mid_selection = (idx_order[start:stop],)
                if (
                    arr.shape[-1] == 1
                ):  # can avoid loading in whole chunk if 1D for array indexed chunks, a bit faster
                    chunk_begin = selection[0]
                    chunk_end = selection[-1] + 1
                else:
                    chunk_begin = chunk_idx * chunks[begin:end]
                    chunk_end = np.minimum((chunk_idx + 1) * chunks[begin:end], shape[begin:end])
                loc_mid_selection = tuple(a for a in (selection - chunk_begin).T)

                # loop over chunks coming from slices before and after array indices
                for cprior_tuple in product(*cprior_slices):
                    out_prior_selection, prior_selection, loc_prior_selection = _get_selection(
                        cprior_tuple, prior_tuple, chunks[:begin]
                    )
                    for cpost_tuple in product(*cpost_slices):
                        out_post_selection, post_selection, loc_post_selection = _get_selection(
                            cpost_tuple, post_tuple, chunks[end:]
                        )
                        locbegin, locend = _get_local_slice(
                            prior_selection, post_selection, (chunk_begin, chunk_end)
                        )
                        to_be_loaded = np.empty(locend - locbegin, dtype=self.dtype)
                        # basically load whole chunk, except for slice part at beginning and end
                        super().get_slice_numpy(to_be_loaded, (locbegin, locend))
                        loc_idx = loc_prior_selection + loc_mid_selection + loc_post_selection
                        out_idx = out_prior_selection + out_mid_selection + out_post_selection
                        out[out_idx] = to_be_loaded[loc_idx]
            return out.reshape(out_shape)  # should have filled in correct order, just need to reshape

        # Default when there are booleans
        # TODO: for boolean indexing could be optimised by avoiding
        # calculating out_shape prior to loop and keeping track on-the-fly (like in LazyExpr machinery)
        out = np.empty(out_shape, dtype=self.dtype)
        return self._get_set_findex_default(_slice, out)

    def _get_set_findex_default(self, _slice, out=None, value=None):
        _get = out is not None
        out = self if out is None else out  # default return for setitem with no intersecting chunks
        if 0 in self.shape:
            return out
        chunk_size = ndindex.ChunkSize(self.chunks)  # only works with nonzero chunks
        # repeated indices are grouped together
        intersecting_chunks = chunk_size.as_subchunks(
            _slice, self.shape
        )  # if _slice is (), returns all chunks
        for c in intersecting_chunks:
            sub_idx = _slice.as_subindex(c).raw
            sel_idx = c.as_subindex(_slice)
            start, stop, step, _ = get_ndarray_start_stop(self.ndim, c.raw, self.shape)
            chunk = np.empty(tuple(sp - st for st, sp in zip(start, stop, strict=True)), dtype=self.dtype)
            super().get_slice_numpy(chunk, (start, stop))
            if _get:
                new_shape = sel_idx.newshape(out.shape)
                out[sel_idx.raw] = chunk[sub_idx].reshape(new_shape)
            else:
                chunk[sub_idx] = value if np.isscalar(value) else value[sel_idx.raw]
                out = super().set_slice((start, stop), chunk)
        return out

    def get_oselection_numpy(self, key: list | np.ndarray) -> np.ndarray:
        """
        Select independently from self along axes specified in key. Key must be same length as self shape.
        See Zarr https://zarr.readthedocs.io/en/stable/user-guide/arrays.html#orthogonal-indexing.
        """
        shape = tuple(len(k) for k in key) + self.shape[len(key) :]
        # Create the array to store the result
        arr = np.empty(shape, dtype=self.dtype)
        return super().get_oindex_numpy(arr, key)

    def set_oselection_numpy(self, key: list | np.ndarray, arr: NDArray) -> np.ndarray:
        """
        Select independently from self along axes specified in key and set to entries in arr.
        Key must be same length as self shape.
        See Zarr https://zarr.readthedocs.io/en/stable/user-guide/arrays.html#orthogonal-indexing.
        """
        return super().set_oindex_numpy(key, arr)

    def _get_set_nonunit_steps(self, _slice, out=None, value=None):
        start, stop, step, mask = _slice
        _get = out is not None
        out = self if out is None else out  # default return for setitem with no intersecting chunks
        if 0 in self.shape:
            return out

        chunks = self.chunks
        _slice = tuple(slice(s, st, stp) for s, st, stp in zip(start, stop, step, strict=True))
        intersecting_chunks = [
            slice_to_chunktuple(s, c) for s, c in zip(_slice, chunks, strict=True)
        ]  # internally handles negative steps
        for c in product(*intersecting_chunks):
            sel_idx, glob_selection, sub_idx = _get_selection(c, _slice, chunks)
            sel_idx = tuple(s for s, m in zip(sel_idx, mask, strict=True) if not m)
            sub_idx = tuple(s if not m else s.start for s, m in zip(sub_idx, mask, strict=True))
            locstart, locstop = _get_local_slice(
                glob_selection,
                (),
                ((), ()),  # switches start and stop for negative steps
            )
            chunk = np.empty(
                tuple(sp - st for st, sp in zip(locstart, locstop, strict=True)), dtype=self.dtype
            )
            # basically load whole chunk, except for slice part at beginning and end
            super().get_slice_numpy(chunk, (locstart, locstop))  # copy relevant slice of chunk
            if _get:
                out[sel_idx] = chunk[sub_idx]  # update relevant parts of chunk
            else:
                chunk[sub_idx] = (
                    value if np.isscalar(value) else value[sel_idx]
                )  # update relevant parts of chunk
                out = super().set_slice((locstart, locstop), chunk)  # load updated partial chunk into array
        return out

    def __getitem__(
        self,
        key: None
        | int
        | slice
        | Sequence[slice | int | np.bool_ | np.ndarray[int | np.bool_] | None]
        | NDArray[int | np.bool_]
        | blosc2.LazyExpr
        | str,
    ) -> np.ndarray | blosc2.LazyExpr:
        """
        Retrieve a (multidimensional) slice as specified by the key.

        Note that this __getitem__ closely matches NumPy fancy indexing behaviour, except in
        some edge cases which are not supported by ndindex.
        Array indices separated by slice object - e.g. arr[0, :10, [0,1]] - are NOT supported.
        See https://www.blosc.org/posts/blosc2-fancy-indexing for more details.

        Parameters
        ----------
        key: int, slice, sequence of (slices, int), array of bools, LazyExpr or str
            The slice(s) to be retrieved. Note that step parameter is not yet honored
            in slices. If a LazyExpr is provided, the expression is expected to be of
            boolean type, and the result will be another LazyExpr returning the values
            of this array where the expression is True.
            When key is a (nd-)array of bools, the result will be the values of ``self``
            where the bool values are True (similar to NumPy).
            If key is an N-dim array of integers, the result will be the values of
            this array at the specified indices with the shape of the index.
            If the key is a string, and it is a field name of self, a :ref:`NDField`
            accessor will be returned; if not, it will be attempted to convert to a
            :ref:`LazyExpr`, and will search for its operands in the fields of ``self``.

        Returns
        -------
        out: np.ndarray | blosc2.LazyExpr
            The requested data as a NumPy array or a :ref:`LazyExpr`.

        Examples
        --------
        >>> import blosc2
        >>> shape = [25, 10]
        >>> # Create an array
        >>> a = blosc2.full(shape, 3.3333)
        >>> # Get slice as a NumPy array
        >>> a[:5, :5]
        array([[3.3333, 3.3333, 3.3333, 3.3333, 3.3333],
               [3.3333, 3.3333, 3.3333, 3.3333, 3.3333],
               [3.3333, 3.3333, 3.3333, 3.3333, 3.3333],
               [3.3333, 3.3333, 3.3333, 3.3333, 3.3333],
               [3.3333, 3.3333, 3.3333, 3.3333, 3.3333]])
        """
        # The more general case (this is quite slow)
        # If the key is a LazyExpr, decorate with ``where`` and return it
        if isinstance(key, blosc2.LazyExpr):
            return key.where(self)
        if isinstance(key, str):
            if self.dtype.fields is None:
                raise ValueError("The array is not structured (its dtype does not have fields)")
            if key in self.fields:
                # A shortcut to access fields
                return self.fields[key]
            # Assume that the key is a boolean expression
            expr = blosc2.LazyExpr._new_expr(key, self.fields, guess=False)
            return expr.where(self)

        key = key[()] if isinstance(key, NDArray) else key  # key not iterable
        key = tuple(k[()] if isinstance(k, NDArray) else k for k in key) if isinstance(key, tuple) else key

        # decompress NDArrays
        key_, mask = process_key(key, self.shape)  # internally handles key an integer
        key = key[()] if hasattr(key, "shape") and key.shape == () else key  # convert to scalar

        # fancy indexing
        if isinstance(key_, (list, np.ndarray)) or builtins.any(
            isinstance(k, (list, np.ndarray)) for k in key_
        ):
            # check scalar booleans, which add 1 dim to beginning
            if np.issubdtype(type(key), bool) and np.isscalar(key):
                if key:
                    _slice = ndindex.ndindex(()).expand(self.shape)  # just get whole array
                    out_shape = _slice.newshape(self.shape)
                    out = np.empty(out_shape, dtype=self.dtype)
                    return np.expand_dims(self._get_set_findex_default(_slice, out=out), 0)
                else:  # do nothing
                    return np.empty((0,) + self.shape, dtype=self.dtype)
            elif (
                hasattr(key, "dtype") and np.issubdtype(key.dtype, np.bool_) and key.shape == self.shape
            ):  # check ORIGINAL key
                # This can be interpreted as a boolean expression but only for key shape same as self shape
                expr = blosc2.LazyExpr._new_expr("key", {"key": key}, guess=False).where(self)
                # Decorate with where and force a getitem operation to return actual values.
                # This behavior is consistent with NumPy, although different from e.g. ['expr']
                # which returns a lazy expression.
                # This is faster than the fancy indexing path
                return expr[:]
            return self.get_fselection_numpy(key)  # fancy index default, can be quite slow

        start, stop, step, none_mask = get_ndarray_start_stop(self.ndim, key_, self.shape)
        shape = np.array(
            [(sp - st - np.sign(stp)) // stp + 1 for st, sp, stp in zip(start, stop, step, strict=True)]
        )
        if mask is not None:  # there are some dummy dims from ints
            # only get mask for not Nones in key to have nm_ same length as shape
            nm_ = [not m for m, n in zip(mask, none_mask, strict=True) if not n]
            # have to make none_mask refer to sliced dims (which will be less if ints present)
            none_mask = [n for m, n in zip(mask, none_mask, strict=True) if not m]
            shape = tuple(shape[nm_])

        # Create the array to store the result
        nparr = np.empty(shape, dtype=self.dtype)
        if step != (1,) * self.ndim:
            nparr = self._get_set_nonunit_steps((start, stop, step, [not i for i in nm_]), out=nparr)
        else:
            nparr = super().get_slice_numpy(nparr, (start, stop))

        if np.any(none_mask):
            nparr = np.expand_dims(nparr, axis=[i for i, n in enumerate(none_mask) if n])

        if self._keep_last_read:
            self._last_read.clear()
            inmutable_key = make_key_hashable(key)
            self._last_read[inmutable_key] = nparr

        return nparr

    def __setitem__(
        self,
        key: None | int | slice | Sequence[slice | int | np.bool_ | np.ndarray[int | np.bool_] | None],
        value: object,
    ):
        """Set a slice of the array.

        Parameters
        ----------
        key: int, slice or sequence of slices
            The index or indices specifying the slice(s) to be updated. Note that the step parameter
            is not yet supported.
        value: Py_Object Supporting the Buffer Protocol
            An object supporting the
            `Buffer Protocol <https://docs.python.org/3/c-api/buffer.html>`_
            which will be used to overwrite the specified slice(s).

        Examples
        --------
        >>> import blosc2
        >>> # Create an array
        >>> a = blosc2.full([8, 8], 3.3333)
        >>> # Set a slice to 0
        >>> a[:5, :5] = 0
        >>> a[:]
        array([[0.    , 0.    , 0.    , 0.    , 0.    , 3.3333, 3.3333, 3.3333],
               [0.    , 0.    , 0.    , 0.    , 0.    , 3.3333, 3.3333, 3.3333],
               [0.    , 0.    , 0.    , 0.    , 0.    , 3.3333, 3.3333, 3.3333],
               [0.    , 0.    , 0.    , 0.    , 0.    , 3.3333, 3.3333, 3.3333],
               [0.    , 0.    , 0.    , 0.    , 0.    , 3.3333, 3.3333, 3.3333],
               [3.3333, 3.3333, 3.3333, 3.3333, 3.3333, 3.3333, 3.3333, 3.3333],
               [3.3333, 3.3333, 3.3333, 3.3333, 3.3333, 3.3333, 3.3333, 3.3333],
               [3.3333, 3.3333, 3.3333, 3.3333, 3.3333, 3.3333, 3.3333, 3.3333]])
        """
        blosc2_ext.check_access_mode(self.schunk.urlpath, self.schunk.mode)

        # key not iterable
        key = key[()] if isinstance(key, NDArray) else key
        key = tuple(k[()] if isinstance(k, NDArray) else k for k in key) if isinstance(key, tuple) else key

        key_, mask = process_key(key, self.shape)  # internally handles key an integer
        if hasattr(value, "shape") and value.shape == ():
            value = value.item()
        value = (
            value if np.isscalar(value) else blosc2.as_simpleproxy(value)
        )  # convert to SimpleProxy for e.g. JAX, Tensorflow, PyTorch

        if builtins.any(isinstance(k, (list, np.ndarray)) for k in key_):  # fancy indexing
            _slice = ndindex.ndindex(key_).expand(
                self.shape
            )  # handles negative indices -> positive internally
            # check scalar booleans, which add 1 dim to beginning but which cause problems for ndindex.as_subindex
            if (
                key.shape == () and hasattr(key, "dtype") and np.issubdtype(key.dtype, np.bool_)
            ):  # check ORIGINAL key after decompression
                if key:
                    _slice = ndindex.ndindex(()).expand(self.shape)  # just get whole array
                else:  # do nothing
                    return self
            return self._get_set_findex_default(_slice, value=value)

        start, stop, step, none_mask = get_ndarray_start_stop(self.ndim, key_, self.shape)

        if step != (1,) * self.ndim:  # handle non-unit or negative steps
            if np.any(none_mask):
                raise ValueError("Cannot mix non-unit steps and None indexing for __setitem__.")
            return self._get_set_nonunit_steps((start, stop, step, mask), value=value)

        shape = [sp - st for sp, st in zip(stop, start, strict=False)]
        if isinstance(value, blosc2.Operand):  # handles SimpleProxy, NDArray, LazyExpr etc.
            value = value[()]  # convert to numpy
        if np.isscalar(value) or value.shape == ():
            value = np.full(shape, value, dtype=self.dtype)
        if value.dtype != self.dtype:  # handles decompressed NDArray too
            try:
                value = value.astype(self.dtype)
            except ComplexWarning:
                # numexpr type inference can lead to unnecessary type promotions
                # when using complex functions (e.g. conj) with real arrays
                value = value.real.astype(self.dtype)

        return super().set_slice((start, stop), value)

    def __iter__(self):
        """Iterate over the (outer) elements of the array.

        Returns
        -------
        out: iterator
        """
        return NDOuterIterator(self)

    def __len__(self) -> int:
        """Returns the length of the first dimension of the array.
        This is equivalent to ``self.shape[0]``.
        """
        if self.shape == ():
            raise TypeError("len() of unsized object")
        return self.shape[0]

    def get_chunk(self, nchunk: int) -> bytes:
        """Shortcut to :meth:`SChunk.get_chunk <blosc2.schunk.SChunk.get_chunk>`. This can be accessed
        through the :attr:`schunk` attribute as well.

        Parameters
        ----------
        nchunk: int
            The index of the chunk to retrieve.

        Returns
        -------
        chunk: bytes
            The chunk data at the specified index.

        See Also
        --------
        :attr:`schunk`
            The attribute that provides access to the underlying `SChunk` object.

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> # Create an SChunk with some data
        >>> array = np.arange(10)
        >>> ndarray = blosc2.asarray(array)
        >>> chunk = ndarray.get_chunk(0)
        >>> # Decompress the chunk to convert it into a numpy array
        >>> decompressed_chunk = blosc2.decompress(chunk)
        >>> np_array_chunk = np.frombuffer(decompressed_chunk, dtype=np.int64)
        >>> # Verify the content of the chunk
        >>> if isinstance(np_array_chunk, np.ndarray):
        >>>         print(np_array_chunk)
        >>>         print(np_array_chunk.shape) # Assuming chunk is a list or numpy array
        [ 0  1  2  3  4  5  6  7  8  9]
        (10,)
        """
        return self.schunk.get_chunk(nchunk)

    def reshape(self, shape: tuple[int], **kwargs: Any) -> NDArray:
        """Return a new array with the specified shape.

        See full documentation in :func:`reshape`.

        See Also
        --------
        :func:`reshape`
        """
        return reshape(self, shape, **kwargs)

    def iterchunks_info(
        self,
    ) -> Iterator[
        NamedTuple(
            "info",
            nchunk=int,
            coords=tuple,
            cratio=float,
            special=blosc2.SpecialValue,
            repeated_value=bytes | None,
            lazychunk=bytes,
        )
    ]:
        """
        Iterate over :paramref:`self` chunks of the array, providing information on index
        and special values.

        Yields
        ------
        info: namedtuple
            A namedtuple with the following fields:

                nchunk: int
                    The index of the chunk.
                coords: tuple
                    The coordinates of the chunk, in chunk units.
                cratio: float
                    The compression ratio of the chunk.
                special: :class:`SpecialValue`
                    The special value enum of the chunk; if 0, the chunk is not special.
                repeated_value: :attr:`self.dtype` or None
                    The repeated value for the chunk; if not SpecialValue.VALUE, it is None.
                lazychunk: bytes
                    A buffer containing the complete lazy chunk.

        Examples
        --------
        >>> import blosc2
        >>> a = blosc2.full(shape=(1000, ) * 3, fill_value=9, chunks=(500, ) * 3, dtype="f4")
        >>> for info in a.iterchunks_info():
        ...     print(info.coords)
        (0, 0, 0)
        (0, 0, 1)
        (0, 1, 0)
        (0, 1, 1)
        (1, 0, 0)
        (1, 0, 1)
        (1, 1, 0)
        (1, 1, 1)
        """
        ChunkInfoNDArray = namedtuple(
            "ChunkInfoNDArray", ["nchunk", "coords", "cratio", "special", "repeated_value", "lazychunk"]
        )
        chunks_idx = np.array(self.ext_shape) // np.array(self.chunks)
        for cinfo in self.schunk.iterchunks_info():
            nchunk, cratio, special, repeated_value, lazychunk = cinfo
            coords = tuple(np.unravel_index(cinfo.nchunk, chunks_idx))
            if cinfo.special == SpecialValue.VALUE:
                repeated_value = np.frombuffer(cinfo.repeated_value, dtype=self.dtype)[0]
            yield ChunkInfoNDArray(nchunk, coords, cratio, special, repeated_value, lazychunk)

    def tobytes(self) -> bytes:
        """Returns a buffer containing the data of the entire array.

        Returns
        -------
        out: bytes
            The buffer with the data of the whole array.

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> dtype = np.dtype("i4")
        >>> shape = [23, 11]
        >>> a = np.arange(0, int(np.prod(shape)), dtype=dtype).reshape(shape)
        >>> # Create an array
        >>> b = blosc2.asarray(a)
        >>> b.tobytes() == bytes(a[...])
        True
        """
        return super().tobytes()

    def to_cframe(self) -> bytes:
        """Get a bytes object containing the serialized :ref:`NDArray` instance.

        Returns
        -------
        out: bytes
            The buffer containing the serialized :ref:`NDArray` instance.

        See Also
        --------
        :func:`~blosc2.ndarray_from_cframe`
            This function can be used to reconstruct a NDArray from the serialized bytes.

        Examples
        --------
        >>> import blosc2
        >>> a = blosc2.full(shape=(1000, 1000), fill_value=9, dtype='i4')
        >>> # Get the bytes object containing the serialized instance
        >>> cframe_bytes = a.to_cframe()
        >>> blosc_array = blosc2.ndarray_from_cframe(cframe_bytes)
        >>> print("Shape of the NDArray:", blosc_array.shape)
        >>> print("Data type of the NDArray:", blosc_array.dtype)
        Shape of the NDArray: (1000, 1000)
        Data type of the NDArray: int32
        """
        return super().to_cframe()

    def copy(self, dtype: np.dtype | str = None, **kwargs: Any) -> NDArray:
        """Create a copy of an array with different parameters.

        Parameters
        ----------
        dtype: np.dtype or list str
            The new array dtype. Default is `self.dtype`.

        Other Parameters
        ----------------
        kwargs: dict, optional
            Additional keyword arguments supported by the :func:`empty` constructor.
            If not specified, the defaults will be taken from the original
            array (except for the urlpath).

        Returns
        -------
        out: :ref:`NDArray`
            A :ref:`NDArray` with a copy of the data.

        See Also
        --------
        :func:`copy`

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> shape = (10, 10)
        >>> blocks = (10, 10)
        >>> dtype = np.bool_
        >>> # Create a NDArray with default chunks
        >>> a = blosc2.zeros(shape, blocks=blocks, dtype=dtype)
        >>> # Get a copy with default chunks and blocks
        >>> b = a.copy(chunks=None, blocks=None)
        >>> np.array_equal(b[...], a[...])
        True
        """
        if dtype is None:
            dtype = self.dtype

        # Add the default parameters
        kwargs["cparams"] = kwargs.get("cparams", self.cparams)
        kwargs["dparams"] = kwargs.get("dparams", self.dparams)
        if "meta" in kwargs:
            # Do not allow to pass meta to copy
            raise ValueError("meta should not be passed to copy")

        kwargs = _check_ndarray_kwargs(**kwargs)
        return super().copy(dtype, **kwargs)

    def save(self, urlpath: str, contiguous=True, **kwargs: Any) -> None:
        """Save the array to a file.

        This is a convenience function that calls the :func:`copy` method with the
        `urlpath` parameter and the additional keyword arguments provided.

        See :func:`save` for more information.

        Parameters
        ----------
        urlpath: str
            The path where the array will be saved.
        contiguous: bool, optional
            Whether to save the array contiguously.

        Other Parameters
        ----------------
        kwargs: dict, optional
            Additional keyword arguments supported by the :func:`save` method.

        Returns
        -------
        out: None

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> shape = (10, 10)
        >>> blocks = (10, 10)
        >>> dtype = np.bool_
        >>> # Create a NDArray with default chunks
        >>> a = blosc2.zeros(shape, blocks=blocks, dtype=dtype)
        >>> # Save the array to a file
        >>> a.save("array.b2frame")
        """
        blosc2_ext.check_access_mode(urlpath, "w")
        # Add urlpath to kwargs
        kwargs["urlpath"] = urlpath
        # Add the contiguous parameter
        kwargs["contiguous"] = contiguous

        super().copy(self.dtype, **kwargs)

    def resize(self, newshape: tuple | list) -> None:
        """Change the shape of the array by growing or shrinking one or more dimensions.

        Parameters
        ----------
        newshape : tuple or list
            The new shape of the array. It should have the same number of dimensions
            as :paramref:`self`, the current shape.

        Returns
        -------
        out: None

        Notes
        -----
        The array values in the newly added positions are not initialized.
        The user is responsible for initializing them.

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> import math
        >>> dtype = np.dtype(np.float32)
        >>> shape = [23, 11]
        >>> a = np.linspace(1, 3, num=math.prod(shape)).reshape(shape)
        >>> # Create an array
        >>> b = blosc2.asarray(a)
        >>> newshape = [50, 10]
        >>> # Extend first dimension, shrink second dimension
        >>> b.resize(newshape)
        >>> b.shape
        (50, 10)
        """
        blosc2_ext.check_access_mode(self.schunk.urlpath, self.schunk.mode)
        super().resize(newshape)

    def slice(self, key: int | slice | Sequence[slice], **kwargs: Any) -> NDArray:
        """Get a (multidimensional) slice as a new :ref:`NDArray`.

        Parameters
        ----------
        key: int, slice or sequence of slices
            The index for the slices to be retrieved. Note that the step parameter is
            not yet supported in slices.

        Other Parameters
        ----------------
        kwargs: dict, optional
            Additional keyword arguments supported by the :func:`empty` constructor.

        Returns
        -------
        out: :ref:`NDArray`
            An array containing the requested data. The dtype will match that of `self`.

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> shape = [23, 11]
        >>> a = np.arange(np.prod(shape)).reshape(shape)
        >>> # Create an array
        >>> b = blosc2.asarray(a)
        >>> slices = (slice(3, 7), slice(1, 11))
        >>> # Get a slice as a new NDArray
        >>> c = b.slice(slices)
        >>> print(c.shape)
        (4, 10)
        >>> print(type(c))
        <class 'blosc2.ndarray.NDArray'>

        Notes
        -----
        There is a fast path for slices that are aligned with underlying chunks.
        Aligned means that the slices are made entirely with complete chunks.
        """
        if "cparams" not in kwargs:
            kwargs["cparams"] = {
                "codec": self.cparams.codec,
                "clevel": self.cparams.clevel,
                "filters": self.cparams.filters,
            }
        kwargs = _check_ndarray_kwargs(**kwargs)  # sets cparams to defaults
        key, mask = process_key(key, self.shape)
        start, stop, step, _ = get_ndarray_start_stop(self.ndim, key, self.shape)

        # Fast path for slices made with aligned chunks
        if step == (1,) * self.ndim:
            aligned_chunks = detect_aligned_chunks(key, self.shape, self.chunks, consecutive=False)
            if aligned_chunks:
                # print("Aligned chunks detected", aligned_chunks)
                # Create a new ndarray for the key slice
                new_shape = [
                    sp - st for sp, st in zip([k.stop for k in key], [k.start for k in key], strict=False)
                ]
                newarr = blosc2.empty(
                    shape=new_shape,
                    dtype=self.dtype,
                    chunks=self.chunks,
                    blocks=self.blocks,
                    **kwargs,
                )
                # Get the chunks from the original array and update the new array
                # No need for chunks to decompress and compress again
                for order, nchunk in enumerate(aligned_chunks):
                    chunk = self.schunk.get_chunk(nchunk)
                    newarr.schunk.update_chunk(order, chunk)
                return newarr.squeeze(axis=np.where(mask)[0])  # remove any dummy dims introduced

        key = (start, stop)
        ndslice = super().get_slice(key, mask, **kwargs)

        # This is memory intensive, but we have not a better way to do it yet
        # TODO: perhaps add a step param in the get_slice method in the future?
        if step != (1,) * self.ndim:
            nparr = ndslice[...]
            if len(step) == 1:
                nparr = nparr[:: step[0]]
            else:
                slice_ = tuple(slice(None, None, st) for st in step)
                nparr = nparr[slice_]
            return asarray(nparr, **kwargs)

        return ndslice

    def squeeze(self, axis: int | Sequence[int]) -> NDArray:
        """Remove single-dimensional entries from the shape of the array.

        This method modifies the array in-place. If mask is None removes any dimensions with size 1.
        If axis is provided, it should be an int or tuple of ints and the corresponding
        dimensions (of size 1) will be removed.

        Returns
        -------
        out: NDArray

        Examples
        --------
        >>> import blosc2
        >>> shape = [1, 23, 1, 11, 1]
        >>> # Create an array
        >>> a = blosc2.full(shape, 2**30)
        >>> a.shape
        (1, 23, 1, 11, 1)
        >>> # Squeeze the array
        >>> a.squeeze()
        >>> a.shape
        (23, 11)
        """
        return blosc2.squeeze(self, axis=axis)

    def indices(self, order: str | list[str] | None = None, **kwargs: Any) -> NDArray:
        """
        Return the indices of a sorted array following the specified order.

        This is only valid for 1-dim structured arrays.

        See full documentation in :func:`indices`.
        """
        return indices(self, order, **kwargs)

    def sort(self, order: str | list[str] | None = None, **kwargs: Any) -> NDArray:
        """
        Return a sorted array following the specified order, or the order of the fields.

        This is only valid for 1-dim structured arrays.

        See full documentation in :func:`sort`.
        """
        return sort(self, order, **kwargs)

    def as_ffi_ptr(self):
        """Returns the pointer to the raw FFI blosc2::b2nd_array_t object.

        This function is useful for passing the array to C functions.
        """
        return super().as_ffi_ptr()

    def __matmul__(self, other):
        return blosc2.linalg.matmul(self, other)


def squeeze(x: Array, axis: int | Sequence[int]) -> NDArray:
    """
    Remove single-dimensional entries from the shape of the array.

    This method modifies the array in-place.

    Parameters
    ----------
    x: Array
        input array.
    axis: int | Sequence[int]
        Axis (or axes) to squeeze.

    Returns
    -------
    out: Array
        An output array having the same data type and elements as x.

    Examples
    --------
    >>> import blosc2
    >>> shape = [1, 23, 1, 11, 1]
    >>> # Create an array
    >>> b = blosc2.full(shape, 2**30)
    >>> b.shape
    (1, 23, 1, 11, 1)
    >>> # Squeeze the array
    >>> blosc2.squeeze(b)
    >>> b.shape
    (23, 11)
    """
    axis = [axis] if isinstance(axis, int) else axis
    mask = [False for i in range(x.ndim)]
    for a in axis:
        if a < 0:
            a += x.ndim  # Adjust axis to be within the array's dimensions
        if mask[a]:
            raise ValueError("Axis values must be unique.")
        mask[a] = True
    return blosc2_ext.squeeze(x, axis_mask=mask)


def array_from_ffi_ptr(array_ptr) -> NDArray:
    """
    Create an NDArray from a raw FFI pointer.

    This function is useful for passing arrays across FFI boundaries.
    This function move the ownership of the underlying `b2nd_array_t*` object to the new NDArray, and it will be freed
    when the object is destroyed.
    """
    return blosc2_ext.array_from_ffi_ptr(array_ptr)


def where(
    condition: blosc2.LazyExpr | NDArray,
    x: blosc2.Array | int | float | complex | bool | str | bytes | None = None,
    y: blosc2.Array | int | float | complex | bool | str | bytes | None = None,
) -> blosc2.LazyExpr:
    """
    Return elements chosen from `x` or `y` depending on `condition`.

    Parameters
    ----------
    condition: :ref:`LazyExpr`
        Where True, yield `x`, otherwise yield `y`.
    x: :ref:`NDArray` or :ref:`NDField` or np.ndarray or scalar or bytes
        Values from which to choose when `condition` is True.
    y: :ref:`NDArray` or :ref:`NDField` or np.ndarray or scalar or bytes
        Values from which to choose when `condition` is False.

    References
    ----------
    `np.where <https://numpy.org/doc/stable/reference/generated/numpy.where.html>`_
    """
    return condition.where(x, y)


def lazywhere(value1=None, value2=None):
    """Decorator to apply a where condition to a LazyExpr."""

    def inner_decorator(func):
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs).where(value1, value2)

        return wrapper

    return inner_decorator


def _check_shape(shape):
    if isinstance(shape, int | np.integer):
        shape = (shape,)
    elif not isinstance(shape, tuple | list):
        raise TypeError("shape should be a tuple or a list!")
    if len(shape) > blosc2.MAX_DIM:
        raise ValueError(f"shape length {len(shape)} is too large (>{blosc2.MAX_DIM})!")
    return shape


def _check_dtype(dtype):
    dtype = np.dtype(dtype)
    if dtype.itemsize > blosc2.MAX_TYPESIZE:
        raise ValueError(f"dtype itemsize {dtype.itemsize} is too large (>{blosc2.MAX_TYPESIZE})!")
    return dtype


def empty(shape: int | tuple | list, dtype: np.dtype | str | None = np.float64, **kwargs: Any) -> NDArray:
    """Create an empty array.

    Parameters
    ----------
    shape: int, tuple or list
        The shape for the final array.
    dtype: np.dtype or list str
        The data type of the array elements in NumPy format. Default is `np.uint8`.
        This will override the `typesize`
        in the compression parameters if they are provided.

    Other Parameters
    ----------------
    kwargs: dict, optional
        Keyword arguments supported:
            chunks: tuple or list
                The chunk shape. If None (default), Blosc2 will compute
                an efficient chunk shape.
            blocks: tuple or list
                The block shape. If None (default), Blosc2 will compute
                an efficient block shape. This will override the `blocksize`
                in the cparams if they are provided.

        The other keyword arguments supported are the same as for the
        :obj:`SChunk.__init__ <blosc2.schunk.SChunk.__init__>` constructor.

    Returns
    -------
    out: :ref:`NDArray`
        A :ref:`NDArray` is returned.

    Examples
    --------
    >>> import blosc2
    >>> import numpy as np
    >>> shape = [20, 20]
    >>> dtype = np.int32
    >>> # Create empty array with default chunks and blocks
    >>> array = blosc2.empty(shape, dtype=dtype)
    >>> array.shape
    (20, 20)
    >>> array.dtype
    dtype('int32')
    """
    dtype = _check_dtype(dtype)
    shape = _check_shape(shape)
    kwargs = _check_ndarray_kwargs(**kwargs)
    chunks = kwargs.pop("chunks", None)
    blocks = kwargs.pop("blocks", None)
    chunks, blocks = compute_chunks_blocks(shape, chunks, blocks, dtype, **kwargs)
    return blosc2_ext.empty(shape, chunks, blocks, dtype, **kwargs)


def uninit(shape: int | tuple | list, dtype: np.dtype | str = np.float64, **kwargs: Any) -> NDArray:
    """Create an array with uninitialized values.

    The parameters and keyword arguments are the same as for the
    :func:`empty` constructor.

    Returns
    -------
    out: :ref:`NDArray`
        A :ref:`NDArray` is returned.

    Examples
    --------
    >>> import blosc2
    >>> shape = [8, 8]
    >>> chunks = [6, 5]
    >>> # Create uninitialized array
    >>> array = blosc2.uninit(shape, dtype='f8', chunks=chunks)
    >>> array.shape
    (8, 8)
    >>> array.chunks
    (6, 5)
    >>> array.dtype
    dtype('float64')
    """
    dtype = _check_dtype(dtype)
    shape = _check_shape(shape)
    kwargs = _check_ndarray_kwargs(**kwargs)
    chunks = kwargs.pop("chunks", None)
    blocks = kwargs.pop("blocks", None)
    chunks, blocks = compute_chunks_blocks(shape, chunks, blocks, dtype, **kwargs)
    return blosc2_ext.uninit(shape, chunks, blocks, dtype, **kwargs)


def nans(shape: int | tuple | list, dtype: np.dtype | str = np.float64, **kwargs: Any) -> NDArray:
    """Create an array with NaNs values.

    The parameters and keyword arguments are the same as for the
    :func:`empty` constructor.

    Returns
    -------
    out: :ref:`NDArray <NDArray>`
        A :ref:`NDArray <NDArray>` is returned.

    Examples
    --------
    >>> import blosc2
    >>> shape = [8, 8]
    >>> chunks = [6, 5]
    >>> # Create an array of NaNs
    >>> array = blosc2.nans(shape, dtype='f8', chunks=chunks)
    >>> array.shape
    (8, 8)
    >>> array.chunks
    (6, 5)
    >>> array.dtype
    dtype('float64')
    """
    dtype = _check_dtype(dtype)
    shape = _check_shape(shape)
    kwargs = _check_ndarray_kwargs(**kwargs)
    chunks = kwargs.pop("chunks", None)
    blocks = kwargs.pop("blocks", None)
    chunks, blocks = compute_chunks_blocks(shape, chunks, blocks, dtype, **kwargs)
    return blosc2_ext.nans(shape, chunks, blocks, dtype, **kwargs)


def zeros(shape: int | tuple | list, dtype: np.dtype | str = np.float64, **kwargs: Any) -> NDArray:
    """Create an array with zero as the default value
    for uninitialized portions of the array.

    The parameters and keyword arguments are the same as for the
    :func:`empty` constructor.

    Returns
    -------
    out: :ref:`NDArray`
        A :ref:`NDArray` is returned.

    Examples
    --------
    >>> import blosc2
    >>> import numpy as np
    >>> shape = [8, 8]
    >>> chunks = [6, 5]
    >>> blocks = [5, 5]
    >>> dtype = np.float64
    >>> # Create zeros array
    >>> array = blosc2.zeros(shape, dtype=dtype, chunks=chunks, blocks=blocks)
    >>> array.shape
    (8, 8)
    >>> array.chunks
    (6, 5)
    >>> array.blocks
    (5, 5)
    >>> array.dtype
    dtype('float64')
    """
    dtype = _check_dtype(dtype)
    shape = _check_shape(shape)
    kwargs = _check_ndarray_kwargs(**kwargs)
    chunks = kwargs.pop("chunks", None)
    blocks = kwargs.pop("blocks", None)
    chunks, blocks = compute_chunks_blocks(shape, chunks, blocks, dtype, **kwargs)
    return blosc2_ext.zeros(shape, chunks, blocks, dtype, **kwargs)


def full(
    shape: int | tuple | list,
    fill_value: bytes | int | float | bool,
    dtype: np.dtype | str = None,
    **kwargs: Any,
) -> NDArray:
    """Create an array, with :paramref:`fill_value` being used as the default value
    for uninitialized portions of the array.

    Parameters
    ----------
    shape: int, tuple or list
        The shape of the final array.
    fill_value: bytes, int, float or bool
        Default value to use for uninitialized portions of the array.
        Its size will override the `typesize`
        in the cparams if they are passed.
    dtype: np.dtype or list str
        The ndarray dtype in NumPy format. By default, this will
        be taken from the :paramref:`fill_value`.
        This will override the `typesize`
        in the cparams if they are passed.

    Other Parameters
    ----------------
    kwargs: dict, optional
        Keyword arguments that are supported by the :func:`empty` constructor.

    Returns
    -------
    out: :ref:`NDArray`
        A :ref:`NDArray` is returned.

    Examples
    --------
    >>> import blosc2
    >>> import numpy as np
    >>> shape = [25, 10]
    >>> # Create array filled with True
    >>> array = blosc2.full(shape, True)
    >>> array.shape
    (25, 10)
    >>> array.dtype
    dtype('bool')
    """
    if isinstance(fill_value, bytes):
        dtype = np.dtype(f"S{len(fill_value)}")
    if dtype is None:
        dtype = np.dtype(type(fill_value))
    else:
        dtype = np.dtype(dtype)
    dtype = _check_dtype(dtype)
    shape = _check_shape(shape)
    kwargs = _check_ndarray_kwargs(**kwargs)
    chunks = kwargs.pop("chunks", None)
    blocks = kwargs.pop("blocks", None)
    chunks, blocks = compute_chunks_blocks(shape, chunks, blocks, dtype, **kwargs)
    return blosc2_ext.full(shape, chunks, blocks, fill_value, dtype, **kwargs)


def ones(shape: int | tuple | list, dtype: np.dtype | str = None, **kwargs: Any) -> NDArray:
    """Create an array with one as values.

    The parameters and keyword arguments are the same as for the
    :func:`empty` constructor.

    Returns
    -------
    out: :ref:`NDArray`
            A :ref:`NDArray` is returned.

    Examples
    --------
    >>> import blosc2
    >>> import numpy as np
    >>> shape = [8, 8]
    >>> chunks = [6, 5]
    >>> blocks = [5, 5]
    >>> dtype = np.float64
    >>> # Create ones array
    >>> array = blosc2.ones(shape, dtype=dtype, chunks=chunks, blocks=blocks)
    >>> array.shape
    (8, 8)
    >>> array.chunks
    (6, 5)
    >>> array.blocks
    (5, 5)
    >>> array.dtype
    dtype('float64')
    """
    if dtype is None:
        dtype = blosc2.DEFAULT_FLOAT
    return full(shape, 1, dtype, **kwargs)


def arange(
    start: int | float,
    stop: int | float | None = None,
    step: int | float | None = 1,
    dtype: np.dtype | str = None,
    shape: int | tuple | list | None = None,
    c_order: bool = True,
    **kwargs: Any,
) -> NDArray:
    """
    Return evenly spaced values within a given interval.
    Due to rounding errors for chunkwise filling, may differ
    from numpy.arange in edge cases.

    Parameters
    ----------
    start: int, float
        The starting value of the sequence.
    stop: int, float
        The end value of the sequence.
    step: int, float or None
        Spacing between values.
    dtype: np.dtype or list str
        The data type of the array elements in NumPy format. Default is
        None.  If dtype is None, inferred from start, stop and step.
        Output type is integer unless one or more have type float.
        This will override the `typesize` in the compression parameters if
        they are provided.
    shape: int, tuple or list
        The shape of the final array. If None, the shape will be computed.
    c_order: bool
        Whether to store the array in C order (row-major) or insertion order.
        Insertion order means that values will be stored in the array
        following the order of chunks in the array; this is more memory
        efficient, as it does not require an intermediate copy of the array.
        Default is C order.

    Other Parameters
    ----------------
    kwargs: dict, optional
        Keyword arguments that are supported by the :func:`empty` constructor.

    Returns
    -------
    out: :ref:`NDArray`
        A :ref:`NDArray` is returned.

    Examples
    --------
    >>> import blosc2
    >>> import numpy as np
    >>> # Create an array with values from 0 to 10
    >>> array = blosc2.arange(0, 10, 1)
    >>> print(array)
    [0 1 2 3 4 5 6 7 8 9]
    """

    def arange_fill(inputs, output, offset):
        lout = len(output)
        start, _, step = inputs
        start += offset[0] * step
        stop = start + lout * step
        if math.ceil((stop - start) / step) == lout:  # USE ARANGE IF POSSIBLE (2X FASTER)
            output[:] = np.arange(start, stop, step, dtype=output.dtype)
        else:  # use linspace to have finer control over exclusion of endpoint for float types
            output[:] = np.linspace(start, stop, lout, endpoint=False, dtype=output.dtype)

    if step is None:  # not array-api compliant but for backwards compatibility
        step = 1
    if stop is None:
        stop = start
        start = 0
    NUM = int((stop - start) / step)
    if shape is None:
        shape = (builtins.max(NUM, 0),)
    else:
        # Check that the shape is consistent with the start, stop and step values
        if math.prod(shape) != NUM:
            raise ValueError("The shape is not consistent with the start, stop and step values")
    if dtype is None:
        dtype = (
            blosc2.DEFAULT_FLOAT
            if np.any([np.issubdtype(type(d), float) for d in (start, stop, step)])
            else blosc2.DEFAULT_INT
        )
    dtype = _check_dtype(dtype)

    if is_inside_new_expr() or NUM < 0:
        # We already have the dtype and shape, so return immediately
        return blosc2.zeros(shape, dtype=dtype, **kwargs)

    lshape = (math.prod(shape),)
    lazyarr = blosc2.lazyudf(arange_fill, (start, stop, step), dtype=dtype, shape=lshape)

    if len(shape) == 1:
        # C order is guaranteed, and no reshape is needed
        return lazyarr.compute(**kwargs)

    return reshape(lazyarr, shape, c_order=c_order, **kwargs)


# Define a numpy linspace-like function
def linspace(
    start: int | float | complex,
    stop: int | float | complex,
    num: int | None = None,
    dtype=None,
    endpoint: bool = True,
    shape=None,
    c_order: bool = True,
    **kwargs: Any,
) -> NDArray:
    """Return evenly spaced numbers over a specified interval.

    This is similar to `numpy.linspace` but it returns a `NDArray`
    instead of a numpy array.  Also, it supports a `shape` parameter
    to return a ndim array.

    Parameters
    ----------
    start: int, float, complex
        The starting value of the sequence.
    stop: int, float, complex
        The end value of the sequence.
    num: int | None
        Number of samples to generate. Default None.
    dtype: np.dtype or list str
        The data type of the array elements in NumPy format. If None, inferred from
        start, stop, step. Default is None.
    endpoint: bool
        If True, `stop` is the last sample. Otherwise, it is not included.
    shape: int, tuple or list
        The shape of the final array. If None, the shape will be guessed from `num`.
    c_order: bool
        Whether to store the array in C order (row-major) or insertion order.
        Insertion order means that values will be stored in the array
        following the order of chunks in the array; this is more memory
        efficient, as it does not require an intermediate copy of the array.
        Default is True.
    **kwargs: Any
        Keyword arguments accepted by the :func:`empty` constructor.


    Returns
    -------
    out: :ref:`NDArray`
        A :ref:`NDArray` is returned.
    """

    def linspace_fill(inputs, output, offset):
        lout = len(output)
        start, stop, num, endpoint = inputs
        # if num = 1 do nothing
        step = (stop - start) / (num - 1) if endpoint and num > 1 else (stop - start) / num
        # Compute proper start and stop values for the current chunk
        # except for 0th iter, have already included start_ in prev iter
        start_ = start + offset[0] * step
        stop_ = start_ + lout * step
        if offset[0] + lout == num:  # reached end, include stop if necessary
            output[:] = np.linspace(start_, stop, lout, endpoint=endpoint, dtype=output.dtype)
        else:
            output[:] = np.linspace(start_, stop_, lout, endpoint=False, dtype=output.dtype)

    if shape is None:
        if num is None:
            raise ValueError("Either `shape` or `num` must be specified.")
        # num is not None
        shape = (num,)
    else:
        num = math.prod(shape) if num is None else num

    # check compatibility of shape and num
    if math.prod(shape) != num or num < 0:
        raise ValueError(
            f"Shape is not consistent with the specified num value {num}." + "num must be nonnegative."
            if num < 0
            else ""
        )

    if dtype is None:
        dtype = (
            blosc2.DEFAULT_COMPLEX
            if np.any([np.issubdtype(type(d), complex) for d in (start, stop)])
            else blosc2.DEFAULT_FLOAT
        )

    dtype = _check_dtype(dtype)

    if is_inside_new_expr() or num == 0:
        # We already have the dtype and shape, so return immediately
        return blosc2.zeros(shape, dtype=dtype, **kwargs)  # will return empty array for num == 0

    inputs = (start, stop, num, endpoint)
    lazyarr = blosc2.lazyudf(linspace_fill, inputs, dtype=dtype, shape=(num,))
    if len(shape) == 1:
        # C order is guaranteed, and no reshape is needed
        return lazyarr.compute(**kwargs)

    return reshape(lazyarr, shape, c_order=c_order, **kwargs)


def eye(N, M=None, k=0, dtype=np.float64, **kwargs: Any) -> NDArray:
    """Return a 2-D array with ones on the diagonal and zeros elsewhere.

    Parameters
    ----------
    N: int
        Number of rows in the output.
    M: int, optional
        Number of columns in the output. If None, defaults to `N`.
    k: int, optional
        Index of the diagonal: 0 (the default) refers to the main diagonal,
        a positive value refers to an upper diagonal, and a negative value
        to a lower diagonal.
    dtype: np.dtype or list str
        The data type of the array elements in NumPy format. Default is `np.float64`.

    Returns
    -------
    out: :ref:`NDArray`
        A :ref:`NDArray` is returned.

    Examples
    --------
    >>> import blosc2
    >>> import numpy as np
    >>> array = blosc2.eye(2, 3, dtype=np.int32)
    >>> print(array[:])
    [[1 0 0]
     [0 1 0]]
    """

    def fill_eye(inputs, output: np.array, offset: tuple):
        out_k = offset[0] - offset[1] + inputs[0]
        output[:] = np.eye(*output.shape, out_k, dtype=output.dtype)

    if M is None:
        M = N
    shape = (N, M)
    dtype = _check_dtype(dtype)

    if is_inside_new_expr():
        # We already have the dtype and shape, so return immediately
        return blosc2.zeros(shape, dtype=dtype)

    lazyarr = blosc2.lazyudf(fill_eye, (k,), dtype=dtype, shape=shape)
    return lazyarr.compute(**kwargs)


def fromiter(iterable, shape, dtype, c_order=True, **kwargs) -> NDArray:
    """Create a new array from an iterable object.

    Parameters
    ----------
    iterable: iterable
        An iterable object providing data for the array.
    shape: int, tuple or list
        The shape of the final array.
    dtype: np.dtype or list str
        The data type of the array elements in NumPy format.
    c_order: bool
        Whether to store the array in C order (row-major) or insertion order.
        Insertion order means that iterable values will be stored in the array
        following the order of chunks in the array; this is more memory
        efficient, as it does not require an intermediate copy of the array.
        Default is C order.

    Other Parameters
    ----------------
    kwargs: dict, optional
        Keyword arguments that are supported by the :func:`empty` constructor.

    Returns
    -------
    out: :ref:`NDArray`
        A :ref:`NDArray` is returned.

    Examples
    --------
    >>> import blosc2
    >>> import numpy as np
    >>> # Create an array from an iterable
    >>> array = blosc2.fromiter(range(10), shape=(10,), dtype=np.int64)
    >>> print(array[:])
    [0 1 2 3 4 5 6 7 8 9]
    """

    def iter_fill(inputs, output, offset):
        nout = math.prod(output.shape)
        (iterable,) = inputs
        output[:] = np.fromiter(iterable, dtype=output.dtype, count=nout).reshape(output.shape)

    dtype = _check_dtype(dtype)

    if is_inside_new_expr():
        # We already have the dtype and shape, so return immediately
        return blosc2.zeros(shape, dtype=dtype)

    lshape = (math.prod(shape),)
    inputs = (iterable,)
    lazyarr = blosc2.lazyudf(iter_fill, inputs, dtype=dtype, shape=lshape)

    if len(shape) == 1:
        # C order is guaranteed, and no reshape is needed
        return lazyarr.compute(**kwargs)

    # TODO: in principle, the next should work, but tests still fail:
    # return reshape(lazyarr, shape, c_order=c_order, **kwargs)
    # Creating a temporary file is a workaround for the issue
    with tempfile.NamedTemporaryFile(suffix=".b2nd", delete=True) as tmp_file:
        larr = lazyarr.compute(urlpath=tmp_file.name, mode="w")  # intermediate array
        return reshape(larr, shape, c_order=c_order, **kwargs)


def frombuffer(
    buffer: bytes, shape: int | tuple | list, dtype: np.dtype | str = np.uint8, **kwargs: Any
) -> NDArray:
    """Create an array out of a buffer.

    Parameters
    ----------
    buffer: bytes
        The buffer of the data to populate the container.
    shape: int, tuple or list
        The shape for the final container.
    dtype: np.dtype or list str
        The ndarray dtype in NumPy format. Default is `np.uint8`.
        This will override the `typesize`
        in the cparams if they are passed.

    Other Parameters
    ----------------
    kwargs: dict, optional
        Keyword arguments that are supported by the :func:`empty` constructor.

    Returns
    -------
    out: :ref:`NDArray`
        A :ref:`NDArray` is returned.

    Examples
    --------
    >>> import blosc2
    >>> import numpy as np
    >>> shape = [25, 10]
    >>> chunks = (49, 49)
    >>> dtype = np.dtype("|S8")
    >>> typesize = dtype.itemsize
    >>> # Create a buffer
    >>> buffer = bytes(np.random.normal(0, 1, np.prod(shape)) * typesize)
    >>> # Create a NDArray from a buffer with default blocks
    >>> a = blosc2.frombuffer(buffer, shape, chunks=chunks, dtype=dtype)
    """
    shape = _check_shape(shape)
    kwargs = _check_ndarray_kwargs(**kwargs)
    chunks = kwargs.pop("chunks", None)
    blocks = kwargs.pop("blocks", None)
    chunks, blocks = compute_chunks_blocks(shape, chunks, blocks, dtype, **kwargs)
    return blosc2_ext.from_buffer(buffer, shape, chunks, blocks, dtype, **kwargs)


def copy(array: NDArray, dtype: np.dtype | str = None, **kwargs: Any) -> NDArray:
    """
    This is equivalent to :meth:`NDArray.copy`

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> # Create an instance of NDArray with some data
    >>> original_array = blosc2.asarray(np.array([[1.1, 2.2, 3.3], [4.4, 5.5, 6.6]]))
    >>> # Create a copy of the array without changing dtype
    >>> copied_array = blosc2.copy(original_array)
    >>> print("Copied array (default dtype):")
    >>> print(copied_array)
    Copied array (default dtype):
    [[1.1 2.2 3.3]
    [4.4 5.5 6.6]]
    """
    return array.copy(dtype, **kwargs)


def concat(arrays: list[NDArray], /, axis=0, **kwargs: Any) -> NDArray:
    """Concatenate a list of arrays along a specified axis.

    Parameters
    ----------
    arrays: list of :ref:`NDArray`
        A list containing two or more NDArray instances to be concatenated.
    axis: int, optional
        The axis along which the arrays will be concatenated. Default is 0.

    Other Parameters
    ----------------
    kwargs: dict, optional
        Keyword arguments that are supported by the :func:`empty` constructor.

    Returns
    -------
    out: :ref:`NDArray`
        A new NDArray containing the concatenated data.

    Examples
    --------
    >>> import blosc2
    >>> import numpy as np
    >>> arr1 = blosc2.arange(0, 5, dtype=np.int32)
    >>> arr2 = blosc2.arange(5, 10, dtype=np.int32)
    >>> result = blosc2.concat([arr1, arr2])
    >>> print(result[:])
    [0 1 2 3 4 5 6 7 8 9]
    """
    if len(arrays) < 2:
        return arrays[0]
    arr1 = arrays[0]
    if not isinstance(arr1, blosc2.NDArray):
        raise TypeError("All inputs must be instances of blosc2.NDArray")
    # Do a first pass for checking array compatibility
    if axis < 0:
        axis += arr1.ndim
    if axis >= arr1.ndim:
        raise ValueError(f"Axis {axis} is out of bounds for array of dimension {arr1.ndim}.")
    for arr2 in arrays[1:]:
        if not isinstance(arr2, blosc2.NDArray):
            raise TypeError("All inputs must be instances of blosc2.NDArray")
        if arr1.ndim != arr2.ndim:
            raise ValueError("Both arrays must have the same number of dimensions for concatenation.")
        if arr1.dtype != arr2.dtype:
            raise ValueError("Both arrays must have the same dtype for concatenation.")
        # Check that the shapes match, except for the concatenation axis
        if arr1.shape[:axis] != arr2.shape[:axis] or arr1.shape[axis + 1 :] != arr2.shape[axis + 1 :]:
            raise ValueError(
                f"Shapes of the arrays do not match along the concatenation axis {axis}: "
                f"{arr1.shape} vs {arr2.shape}"
            )

    kwargs = _check_ndarray_kwargs(**kwargs)
    # Proceed with the actual concatenation
    copy = True
    # When provided urlpath coincides with an array
    mode = kwargs.pop("mode", "a")  # default mode for blosc2 is "a"
    for arr2 in arrays[1:]:
        arr1 = blosc2_ext.concat(arr1, arr2, axis, copy=copy, mode=mode, **kwargs)
        # Have now overwritten existing file (if mode ='w'), need to change mode
        # for concatenating to the same file
        mode = "r" if mode == "r" else "a"
        # arr1 is now the result of the concatenation, so we can now just enlarge it
        copy = False

    return arr1


def expand_dims(array: NDArray, axis=0) -> NDArray:
    """
    Expand the shape of an array by adding new axes at the specified positions.

    Parameters
    ----------
    array: :ref:`NDArray`
        The array to be expanded.
    axis: int or list of int, optional
        Position in the expanded axes where the new axis (or axes) is placed. Default is 0.

    Returns
    -------
    out: :ref:`NDArray`
        A new NDArray with the expanded shape.
    """
    array = blosc2.asarray(array)
    if not isinstance(array, blosc2.NDArray):
        raise TypeError("Argument array must be instance of blosc2.NDArray")
    axis = [axis] if isinstance(axis, int) else axis
    final_dims = array.ndim + len(axis)
    mask = [False for i in range(final_dims)]
    for a in axis:
        if a < 0:
            a += final_dims  # Adjust axis to be within the new stacked array's dimensions
        if mask[a]:
            raise ValueError("Axis values must be unique.")
        mask[a] = True
    return blosc2_ext.expand_dims(array, axis_mask=mask, final_dims=final_dims)


def stack(arrays: list[NDArray], axis=0, **kwargs: Any) -> NDArray:
    """Stack multiple arrays, creating a new axis.

    Parameters
    ----------
    arrays: list of :ref:`NDArray`
        A list containing two or more NDArray instances to be stacked.
    axis: int, optional
        The new axis along which the arrays will be stacked. Default is 0.

    Other Parameters
    ----------------
    kwargs: dict, optional
        Keyword arguments that are supported by the :func:`empty` constructor.

    Returns
    -------
    out: :ref:`NDArray`
        A new NDArray containing the stacked data.

    Examples
    --------
    >>> import blosc2
    >>> import numpy as np
    >>> arr1 = blosc2.arange(0, 6, dtype=np.int32, shape=(2,3))
    >>> arr2 = blosc2.arange(6, 12, dtype=np.int32, shape=(2,3))
    >>> result = blosc2.stack([arr1, arr2])
    >>> print(result.shape)
    (2, 2, 3)
    """
    if axis < 0:
        axis += arrays[0].ndim + 1  # Adjust axis to be within the new stacked array's dimensions
    newarrays = []
    for arr in arrays:
        newarrays += [blosc2.expand_dims(arr, axis=axis)]
    return blosc2.concat(newarrays, axis, **kwargs)


def save(array: NDArray, urlpath: str, contiguous=True, **kwargs: Any) -> None:
    """Save an array to a file.

    Parameters
    ----------
    array: :ref:`NDArray`
        The array to be saved.
    urlpath: str
        The path to the file where the array will be saved.
    contiguous: bool, optional
        Whether to store the array contiguously.

    Other Parameters
    ----------------
    kwargs: dict, optional
        Keyword arguments that are supported by the :func:`save` method.

    Examples
    --------
    >>> import blosc2
    >>> import numpy as np
    >>> # Create an array
    >>> array = blosc2.arange(0, 100, dtype=np.int64, shape=(10, 10))
    >>> # Save the array to a file
    >>> blosc2.save(array, "array.b2", mode="w")
    """
    array.save(urlpath, contiguous, **kwargs)


def asarray(array: Sequence | blosc2.Array, copy: bool | None = None, **kwargs: Any) -> NDArray:
    """Convert the `array` to an `NDArray`.

    Parameters
    ----------
    array: array_like
        An array supporting numpy array interface.

    copy: bool | None, optional
        Whether to copy the input. If True, the function copies.
        If False, raise a ValueError if copy is necessary. If None and
        input is NDArray, avoid copy by returning lazyexpr.
        Default: None.

    kwargs: dict, optional
        Keyword arguments that are supported by the :func:`empty` constructor.

    Returns
    -------
    out: :ref:`NDArray` or :ref:`LazyExpr`
        An new NDArray or LazyExpr made of :paramref:`array`.

    Notes
    -----
    This will create the NDArray chunk-by-chunk directly from the input array,
    without the need to create a contiguous NumPy array internally.  This can
    be used for ingesting e.g. disk or network based arrays very effectively
    and without consuming lots of memory.

    Examples
    --------
    >>> import blosc2
    >>> import numpy as np
    >>> # Create some data
    >>> shape = [25, 10]
    >>> a = np.arange(0, np.prod(shape), dtype=np.int64).reshape(shape)
    >>> # Create a NDArray from a NumPy array
    >>> nda = blosc2.asarray(a)
    """
    # Convert scalars to numpy array
    casting = kwargs.pop("casting", "unsafe")
    if casting != "unsafe":
        raise ValueError("Only unsafe casting is supported at the moment.")
    if not hasattr(array, "shape"):
        array = np.asarray(array)  # defaults if dtype=None
    dtype_ = blosc2.proxy._convert_dtype(array.dtype)
    dtype = kwargs.pop("dtype", dtype_)  # check if dtype provided
    kwargs = _check_ndarray_kwargs(**kwargs)
    chunks = kwargs.pop("chunks", None)
    blocks = kwargs.pop("blocks", None)
    # Use the chunks and blocks from the array if they are not passed
    if chunks is None and hasattr(array, "chunks"):
        chunks = array.chunks
    # Zarr adds a .blocks property that maps to a zarr.indexing.BlockIndex object
    # Let's avoid this
    if blocks is None and hasattr(array, "blocks") and isinstance(array.blocks, (tuple, list)):
        blocks = array.blocks

    copy = True if copy is None and not isinstance(array, NDArray) else copy
    if copy:
        chunks, blocks = compute_chunks_blocks(array.shape, chunks, blocks, dtype_, **kwargs)
        # Fast path for small arrays. This is not too expensive in terms of memory consumption.
        shape = array.shape
        small_size = 2**24  # 16 MB
        array_nbytes = math.prod(shape) * dtype_.itemsize
        if array_nbytes < small_size:
            if not isinstance(array, np.ndarray) and hasattr(array, "chunks"):
                # A getitem operation should be enough to get a numpy array
                array = array[()]

            array = np.require(array, dtype=dtype, requirements="C")  # require contiguous array

            return blosc2_ext.asarray(array, chunks, blocks, **kwargs)

        # Create the empty array
        ndarr = empty(shape, dtype_, chunks=chunks, blocks=blocks, **kwargs)
        behaved = are_partitions_behaved(shape, chunks, blocks)

        # Get the coordinates of the chunks
        chunks_idx, nchunks = get_chunks_idx(shape, chunks)

        # Iterate over the chunks and update the empty array
        for nchunk in range(nchunks):
            # Compute current slice coordinates
            coords = tuple(np.unravel_index(nchunk, chunks_idx))
            slice_ = tuple(
                slice(c * s, builtins.min((c + 1) * s, shape[i]))
                for i, (c, s) in enumerate(zip(coords, chunks, strict=True))
            )
            # Ensure the array slice is contiguous and of correct dtype
            array_slice = np.require(array[slice_], dtype=dtype, requirements="C")
            if behaved:
                # The whole chunk is to be updated, so this fastpath is safe
                ndarr.schunk.update_data(nchunk, array_slice, copy=False)
            else:
                ndarr[slice_] = array_slice
    else:
        if not isinstance(array, NDArray):
            raise ValueError("Must always do a copy for asarray unless NDArray provided.")
        # TODO: make a direct view possible
        return array

    return ndarr


def astype(
    array: Sequence | blosc2.Array,
    dtype,
    casting: str = "unsafe",
    copy: bool = True,
    **kwargs: Any,
) -> NDArray:
    """
    Copy of the array, cast to a specified type. Does not support copy = False.

    Parameters
    ----------
    array: Sequence | blosc2.Array
        The array to be cast to a different type.
    dtype: DType-like
        The desired data type to cast to.
    casting: str = 'unsafe'
        Controls what kind of data casting may occur. Defaults to 'unsafe' for backwards compatibility.
        * 'no' means the data types should not be cast at all.
        * 'equiv' means only byte-order changes are allowed.
        * 'safe' means only casts which can preserve values are allowed.
        * 'same_kind' means only safe casts or casts within a kind, like float64 to float32, are allowed.
        * 'unsafe' means any data conversions may be done.
    copy: bool = True
        Must always be True as copy is made by default. Will be changed in a future version

    Returns
    -------
    out: NDArray
        New array with specified data type.
    """
    return asarray(array, dtype=dtype, casting=casting, copy=copy, **kwargs)


def _check_ndarray_kwargs(**kwargs):
    storage = kwargs.get("storage")
    if storage is not None:
        for key in kwargs:
            if key in list(blosc2.Storage.__annotations__):
                raise AttributeError(
                    "Cannot pass both `storage` and other kwargs already included in Storage"
                )
        if isinstance(storage, blosc2.Storage):
            kwargs = {**kwargs, **asdict(storage)}
        else:
            kwargs = {**kwargs, **storage}
    else:
        # Add the default storage values as long as they are not already passed
        storage_dflts = asdict(blosc2.Storage(urlpath=kwargs.get("urlpath")))  # urlpath can affect defaults
        # If a key appears in both operands, the one from the right-hand operand wins
        kwargs = storage_dflts | kwargs

    supported_keys = [
        "chunks",
        "blocks",
        "cparams",
        "dparams",
        "meta",
        "urlpath",
        "contiguous",
        "mode",
        "mmap_mode",
        "initial_mapping_size",
        "storage",
        "out",
        "_chunksize_reduc_factor",
    ]
    _ = kwargs.pop("device", None)  # pop device (not used, but needs to be discarded)
    for key in kwargs:
        if key not in supported_keys:
            raise KeyError(
                f"Only {supported_keys} are supported as keyword arguments, and you passed '{key}'"
            )

    if "cparams" in kwargs:
        cparams = kwargs["cparams"]
        if cparams is None:
            kwargs["cparams"] = blosc2.cparams_dflts
        if isinstance(cparams, blosc2.CParams):
            kwargs["cparams"] = asdict(kwargs["cparams"])
        else:
            if "chunks" in kwargs["cparams"]:
                raise ValueError("You cannot pass chunks in cparams, use `chunks` argument instead")
            if "blocks" in kwargs["cparams"]:
                raise ValueError("You cannot pass chunks in cparams, use `blocks` argument instead")
    if "dparams" in kwargs and isinstance(kwargs["dparams"], blosc2.DParams):
        kwargs["dparams"] = asdict(kwargs["dparams"])

    return kwargs


def get_slice_nchunks(
    schunk: blosc2.SChunk, key: tuple[(int, int)] | int | slice | Sequence[slice]
) -> np.ndarray:
    """
    Get the unidimensional chunk indexes needed to obtain a
    slice of a :ref:`SChunk <SChunk>` or a :ref:`NDArray`.

    Parameters
    ----------
    schunk: :ref:`SChunk <SChunk>` or :ref:`NDArray`
        The super-chunk or ndarray container.
    key: tuple(int, int), int, slice or sequence of slices
        For a SChunk: a tuple with the start and stop of the slice, an integer,
        or a single slice. For a ndarray, sequences of slices (one per dimension) are accepted.

    Returns
    -------
    out: np.ndarray
        An array with the unidimensional chunk indexes.
    """
    if isinstance(schunk, NDArray):
        array = schunk
        key, _ = process_key(key, array.shape)
        start, stop, step, _ = get_ndarray_start_stop(array.ndim, key, array.shape)
        if step != (1,) * array.ndim:
            raise IndexError("Step parameter is not supported yet")
        key = (start, stop)
        return blosc2_ext.array_get_slice_nchunks(array, key)
    else:
        if isinstance(key, int):
            key = (key, key + 1)
        elif isinstance(key, slice):
            if key.step not in (1, None):
                raise IndexError("Only step=1 is supported")
            key = (key.start, key.stop)
        return blosc2_ext.schunk_get_slice_nchunks(schunk, key)


def indices(array: blosc2.Array, order: str | list[str] | None = None, **kwargs: Any) -> NDArray:
    """
    Return the indices of a sorted array following the specified order.

    This is only valid for 1-dim structured arrays.

    Parameters
    ----------
    array: :ref:`blosc2.Array`
        The (structured) array to be sorted.
    order: str, list of str, optional
        Specifies which fields to compare first, second, etc. A single
        field can be specified as a string. Not all fields need to be
        specified, only the ones by which the array is to be sorted.
        If None, the array is not sorted.
    kwargs: Any, optional
        Keyword arguments that are supported by the :func:`empty` constructor.

    Returns
    -------
    out: :ref:`NDArray`
        The sorted array.
    """
    if not order:
        # Shortcut for this relatively rare case
        return arange(array.shape[0], dtype=np.int64)

    # Create a lazy array to access the sort machinery there
    # This is a bit of a hack, but it is the simplest way to do it
    # (the sorting mechanism in LazyExpr should be improved to avoid this)
    lbool = blosc2.lazyexpr(blosc2.ones(array.shape, dtype=np.bool_))
    larr = array[lbool]
    return larr.indices(order).compute(**kwargs)


def sort(array: blosc2.Array, order: str | list[str] | None = None, **kwargs: Any) -> NDArray:
    """
    Return a sorted array following the specified order.

    This is only valid for 1-dim structured arrays.

    Parameters
    ----------
    array: :ref:`blosc2.Array`
        The (structured) array to be sorted.
    order: str, list of str, optional
        Specifies which fields to compare first, second, etc. A single
        field can be specified as a string. Not all fields need to be
        specified, only the ones by which the array is to be sorted.
    kwargs: Any, optional
        Keyword arguments that are supported by the :func:`empty` constructor.

    Returns
    -------
    out: :ref:`NDArray`
        The sorted array.
    """
    if not order:
        return array

    # Create a lazy array to access the sort machinery there
    # This is a bit of a hack, but it is the simplest way to do it
    # (the sorting mechanism in LazyExpr should be improved to avoid this)
    lbool = blosc2.lazyexpr(blosc2.ones(array.shape, dtype=np.bool_))
    larr = array[lbool]
    return larr.sort(order).compute(**kwargs)


# Class for dealing with fields in an NDArray
# This will allow to access fields by name in the dtype of the NDArray
class NDField(Operand):
    def __init__(self, ndarr: NDArray, field: str):
        """
        Create a new NDField.

        Parameters
        ----------
        ndarr: :ref:`NDArray`
            The NDArray to which assign the field.
        field: str
            The field's name.

        Returns
        -------
        out: :ref:`NDField`
            The corresponding :ref:`NDField`.
        """
        if not isinstance(ndarr, NDArray):
            raise TypeError("ndarr should be a NDArray!")
        if not isinstance(field, str):
            raise TypeError("field should be a string!")
        if ndarr.dtype.fields is None:
            raise TypeError("NDArray does not have a structured dtype!")
        if field not in ndarr.dtype.fields:
            raise TypeError(f"Field {field} not found in the dtype of the NDArray")
        # Store immutable properties
        self.ndarr = ndarr
        self.chunks = ndarr.chunks
        self.blocks = ndarr.blocks
        self.field = field
        self._dtype = ndarr.dtype.fields[field][0]
        self.offset = ndarr.dtype.fields[field][1]

    def __repr__(self):
        """
        Get a string as a representation.

        Returns
        -------
        out: str
        """
        return f"NDField({self.ndarr}, {self.field})"

    @property
    def shape(self) -> tuple[int]:
        """The shape of the associated :ref:`NDArray`."""
        return self.ndarr.shape

    @property
    def dtype(self) -> np.dtype:
        """The dtype of the field of associated :ref:`NDArray`."""
        return self._dtype

    @property
    def schunk(self) -> blosc2.SChunk:
        """The associated :ref:`SChunk <SChunk>`."""
        return self.ndarr.schunk

    def __getitem__(self, key: int | slice | Sequence[slice]) -> np.ndarray:
        """
        Get a slice of :paramref:`self`.

        Parameters
        ----------
        key: int or slice or Sequence[slice]
            The slice to be retrieved.

        Returns
        -------
        out: NumPy.ndarray
            A NumPy array with the data slice.

        """
        # If key is a LazyExpr, decorate it with ``where`` and return it
        if isinstance(key, blosc2.LazyExpr):
            return key.where(self)

        if isinstance(key, str):
            # Try to compute the key as a boolean expression
            # Operands will be a dict with all the fields in the NDArray
            operands = {field: NDField(self.ndarr, field) for field in self.ndarr.dtype.names}
            expr = blosc2.lazyexpr(key, operands)
            if expr.dtype != np.bool_:
                raise TypeError("The expression should return a boolean array")
            return expr.where(self)
            # raise TypeError("This array is a NDField; use a structured NDArray for bool expressions")

        # Check if the key is in the last read cache
        inmutable_key = make_key_hashable(key)
        if inmutable_key in self.ndarr._last_read:
            return self.ndarr._last_read[inmutable_key][self.field]

        # Do the actual read in the parent NDArray
        nparr = self.ndarr[key]
        # And return the field
        return nparr[self.field]

    def __setitem__(self, key: int | slice | Sequence[slice], value: blosc2.Array) -> None:
        """
        Set a slice of :paramref:`self` to a value.

        Parameters
        ----------
        key: int or slice or Sequence[slice]
            The slice to be set.
        value: blosc2.Array
            The value to be set.
        """
        if isinstance(key, str):
            raise TypeError("This array is a NDField; use a structured NDArray for bool expressions")
        if not isinstance(value, np.ndarray):
            value = value[:]
        # Get the values in the parent NDArray
        nparr = self.ndarr[key]
        # Set the field
        nparr[self.field] = value
        # Save the values in the parent NDArray
        self.ndarr[key] = nparr

    def __iter__(self):
        """
        Iterate over the elements in the field.

        Returns
        -------
        out: iterator
        """
        return NDOuterIterator(self)

    def __len__(self) -> int:
        """
        Returns the length of the first dimension of the field.
        """
        return self.shape[0]


class OIndex:
    def __init__(self, array: NDArray):
        self.array = array

    def __getitem__(self, selection) -> np.ndarray:
        return self.array.get_oselection_numpy(selection)

    def __setitem__(self, selection, input) -> np.ndarray:
        return self.array.set_oselection_numpy(selection, input)


# class VIndex:
#     def __init__(self, array: NDArray):
#         self.array = array

#     # TODO: all this
#     def __getitem__(self, selection) -> np.ndarray:
#         return NotImplementedError

#     def __setitem__(self, selection, input) -> np.ndarray:
#         return NotImplementedError


def empty_like(x: blosc2.Array, dtype=None, **kwargs) -> NDArray:
    """
    Returns an uninitialized array with the same shape as an input array x.

    Parameters
    ----------
    x : blosc2.Array
        Input array from which to derive the output array shape.

    dtype (Optional):
        Output array data type. If dtype is None, the output array data type
        is inferred from x. Default: None.

    kwargs: Any, optional
        Keyword arguments that are supported by the :func:`empty` constructor.
        These arguments will be set in the resulting :ref:`NDArray`.

    Returns
    ------
    out : NDArray
        An array having the same shape as x and containing uninitialized data.
    """
    if dtype is None:
        dtype = x.dtype
    return blosc2.empty(shape=x.shape, dtype=dtype, **kwargs)


def ones_like(x: blosc2.Array, dtype=None, **kwargs) -> NDArray:
    """
    Returns an array of ones with the same shape as an input array x.

    Parameters
    ----------
    x : blosc2.Array
        Input array from which to derive the output array shape.

    dtype (Optional):
        Output array data type. If dtype is None, the output array data type
        is inferred from x. Default: None.

    kwargs: Any, optional
        Keyword arguments that are supported by the :func:`empty` constructor.
        These arguments will be set in the resulting :ref:`NDArray`.

    Returns
    ------
    out : NDArray
        An array having the same shape as x and containing ones.
    """
    if dtype is None:
        dtype = x.dtype
    return blosc2.ones(shape=x.shape, dtype=dtype, **kwargs)


def zeros_like(x: blosc2.Array, dtype=None, **kwargs) -> NDArray:
    """
    Returns an array of zeros with the same shape as an input array x.

    Parameters
    ----------
    x : blosc2.Array
        Input array from which to derive the output array shape.

    dtype (Optional):
        Output array data type. If dtype is None, the output array data type
        is inferred from x. Default: None.

    kwargs: Any, optional
        Keyword arguments that are supported by the :func:`empty` constructor.
        These arguments will be set in the resulting :ref:`NDArray`.

    Returns
    ------
    out : NDArray
        An array having the same shape as x and containing zeros.
    """
    if dtype is None:
        dtype = x.dtype
    return blosc2.zeros(shape=x.shape, dtype=dtype, **kwargs)


def full_like(x: blosc2.Array, fill_value: bool | int | float | complex, dtype=None, **kwargs) -> NDArray:
    """
    Returns an array filled with a value with the same shape as an input array x.

    Parameters
    ----------
    x : blosc2.Array
        Input array from which to derive the output array shape.

    fill_value: bool | int | float | complex
        The fill value.

    dtype (Optional):
        Output array data type. If dtype is None, the output array data type
        is inferred from x. Default: None.

    kwargs: Any, optional
        Keyword arguments that are supported by the :func:`empty` constructor.
        These arguments will be set in the resulting :ref:`NDArray`.

    Returns
    ------
    out : NDArray
        An array having the same shape as x and containing the fill value.
    """
    if dtype is None:
        dtype = x.dtype
    return blosc2.full(shape=x.shape, fill_value=fill_value, dtype=dtype, **kwargs)


def take(x: blosc2.Array, indices: blosc2.Array, axis: int | None = None) -> NDArray:
    """
    Returns elements of an array along an axis.

    Parameters
    ----------
    x: blosc2.Array
        Input array. Should have one or more dimensions (axes).

    indices: array-like
        Array indices. The array must be one-dimensional and have an integer data type.

    axis: int | None
        Axis over which to select values.
        If x is a one-dimensional array, providing an axis is optional; however, if x
        has more than one dimension, providing an axis is required. Default: None.

    Returns
    -------
    out: NDArray
        Selected indices of x.
    """
    if axis is None:
        axis = 0
        if x.ndim != 1:
            raise ValueError("Must specify axis parameter if x is not 1D.")
    if axis < 0:
        axis += x.ndim
    if not isinstance(axis, (int, np.integer)):
        raise ValueError("Axis must be integer.")
    if isinstance(indices, list):
        indices = np.asarray(indices)
    if indices.ndim != 1:
        raise ValueError("Indices must be 1D array.")
    key = tuple(indices if i == axis else slice(None, None, 1) for i in range(x.ndim))
    # TODO: Implement fancy indexing in .slice so that this is more efficient
    return blosc2.asarray(x[key])


def take_along_axis(x: blosc2.Array, indices: blosc2.Array, axis: int = -1) -> NDArray:
    """
    Returns elements of an array along an axis.

    Parameters
    ----------
    x: blosc2.Array
        Input array. Should have one or more dimensions (axes).

    indices: array-like
        Array indices. The array must have same number of dimensions as x and
        have an integer data type.

    axis: int
        Axis over which to select values. Default: -1.

    Returns
    -------
    out: NDArray
        Selected indices of x.
    """
    if not isinstance(axis, (int, np.integer)):
        raise ValueError("Axis must be integer.")
    if indices.ndim != x.ndim:
        raise ValueError("Indices must have same dimensions as x.")
    if axis < 0:
        axis += x.ndim
    if indices.shape[axis] == 0:
        return blosc2.empty(x.shape[:axis] + (0,) + x.shape[axis + 1 :], dtype=x.dtype)
    ones = (1,) * x.ndim
    # TODO: Implement fancy indexing in .slice so that this is more efficient and possibly use oindex(?)
    key = tuple(
        indices if i == axis else np.arange(x.shape[i]).reshape(ones[:i] + (-1,) + ones[i + 1 :])
        for i in range(x.ndim)
    )
    return blosc2.asarray(x[key])


def broadcast_to(arr: blosc2.Array, shape: tuple[int, ...]) -> NDArray:
    """
    Broadcast an array to a new shape.
    Warning: Computes a lazyexpr, so probably a bit suboptimal

    Parameters
    ----------
    arr: blosc2.Array
        The array to broadcast.

    shape: tuple
        The shape of the desired array.

    Returns
    -------
    broadcast: NDArray
    A new array with the given shape.
    """
    return (arr + blosc2.zeros(shape, dtype=arr.dtype)).compute()


def meshgrid(*arrays: blosc2.Array, indexing: str = "xy") -> Sequence[NDArray]:
    """
    Returns coordinate matrices from coordinate vectors.

    Parameters
    ----------
    *arrays: blosc2.Array
        An arbitrary number of one-dimensional arrays representing grid coordinates. Each array should have the same numeric data type.

    indexing: str
        Cartesian 'xy' or matrix 'ij' indexing of output. If provided zero or one one-dimensional vector(s) the indexing keyword is ignored.
        Default: 'xy'.

    Returns
    -------
    out: (List[NDArray])
        List of N arrays, where N is the number of provided one-dimensional input arrays, with same dtype.
        For N one-dimensional arrays having lengths Ni = len(xi),

        * if matrix indexing ij, then each returned array has shape (N1, N2, N3, ..., Nn).
        * if Cartesian indexing xy, then each returned array has shape (N2, N1, N3, ..., Nn).
    """
    out = ()
    shape = np.ones(len(arrays))
    first_arr = arrays[0]
    myarrs = ()
    if indexing == "xy" and len(shape) > 1:
        # switch 0th and 1st shapes around
        def mygen(i):
            if i not in (0, 1):
                return (j for j in range(len(arrays)) if j != i)
            else:
                return (j for j in range(len(arrays)) if j != builtins.abs(i - 1))
    else:
        mygen = lambda i: (j for j in range(len(arrays)) if j != i)  # noqa : E731

    for i, a in enumerate(arrays):
        if len(a.shape) != 1 or a.dtype != first_arr.dtype:
            raise ValueError("All arrays must be 1D and of same dtype.")
        shape[i] = a.shape[0]
        myarrs += (blosc2.expand_dims(a, tuple(mygen(i))),)  # cheap, creates a view

    # handle Cartesian indexing
    shape = tuple(shape)
    if indexing == "xy" and len(shape) > 1:
        shape = (shape[1], shape[0]) + shape[2:]

    # do broadcast
    for a in myarrs:
        out += (broadcast_to(a, shape),)
    return out
