#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# This source code is licensed under a BSD-style license (found in the
# LICENSE file in the root directory of this source tree)
#######################################################################

from __future__ import annotations

import builtins
import inspect
import math
from collections import OrderedDict, namedtuple
from functools import reduce
from typing import TYPE_CHECKING, Any, NamedTuple

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


def process_key(key, shape):
    if key is None:
        key = tuple(slice(None) for _ in range(len(shape)))
    key = ndindex.ndindex(key).expand(shape).raw
    mask = tuple(isinstance(k, int) for k in key)
    key = tuple(k if isinstance(k, slice) else slice(k, k + 1, None) for k in key)
    return key, mask


def get_ndarray_start_stop(ndim, key, shape):
    start = [s.start if s.start is not None else 0 for s in key]
    stop = [s.stop if s.stop is not None else sh for s, sh in zip(key, shape, strict=False)]
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
    step = tuple(s.step if s.step is not None else 1 for s in key)
    return start, stop, step


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


def get_chunks_idx(shape, chunks):
    chunks_idx = tuple(math.ceil(s / c) for s, c in zip(shape, chunks, strict=True))
    nchunks = math.prod(chunks_idx)
    return chunks_idx, nchunks


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
    src: NDArray | NDField | blosc2.LazyArray | blosc2.C2Array,
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

    if is_inside_new_expr():
        # We already have the dtype and shape, so return immediately
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
            dst_buf[dst_buf_slice] = src[src_slice]
        # Compute the shape of dst_slice
        dst_slice_shape = tuple(s.stop - s.start for s in dst_slice)
        # ... and assign the buffer to the destination array
        dst[dst_slice] = dst_buf.reshape(dst_slice_shape)

    return dst


def _check_allowed_dtypes(
    value: bool | int | float | str | blosc2.NDArray | blosc2.NDField | blosc2.C2Array | blosc2.Proxy,
):
    if not (
        isinstance(
            value,
            blosc2.LazyExpr
            | blosc2.NDArray
            | blosc2.NDField
            | blosc2.C2Array
            | blosc2.Proxy
            | blosc2.ProxyNDField
            | blosc2.SimpleProxy
            | np.ndarray,
        )
        or np.isscalar(value)
    ):
        raise RuntimeError(
            "Expected LazyExpr, NDArray, NDField, C2Array, Proxy, np.ndarray or scalar instances"
            f" and you provided a '{type(value)}' instance"
        )


def sum(
    ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr,
    axis: int | tuple[int] | None = None,
    dtype: np.dtype | str = None,
    keepdims: bool = False,
    **kwargs: Any,
) -> np.ndarray | NDArray | int | float | complex | bool:
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
    ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr,
    axis: int | tuple[int] | None = None,
    dtype: np.dtype | str = None,
    keepdims: bool = False,
    **kwargs: Any,
) -> np.ndarray | NDArray | int | float | complex | bool:
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
    ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr,
    axis: int | tuple[int] | None = None,
    dtype: np.dtype | str = None,
    ddof: int = 0,
    keepdims: bool = False,
    **kwargs: Any,
) -> np.ndarray | NDArray | int | float | bool:
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
    ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr,
    axis: int | tuple[int] | None = None,
    dtype: np.dtype | str = None,
    ddof: int = 0,
    keepdims: bool = False,
    **kwargs: Any,
) -> np.ndarray | NDArray | int | float | bool:
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
    ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr,
    axis: int | tuple[int] | None = None,
    dtype: np.dtype | str = None,
    keepdims: bool = False,
    **kwargs: Any,
) -> np.ndarray | NDArray | int | float | complex | bool:
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
    ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr,
    axis: int | tuple[int] | None = None,
    keepdims: bool = False,
    **kwargs: Any,
) -> np.ndarray | NDArray | int | float | complex | bool:
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
    ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr,
    axis: int | tuple[int] | None = None,
    keepdims: bool = False,
    **kwargs: Any,
) -> np.ndarray | NDArray | int | float | complex | bool:
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
    ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr,
    axis: int | tuple[int] | None = None,
    keepdims: bool = False,
    **kwargs: Any,
) -> np.ndarray | NDArray | bool:
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


def all(
    ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr,
    axis: int | tuple[int] | None = None,
    keepdims: bool = False,
    **kwargs: Any,
) -> np.ndarray | NDArray | bool:
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


class Operand:
    """Base class for all operands in expressions."""

    def __neg__(self) -> blosc2.LazyExpr:
        return blosc2.LazyExpr(new_op=(0, "-", self))

    def __and__(self, value: int | float | NDArray | NDField | blosc2.C2Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "&", value))

    def __add__(self, value: int | float | NDArray | NDField | blosc2.C2Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "+", value))

    def __iadd__(self, value: int | float | NDArray | NDField | blosc2.C2Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "+", value))

    # Provide minimal __array_interface__ to allow NumPy to work with this object
    @property
    def __array_interface__(self):
        return {
            "shape": self.shape,
            "typestr": self.dtype.str,
            "data": self[()],
            "version": 3,
        }

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        # Handle operations at the array level
        if method != "__call__":
            return NotImplemented

        ufunc_map = {
            np.add: "+",
            np.subtract: "-",
            np.multiply: "*",
            np.divide: "/",
            np.true_divide: "/",
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
        }

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
            np.abs: "abs",
            np.conj: "conj",
            np.real: "real",
            np.imag: "imag",
            np.bitwise_not: "~",
        }

        if ufunc in ufunc_map:
            value = inputs[0] if inputs[1] is self else inputs[1]
            _check_allowed_dtypes(value)
            return blosc2.LazyExpr(new_op=(value, ufunc_map[ufunc], self))

        if ufunc in ufunc_map_1param:
            value = inputs[0]
            _check_allowed_dtypes(value)
            return blosc2.LazyExpr(new_op=(value, ufunc_map_1param[ufunc], None))

        return NotImplemented

    def __radd__(self, value: int | float | NDArray | NDField | blosc2.C2Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(value, "+", self))

    def __sub__(self, value: int | float | NDArray | NDField | blosc2.C2Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "-", value))

    def __isub__(self, value: int | float | NDArray | NDField | blosc2.C2Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "-", value))

    def __rsub__(self, value: int | float | NDArray | NDField | blosc2.C2Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(value, "-", self))

    def __mul__(self, value: int | float | NDArray | NDField | blosc2.C2Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "*", value))

    def __imul__(self, value: int | float | NDArray | NDField | blosc2.C2Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "*", value))

    def __rmul__(self, value: int | float | NDArray | NDField | blosc2.C2Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(value, "*", self))

    def __truediv__(self, value: int | float | NDArray | NDField | blosc2.C2Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "/", value))

    def __itruediv__(self, value: int | float | NDArray | NDField | blosc2.C2Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "/", value))

    def __rtruediv__(self, value: int | float | NDArray | NDField | blosc2.C2Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(value, "/", self))

    def __lt__(self, value: int | float | NDArray | NDField | blosc2.C2Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "<", value))

    def __le__(self, value: int | float | NDArray | NDField | blosc2.C2Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "<=", value))

    def __gt__(self, value: int | float | NDArray | NDField | blosc2.C2Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, ">", value))

    def __ge__(self, value: int | float | NDArray | NDField | blosc2.C2Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, ">=", value))

    def __eq__(self, value: int | float | NDArray | NDField | blosc2.C2Array, /):
        _check_allowed_dtypes(value)
        if blosc2._disable_overloaded_equal:
            return self is value
        return blosc2.LazyExpr(new_op=(self, "==", value))

    def __ne__(self, value: int | float | NDArray | NDField | blosc2.C2Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "!=", value))

    def __pow__(self, value: int | float | NDArray | NDField | blosc2.C2Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "**", value))

    def __ipow__(self, value: int | float | NDArray | NDField | blosc2.C2Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(self, "**", value))

    def __rpow__(self, value: int | float | NDArray | NDField | blosc2.C2Array, /) -> blosc2.LazyExpr:
        _check_allowed_dtypes(value)
        return blosc2.LazyExpr(new_op=(value, "**", self))

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


def extract_values(arr, indices: np.ndarray[np.int_], max_cache_size: int = 10) -> np.ndarray:
    """
    Extract values from a chunked and compressed array using an array of indices.

    Parameters
    ----------
    arr : blosc2.NDArray
        The chunked and compressed array.
    indices : np.ndarray
        The array of indices to extract values from.
    max_cache_size : int
        The maximum number of chunks to cache.

    Returns
    -------
    extracted_values : np.ndarray
        The extracted values.
    """
    # Initialize the result array
    extracted_values = np.empty(len(indices), dtype=arr.dtype)

    # Limited size dictionary to store decompressed chunks
    chunk_cache = LimitedSizeDict(max_cache_size)

    # Iterate through the indices and extract values
    chunk_size = int(arr.chunks[0])
    for i, idx in enumerate(indices):
        chunk_idx = idx // chunk_size
        if chunk_idx not in chunk_cache:
            # Compute the bounds for this chunk
            start = chunk_idx * chunk_size
            end = start + chunk_size
            chunk_cache[chunk_idx] = arr[start:end]

        local_idx = idx % chunk_size
        extracted_values[i] = chunk_cache[chunk_idx][local_idx]

    return extracted_values


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
        super().__init__(kwargs["_array"])
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
    def storage(self) -> blosc2.Storage:
        """The storage of the array."""
        return self.schunk.storage

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
        items = []
        items += [("type", f"{self.__class__.__name__}")]
        items += [("shape", self.shape)]
        items += [("chunks", self.chunks)]
        items += [("blocks", self.blocks)]
        items += [("dtype", self.dtype)]
        items += [("cratio", f"{self.schunk.cratio:.2f}")]
        items += [("cparams", self.schunk.cparams)]
        items += [("dparams", self.schunk.dparams)]
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
        """The size (in bytes) for this container."""
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

    def __getitem__(  # noqa: C901
        self,
        key: int | slice | Sequence[slice | int] | np.ndarray[np.bool_] | NDArray | blosc2.LazyExpr | str,
    ) -> np.ndarray | blosc2.LazyExpr:
        """Retrieve a (multidimensional) slice as specified by the key.

        Parameters
        ----------
        key: int, slice, sequence of (slices, int), array of bools, LazyExpr or str
            The slice(s) to be retrieved. Note that step parameter is not yet honored
            in slices. If a LazyExpr is provided, the expression is expected to be of
            boolean type, and the result will be another LazyExpr returning the values
            of this array where the expression is True.
            When key is a (nd-)array of bools, the result will be the values of ``self``
            where the bool values are True (similar to NumPy).
            If key is a 1-dim sequence of integers, the result will be the values of
            this array at the specified indices. N-dim indices are not yet supported.
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
        # First try some fast paths for common cases
        if isinstance(key, np.integer):
            # Massage the key to a tuple and go the fast path
            key_ = (slice(key, key + 1), *(slice(None),) * (self.ndim - 1))
            start, stop, step = get_ndarray_start_stop(self.ndim, key_, self.shape)
            shape = tuple(sp - st for st, sp in zip(start, stop, strict=True))
        elif isinstance(key, tuple) and (
            builtins.sum(isinstance(k, builtins.slice) for k in key) == self.ndim
        ):
            # This can be processed in a fast way already
            start, stop, step = get_ndarray_start_stop(self.ndim, key, self.shape)
            shape = tuple(sp - st for st, sp in zip(start, stop, strict=True))
        elif isinstance(key, (list, np.ndarray, NDArray)):
            if isinstance(key, list):
                key = np.array(key, dtype=np.int64)
            if np.issubdtype(key.dtype, np.bool_):
                # This can be interpreted as a boolean expression
                if key.shape != self.shape:
                    raise ValueError("The shape of the boolean expression should match the array shape")
                # expr = blosc2.lazyexpr(f"(key)")
                # The next should be a bit faster
                expr = blosc2.LazyExpr._new_expr("key", {"key": key}, guess=False)
                # Decorate with where and force a getitem operation to return actual values.
                # This behavior is consistent with NumPy, although different from e.g. ['expr']
                # which returns a lazy expression.
                return expr.where(self)[:]
            if key.dtype != np.int64 or key.ndim != 1:
                raise ValueError("Only 1-dim sequences of ints are supported for now")
            if isinstance(key, NDArray):
                key = key[:]
            # This is a fast path for the case where the key is a short list of 1-d indices
            # For the rest, use an array of booleans.
            return extract_values(self, key)
        else:
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
            key_, mask = process_key(key, self.shape)
            start, stop, step = get_ndarray_start_stop(self.ndim, key_, self.shape)
            shape = np.array([sp - st for st, sp in zip(start, stop, strict=True)])
            shape = tuple(shape[[not m for m in mask]])

        # Create the array to store the result
        arr = np.empty(shape, dtype=self.dtype)
        nparr = super().get_slice_numpy(arr, (start, stop))
        if step != (1,) * self.ndim:
            if len(step) == 1:
                return nparr[:: step[0]]
            slice_ = tuple(slice(None, None, st) for st in step)
            return nparr[slice_]

        if self._keep_last_read:
            self._last_read.clear()
            inmutable_key = make_key_hashable(key)
            self._last_read[inmutable_key] = nparr

        return nparr

    def __setitem__(self, key: int | slice | Sequence[slice], value: object):
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
        key, _ = process_key(key, self.shape)
        start, stop, step = get_ndarray_start_stop(self.ndim, key, self.shape)
        if step != (1,) * self.ndim:
            raise ValueError("Step parameter is not supported yet")
        key = (start, stop)

        shape = [sp - st for sp, st in zip(stop, start, strict=False)]
        if isinstance(value, int | float | bool):
            value = np.full(shape, value, dtype=self.dtype)
        elif isinstance(value, np.ndarray):
            if value.dtype != self.dtype:
                try:
                    value = value.astype(self.dtype)
                except ComplexWarning:
                    # numexpr type inference can lead to unnecessary type promotions
                    # when using complex functions (e.g. conj) with real arrays
                    value = value.real.astype(self.dtype)
            if value.shape == ():
                value = np.full(shape, value, dtype=self.dtype)
        elif isinstance(value, NDArray):
            value = value[...]

        return super().set_slice(key, value)

    def __iter__(self):
        """Iterate over the (outer) elements of the array.

        Returns
        -------
        out: iterator
        """
        return NDOuterIterator(self)

    def __len__(self):
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
        """
        kwargs = _check_ndarray_kwargs(**kwargs)
        key, mask = process_key(key, self.shape)
        start, stop, step = get_ndarray_start_stop(self.ndim, key, self.shape)
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

    def squeeze(self) -> NDArray:
        """Remove single-dimensional entries from the shape of the array.

        This method modifies the array in-place, removing any dimensions with size 1.

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
        super().squeeze()
        return self

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


def sin(ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr, /) -> blosc2.LazyExpr:
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


def cos(ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr, /) -> blosc2.LazyExpr:
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


def tan(ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr, /) -> blosc2.LazyExpr:
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


def sqrt(ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr, /) -> blosc2.LazyExpr:
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


def sinh(ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr, /) -> blosc2.LazyExpr:
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


def cosh(ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr, /) -> blosc2.LazyExpr:
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


def tanh(ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr, /) -> blosc2.LazyExpr:
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


def arcsin(ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr, /) -> blosc2.LazyExpr:
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


def arccos(ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr, /) -> blosc2.LazyExpr:
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


def arctan(ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr, /) -> blosc2.LazyExpr:
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


def arctan2(
    ndarr1: NDArray | NDField | blosc2.C2Array, ndarr2: NDArray | NDField | blosc2.C2Array, /
) -> blosc2.LazyExpr:
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


def arcsinh(ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr, /) -> blosc2.LazyExpr:
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


def arccosh(ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr, /) -> blosc2.LazyExpr:
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


def arctanh(ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr, /) -> blosc2.LazyExpr:
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


def exp(ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr, /) -> blosc2.LazyExpr:
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


def expm1(ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr, /) -> blosc2.LazyExpr:
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


def log(ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr, /) -> blosc2.LazyExpr:
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


def log10(ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr, /) -> blosc2.LazyExpr:
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


def log1p(ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr, /) -> blosc2.LazyExpr:
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


def conj(ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr, /) -> blosc2.LazyExpr:
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


def real(ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr, /) -> blosc2.LazyExpr:
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


def imag(ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr, /) -> blosc2.LazyExpr:
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


def contains(
    ndarr: NDArray | NDField | blosc2.C2Array, value: str | bytes | NDArray | NDField | blosc2.C2Array, /
) -> blosc2.LazyExpr:
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


def abs(ndarr: NDArray | NDField | blosc2.C2Array | blosc2.LazyExpr, /) -> blosc2.LazyExpr:
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


def where(
    condition: blosc2.LazyExpr,
    x: NDArray | NDField | np.ndarray | int | float | complex | bool | str | bytes | None = None,
    y: NDArray | NDField | np.ndarray | int | float | complex | bool | str | bytes | None = None,
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


def ones(shape: int | tuple | list, dtype: np.dtype | str = np.float64, **kwargs: Any) -> NDArray:
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
    return full(shape, 1, dtype, **kwargs)


def arange(
    start: int | float = 0,
    stop: int | float | None = None,
    step: int | float | None = 1,
    dtype: np.dtype | str = np.int64,
    shape: int | tuple | list | None = None,
    c_order: bool = True,
    **kwargs: Any,
) -> NDArray:
    """Return evenly spaced values within a given interval.

    Parameters
    ----------
    start: int, float, complex or np.number
        The starting value of the sequence.
    stop: int, float, complex or np.number
        The end value of the sequence.
    step: int, float, complex or np.number
        Spacing between values.
    dtype: np.dtype or list str
        The data type of the array elements in NumPy format. Default is `np.uint8`.
        This will override the `typesize`
        in the compression parameters if they are provided.
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
        output[:] = np.arange(start, stop, step, dtype=output.dtype)

    if stop is None:
        stop = start
        start = 0
    if step is None:
        step = 1
    if not shape:
        shape = (int((stop - start) / step),)
    else:
        # Check that the shape is consistent with the start, stop and step values
        if math.prod(shape) != int((stop - start) / step):
            raise ValueError("The shape is not consistent with the start, stop and step values")
    dtype = _check_dtype(dtype)

    if is_inside_new_expr():
        # We already have the dtype and shape, so return immediately
        return blosc2.zeros(shape, dtype=dtype)

    lshape = (math.prod(shape),)
    lazyarr = blosc2.lazyudf(arange_fill, (start, stop, step), dtype=dtype, shape=lshape)

    if len(shape) == 1:
        # C order is guaranteed, and no reshape is needed
        return lazyarr.compute(**kwargs)

    # In principle, when c_order is False, this would be enough:
    # return reshape(lazyarr, shape, c_order=c_order, **kwargs)
    # so that an intermediate NDArray wouldn't be needed, which is more memory efficient.
    # However, benchmarks show that performance is better with the approach below.
    # Incidentally, not requiring C order can be quite illustrative for the user to
    # understand how the process of computing lazy arrays (and chunking) works.
    larr = lazyarr.compute()  # intermediate array
    return reshape(larr, shape, c_order=c_order, **kwargs)


# Define a numpy linspace-like function
def linspace(start, stop, num=50, endpoint=True, dtype=np.float64, shape=None, c_order=True, **kwargs: Any):
    """Return evenly spaced numbers over a specified interval.

    This is similar to `numpy.linspace` but it returns a `NDArray`
    instead of a numpy array.  Also, it supports a `shape` parameter
    to return a ndim array.

    Parameters
    ----------
    start: int, float, complex or np.number
        The starting value of the sequence.
    stop: int, float, complex or np.number
        The end value of the sequence.
    num: int
        Number of samples to generate.
    endpoint: bool
        If True, `stop` is the last sample. Otherwise, it is not included.
    dtype: np.dtype or list str
        The data type of the array elements in NumPy format. Default is `np.float64`.
    shape: int, tuple or list
        The shape of the final array. If None, the shape will be guessed from `num`.
    c_order: bool
        Whether to store the array in C order (row-major) or insertion order.
        Insertion order means that values will be stored in the array
        following the order of chunks in the array; this is more memory
        efficient, as it does not require an intermediate copy of the array.
        Default is C order.

    Returns
    -------
    out: :ref:`NDArray`
        A :ref:`NDArray` is returned.
    """

    def linspace_fill(inputs, output, offset):
        lout = len(output)
        start, stop, num = inputs
        # Compute proper start and stop values for the current chunk
        start_ = start + offset[0] / num * (stop - start)
        stop_ = start_ + lout / num * (stop - start)
        output[:] = np.linspace(start_, stop_, lout, endpoint=False, dtype=output.dtype)

    if not shape:
        shape = (num,)
    dtype = _check_dtype(dtype)

    if is_inside_new_expr():
        # We already have the dtype and shape, so return immediately
        return blosc2.zeros(shape, dtype=dtype)

    lshape = (math.prod(shape),)
    if endpoint:
        stop += (stop - start) / (num - 1)
    inputs = (start, stop, num)
    lazyarr = blosc2.lazyudf(linspace_fill, inputs, dtype=dtype, shape=lshape)
    if len(shape) == 1:
        # C order is guaranteed, and no reshape is needed
        return lazyarr.compute(**kwargs)

    # In principle, when c_order is False, the intermediate array wouldn't be needed,
    # but this is faster; see arange() for more details.
    larr = lazyarr.compute()  # intermediate array
    return reshape(larr, shape, c_order=c_order, **kwargs)


def eye(N, M=None, k=0, dtype=np.float64, **kwargs: Any):
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

    # In principle, when c_order is False, the intermediate array wouldn't be needed,
    # but this is faster; see arange() for more details.
    larr = lazyarr.compute()  # intermediate array
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
    >>> array = np.arange(0, 100, dtype=np.int64).reshape(10, 10)
    >>> # Save the array to a file
    >>> blosc2.save(array, "array.b2")
    """
    array.save(urlpath, contiguous, **kwargs)


def asarray(array: np.ndarray | blosc2.C2Array, **kwargs: Any) -> NDArray:
    """Convert the `array` to an `NDArray`.

    Parameters
    ----------
    array: array_like
        An array supporting numpy array interface.

    Other Parameters
    ----------------
    kwargs: dict, optional
        Keyword arguments that are supported by the :func:`empty` constructor.

    Returns
    -------
    out: :ref:`NDArray`
        An new NDArray made of :paramref:`array`.

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
    kwargs = _check_ndarray_kwargs(**kwargs)
    chunks = kwargs.pop("chunks", None)
    blocks = kwargs.pop("blocks", None)
    # Use the chunks and blocks from the array if they are not passed
    if chunks is None and hasattr(array, "chunks"):
        chunks = array.chunks
    if blocks is None and hasattr(array, "blocks"):
        blocks = array.blocks
    chunks, blocks = compute_chunks_blocks(array.shape, chunks, blocks, array.dtype, **kwargs)

    # Fast path for small arrays. This is not too expensive in terms of memory consumption.
    shape = array.shape
    small_size = 2**24  # 16 MB
    array_nbytes = math.prod(shape) * array.dtype.itemsize
    if array_nbytes < small_size:
        if not isinstance(array, np.ndarray):
            if hasattr(array, "chunks"):
                # A getitem operation should be enough to get a numpy array
                array = array[:]
        else:
            if not array.flags.contiguous:
                array = np.ascontiguousarray(array)
        return blosc2_ext.asarray(array, chunks, blocks, **kwargs)

    # Create the empty array
    ndarr = empty(shape, array.dtype, chunks=chunks, blocks=blocks, **kwargs)
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
        # Ensure the array slice is contiguous
        array_slice = np.ascontiguousarray(array[slice_])
        if behaved:
            # The whole chunk is to be updated, so this fastpath is safe
            ndarr.schunk.update_data(nchunk, array_slice, copy=False)
        else:
            ndarr[slice_] = array_slice

    return ndarr


def _check_ndarray_kwargs(**kwargs):  # noqa: C901
    if "storage" in kwargs:
        for key in kwargs:
            if key in list(blosc2.Storage.__annotations__):
                raise AttributeError(
                    "Cannot pass both `storage` and other kwargs already included in Storage"
                )
        storage = kwargs.get("storage")
        if isinstance(storage, blosc2.Storage):
            kwargs = {**kwargs, **asdict(storage)}
        else:
            kwargs = {**kwargs, **storage}

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
    ]
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
        start, stop, step = get_ndarray_start_stop(array.ndim, key, array.shape)
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


def indices(array: NDArray, order: str | list[str] | None = None, **kwargs: Any) -> NDArray:
    """
    Return the indices of a sorted array following the specified order.

    This is only valid for 1-dim structured arrays.

    Parameters
    ----------
    array: :ref:`NDArray`
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


def sort(array: NDArray, order: str | list[str] | None = None, **kwargs: Any) -> NDArray:
    """
    Return a sorted array following the specified order.

    This is only valid for 1-dim structured arrays.

    Parameters
    ----------
    array: :ref:`NDArray`
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


def matmul(x1: NDArray, x2: NDArray, **kwargs: Any) -> NDArray:
    """
    Computes the matrix product between two Blosc2 NDArrays.

    Parameters
    ----------
    x1: `NDArray`
        The first input array.
    x2: `NDArray`
        The second input array.
    kwargs: Any, optional
        Keyword arguments that are supported by the :func:`empty` constructor.

    Returns
    -------
    out: :ref:`NDArray`
        The matrix product of the inputs. This is a scalar only when both x1,
        x2 are 1-d vectors.

    Raises
    ------
    ValueError
        If the last dimension of ``x1`` is not the same size as
        the second-to-last dimension of ``x2``.

        If a scalar value is passed in.

    References
    ----------
    `numpy.matmul <https://numpy.org/doc/stable/reference/generated/numpy.matmul.html>`_

    Examples
    --------
    For 2-D arrays it is the matrix product:

    >>> import numpy as np
    >>> import blosc2
    >>> a = np.array([[1, 2],
    ...               [3, 4]])
    >>> nd_a = blosc2.asarray(a)
    >>> b = np.array([[2, 3],
    ...               [2, 1]])
    >>> nd_b = blosc2.asarray(b)
    >>> blosc2.matmul(nd_a, nd_b)
    array([[ 6,  5],
           [14, 13]])

    For 2-D mixed with 1-D, the result is the usual.

    >>> a = np.array([[1, 3],
    ...               [0, 1]])
    >>> nd_a = blosc2.asarray(a)
    >>> v = np.array([1, 2])
    >>> nd_v = blosc2.asarray(v)
    >>> blosc2.matmul(nd_a, nd_v)
    array([7, 2])
    >>> blosc2.matmul(nd_v, nd_a)
    array([1, 5])

    """

    # Validate arguments are not scalars
    if np.isscalar(x1) or np.isscalar(x2):
        raise ValueError("Arguments can't be scalars.")

    # Validate arguments are dimension 1 or 2
    if x1.ndim > 2 or x2.ndim > 2:
        raise ValueError("Multiplication of arrays with dimension greater than 2 is not supported yet.")

    # Promote 1D arrays to 2D if necessary
    x1_is_vector = False
    x2_is_vector = False
    if x1.ndim == 1:
        x1 = x1.reshape((1, x1.shape[0]))  # (N,) -> (1, N)
        x1_is_vector = True
    if x2.ndim == 1:
        x2 = x2.reshape((x2.shape[0], 1))  # (M,) -> (M, 1)
        x2_is_vector = True

    # Validate matrix multiplication compatibility
    if x1.shape[-1] != x2.shape[-2]:
        raise ValueError("Shapes are not aligned for matrix multiplication.")

    n, k = x1.shape[-2:]
    m = x2.shape[-1]

    result = blosc2.zeros((n, m), dtype=np.result_type(x1, x2), **kwargs)

    p, q = result.chunks[-2:]
    r = x2.chunks[-1]

    for row in range(0, n, p):
        row_end = (row + p) if (row + p) < n else n
        for col in range(0, m, q):
            col_end = (col + q) if (col + q) < m else m
            for aux in range(0, k, r):
                aux_end = (aux + r) if (aux + r) < k else k
                bx1 = x1[row:row_end, aux:aux_end]
                bx2 = x2[aux:aux_end, col:col_end]
                result[row:row_end, col:col_end] += np.matmul(bx1, bx2)

    if x1_is_vector and x2_is_vector:
        return result[0][0]

    return result.squeeze()


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
        self.dtype = ndarr.dtype.fields[field][0]
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

    def __setitem__(self, key: int | slice | Sequence[slice], value: np.ndarray | NDArray | NDField) -> None:
        """
        Set a slice of :paramref:`self` to a value.

        Parameters
        ----------
        key: int or slice or Sequence[slice]
            The slice to be set.
        value: np.ndarray or NDArray or NDField
            The value to be set.
        """
        if isinstance(key, str):
            raise TypeError("This array is a NDField; use a structured NDArray for bool expressions")
        if isinstance(value, (NDField, NDArray)):
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

    def __len__(self):
        return self.shape[0]
