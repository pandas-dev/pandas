from __future__ import annotations

import functools
from functools import partial
from numbers import Integral

import numpy as np
from tlz import sliding_window

from dask._collections import new_collection
from dask._task_spec import Task
from dask.array._array_expr._collection import asarray, concatenate
from dask.array._array_expr._expr import ArrayExpr
from dask.array.chunk import arange as _arange
from dask.array.chunk import linspace as _linspace
from dask.array.core import normalize_chunks
from dask.array.creation import _get_like_function_shapes_chunks
from dask.array.utils import meta_from_array
from dask.array.wrap import _parse_wrap_args, broadcast_trick
from dask.blockwise import blockwise as core_blockwise
from dask.layers import ArrayChunkShapeDep
from dask.utils import cached_cumsum, derived_from


class Arange(ArrayExpr):
    _parameters = ["start", "stop", "step", "chunks", "like", "dtype", "kwargs"]
    _defaults = {"chunks": "auto", "like": None, "dtype": None, "kwargs": None}

    @functools.cached_property
    def num_rows(self):
        return int(max(np.ceil((self.stop - self.start) / self.step), 0))

    @functools.cached_property
    def dtype(self):
        return (
            self.operand("dtype")
            or np.arange(
                self.start,
                self.stop,
                self.step * self.num_rows if self.num_rows else self.step,
            ).dtype
        )

    @functools.cached_property
    def _meta(self):
        return meta_from_array(self.like, ndim=1, dtype=self.dtype)

    @functools.cached_property
    def chunks(self):
        return normalize_chunks(
            self.operand("chunks"), (self.num_rows,), dtype=self.dtype
        )

    def _layer(self) -> dict:
        dsk = {}
        elem_count = 0
        start, step = self.start, self.step
        like = self.like
        func = partial(_arange, like=like)

        for i, bs in enumerate(self.chunks[0]):
            blockstart = start + (elem_count * step)
            blockstop = start + ((elem_count + bs) * step)
            task = Task(
                (self._name, i),
                func,
                blockstart,
                blockstop,
                step,
                bs,
                self.dtype,
            )
            dsk[(self._name, i)] = task
            elem_count += bs
        return dsk


class Linspace(Arange):
    _parameters = ["start", "stop", "num", "endpoint", "chunks", "dtype"]
    _defaults = {"num": 50, "endpoint": True, "chunks": "auto", "dtype": None}
    like = None

    @functools.cached_property
    def num_rows(self):
        return self.operand("num")

    @functools.cached_property
    def dtype(self):
        return self.operand("dtype") or np.linspace(0, 1, 1).dtype

    @functools.cached_property
    def step(self):
        range_ = self.stop - self.start

        div = (self.num_rows - 1) if self.endpoint else self.num_rows
        if div == 0:
            div = 1

        return float(range_) / div

    def _layer(self) -> dict:
        dsk = {}
        blockstart = self.start
        func = partial(_linspace, endpoint=self.endpoint, dtype=self.dtype)

        for i, bs in enumerate(self.chunks[0]):
            bs_space = bs - 1 if self.endpoint else bs
            blockstop = blockstart + (bs_space * self.step)
            task = Task(
                (self._name, i),
                func,
                blockstart,
                blockstop,
                bs,
            )
            blockstart = blockstart + (self.step * bs)
            dsk[task.key] = task
        return dsk


class BroadcastTrick(ArrayExpr):
    _parameters = ["shape", "dtype", "chunks", "meta", "kwargs"]
    _defaults = {"meta": None}

    @functools.cached_property
    def _meta(self):
        return meta_from_array(
            self.operand("meta"), ndim=self.ndim, dtype=self.operand("dtype")
        )

    def _layer(self) -> dict:
        func = broadcast_trick(self.func)
        k = self.kwargs.copy()
        k.pop("meta", None)
        func = partial(func, meta=self._meta, dtype=self.dtype, **k)
        out_ind = dep_ind = tuple(range(len(self.shape)))
        return core_blockwise(
            func,
            self._name,
            out_ind,
            ArrayChunkShapeDep(self.chunks),
            dep_ind,
            numblocks={},
        )


class Ones(BroadcastTrick):
    func = staticmethod(np.ones_like)


class Zeros(BroadcastTrick):
    func = staticmethod(np.zeros_like)


class Empty(BroadcastTrick):
    func = staticmethod(np.empty_like)


class Full(BroadcastTrick):
    func = staticmethod(np.full_like)


def wrap_func_shape_as_first_arg(*args, klass, **kwargs):
    """
    Transform np creation function into blocked version
    """
    if "shape" not in kwargs:
        shape, args = args[0], args[1:]
    else:
        shape = kwargs.pop("shape")

    if isinstance(shape, ArrayExpr):
        raise TypeError(
            "Dask array input not supported. "
            "Please use tuple, list, or a 1D numpy array instead."
        )

    parsed = _parse_wrap_args(klass.func, args, kwargs, shape)
    return new_collection(
        klass(
            parsed["shape"],
            parsed["dtype"],
            parsed["chunks"],
            kwargs.get("meta"),
            kwargs,
        )
    )


def wrap(func, **kwargs):
    return partial(func, **kwargs)


ones = wrap(wrap_func_shape_as_first_arg, klass=Ones, dtype="f8")
zeros = wrap(wrap_func_shape_as_first_arg, klass=Zeros, dtype="f8")
empty = wrap(wrap_func_shape_as_first_arg, klass=Empty, dtype="f8")
_full = wrap(wrap_func_shape_as_first_arg, klass=Full, dtype="f8")


def arange(start=0, stop=None, step=1, *, chunks="auto", like=None, dtype=None):
    """
    Return evenly spaced values from `start` to `stop` with step size `step`.

    The values are half-open [start, stop), so including start and excluding
    stop. This is basically the same as python's range function but for dask
    arrays.

    When using a non-integer step, such as 0.1, the results will often not be
    consistent. It is better to use linspace for these cases.

    Parameters
    ----------
    start : int, optional
        The starting value of the sequence. The default is 0.
    stop : int
        The end of the interval, this value is excluded from the interval.
    step : int, optional
        The spacing between the values. The default is 1 when not specified.
    chunks :  int
        The number of samples on each block. Note that the last block will have
        fewer samples if ``len(array) % chunks != 0``.
        Defaults to "auto" which will automatically determine chunk sizes.
    dtype : numpy.dtype
        Output dtype. Omit to infer it from start, stop, step
        Defaults to ``None``.
    like : array type or ``None``
        Array to extract meta from. Defaults to ``None``.

    Returns
    -------
    samples : dask array

    See Also
    --------
    dask.array.linspace
    """
    if stop is None:
        stop = start
        start = 0
    return new_collection(Arange(start, stop, step, chunks, like, dtype))


def linspace(
    start, stop, num=50, endpoint=True, retstep=False, chunks="auto", dtype=None
):
    """
    Return `num` evenly spaced values over the closed interval [`start`,
    `stop`].

    Parameters
    ----------
    start : scalar
        The starting value of the sequence.
    stop : scalar
        The last value of the sequence.
    num : int, optional
        Number of samples to include in the returned dask array, including the
        endpoints. Default is 50.
    endpoint : bool, optional
        If True, ``stop`` is the last sample. Otherwise, it is not included.
        Default is True.
    retstep : bool, optional
        If True, return (samples, step), where step is the spacing between
        samples. Default is False.
    chunks :  int
        The number of samples on each block. Note that the last block will have
        fewer samples if `num % blocksize != 0`
    dtype : dtype, optional
        The type of the output array.

    Returns
    -------
    samples : dask array
    step : float, optional
        Only returned if ``retstep`` is True. Size of spacing between samples.

    See Also
    --------
    dask.array.arange
    """
    num = int(num)
    result = new_collection(Linspace(start, stop, num, endpoint, chunks, dtype))
    if retstep:
        return result, result.expr.step
    else:
        return result


def empty_like(a, dtype=None, order="C", chunks=None, name=None, shape=None):
    """
    Return a new array with the same shape and type as a given array.

    Parameters
    ----------
    a : array_like
        The shape and data-type of `a` define these same attributes of the
        returned array.
    dtype : data-type, optional
        Overrides the data type of the result.
    order : {'C', 'F'}, optional
        Whether to store multidimensional data in C- or Fortran-contiguous
        (row- or column-wise) order in memory.
    chunks : sequence of ints
        The number of samples on each block. Note that the last block will have
        fewer samples if ``len(array) % chunks != 0``.
    name : str, optional
        An optional keyname for the array. Defaults to hashing the input
        keyword arguments.
    shape : int or sequence of ints, optional.
        Overrides the shape of the result.

    Returns
    -------
    out : ndarray
        Array of uninitialized (arbitrary) data with the same
        shape and type as `a`.

    See Also
    --------
    ones_like : Return an array of ones with shape and type of input.
    zeros_like : Return an array of zeros with shape and type of input.
    empty : Return a new uninitialized array.
    ones : Return a new array setting values to one.
    zeros : Return a new array setting values to zero.

    Notes
    -----
    This function does *not* initialize the returned array; to do that use
    `zeros_like` or `ones_like` instead.  It may be marginally faster than
    the functions that do set the array values.
    """

    a = asarray(a, name=False)
    shape, chunks = _get_like_function_shapes_chunks(a, chunks, shape)

    # if shape is nan we cannot rely on regular empty function, we use
    # generic map_blocks.
    if np.isnan(shape).any():
        return a.map_blocks(partial(np.empty_like, dtype=(dtype or a.dtype)))

    return empty(
        shape,
        dtype=(dtype or a.dtype),
        order=order,
        chunks=chunks,
        name=name,
        meta=a._meta,
    )


def ones_like(a, dtype=None, order="C", chunks=None, name=None, shape=None):
    """
    Return an array of ones with the same shape and type as a given array.

    Parameters
    ----------
    a : array_like
        The shape and data-type of `a` define these same attributes of
        the returned array.
    dtype : data-type, optional
        Overrides the data type of the result.
    order : {'C', 'F'}, optional
        Whether to store multidimensional data in C- or Fortran-contiguous
        (row- or column-wise) order in memory.
    chunks : sequence of ints
        The number of samples on each block. Note that the last block will have
        fewer samples if ``len(array) % chunks != 0``.
    name : str, optional
        An optional keyname for the array. Defaults to hashing the input
        keyword arguments.
    shape : int or sequence of ints, optional.
        Overrides the shape of the result.

    Returns
    -------
    out : ndarray
        Array of ones with the same shape and type as `a`.

    See Also
    --------
    zeros_like : Return an array of zeros with shape and type of input.
    empty_like : Return an empty array with shape and type of input.
    zeros : Return a new array setting values to zero.
    ones : Return a new array setting values to one.
    empty : Return a new uninitialized array.
    """

    a = asarray(a, name=False)
    shape, chunks = _get_like_function_shapes_chunks(a, chunks, shape)

    # if shape is nan we cannot rely on regular ones function, we use
    # generic map_blocks.
    if np.isnan(shape).any():
        return a.map_blocks(partial(np.ones_like, dtype=(dtype or a.dtype)))

    return ones(
        shape,
        dtype=(dtype or a.dtype),
        order=order,
        chunks=chunks,
        name=name,
        meta=a._meta,
    )


def zeros_like(a, dtype=None, order="C", chunks=None, name=None, shape=None):
    """
    Return an array of zeros with the same shape and type as a given array.

    Parameters
    ----------
    a : array_like
        The shape and data-type of `a` define these same attributes of
        the returned array.
    dtype : data-type, optional
        Overrides the data type of the result.
    order : {'C', 'F'}, optional
        Whether to store multidimensional data in C- or Fortran-contiguous
        (row- or column-wise) order in memory.
    chunks : sequence of ints
        The number of samples on each block. Note that the last block will have
        fewer samples if ``len(array) % chunks != 0``.
    name : str, optional
        An optional keyname for the array. Defaults to hashing the input
        keyword arguments.
    shape : int or sequence of ints, optional.
        Overrides the shape of the result.

    Returns
    -------
    out : ndarray
        Array of zeros with the same shape and type as `a`.

    See Also
    --------
    ones_like : Return an array of ones with shape and type of input.
    empty_like : Return an empty array with shape and type of input.
    zeros : Return a new array setting values to zero.
    ones : Return a new array setting values to one.
    empty : Return a new uninitialized array.
    """

    a = asarray(a, name=False)
    shape, chunks = _get_like_function_shapes_chunks(a, chunks, shape)

    # if shape is nan we cannot rely on regular zeros function, we use
    # generic map_blocks.
    if np.isnan(shape).any():
        return a.map_blocks(partial(np.zeros_like, dtype=(dtype or a.dtype)))

    return zeros(
        shape,
        dtype=(dtype or a.dtype),
        order=order,
        chunks=chunks,
        name=name,
        meta=a._meta,
    )


def full(shape, fill_value, *args, **kwargs):
    # np.isscalar has somewhat strange behavior:
    # https://docs.scipy.org/doc/numpy/reference/generated/numpy.isscalar.html
    if np.ndim(fill_value) != 0:
        raise ValueError(
            f"fill_value must be scalar. Received {type(fill_value).__name__} instead."
        )
    if kwargs.get("dtype") is None:
        if hasattr(fill_value, "dtype"):
            kwargs["dtype"] = fill_value.dtype
        else:
            kwargs["dtype"] = type(fill_value)
    return _full(*args, shape=shape, fill_value=fill_value, **kwargs)


def full_like(a, fill_value, order="C", dtype=None, chunks=None, name=None, shape=None):
    """
    Return a full array with the same shape and type as a given array.

    Parameters
    ----------
    a : array_like
        The shape and data-type of `a` define these same attributes of
        the returned array.
    fill_value : scalar
        Fill value.
    dtype : data-type, optional
        Overrides the data type of the result.
    order : {'C', 'F'}, optional
        Whether to store multidimensional data in C- or Fortran-contiguous
        (row- or column-wise) order in memory.
    chunks : sequence of ints
        The number of samples on each block. Note that the last block will have
        fewer samples if ``len(array) % chunks != 0``.
    name : str, optional
        An optional keyname for the array. Defaults to hashing the input
        keyword arguments.
    shape : int or sequence of ints, optional.
        Overrides the shape of the result.

    Returns
    -------
    out : ndarray
        Array of `fill_value` with the same shape and type as `a`.

    See Also
    --------
    zeros_like : Return an array of zeros with shape and type of input.
    ones_like : Return an array of ones with shape and type of input.
    empty_like : Return an empty array with shape and type of input.
    zeros : Return a new array setting values to zero.
    ones : Return a new array setting values to one.
    empty : Return a new uninitialized array.
    full : Fill a new array.
    """

    a = asarray(a, name=False)
    shape, chunks = _get_like_function_shapes_chunks(a, chunks, shape)

    # if shape is nan we cannot rely on regular full function, we use
    # generic map_blocks.
    if np.isnan(shape).any():
        return a.map_blocks(partial(np.full_like, dtype=(dtype or a.dtype)), fill_value)

    return full(
        shape,
        fill_value,
        dtype=(dtype or a.dtype),
        order=order,
        chunks=chunks,
        name=name,
        meta=a._meta,
    )


@derived_from(np)
def repeat(a, repeats, axis=None):
    if axis is None:
        if a.ndim == 1:
            axis = 0
        else:
            raise NotImplementedError("Must supply an integer axis value")

    if not isinstance(repeats, Integral):
        raise NotImplementedError("Only integer valued repeats supported")

    if -a.ndim <= axis < 0:
        axis += a.ndim
    elif not 0 <= axis <= a.ndim - 1:
        raise ValueError("axis(=%d) out of bounds" % axis)

    if repeats == 0:
        return a[tuple(slice(None) if d != axis else slice(0) for d in range(a.ndim))]
    elif repeats == 1:
        return a

    cchunks = cached_cumsum(a.chunks[axis], initial_zero=True)
    slices = []
    for c_start, c_stop in sliding_window(2, cchunks):
        ls = np.linspace(c_start, c_stop, repeats).round(0)
        for ls_start, ls_stop in sliding_window(2, ls):
            if ls_start != ls_stop:
                slices.append(slice(ls_start, ls_stop))

    all_slice = slice(None, None, None)
    slices = [
        (all_slice,) * axis + (s,) + (all_slice,) * (a.ndim - axis - 1) for s in slices
    ]

    slabs = [a[slc] for slc in slices]

    out = []
    for slab in slabs:
        chunks = list(slab.chunks)
        assert len(chunks[axis]) == 1
        chunks[axis] = (chunks[axis][0] * repeats,)
        chunks = tuple(chunks)
        result = slab.map_blocks(
            np.repeat, repeats, axis=axis, chunks=chunks, dtype=slab.dtype
        )
        out.append(result)

    return concatenate(out, axis=axis)
