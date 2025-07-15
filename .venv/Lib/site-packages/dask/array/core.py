from __future__ import annotations

import contextlib
import functools
import math
import operator
import os
import pickle
import re
import sys
import traceback
import uuid
import warnings
from bisect import bisect
from collections import defaultdict
from collections.abc import Collection, Iterable, Iterator, Mapping, Sequence
from functools import lru_cache, partial, reduce, wraps
from itertools import product, zip_longest
from numbers import Integral, Number
from operator import add, mul
from threading import Lock
from typing import Any, Literal, TypeVar, Union, cast

import numpy as np
from numpy.typing import ArrayLike
from packaging.version import Version
from tlz import accumulate, concat, first, partition
from toolz import frequencies

from dask._compatibility import import_optional_dependency
from dask.core import flatten

xr = import_optional_dependency("xarray", errors="ignore")

from dask import compute, config, core
from dask._task_spec import Alias, List, Task, TaskRef
from dask.array import chunk
from dask.array.chunk import getitem
from dask.array.chunk_types import is_valid_array_chunk, is_valid_chunk_type

# Keep einsum_lookup and tensordot_lookup here for backwards compatibility
from dask.array.dispatch import (  # noqa: F401
    concatenate_lookup,
    einsum_lookup,
    tensordot_lookup,
)
from dask.array.numpy_compat import NUMPY_GE_200, _Recurser
from dask.array.slicing import SlicingNoop, replace_ellipsis, setitem_array, slice_array
from dask.array.utils import (
    asanyarray_safe,
    asarray_safe,
    compute_meta,
    meta_from_array,
)
from dask.base import (
    DaskMethodsMixin,
    compute_as_if_collection,
    dont_optimize,
    is_dask_collection,
    named_schedulers,
    persist,
    tokenize,
)
from dask.blockwise import BlockwiseDep
from dask.blockwise import blockwise as core_blockwise
from dask.blockwise import broadcast_dimensions
from dask.context import globalmethod
from dask.core import quote, reshapelist
from dask.delayed import Delayed, delayed
from dask.highlevelgraph import HighLevelGraph, MaterializedLayer
from dask.layers import ArrayBlockIdDep, ArraySliceDep, ArrayValuesDep
from dask.sizeof import sizeof
from dask.typing import Graph, NestedKeys
from dask.utils import (
    IndexCallable,
    SerializableLock,
    cached_cumsum,
    cached_max,
    cached_property,
    concrete,
    derived_from,
    format_bytes,
    funcname,
    get_scheduler_lock,
    has_keyword,
    is_arraylike,
    is_dataframe_like,
    is_index_like,
    is_integer,
    is_series_like,
    maybe_pluralize,
    ndeepmap,
    ndimlist,
    parse_bytes,
    typename,
)
from dask.widgets import get_template

try:
    ARRAY_TEMPLATE = get_template("array.html.j2")
except ImportError:
    ARRAY_TEMPLATE = None

T_IntOrNaN = Union[int, float]  # Should be Union[int, Literal[np.nan]]

DEFAULT_GET = named_schedulers.get("threads", named_schedulers["sync"])

unknown_chunk_message = (
    "\n\n"
    "A possible solution: "
    "https://docs.dask.org/en/latest/array-chunks.html#unknown-chunks\n"
    "Summary: to compute chunks sizes, use\n\n"
    "   x.compute_chunk_sizes()  # for Dask Array `x`\n"
    "   ddf.to_dask_array(lengths=True)  # for Dask DataFrame `ddf`"
)


class PerformanceWarning(Warning):
    """A warning given when bad chunking may cause poor performance"""


def getter(a, b, asarray=True, lock=None):
    if isinstance(b, tuple) and any(x is None for x in b):
        b2 = tuple(x for x in b if x is not None)
        b3 = tuple(
            None if x is None else slice(None, None)
            for x in b
            if not isinstance(x, Integral)
        )
        return getter(a, b2, asarray=asarray, lock=lock)[b3]

    if lock:
        lock.acquire()
    try:
        c = a[b]
        # Below we special-case `np.matrix` to force a conversion to
        # `np.ndarray` and preserve original Dask behavior for `getter`,
        # as for all purposes `np.matrix` is array-like and thus
        # `is_arraylike` evaluates to `True` in that case.
        if asarray and (not is_arraylike(c) or isinstance(c, np.matrix)):
            c = np.asarray(c)
    finally:
        if lock:
            lock.release()
    return c


def getter_nofancy(a, b, asarray=True, lock=None):
    """A simple wrapper around ``getter``.

    Used to indicate to the optimization passes that the backend doesn't
    support fancy indexing.
    """
    return getter(a, b, asarray=asarray, lock=lock)


def getter_inline(a, b, asarray=True, lock=None):
    """A getter function that optimizations feel comfortable inlining

    Slicing operations with this function may be inlined into a graph, such as
    in the following rewrite

    **Before**

    >>> a = x[:10]  # doctest: +SKIP
    >>> b = a + 1  # doctest: +SKIP
    >>> c = a * 2  # doctest: +SKIP

    **After**

    >>> b = x[:10] + 1  # doctest: +SKIP
    >>> c = x[:10] * 2  # doctest: +SKIP

    This inlining can be relevant to operations when running off of disk.
    """
    return getter(a, b, asarray=asarray, lock=lock)


from dask.array.optimization import fuse_slice, optimize

# __array_function__ dict for mapping aliases and mismatching names
_HANDLED_FUNCTIONS = {}


def implements(*numpy_functions):
    """Register an __array_function__ implementation for dask.array.Array

    Register that a function implements the API of a NumPy function (or several
    NumPy functions in case of aliases) which is handled with
    ``__array_function__``.

    Parameters
    ----------
    \\*numpy_functions : callables
        One or more NumPy functions that are handled by ``__array_function__``
        and will be mapped by `implements` to a `dask.array` function.
    """

    def decorator(dask_func):
        for numpy_function in numpy_functions:
            _HANDLED_FUNCTIONS[numpy_function] = dask_func

        return dask_func

    return decorator


def _should_delegate(self, other) -> bool:
    """Check whether Dask should delegate to the other.
    This implementation follows NEP-13:
    https://numpy.org/neps/nep-0013-ufunc-overrides.html#behavior-in-combination-with-python-s-binary-operations
    """
    if hasattr(other, "__array_ufunc__") and other.__array_ufunc__ is None:
        return True
    elif (
        hasattr(other, "__array_ufunc__")
        and not is_valid_array_chunk(other)
        # don't delegate to our own parent classes
        and not isinstance(self, type(other))
        and type(self) is not type(other)
    ):
        return True
    elif (
        not hasattr(other, "__array_ufunc__")
        and hasattr(other, "__array_priority__")
        and other.__array_priority__ > self.__array_priority__
    ):
        return True
    return False


def check_if_handled_given_other(f):
    """Check if method is handled by Dask given type of other

    Ensures proper deferral to upcast types in dunder operations without
    assuming unknown types are automatically downcast types.
    """

    @wraps(f)
    def wrapper(self, other):
        if _should_delegate(self, other):
            return NotImplemented
        else:
            return f(self, other)

    return wrapper


def slices_from_chunks(chunks):
    """Translate chunks tuple to a set of slices in product order

    >>> slices_from_chunks(((2, 2), (3, 3, 3)))  # doctest: +NORMALIZE_WHITESPACE
     [(slice(0, 2, None), slice(0, 3, None)),
      (slice(0, 2, None), slice(3, 6, None)),
      (slice(0, 2, None), slice(6, 9, None)),
      (slice(2, 4, None), slice(0, 3, None)),
      (slice(2, 4, None), slice(3, 6, None)),
      (slice(2, 4, None), slice(6, 9, None))]
    """
    cumdims = [cached_cumsum(bds, initial_zero=True) for bds in chunks]
    slices = [
        [slice(s, s + dim) for s, dim in zip(starts, shapes)]
        for starts, shapes in zip(cumdims, chunks)
    ]
    return list(product(*slices))


def graph_from_arraylike(
    arr,  # Any array-like which supports slicing
    chunks,
    shape,
    name,
    getitem=getter,
    lock=False,
    asarray=True,
    dtype=None,
    inline_array=False,
) -> HighLevelGraph:
    """
    HighLevelGraph for slicing chunks from an array-like according to a chunk pattern.

    If ``inline_array`` is True, this make a Blockwise layer of slicing tasks where the
    array-like is embedded into every task.,

    If ``inline_array`` is False, this inserts the array-like as a standalone value in
    a MaterializedLayer, then generates a Blockwise layer of slicing tasks that refer
    to it.

    >>> dict(graph_from_arraylike(arr, chunks=(2, 3), shape=(4, 6), name="X", inline_array=True))  # doctest: +SKIP
    {(arr, 0, 0): (getter, arr, (slice(0, 2), slice(0, 3))),
     (arr, 1, 0): (getter, arr, (slice(2, 4), slice(0, 3))),
     (arr, 1, 1): (getter, arr, (slice(2, 4), slice(3, 6))),
     (arr, 0, 1): (getter, arr, (slice(0, 2), slice(3, 6)))}

    >>> dict(  # doctest: +SKIP
            graph_from_arraylike(arr, chunks=((2, 2), (3, 3)), shape=(4,6), name="X", inline_array=False)
        )
    {"original-X": arr,
     ('X', 0, 0): (getter, 'original-X', (slice(0, 2), slice(0, 3))),
     ('X', 1, 0): (getter, 'original-X', (slice(2, 4), slice(0, 3))),
     ('X', 1, 1): (getter, 'original-X', (slice(2, 4), slice(3, 6))),
     ('X', 0, 1): (getter, 'original-X', (slice(0, 2), slice(3, 6)))}
    """
    chunks = normalize_chunks(chunks, shape, dtype=dtype)
    out_ind = tuple(range(len(shape)))

    if (
        has_keyword(getitem, "asarray")
        and has_keyword(getitem, "lock")
        and (not asarray or lock)
    ):
        kwargs = {"asarray": asarray, "lock": lock}
    else:
        # Common case, drop extra parameters
        kwargs = {}

    if inline_array:
        layer = core_blockwise(
            getitem,
            name,
            out_ind,
            arr,
            None,
            ArraySliceDep(chunks),
            out_ind,
            numblocks={},
            _data_producer=True,
            **kwargs,
        )
        return HighLevelGraph.from_collections(name, layer)
    else:
        original_name = "original-" + name

        layers = {}
        layers[original_name] = MaterializedLayer({original_name: arr})
        layers[name] = core_blockwise(
            getitem,
            name,
            out_ind,
            TaskRef(original_name),
            None,
            ArraySliceDep(chunks),
            out_ind,
            numblocks={},
            _data_producer=True,
            **kwargs,
        )

        deps = {
            original_name: set(),
            name: {original_name},
        }
        return HighLevelGraph(layers, deps)


def dotmany(A, B, leftfunc=None, rightfunc=None, **kwargs):
    """Dot product of many aligned chunks

    >>> x = np.array([[1, 2], [1, 2]])
    >>> y = np.array([[10, 20], [10, 20]])
    >>> dotmany([x, x, x], [y, y, y])
    array([[ 90, 180],
           [ 90, 180]])

    Optionally pass in functions to apply to the left and right chunks

    >>> dotmany([x, x, x], [y, y, y], rightfunc=np.transpose)
    array([[150, 150],
           [150, 150]])
    """
    if leftfunc:
        A = map(leftfunc, A)
    if rightfunc:
        B = map(rightfunc, B)
    return sum(map(partial(np.dot, **kwargs), A, B))


def _concatenate2(arrays, axes=None):
    """Recursively concatenate nested lists of arrays along axes

    Each entry in axes corresponds to each level of the nested list.  The
    length of axes should correspond to the level of nesting of arrays.
    If axes is an empty list or tuple, return arrays, or arrays[0] if
    arrays is a list.

    >>> x = np.array([[1, 2], [3, 4]])
    >>> _concatenate2([x, x], axes=[0])
    array([[1, 2],
           [3, 4],
           [1, 2],
           [3, 4]])

    >>> _concatenate2([x, x], axes=[1])
    array([[1, 2, 1, 2],
           [3, 4, 3, 4]])

    >>> _concatenate2([[x, x], [x, x]], axes=[0, 1])
    array([[1, 2, 1, 2],
           [3, 4, 3, 4],
           [1, 2, 1, 2],
           [3, 4, 3, 4]])

    Supports Iterators
    >>> _concatenate2(iter([x, x]), axes=[1])
    array([[1, 2, 1, 2],
           [3, 4, 3, 4]])

    Special Case
    >>> _concatenate2([x, x], axes=())
    array([[1, 2],
           [3, 4]])
    """
    if axes is None:
        axes = []

    if axes == ():
        if isinstance(arrays, list):
            return arrays[0]
        else:
            return arrays

    if isinstance(arrays, Iterator):
        arrays = list(arrays)
    if not isinstance(arrays, (list, tuple)):
        return arrays
    if len(axes) > 1:
        arrays = [_concatenate2(a, axes=axes[1:]) for a in arrays]
    concatenate = concatenate_lookup.dispatch(
        type(max(arrays, key=lambda x: getattr(x, "__array_priority__", 0)))
    )
    if isinstance(arrays[0], dict):
        # Handle concatenation of `dict`s, used as a replacement for structured
        # arrays when that's not supported by the array library (e.g., CuPy).
        keys = list(arrays[0].keys())
        assert all(list(a.keys()) == keys for a in arrays)
        ret = dict()
        for k in keys:
            ret[k] = concatenate(list(a[k] for a in arrays), axis=axes[0])
        return ret
    else:
        return concatenate(arrays, axis=axes[0])


def apply_infer_dtype(func, args, kwargs, funcname, suggest_dtype="dtype", nout=None):
    """
    Tries to infer output dtype of ``func`` for a small set of input arguments.

    Parameters
    ----------
    func: Callable
        Function for which output dtype is to be determined

    args: List of array like
        Arguments to the function, which would usually be used. Only attributes
        ``ndim`` and ``dtype`` are used.

    kwargs: dict
        Additional ``kwargs`` to the ``func``

    funcname: String
        Name of calling function to improve potential error messages

    suggest_dtype: None/False or String
        If not ``None`` adds suggestion to potential error message to specify a dtype
        via the specified kwarg. Defaults to ``'dtype'``.

    nout: None or Int
        ``None`` if function returns single output, integer if many.
        Defaults to ``None``.

    Returns
    -------
    : dtype or List of dtype
        One or many dtypes (depending on ``nout``)
    """
    from dask.array.utils import meta_from_array

    # make sure that every arg is an evaluated array
    args = [
        (
            np.ones_like(meta_from_array(x), shape=((1,) * x.ndim), dtype=x.dtype)
            if is_arraylike(x)
            else x
        )
        for x in args
    ]
    try:
        with np.errstate(all="ignore"):
            o = func(*args, **kwargs)
    except Exception as e:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        tb = "".join(traceback.format_tb(exc_traceback))
        suggest = (
            (
                "Please specify the dtype explicitly using the "
                "`{dtype}` kwarg.\n\n".format(dtype=suggest_dtype)
            )
            if suggest_dtype
            else ""
        )
        msg = (
            f"`dtype` inference failed in `{funcname}`.\n\n"
            f"{suggest}"
            "Original error is below:\n"
            "------------------------\n"
            f"{e!r}\n\n"
            "Traceback:\n"
            "---------\n"
            f"{tb}"
        )
    else:
        msg = None
    if msg is not None:
        raise ValueError(msg)
    return getattr(o, "dtype", type(o)) if nout is None else tuple(e.dtype for e in o)


def normalize_arg(x):
    """Normalize user provided arguments to blockwise or map_blocks

    We do a few things:

    1.  If they are string literals that might collide with blockwise_token then we
        quote them
    2.  IF they are large (as defined by sizeof) then we put them into the
        graph on their own by using dask.delayed
    """
    if is_dask_collection(x):
        return x
    elif isinstance(x, str) and re.match(r"_\d+", x):
        return delayed(x)
    elif isinstance(x, list) and len(x) >= 10:
        return delayed(x)
    elif sizeof(x) > 1e6:
        return delayed(x)
    else:
        return x


def _pass_extra_kwargs(func, keys, *args, **kwargs):
    """Helper for :func:`dask.array.map_blocks` to pass `block_info` or `block_id`.

    For each element of `keys`, a corresponding element of args is changed
    to a keyword argument with that key, before all arguments re passed on
    to `func`.
    """
    kwargs.update(zip(keys, args))
    return func(*args[len(keys) :], **kwargs)


def map_blocks(
    func,
    *args,
    name=None,
    token=None,
    dtype=None,
    chunks=None,
    drop_axis=None,
    new_axis=None,
    enforce_ndim=False,
    meta=None,
    **kwargs,
):
    """Map a function across all blocks of a dask array.

    Note that ``map_blocks`` will attempt to automatically determine the output
    array type by calling ``func`` on 0-d versions of the inputs. Please refer to
    the ``meta`` keyword argument below if you expect that the function will not
    succeed when operating on 0-d arrays.

    Parameters
    ----------
    func : callable
        Function to apply to every block in the array.
        If ``func`` accepts ``block_info=`` or ``block_id=``
        as keyword arguments, these will be passed dictionaries
        containing information about input and output chunks/arrays
        during computation. See examples for details.
    args : dask arrays or other objects
    dtype : np.dtype, optional
        The ``dtype`` of the output array. It is recommended to provide this.
        If not provided, will be inferred by applying the function to a small
        set of fake data.
    chunks : tuple, optional
        Chunk shape of resulting blocks if the function does not preserve
        shape. If not provided, the resulting array is assumed to have the same
        block structure as the first input array.
    drop_axis : number or iterable, optional
        Dimensions lost by the function.
    new_axis : number or iterable, optional
        New dimensions created by the function. Note that these are applied
        after ``drop_axis`` (if present). The size of each chunk along this
        dimension will be set to 1. Please specify ``chunks`` if the individual
        chunks have a different size.
    enforce_ndim : bool, default False
        Whether to enforce at runtime that the dimensionality of the array
        produced by ``func`` actually matches that of the array returned by
        ``map_blocks``.
        If True, this will raise an error when there is a mismatch.
    token : string, optional
        The key prefix to use for the output array. If not provided, will be
        determined from the function name.
    name : string, optional
        The key name to use for the output array. Note that this fully
        specifies the output key name, and must be unique. If not provided,
        will be determined by a hash of the arguments.
    meta : array-like, optional
        The ``meta`` of the output array, when specified is expected to be an
        array of the same type and dtype of that returned when calling ``.compute()``
        on the array returned by this function. When not provided, ``meta`` will be
        inferred by applying the function to a small set of fake data, usually a
        0-d array. It's important to ensure that ``func`` can successfully complete
        computation without raising exceptions when 0-d is passed to it, providing
        ``meta`` will be required otherwise. If the output type is known beforehand
        (e.g., ``np.ndarray``, ``cupy.ndarray``), an empty array of such type dtype
        can be passed, for example: ``meta=np.array((), dtype=np.int32)``.
    **kwargs :
        Other keyword arguments to pass to function. Values must be constants
        (not dask.arrays)

    See Also
    --------
    dask.array.map_overlap : Generalized operation with overlap between neighbors.
    dask.array.blockwise : Generalized operation with control over block alignment.

    Examples
    --------
    >>> import dask.array as da
    >>> x = da.arange(6, chunks=3)

    >>> x.map_blocks(lambda x: x * 2).compute()
    array([ 0,  2,  4,  6,  8, 10])

    The ``da.map_blocks`` function can also accept multiple arrays.

    >>> d = da.arange(5, chunks=2)
    >>> e = da.arange(5, chunks=2)

    >>> f = da.map_blocks(lambda a, b: a + b**2, d, e)
    >>> f.compute()
    array([ 0,  2,  6, 12, 20])

    If the function changes shape of the blocks then you must provide chunks
    explicitly.

    >>> y = x.map_blocks(lambda x: x[::2], chunks=((2, 2),))

    You have a bit of freedom in specifying chunks.  If all of the output chunk
    sizes are the same, you can provide just that chunk size as a single tuple.

    >>> a = da.arange(18, chunks=(6,))
    >>> b = a.map_blocks(lambda x: x[:3], chunks=(3,))

    If the function changes the dimension of the blocks you must specify the
    created or destroyed dimensions.

    >>> b = a.map_blocks(lambda x: x[None, :, None], chunks=(1, 6, 1),
    ...                  new_axis=[0, 2])

    If ``chunks`` is specified but ``new_axis`` is not, then it is inferred to
    add the necessary number of axes on the left.

    Note that ``map_blocks()`` will concatenate chunks along axes specified by
    the keyword parameter ``drop_axis`` prior to applying the function.
    This is illustrated in the figure below:

    .. image:: /images/map_blocks_drop_axis.png

    Due to memory-size-constraints, it is often not advisable to use ``drop_axis``
    on an axis that is chunked.  In that case, it is better not to use
    ``map_blocks`` but rather
    ``dask.array.reduction(..., axis=dropped_axes, concatenate=False)`` which
    maintains a leaner memory footprint while it drops any axis.

    Map_blocks aligns blocks by block positions without regard to shape. In the
    following example we have two arrays with the same number of blocks but
    with different shape and chunk sizes.

    >>> x = da.arange(1000, chunks=(100,))
    >>> y = da.arange(100, chunks=(10,))

    The relevant attribute to match is numblocks.

    >>> x.numblocks
    (10,)
    >>> y.numblocks
    (10,)

    If these match (up to broadcasting rules) then we can map arbitrary
    functions across blocks

    >>> def func(a, b):
    ...     return np.array([a.max(), b.max()])

    >>> da.map_blocks(func, x, y, chunks=(2,), dtype='i8')
    dask.array<func, shape=(20,), dtype=int64, chunksize=(2,), chunktype=numpy.ndarray>

    >>> _.compute()
    array([ 99,   9, 199,  19, 299,  29, 399,  39, 499,  49, 599,  59, 699,
            69, 799,  79, 899,  89, 999,  99])

    Your block function can get information about where it is in the array by
    accepting a special ``block_info`` or ``block_id`` keyword argument.
    During computation, they will contain information about each of the input
    and output chunks (and dask arrays) relevant to each call of ``func``.

    >>> def func(block_info=None):
    ...     pass

    This will receive the following information:

    >>> block_info  # doctest: +SKIP
    {0: {'shape': (1000,),
         'num-chunks': (10,),
         'chunk-location': (4,),
         'array-location': [(400, 500)]},
     None: {'shape': (1000,),
            'num-chunks': (10,),
            'chunk-location': (4,),
            'array-location': [(400, 500)],
            'chunk-shape': (100,),
            'dtype': dtype('float64')}}

    The keys to the ``block_info`` dictionary indicate which is the input and
    output Dask array:

    - **Input Dask array(s):** ``block_info[0]`` refers to the first input Dask array.
      The dictionary key is ``0`` because that is the argument index corresponding
      to the first input Dask array.
      In cases where multiple Dask arrays have been passed as input to the function,
      you can access them with the number corresponding to the input argument,
      eg: ``block_info[1]``, ``block_info[2]``, etc.
      (Note that if you pass multiple Dask arrays as input to map_blocks,
      the arrays must match each other by having matching numbers of chunks,
      along corresponding dimensions up to broadcasting rules.)
    - **Output Dask array:** ``block_info[None]`` refers to the output Dask array,
      and contains information about the output chunks.
      The output chunk shape and dtype may may be different than the input chunks.

    For each dask array, ``block_info`` describes:

    - ``shape``: the shape of the full Dask array,
    - ``num-chunks``: the number of chunks of the full array in each dimension,
    - ``chunk-location``: the chunk location (for example the fourth chunk over
      in the first dimension), and
    - ``array-location``: the array location within the full Dask array
      (for example the slice corresponding to ``40:50``).

    In addition to these, there are two extra parameters described by
    ``block_info`` for the output array (in ``block_info[None]``):

    - ``chunk-shape``: the output chunk shape, and
    - ``dtype``: the output dtype.

    These features can be combined to synthesize an array from scratch, for
    example:

    >>> def func(block_info=None):
    ...     loc = block_info[None]['array-location'][0]
    ...     return np.arange(loc[0], loc[1])

    >>> da.map_blocks(func, chunks=((4, 4),), dtype=np.float64)
    dask.array<func, shape=(8,), dtype=float64, chunksize=(4,), chunktype=numpy.ndarray>

    >>> _.compute()
    array([0, 1, 2, 3, 4, 5, 6, 7])

    ``block_id`` is similar to ``block_info`` but contains only the ``chunk_location``:

    >>> def func(block_id=None):
    ...     pass

    This will receive the following information:

    >>> block_id  # doctest: +SKIP
    (4, 3)

    You may specify the key name prefix of the resulting task in the graph with
    the optional ``token`` keyword argument.

    >>> x.map_blocks(lambda x: x + 1, name='increment')
    dask.array<increment, shape=(1000,), dtype=int64, chunksize=(100,), chunktype=numpy.ndarray>

    For functions that may not handle 0-d arrays, it's also possible to specify
    ``meta`` with an empty array matching the type of the expected result. In
    the example below, ``func`` will result in an ``IndexError`` when computing
    ``meta``:

    >>> rng = da.random.default_rng()
    >>> da.map_blocks(lambda x: x[2], rng.random(5), meta=np.array(()))
    dask.array<lambda, shape=(5,), dtype=float64, chunksize=(5,), chunktype=numpy.ndarray>

    Similarly, it's possible to specify a non-NumPy array to ``meta``, and provide
    a ``dtype``:

    >>> import cupy  # doctest: +SKIP
    >>> rng = da.random.default_rng(cupy.random.default_rng())  # doctest: +SKIP
    >>> dt = np.float32
    >>> da.map_blocks(lambda x: x[2], rng.random(5, dtype=dt), meta=cupy.array((), dtype=dt))  # doctest: +SKIP
    dask.array<lambda, shape=(5,), dtype=float32, chunksize=(5,), chunktype=cupy.ndarray>
    """
    if drop_axis is None:
        drop_axis = []

    if not callable(func):
        msg = (
            "First argument must be callable function, not %s\n"
            "Usage:   da.map_blocks(function, x)\n"
            "   or:   da.map_blocks(function, x, y, z)"
        )
        raise TypeError(msg % type(func).__name__)
    if token:
        warnings.warn(
            "The `token=` keyword to `map_blocks` has been moved to `name=`. "
            "Please use `name=` instead as the `token=` keyword will be removed "
            "in a future release.",
            category=FutureWarning,
        )
        name = token

    token = f"{name or funcname(func)}"
    new_axes = {}

    if isinstance(drop_axis, Number):
        drop_axis = [drop_axis]
    if isinstance(new_axis, Number):
        new_axis = [new_axis]  # TODO: handle new_axis

    arrs = [a for a in args if isinstance(a, Array)]
    argpairs = []
    for a in args:
        if isinstance(a, Array):
            argpairs.append((a, tuple(range(a.ndim))[::-1]))
        elif isinstance(a, BlockwiseDep):
            argpairs.append((a, tuple(range(args[0].ndim))[::-1]))
        else:
            argpairs.append((a, None))

    if arrs:
        out_ind = tuple(range(max(a.ndim for a in arrs)))[::-1]
    else:
        out_ind = ()

    original_kwargs = kwargs

    if dtype is None and meta is None:
        try:
            meta = compute_meta(func, dtype, *args, **kwargs)
        except Exception:
            pass

        dtype = apply_infer_dtype(func, args, original_kwargs, "map_blocks")

    if drop_axis:
        ndim_out = len(out_ind)
        if any(i < -ndim_out or i >= ndim_out for i in drop_axis):
            raise ValueError(
                f"drop_axis out of range (drop_axis={drop_axis}, "
                f"but output is {ndim_out}d)."
            )
        drop_axis = [i % ndim_out for i in drop_axis]
        out_ind = tuple(x for i, x in enumerate(out_ind) if i not in drop_axis)
    if new_axis is None and chunks is not None and len(out_ind) < len(chunks):
        new_axis = range(len(chunks) - len(out_ind))
    if new_axis:
        # new_axis = [x + len(drop_axis) for x in new_axis]
        out_ind = list(out_ind)
        for ax in sorted(new_axis):
            n = len(out_ind) + len(drop_axis)
            out_ind.insert(ax, n)
            if chunks is not None:
                new_axes[n] = chunks[ax]
            else:
                new_axes[n] = 1
        out_ind = tuple(out_ind)
        if max(new_axis) > max(out_ind):
            raise ValueError("New_axis values do not fill in all dimensions")

    if chunks is not None:
        if len(chunks) != len(out_ind):
            raise ValueError(
                f"Provided chunks have {len(chunks)} dims; expected {len(out_ind)} dims"
            )
        adjust_chunks = dict(zip(out_ind, chunks))
    else:
        adjust_chunks = None

    if enforce_ndim:
        out = blockwise(
            apply_and_enforce,
            out_ind,
            *concat(argpairs),
            expected_ndim=len(out_ind),
            _func=func,
            token=token,
            new_axes=new_axes,
            dtype=dtype,
            concatenate=True,
            align_arrays=False,
            adjust_chunks=adjust_chunks,
            meta=meta,
            **kwargs,
        )
    else:
        out = blockwise(
            func,
            out_ind,
            *concat(argpairs),
            token=token,
            new_axes=new_axes,
            dtype=dtype,
            concatenate=True,
            align_arrays=False,
            adjust_chunks=adjust_chunks,
            meta=meta,
            **kwargs,
        )

    extra_argpairs = []
    extra_names = []

    # If func has block_id as an argument, construct an object to inject it.
    if has_keyword(func, "block_id"):
        extra_argpairs.append((ArrayBlockIdDep(out.chunks), out_ind))
        extra_names.append("block_id")

    if has_keyword(func, "_overlap_trim_info"):
        # Internal for map overlap to reduce size of graph
        num_chunks = out.numblocks
        block_id_dict = {
            block_id: (block_id, num_chunks)
            for block_id in product(*(range(len(c)) for c in out.chunks))
        }
        extra_argpairs.append((ArrayValuesDep(out.chunks, block_id_dict), out_ind))
        extra_names.append("_overlap_trim_info")

    # If func has block_info as an argument, construct a dict of block info
    # objects and prepare to inject it.
    if has_keyword(func, "block_info"):
        starts = {}
        num_chunks = {}
        shapes = {}

        for i, (arg, in_ind) in enumerate(argpairs):
            if in_ind is not None:
                shapes[i] = arg.shape
                if drop_axis:
                    # We concatenate along dropped axes, so we need to treat them
                    # as if there is only a single chunk.
                    starts[i] = [
                        (
                            cached_cumsum(arg.chunks[j], initial_zero=True)
                            if ind in out_ind
                            else [0, arg.shape[j]]
                        )
                        for j, ind in enumerate(in_ind)
                    ]
                    num_chunks[i] = tuple(len(s) - 1 for s in starts[i])
                else:
                    starts[i] = [
                        cached_cumsum(c, initial_zero=True) for c in arg.chunks
                    ]
                    num_chunks[i] = arg.numblocks
        out_starts = [cached_cumsum(c, initial_zero=True) for c in out.chunks]

        block_info_dict = {}
        for block_id in product(*(range(len(c)) for c in out.chunks)):
            # Get position of chunk, indexed by axis labels
            location = {out_ind[i]: loc for i, loc in enumerate(block_id)}
            info = {}
            for i, shape in shapes.items():
                # Compute chunk key in the array, taking broadcasting into
                # account. We don't directly know which dimensions are
                # broadcast, but any dimension with only one chunk can be
                # treated as broadcast.
                arr_k = tuple(
                    location.get(ind, 0) if num_chunks[i][j] > 1 else 0
                    for j, ind in enumerate(argpairs[i][1])
                )
                info[i] = {
                    "shape": shape,
                    "num-chunks": num_chunks[i],
                    "array-location": [
                        (starts[i][ij][j], starts[i][ij][j + 1])
                        for ij, j in enumerate(arr_k)
                    ],
                    "chunk-location": arr_k,
                }

            info[None] = {
                "shape": out.shape,
                "num-chunks": out.numblocks,
                "array-location": [
                    (out_starts[ij][j], out_starts[ij][j + 1])
                    for ij, j in enumerate(block_id)
                ],
                "chunk-location": block_id,
                "chunk-shape": tuple(
                    out.chunks[ij][j] for ij, j in enumerate(block_id)
                ),
                "dtype": dtype,
            }
            block_info_dict[block_id] = info

        extra_argpairs.append((ArrayValuesDep(out.chunks, block_info_dict), out_ind))
        extra_names.append("block_info")

    if extra_argpairs:
        # Rewrite the Blockwise layer. It would be nice to find a way to
        # avoid doing it twice, but it's currently needed to determine
        # out.chunks from the first pass. Since it constructs a Blockwise
        # rather than an expanded graph, it shouldn't be too expensive.
        out = blockwise(
            _pass_extra_kwargs,
            out_ind,
            func,
            None,
            tuple(extra_names),
            None,
            *concat(extra_argpairs),
            *concat(argpairs),
            name=out.name,
            dtype=out.dtype,
            concatenate=True,
            align_arrays=False,
            adjust_chunks=dict(zip(out_ind, out.chunks)),
            meta=meta,
            **kwargs,
        )

    return out


def apply_and_enforce(*args, **kwargs):
    """Apply a function, and enforce the output.ndim to match expected_ndim

    Ensures the output has the expected dimensionality."""
    func = kwargs.pop("_func")
    expected_ndim = kwargs.pop("expected_ndim")
    out = func(*args, **kwargs)
    if getattr(out, "ndim", 0) != expected_ndim:
        out_ndim = getattr(out, "ndim", 0)
        raise ValueError(
            f"Dimension mismatch: expected output of {func} "
            f"to have dims = {expected_ndim}.  Got {out_ndim} instead."
        )
    return out


def broadcast_chunks(*chunkss):
    """Construct a chunks tuple that broadcasts many chunks tuples

    >>> a = ((5, 5),)
    >>> b = ((5, 5),)
    >>> broadcast_chunks(a, b)
    ((5, 5),)

    >>> a = ((10, 10, 10), (5, 5),)
    >>> b = ((5, 5),)
    >>> broadcast_chunks(a, b)
    ((10, 10, 10), (5, 5))

    >>> a = ((10, 10, 10), (5, 5),)
    >>> b = ((1,), (5, 5),)
    >>> broadcast_chunks(a, b)
    ((10, 10, 10), (5, 5))

    >>> a = ((10, 10, 10), (5, 5),)
    >>> b = ((3, 3,), (5, 5),)
    >>> broadcast_chunks(a, b)
    Traceback (most recent call last):
        ...
    ValueError: Chunks do not align: [(10, 10, 10), (3, 3)]
    """
    if not chunkss:
        return ()
    elif len(chunkss) == 1:
        return chunkss[0]
    n = max(map(len, chunkss))
    chunkss2 = [((1,),) * (n - len(c)) + c for c in chunkss]
    result = []
    for i in range(n):
        step1 = [c[i] for c in chunkss2]
        if all(c == (1,) for c in step1):
            step2 = step1
        else:
            step2 = [c for c in step1 if c != (1,)]
        if len(set(step2)) != 1:
            raise ValueError("Chunks do not align: %s" % str(step2))
        result.append(step2[0])
    return tuple(result)


def store(
    sources: Array | Collection[Array],
    targets: ArrayLike | Delayed | Collection[ArrayLike | Delayed],
    lock: bool | Lock = True,
    regions: tuple[slice, ...] | Collection[tuple[slice, ...]] | None = None,
    compute: bool = True,
    return_stored: bool = False,
    load_stored: bool | None = None,
    **kwargs,
):
    """Store dask arrays in array-like objects, overwrite data in target

    This stores dask arrays into object that supports numpy-style setitem
    indexing.  It stores values chunk by chunk so that it does not have to
    fill up memory.  For best performance you can align the block size of
    the storage target with the block size of your array.

    If your data fits in memory then you may prefer calling
    ``np.array(myarray)`` instead.

    Parameters
    ----------

    sources: Array or collection of Arrays
    targets: array-like or Delayed or collection of array-likes and/or Delayeds
        These should support setitem syntax ``target[10:20] = ...``.
        If sources is a single item, targets must be a single item; if sources is a
        collection of arrays, targets must be a matching collection.
    lock: boolean or threading.Lock, optional
        Whether or not to lock the data stores while storing.
        Pass True (lock each file individually), False (don't lock) or a
        particular :class:`threading.Lock` object to be shared among all writes.
    regions: tuple of slices or collection of tuples of slices, optional
        Each ``region`` tuple in ``regions`` should be such that
        ``target[region].shape = source.shape``
        for the corresponding source and target in sources and targets,
        respectively. If this is a tuple, the contents will be assumed to be
        slices, so do not provide a tuple of tuples.
    compute: boolean, optional
        If true compute immediately; return :class:`dask.delayed.Delayed` otherwise.
    return_stored: boolean, optional
        Optionally return the stored result (default False).
    load_stored: boolean, optional
        Optionally return the stored result, loaded in to memory (default None).
        If None, ``load_stored`` is True if ``return_stored`` is True and
        ``compute`` is False. *This is an advanced option.*
        When False, store will return the appropriate ``target`` for each chunk that is stored.
        Directly computing this result is not what you want.
        Instead, you can use the returned ``target`` to execute followup operations to the store.
    kwargs:
        Parameters passed to compute/persist (only used if compute=True)

    Returns
    -------

    If return_stored=True
        tuple of Arrays
    If return_stored=False and compute=True
        None
    If return_stored=False and compute=False
        Delayed

    Examples
    --------

    >>> import h5py  # doctest: +SKIP
    >>> f = h5py.File('myfile.hdf5', mode='a')  # doctest: +SKIP
    >>> dset = f.create_dataset('/data', shape=x.shape,
    ...                                  chunks=x.chunks,
    ...                                  dtype='f8')  # doctest: +SKIP

    >>> store(x, dset)  # doctest: +SKIP

    Alternatively store many arrays at the same time

    >>> store([x, y, z], [dset1, dset2, dset3])  # doctest: +SKIP
    """

    if isinstance(sources, Array):
        sources = [sources]
        # There's no way to test that targets is a single array-like.
        # We need to trust the user.
        targets = [targets]  # type: ignore
    targets = cast("Collection[ArrayLike | Delayed]", targets)

    if any(not isinstance(s, Array) for s in sources):
        raise ValueError("All sources must be dask array objects")

    if len(sources) != len(targets):
        raise ValueError(
            "Different number of sources [%d] and targets [%d]"
            % (len(sources), len(targets))
        )

    if isinstance(regions, tuple) or regions is None:
        regions_list = [regions] * len(sources)
    else:
        regions_list = list(regions)
        if len(sources) != len(regions_list):
            raise ValueError(
                f"Different number of sources [{len(sources)}] and "
                f"targets [{len(targets)}] than regions [{len(regions_list)}]"
            )
    del regions

    if load_stored is None:
        load_stored = return_stored and not compute

    if lock is True:
        lock = get_scheduler_lock(collection=Array, scheduler=kwargs.get("scheduler"))

    arrays = []
    for s, t, r in zip(sources, targets, regions_list):
        slices = ArraySliceDep(s.chunks)
        arrays.append(
            s.map_blocks(
                load_store_chunk,  # type: ignore[arg-type]
                t,
                # Note: slices / BlockwiseDep have to be passed by arg, not by kwarg
                slices,
                region=r,
                lock=lock,
                return_stored=return_stored,
                load_stored=load_stored,
                name="store-map",
                meta=s._meta,
            )
        )

    if compute:
        if not return_stored:
            import dask

            dask.compute(arrays, **kwargs)
            return None
        else:
            stored_persisted = persist(*arrays, **kwargs)
            arrays = []
            for s, r in zip(stored_persisted, regions_list):
                slices = ArraySliceDep(s.chunks)
                arrays.append(
                    s.map_blocks(
                        load_chunk,  # type: ignore[arg-type]
                        # Note: slices / BlockwiseDep have to be passed by arg, not by kwarg
                        slices,
                        lock=lock,
                        region=r,
                        meta=s._meta,
                    )
                )
    if len(arrays) == 1:
        return arrays[0]
    return tuple(arrays)


def blockdims_from_blockshape(shape, chunks):
    """

    >>> blockdims_from_blockshape((10, 10), (4, 3))
    ((4, 4, 2), (3, 3, 3, 1))
    >>> blockdims_from_blockshape((10, 0), (4, 0))
    ((4, 4, 2), (0,))
    """
    if chunks is None:
        raise TypeError("Must supply chunks= keyword argument")
    if shape is None:
        raise TypeError("Must supply shape= keyword argument")
    if np.isnan(sum(shape)) or np.isnan(sum(chunks)):
        raise ValueError(
            "Array chunk sizes are unknown. shape: %s, chunks: %s%s"
            % (shape, chunks, unknown_chunk_message)
        )
    if not all(map(is_integer, chunks)):
        raise ValueError("chunks can only contain integers.")
    if not all(map(is_integer, shape)):
        raise ValueError("shape can only contain integers.")
    shape = tuple(map(int, shape))
    chunks = tuple(map(int, chunks))
    return tuple(
        ((bd,) * (d // bd) + ((d % bd,) if d % bd else ()) if d else (0,))
        for d, bd in zip(shape, chunks)
    )


def finalize(results):
    if not results:
        return concatenate3(results)
    results2 = results
    while isinstance(results2, (tuple, list)):
        if len(results2) > 1:
            return concatenate3(results)
        else:
            results2 = results2[0]

    results = unpack_singleton(results)
    # Single chunk. There is a risk that the result holds a buffer stored in the
    # graph or on a process-local Worker. Deep copy to make sure that nothing can
    # accidentally write back to it.
    try:
        return results.copy()  # numpy, sparse, scipy.sparse (any version)
    except AttributeError:
        # Not an Array API object
        return results


CHUNKS_NONE_ERROR_MESSAGE = """
You must specify a chunks= keyword argument.
This specifies the chunksize of your array blocks.

See the following documentation page for details:
  https://docs.dask.org/en/latest/array-creation.html#chunks
""".strip()


class Array(DaskMethodsMixin):
    """Parallel Dask Array

    A parallel nd-array comprised of many numpy arrays arranged in a grid.

    This constructor is for advanced uses only.  For normal use see the
    :func:`dask.array.from_array` function.

    Parameters
    ----------
    dask : dict
        Task dependency graph
    name : string
        Name of array in dask
    shape : tuple of ints
        Shape of the entire array
    chunks: iterable of tuples
        block sizes along each dimension
    dtype : str or dtype
        Typecode or data-type for the new Dask Array
    meta : empty ndarray
        empty ndarray created with same NumPy backend, ndim and dtype as the
        Dask Array being created (overrides dtype)

    See Also
    --------
    dask.array.from_array
    """

    __slots__ = "dask", "__name", "_cached_keys", "__chunks", "_meta", "__dict__"

    def __new__(cls, dask, name, chunks, dtype=None, meta=None, shape=None):
        self = super().__new__(cls)
        assert isinstance(dask, Mapping)
        if not isinstance(dask, HighLevelGraph):
            dask = HighLevelGraph.from_collections(name, dask, dependencies=())
        self.dask = dask
        self._name = str(name)
        meta = meta_from_array(meta, dtype=dtype)

        if (
            isinstance(chunks, str)
            or isinstance(chunks, tuple)
            and chunks
            and any(isinstance(c, str) for c in chunks)
        ):
            dt = meta.dtype
        else:
            dt = None
        self._chunks = normalize_chunks(chunks, shape, dtype=dt)
        if self.chunks is None:
            raise ValueError(CHUNKS_NONE_ERROR_MESSAGE)
        self._meta = meta_from_array(meta, ndim=self.ndim, dtype=dtype)

        for plugin in config.get("array_plugins", ()):
            result = plugin(self)
            if result is not None:
                self = result

        try:
            layer = self.dask.layers[name]
        except (AttributeError, KeyError):
            # self is no longer an Array after applying the plugins, OR
            # a plugin replaced the HighLevelGraph with a plain dict, OR
            # name is not the top layer's name (this can happen after the layer is
            # manipulated, to avoid a collision)
            pass
        else:
            if layer.collection_annotations is None:
                layer.collection_annotations = {
                    "shape": self.shape,
                    "dtype": self.dtype,
                    "chunksize": self.chunksize,
                    "chunks": self.chunks,
                    "type": typename(type(self)),
                    "chunk_type": typename(type(self._meta)),
                }
            else:
                layer.collection_annotations.update(
                    {
                        "shape": self.shape,
                        "dtype": self.dtype,
                        "chunksize": self.chunksize,
                        "chunks": self.chunks,
                        "type": typename(type(self)),
                        "chunk_type": typename(type(self._meta)),
                    }
                )

        return self

    def __reduce__(self):
        return (Array, (self.dask, self.name, self.chunks, self.dtype, self._meta))

    def __dask_graph__(self) -> Graph:
        return self.dask

    def __dask_layers__(self) -> Sequence[str]:
        return (self.name,)

    def __dask_keys__(self) -> NestedKeys:
        if self._cached_keys is not None:
            return self._cached_keys

        name, chunks, numblocks = self.name, self.chunks, self.numblocks

        def keys(*args):
            if not chunks:
                return [(name,)]
            ind = len(args)
            if ind + 1 == len(numblocks):
                result = [(name,) + args + (i,) for i in range(numblocks[ind])]
            else:
                result = [keys(*(args + (i,))) for i in range(numblocks[ind])]
            return result

        self._cached_keys = result = keys()
        return result

    def __dask_tokenize__(self):
        return self.name

    __dask_optimize__ = globalmethod(
        optimize, key="array_optimize", falsey=dont_optimize
    )
    __dask_scheduler__ = staticmethod(DEFAULT_GET)

    def __dask_postcompute__(self):
        return finalize, ()

    def __dask_postpersist__(self):
        return self._rebuild, ()

    def _rebuild(self, dsk, *, rename=None):
        name = self._name
        if rename:
            name = rename.get(name, name)
        return Array(dsk, name, self.chunks, self.dtype, self._meta)

    def _reset_cache(self, key=None):
        """
        Reset cached properties.

        Parameters
        ----------
        key : str, optional
            Remove specified key. The default removes all items.
        """
        if key is None:
            self.__dict__.clear()
        else:
            self.__dict__.pop(key, None)

    @cached_property
    def _key_array(self):
        return np.array(self.__dask_keys__(), dtype=object)

    @cached_property
    def numblocks(self):
        return tuple(map(len, self.chunks))

    @cached_property
    def npartitions(self):
        return reduce(mul, self.numblocks, 1)

    def compute_chunk_sizes(self):
        """
        Compute the chunk sizes for a Dask array. This is especially useful
        when the chunk sizes are unknown (e.g., when indexing one Dask array
        with another).

        Notes
        -----
        This function modifies the Dask array in-place.

        Examples
        --------
        >>> import dask.array as da
        >>> import numpy as np
        >>> x = da.from_array([-2, -1, 0, 1, 2], chunks=2)
        >>> x.chunks
        ((2, 2, 1),)
        >>> y = x[x <= 0]
        >>> y.chunks
        ((nan, nan, nan),)
        >>> y.compute_chunk_sizes()  # in-place computation
        dask.array<getitem, shape=(3,), dtype=int64, chunksize=(2,), chunktype=numpy.ndarray>
        >>> y.chunks
        ((2, 1, 0),)

        """
        x = self
        chunk_shapes = x.map_blocks(
            _get_chunk_shape,
            dtype=int,
            chunks=tuple(len(c) * (1,) for c in x.chunks) + ((x.ndim,),),
            new_axis=x.ndim,
        )

        c = []
        for i in range(x.ndim):
            s = x.ndim * [0] + [i]
            s[i] = slice(None)
            s = tuple(s)

            c.append(tuple(chunk_shapes[s]))

        # `map_blocks` assigns numpy dtypes
        # cast chunk dimensions back to python int before returning
        x._chunks = tuple(
            tuple(int(chunk) for chunk in chunks) for chunks in compute(tuple(c))[0]
        )
        return x

    @cached_property
    def shape(self) -> tuple[T_IntOrNaN, ...]:
        return tuple(cached_cumsum(c, initial_zero=True)[-1] for c in self.chunks)

    @property
    def chunksize(self) -> tuple[T_IntOrNaN, ...]:
        return tuple(cached_max(c) for c in self.chunks)

    @property
    def dtype(self):
        if isinstance(self._meta, tuple):
            dtype = self._meta[0].dtype
        else:
            dtype = self._meta.dtype
        return dtype

    @property
    def _chunks(self):
        """Non-public chunks property. Allows setting a chunk value."""
        return self.__chunks

    @_chunks.setter
    def _chunks(self, chunks):
        self.__chunks = chunks

        # When the chunks changes the cached properties that was
        # dependent on it needs to be deleted:
        for key in ["numblocks", "npartitions", "shape", "ndim", "size", "_key_array"]:
            self._reset_cache(key)

    @property
    def chunks(self):
        """Chunks property."""
        return self.__chunks

    @chunks.setter
    def chunks(self, chunks):
        raise TypeError(
            "Can not set chunks directly\n\n"
            "Please use the rechunk method instead:\n"
            f"  x.rechunk({chunks})\n\n"
            "If trying to avoid unknown chunks, use\n"
            "  x.compute_chunk_sizes()"
        )

    def __len__(self):
        if not self.chunks:
            raise TypeError("len() of unsized object")
        if np.isnan(self.chunks[0]).any():
            msg = (
                "Cannot call len() on object with unknown chunk size."
                f"{unknown_chunk_message}"
            )
            raise ValueError(msg)
        return int(sum(self.chunks[0]))

    def __array_ufunc__(self, numpy_ufunc, method, *inputs, **kwargs):
        out = kwargs.get("out", ())
        for x in inputs + out:
            if _should_delegate(self, x):
                return NotImplemented

        if method == "__call__":
            if numpy_ufunc is np.matmul:
                from dask.array.routines import matmul

                # special case until apply_gufunc handles optional dimensions
                return matmul(*inputs, **kwargs)
            if numpy_ufunc.signature is not None:
                from dask.array.gufunc import apply_gufunc

                return apply_gufunc(
                    numpy_ufunc, numpy_ufunc.signature, *inputs, **kwargs
                )
            if numpy_ufunc.nout > 1:
                from dask.array import ufunc

                try:
                    da_ufunc = getattr(ufunc, numpy_ufunc.__name__)
                except AttributeError:
                    return NotImplemented
                return da_ufunc(*inputs, **kwargs)
            else:
                return elemwise(numpy_ufunc, *inputs, **kwargs)
        elif method == "outer":
            from dask.array import ufunc

            try:
                da_ufunc = getattr(ufunc, numpy_ufunc.__name__)
            except AttributeError:
                return NotImplemented
            return da_ufunc.outer(*inputs, **kwargs)
        else:
            return NotImplemented

    def __repr__(self):
        """

        >>> import dask.array as da
        >>> da.ones((10, 10), chunks=(5, 5), dtype='i4')
        dask.array<..., shape=(10, 10), dtype=int32, chunksize=(5, 5), chunktype=numpy.ndarray>
        """
        chunksize = str(self.chunksize)
        name = self.name.rsplit("-", 1)[0]
        return (
            "dask.array<{}, shape={}, dtype={}, chunksize={}, chunktype={}.{}>".format(
                name,
                self.shape,
                self.dtype,
                chunksize,
                type(self._meta).__module__.split(".")[0],
                type(self._meta).__name__,
            )
        )

    def _repr_html_(self):
        if ARRAY_TEMPLATE is None:
            # if the jinja template is not available, (e.g. because jinja2 is not installed)
            # fall back to the textual representation
            return repr(self)

        try:
            grid = self.to_svg(size=config.get("array.svg.size", 120))
        except NotImplementedError:
            grid = ""

        if "sparse" in typename(type(self._meta)):
            nbytes = None
            cbytes = None
        elif not math.isnan(self.nbytes):
            nbytes = format_bytes(self.nbytes)
            cbytes = format_bytes(math.prod(self.chunksize) * self.dtype.itemsize)
        else:
            nbytes = "unknown"
            cbytes = "unknown"

        return ARRAY_TEMPLATE.render(
            array=self,
            grid=grid,
            nbytes=nbytes,
            cbytes=cbytes,
            layers=maybe_pluralize(len(self.dask.layers), "graph layer"),
        )

    @cached_property
    def ndim(self) -> int:
        return len(self.shape)

    @cached_property
    def size(self) -> T_IntOrNaN:
        """Number of elements in array"""
        return reduce(mul, self.shape, 1)

    @property
    def nbytes(self) -> T_IntOrNaN:
        """Number of bytes in array"""
        return self.size * self.dtype.itemsize

    @property
    def itemsize(self) -> int:
        """Length of one array element in bytes"""
        return self.dtype.itemsize

    @property
    def _name(self):
        return self.__name

    @_name.setter
    def _name(self, val):
        self.__name = val
        # Clear the key cache when the name is reset
        self._cached_keys = None
        self._reset_cache("_key_array")

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self, val):
        raise TypeError(
            "Cannot set name directly\n\n"
            "Name is used to relate the array to the task graph.\n"
            "It is uncommon to need to change it, but if you do\n"
            "please set ``._name``"
        )

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    __array_priority__ = 11  # higher than numpy.ndarray and numpy.matrix

    def __array__(self, dtype=None, copy=None, **kwargs):
        if kwargs:
            warnings.warn(
                f"Extra keyword arguments {kwargs} are ignored and won't be "
                "accepted in the future",
                FutureWarning,
            )
        if copy is False:
            warnings.warn(
                "Can't acquire a memory view of a Dask array. "
                "This will raise in the future.",
                FutureWarning,
            )

        x = self.compute()

        # Apply requested dtype and convert non-numpy backends to numpy.
        # If copy is True, numpy is going to perform its own deep copy
        # after this method returns.
        # If copy is None, finalize() ensures that the returned object
        # does not share memory with an object stored in the graph or on a
        # process-local Worker.
        return np.asarray(x, dtype=dtype)

    def __array_function__(self, func, types, args, kwargs):
        import dask.array as module

        def handle_nonmatching_names(func, args, kwargs):
            if func not in _HANDLED_FUNCTIONS:
                warnings.warn(
                    "The `{}` function is not implemented by Dask array. "
                    "You may want to use the da.map_blocks function "
                    "or something similar to silence this warning. "
                    "Your code may stop working in a future release.".format(
                        func.__module__ + "." + func.__name__
                    ),
                    FutureWarning,
                )
                # Need to convert to array object (e.g. numpy.ndarray or
                # cupy.ndarray) as needed, so we can call the NumPy function
                # again and it gets the chance to dispatch to the right
                # implementation.
                args, kwargs = compute(args, kwargs)
                return func(*args, **kwargs)

            return _HANDLED_FUNCTIONS[func](*args, **kwargs)

        # First, verify that all types are handled by Dask. Otherwise, return NotImplemented.
        if not all(
            # Accept our own superclasses as recommended by NEP-13
            # (https://numpy.org/neps/nep-0013-ufunc-overrides.html#subclass-hierarchies)
            issubclass(type(self), type_) or is_valid_chunk_type(type_)
            for type_ in types
        ):
            return NotImplemented

        # Now try to find a matching function name.  If that doesn't work, we may
        # be dealing with an alias or a function that's simply not in the Dask API.
        # Handle aliases via the _HANDLED_FUNCTIONS dict mapping, and warn otherwise.
        for submodule in func.__module__.split(".")[1:]:
            try:
                module = getattr(module, submodule)
            except AttributeError:
                return handle_nonmatching_names(func, args, kwargs)

        if not hasattr(module, func.__name__):
            return handle_nonmatching_names(func, args, kwargs)

        da_func = getattr(module, func.__name__)
        if da_func is func:
            return handle_nonmatching_names(func, args, kwargs)

        # If ``like`` is contained in ``da_func``'s signature, add ``like=self``
        # to the kwargs dictionary.
        if has_keyword(da_func, "like"):
            kwargs["like"] = self

        return da_func(*args, **kwargs)

    @property
    def _elemwise(self):
        return elemwise

    @wraps(store)
    def store(self, target, **kwargs):
        return store([self], [target], **kwargs)

    def to_svg(self, size=500):
        """Convert chunks from Dask Array into an SVG Image

        Parameters
        ----------
        chunks: tuple
        size: int
            Rough size of the image

        Examples
        --------
        >>> x.to_svg(size=500)  # doctest: +SKIP

        Returns
        -------
        text: An svg string depicting the array as a grid of chunks
        """
        from dask.array.svg import svg

        return svg(self.chunks, size=size)

    def to_hdf5(self, filename, datapath, **kwargs):
        """Store array in HDF5 file

        >>> x.to_hdf5('myfile.hdf5', '/x')  # doctest: +SKIP

        Optionally provide arguments as though to ``h5py.File.create_dataset``

        >>> x.to_hdf5('myfile.hdf5', '/x', compression='lzf', shuffle=True)  # doctest: +SKIP

        See Also
        --------
        dask.array.store
        h5py.File.create_dataset
        """
        return to_hdf5(filename, datapath, self, **kwargs)

    def to_dask_dataframe(self, columns=None, index=None, meta=None):
        """Convert dask Array to dask Dataframe

        Parameters
        ----------
        columns: list or string
            list of column names if DataFrame, single string if Series
        index : dask.dataframe.Index, optional
            An optional *dask* Index to use for the output Series or DataFrame.

            The default output index depends on whether the array has any unknown
            chunks. If there are any unknown chunks, the output has ``None``
            for all the divisions (one per chunk). If all the chunks are known,
            a default index with known divisions is created.

            Specifying ``index`` can be useful if you're conforming a Dask Array
            to an existing dask Series or DataFrame, and you would like the
            indices to match.
        meta : object, optional
            An optional `meta` parameter can be passed for dask
            to specify the concrete dataframe type to use for partitions of
            the Dask dataframe. By default, pandas DataFrame is used.

        See Also
        --------
        dask.dataframe.from_dask_array
        """
        from dask.dataframe import from_dask_array

        return from_dask_array(self, columns=columns, index=index, meta=meta)

    def to_backend(self, backend: str | None = None, **kwargs):
        """Move to a new Array backend

        Parameters
        ----------
        backend : str, Optional
            The name of the new backend to move to. The default
            is the current "array.backend" configuration.

        Returns
        -------
        Array
        """
        from dask.array.creation import to_backend

        return to_backend(self, backend=backend, **kwargs)

    def __bool__(self):
        if self.size > 1:
            raise ValueError(
                f"The truth value of a {self.__class__.__name__} is ambiguous. "
                "Use a.any() or a.all()."
            )
        else:
            return bool(self.compute())

    __nonzero__ = __bool__  # python 2

    def _scalarfunc(self, cast_type):
        if self.size > 1:
            raise TypeError("Only length-1 arrays can be converted to Python scalars")
        else:
            return cast_type(self.compute().item())

    def __int__(self):
        return self._scalarfunc(int)

    __long__ = __int__  # python 2

    def __float__(self):
        return self._scalarfunc(float)

    def __complex__(self):
        return self._scalarfunc(complex)

    def __index__(self):
        return self._scalarfunc(operator.index)

    def __setitem__(self, key, value):
        if value is np.ma.masked:
            value = np.ma.masked_all((), dtype=self.dtype)

        if not is_dask_collection(value) and self.dtype.kind in "iu":
            if np.isnan(value).any():
                raise ValueError("cannot convert float NaN to integer")
            if np.isinf(value).any():
                raise ValueError("cannot convert float infinity to integer")

        # Suppress dtype broadcasting; __setitem__ can't change the dtype.
        # Use asanyarray to retain e.g. np.ma objects.
        value = asanyarray(value, dtype=self.dtype, like=self)

        if isinstance(key, Array) and (
            key.dtype.kind in "iu"
            or (key.dtype == bool and key.ndim == 1 and self.ndim > 1)
        ):
            key = (key,)

        ## Use the "where" method for cases when key is an Array of bools
        if isinstance(key, Array):
            from dask.array.routines import where

            left_shape = np.array(key.shape)
            right_shape = np.array(self.shape)

            # We want to treat unknown shape on *either* sides as a match
            match = left_shape == right_shape
            match |= np.isnan(left_shape) | np.isnan(right_shape)

            if not match.all():
                raise IndexError(
                    f"boolean index shape {key.shape} must match indexed array's "
                    f"{self.shape}."
                )

            # If value has ndim > 0, they must be broadcastable to self.shape[idx].
            # This raises when the bool mask causes the size to become unknown,
            # e.g. this is valid in numpy but raises here:
            # x = da.array([1,2,3])
            # x[da.array([True, True, False])] = [4, 5]
            if value.ndim:
                value = broadcast_to(value, self[key].shape)

            y = where(key, value, self)
            # FIXME does any backend allow mixed ops vs. numpy?
            # If yes, is it wise to let them change the meta?
            self._meta = y._meta
            self.dask = y.dask
            self._name = y.name
            self._chunks = y.chunks
            return

        if np.isnan(self.shape).any():
            raise ValueError(f"Arrays chunk sizes are unknown. {unknown_chunk_message}")

        # Still here? Then apply the assignment to other type of
        # indices via the `setitem_array` function.

        out = "setitem-" + tokenize(self, key, value)
        dsk = setitem_array(out, self, key, value)

        meta = meta_from_array(self._meta)
        if np.isscalar(meta):
            meta = np.array(meta)

        graph = HighLevelGraph.from_collections(out, dsk, dependencies=[self])
        y = Array(graph, out, chunks=self.chunks, dtype=self.dtype, meta=meta)

        self._meta = y._meta
        self.dask = y.dask
        self._name = y.name
        self._chunks = y.chunks

    def __getitem__(self, index):
        # Field access, e.g. x['a'] or x[['a', 'b']]
        if isinstance(index, str) or (
            isinstance(index, list) and index and all(isinstance(i, str) for i in index)
        ):
            if isinstance(index, str):
                dt = self.dtype[index]
            else:
                dt = np.dtype(
                    {
                        "names": index,
                        "formats": [self.dtype.fields[name][0] for name in index],
                        "offsets": [self.dtype.fields[name][1] for name in index],
                        "itemsize": self.dtype.itemsize,
                    }
                )

            if dt.shape:
                new_axis = list(range(self.ndim, self.ndim + len(dt.shape)))
                chunks = self.chunks + tuple((i,) for i in dt.shape)
                return self.map_blocks(
                    getitem, index, dtype=dt.base, chunks=chunks, new_axis=new_axis
                )
            else:
                return self.map_blocks(getitem, index, dtype=dt)

        if not isinstance(index, tuple):
            index = (index,)

        from dask.array.slicing import (
            normalize_index,
            slice_with_bool_dask_array,
            slice_with_int_dask_array,
        )

        index2 = normalize_index(index, self.shape)
        dependencies = {self.name}
        for i in index2:
            if isinstance(i, Array):
                dependencies.add(i.name)

        if any(isinstance(i, Array) and i.dtype.kind in "iu" for i in index2):
            self, index2 = slice_with_int_dask_array(self, index2)
        if any(isinstance(i, Array) and i.dtype == bool for i in index2):
            self, index2 = slice_with_bool_dask_array(self, index2)

        if all(isinstance(i, slice) and i == slice(None) for i in index2):
            return self

        out = "getitem-" + tokenize(self, index2)
        try:
            dsk, chunks = slice_array(out, self.name, self.chunks, index2)
        except SlicingNoop:
            return self

        graph = HighLevelGraph.from_collections(out, dsk, dependencies=[self])

        meta = meta_from_array(self._meta, ndim=len(chunks))
        if np.isscalar(meta):
            meta = np.array(meta)

        return Array(graph, out, chunks, meta=meta)

    def _vindex(self, key):
        if not isinstance(key, tuple):
            key = (key,)
        if any(k is None for k in key):
            raise IndexError(
                "vindex does not support indexing with None (np.newaxis), "
                "got {}".format(key)
            )
        if all(isinstance(k, slice) for k in key):
            if all(
                k.indices(d) == slice(0, d).indices(d) for k, d in zip(key, self.shape)
            ):
                return self
            raise IndexError(
                "vindex requires at least one non-slice to vectorize over "
                "when the slices are not over the entire array (i.e, x[:]). "
                "Use normal slicing instead when only using slices. Got: {}".format(key)
            )
        elif any(is_dask_collection(k) for k in key):
            if math.prod(self.numblocks) == 1 and len(key) == 1 and self.ndim == 1:
                idxr = key[0]
                # we can broadcast in this case
                return idxr.map_blocks(
                    _numpy_vindex, self, dtype=self.dtype, chunks=idxr.chunks
                )
            else:
                raise IndexError(
                    "vindex does not support indexing with dask objects. Call compute "
                    "on the indexer first to get an evalurated array. Got: {}".format(
                        key
                    )
                )
        return _vindex(self, *key)

    @property
    def vindex(self):
        """Vectorized indexing with broadcasting.

        This is equivalent to numpy's advanced indexing, using arrays that are
        broadcast against each other. This allows for pointwise indexing:

        >>> import dask.array as da
        >>> x = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> x = da.from_array(x, chunks=2)
        >>> x.vindex[[0, 1, 2], [0, 1, 2]].compute()
        array([1, 5, 9])

        Mixed basic/advanced indexing with slices/arrays is also supported. The
        order of dimensions in the result follows those proposed for
        `ndarray.vindex <https://github.com/numpy/numpy/pull/6256>`_:
        the subspace spanned by arrays is followed by all slices.

        Note: ``vindex`` provides more general functionality than standard
        indexing, but it also has fewer optimizations and can be significantly
        slower.
        """
        return IndexCallable(self._vindex)

    @property
    def blocks(self):
        """An array-like interface to the blocks of an array.

        This returns a ``Blockview`` object that provides an array-like interface
        to the blocks of a dask array.  Numpy-style indexing of a ``Blockview`` object
        returns a selection of blocks as a new dask array.

        You can index ``array.blocks`` like a numpy array of shape
        equal to the number of blocks in each dimension, (available as
        array.blocks.size).  The dimensionality of the output array matches
        the dimension of this array, even if integer indices are passed.
        Slicing with ``np.newaxis`` or multiple lists is not supported.

        Examples
        --------
        >>> import dask.array as da
        >>> x = da.arange(8, chunks=2)
        >>> x.blocks.shape # aliases x.numblocks
        (4,)
        >>> x.blocks[0].compute()
        array([0, 1])
        >>> x.blocks[:3].compute()
        array([0, 1, 2, 3, 4, 5])
        >>> x.blocks[::2].compute()
        array([0, 1, 4, 5])
        >>> x.blocks[[-1, 0]].compute()
        array([6, 7, 0, 1])
        >>> x.blocks.ravel() # doctest: +NORMALIZE_WHITESPACE
        [dask.array<blocks, shape=(2,), dtype=int64, chunksize=(2,), chunktype=numpy.ndarray>,
         dask.array<blocks, shape=(2,), dtype=int64, chunksize=(2,), chunktype=numpy.ndarray>,
         dask.array<blocks, shape=(2,), dtype=int64, chunksize=(2,), chunktype=numpy.ndarray>,
         dask.array<blocks, shape=(2,), dtype=int64, chunksize=(2,), chunktype=numpy.ndarray>]

        Returns
        -------
        An instance of ``dask.array.Blockview``
        """
        return BlockView(self)

    @property
    def partitions(self):
        """Slice an array by partitions. Alias of dask array .blocks attribute.

        This alias allows you to write agnostic code that works with both
        dask arrays and dask dataframes.

        This returns a ``Blockview`` object that provides an array-like interface
        to the blocks of a dask array.  Numpy-style indexing of a ``Blockview`` object
        returns a selection of blocks as a new dask array.

        You can index ``array.blocks`` like a numpy array of shape
        equal to the number of blocks in each dimension, (available as
        array.blocks.size).  The dimensionality of the output array matches
        the dimension of this array, even if integer indices are passed.
        Slicing with ``np.newaxis`` or multiple lists is not supported.

        Examples
        --------
        >>> import dask.array as da
        >>> x = da.arange(8, chunks=2)
        >>> x.partitions.shape # aliases x.numblocks
        (4,)
        >>> x.partitions[0].compute()
        array([0, 1])
        >>> x.partitions[:3].compute()
        array([0, 1, 2, 3, 4, 5])
        >>> x.partitions[::2].compute()
        array([0, 1, 4, 5])
        >>> x.partitions[[-1, 0]].compute()
        array([6, 7, 0, 1])
        >>> x.partitions.ravel() # doctest: +NORMALIZE_WHITESPACE
        [dask.array<blocks, shape=(2,), dtype=int64, chunksize=(2,), chunktype=numpy.ndarray>,
         dask.array<blocks, shape=(2,), dtype=int64, chunksize=(2,), chunktype=numpy.ndarray>,
         dask.array<blocks, shape=(2,), dtype=int64, chunksize=(2,), chunktype=numpy.ndarray>,
         dask.array<blocks, shape=(2,), dtype=int64, chunksize=(2,), chunktype=numpy.ndarray>]

        Returns
        -------
        An instance of ``da.array.Blockview``
        """
        return self.blocks

    def dot(self, other):
        """Dot product of self and other.

        Refer to :func:`dask.array.tensordot` for full documentation.

        See Also
        --------
        dask.array.dot : equivalent function
        """
        from dask.array.routines import tensordot

        return tensordot(self, other, axes=((self.ndim - 1,), (other.ndim - 2,)))

    @property
    def A(self):
        return self

    @property
    def T(self):
        return self.transpose()

    def transpose(self, *axes):
        """Reverse or permute the axes of an array. Return the modified array.

        Refer to :func:`dask.array.transpose` for full documentation.

        See Also
        --------
        dask.array.transpose : equivalent function
        """
        from dask.array.routines import transpose

        if not axes:
            axes = None
        elif len(axes) == 1 and isinstance(axes[0], Iterable):
            axes = axes[0]
        if (axes == tuple(range(self.ndim))) or (axes == tuple(range(-self.ndim, 0))):
            # no transpose necessary
            return self
        else:
            return transpose(self, axes=axes)

    def ravel(self):
        """Return a flattened array.

        Refer to :func:`dask.array.ravel` for full documentation.

        See Also
        --------
        dask.array.ravel : equivalent function
        """
        from dask.array.routines import ravel

        return ravel(self)

    flatten = ravel

    def choose(self, choices):
        """Use an index array to construct a new array from a set of choices.

        Refer to :func:`dask.array.choose` for full documentation.

        See Also
        --------
        dask.array.choose : equivalent function
        """
        from dask.array.routines import choose

        return choose(self, choices)

    def reshape(self, *shape, merge_chunks=True, limit=None):
        """Reshape array to new shape

        Refer to :func:`dask.array.reshape` for full documentation.

        See Also
        --------
        dask.array.reshape : equivalent function
        """
        from dask.array.reshape import reshape

        if len(shape) == 1 and not isinstance(shape[0], Number):
            shape = shape[0]
        return reshape(self, shape, merge_chunks=merge_chunks, limit=limit)

    def topk(self, k, axis=-1, split_every=None):
        """The top k elements of an array.

        Refer to :func:`dask.array.topk` for full documentation.

        See Also
        --------
        dask.array.topk : equivalent function
        """
        from dask.array.reductions import topk

        return topk(self, k, axis=axis, split_every=split_every)

    def argtopk(self, k, axis=-1, split_every=None):
        """The indices of the top k elements of an array.

        Refer to :func:`dask.array.argtopk` for full documentation.

        See Also
        --------
        dask.array.argtopk : equivalent function
        """
        from dask.array.reductions import argtopk

        return argtopk(self, k, axis=axis, split_every=split_every)

    def astype(self, dtype, **kwargs):
        """Copy of the array, cast to a specified type.

        Parameters
        ----------
        dtype : str or dtype
            Typecode or data-type to which the array is cast.
        casting : {'no', 'equiv', 'safe', 'same_kind', 'unsafe'}, optional
            Controls what kind of data casting may occur. Defaults to 'unsafe'
            for backwards compatibility.

            * 'no' means the data types should not be cast at all.
            * 'equiv' means only byte-order changes are allowed.
            * 'safe' means only casts which can preserve values are allowed.
            * 'same_kind' means only safe casts or casts within a kind,
                like float64 to float32, are allowed.
            * 'unsafe' means any data conversions may be done.
        copy : bool, optional
            By default, astype always returns a newly allocated array. If this
            is set to False and the `dtype` requirement is satisfied, the input
            array is returned instead of a copy.

            .. note::

                Dask does not respect the contiguous memory layout of the array,
                and will ignore the ``order`` keyword argument.
                The default order is 'C' contiguous.
        """
        kwargs.pop("order", None)  # `order` is not respected, so we remove this kwarg
        # Scalars don't take `casting` or `copy` kwargs - as such we only pass
        # them to `map_blocks` if specified by user (different than defaults).
        extra = set(kwargs) - {"casting", "copy"}
        if extra:
            raise TypeError(
                f"astype does not take the following keyword arguments: {list(extra)}"
            )
        casting = kwargs.get("casting", "unsafe")
        dtype = np.dtype(dtype)
        if self.dtype == dtype:
            return self
        elif not np.can_cast(self.dtype, dtype, casting=casting):
            raise TypeError(
                f"Cannot cast array from {self.dtype!r} to {dtype!r} "
                f"according to the rule {casting!r}"
            )
        return self.map_blocks(chunk.astype, dtype=dtype, astype_dtype=dtype, **kwargs)

    def __abs__(self):
        return elemwise(operator.abs, self)

    @check_if_handled_given_other
    def __add__(self, other):
        return elemwise(operator.add, self, other)

    @check_if_handled_given_other
    def __radd__(self, other):
        return elemwise(operator.add, other, self)

    @check_if_handled_given_other
    def __and__(self, other):
        return elemwise(operator.and_, self, other)

    @check_if_handled_given_other
    def __rand__(self, other):
        return elemwise(operator.and_, other, self)

    @check_if_handled_given_other
    def __div__(self, other):
        return elemwise(operator.div, self, other)

    @check_if_handled_given_other
    def __rdiv__(self, other):
        return elemwise(operator.div, other, self)

    @check_if_handled_given_other
    def __eq__(self, other):
        return elemwise(operator.eq, self, other)

    @check_if_handled_given_other
    def __gt__(self, other):
        return elemwise(operator.gt, self, other)

    @check_if_handled_given_other
    def __ge__(self, other):
        return elemwise(operator.ge, self, other)

    def __invert__(self):
        return elemwise(operator.invert, self)

    @check_if_handled_given_other
    def __lshift__(self, other):
        return elemwise(operator.lshift, self, other)

    @check_if_handled_given_other
    def __rlshift__(self, other):
        return elemwise(operator.lshift, other, self)

    @check_if_handled_given_other
    def __lt__(self, other):
        return elemwise(operator.lt, self, other)

    @check_if_handled_given_other
    def __le__(self, other):
        return elemwise(operator.le, self, other)

    @check_if_handled_given_other
    def __mod__(self, other):
        return elemwise(operator.mod, self, other)

    @check_if_handled_given_other
    def __rmod__(self, other):
        return elemwise(operator.mod, other, self)

    @check_if_handled_given_other
    def __mul__(self, other):
        return elemwise(operator.mul, self, other)

    @check_if_handled_given_other
    def __rmul__(self, other):
        return elemwise(operator.mul, other, self)

    @check_if_handled_given_other
    def __ne__(self, other):
        return elemwise(operator.ne, self, other)

    def __neg__(self):
        return elemwise(operator.neg, self)

    @check_if_handled_given_other
    def __or__(self, other):
        return elemwise(operator.or_, self, other)

    def __pos__(self):
        return self

    @check_if_handled_given_other
    def __ror__(self, other):
        return elemwise(operator.or_, other, self)

    @check_if_handled_given_other
    def __pow__(self, other):
        return elemwise(operator.pow, self, other)

    @check_if_handled_given_other
    def __rpow__(self, other):
        return elemwise(operator.pow, other, self)

    @check_if_handled_given_other
    def __rshift__(self, other):
        return elemwise(operator.rshift, self, other)

    @check_if_handled_given_other
    def __rrshift__(self, other):
        return elemwise(operator.rshift, other, self)

    @check_if_handled_given_other
    def __sub__(self, other):
        return elemwise(operator.sub, self, other)

    @check_if_handled_given_other
    def __rsub__(self, other):
        return elemwise(operator.sub, other, self)

    @check_if_handled_given_other
    def __truediv__(self, other):
        return elemwise(operator.truediv, self, other)

    @check_if_handled_given_other
    def __rtruediv__(self, other):
        return elemwise(operator.truediv, other, self)

    @check_if_handled_given_other
    def __floordiv__(self, other):
        return elemwise(operator.floordiv, self, other)

    @check_if_handled_given_other
    def __rfloordiv__(self, other):
        return elemwise(operator.floordiv, other, self)

    @check_if_handled_given_other
    def __xor__(self, other):
        return elemwise(operator.xor, self, other)

    @check_if_handled_given_other
    def __rxor__(self, other):
        return elemwise(operator.xor, other, self)

    @check_if_handled_given_other
    def __matmul__(self, other):
        from dask.array.routines import matmul

        return matmul(self, other)

    @check_if_handled_given_other
    def __rmatmul__(self, other):
        from dask.array.routines import matmul

        return matmul(other, self)

    @check_if_handled_given_other
    def __divmod__(self, other):
        from dask.array.ufunc import divmod

        return divmod(self, other)

    @check_if_handled_given_other
    def __rdivmod__(self, other):
        from dask.array.ufunc import divmod

        return divmod(other, self)

    def any(self, axis=None, keepdims=False, split_every=None, out=None):
        """Returns True if any of the elements evaluate to True.

        Refer to :func:`dask.array.any` for full documentation.

        See Also
        --------
        dask.array.any : equivalent function
        """
        from dask.array.reductions import any

        return any(self, axis=axis, keepdims=keepdims, split_every=split_every, out=out)

    def all(self, axis=None, keepdims=False, split_every=None, out=None):
        """Returns True if all elements evaluate to True.

        Refer to :func:`dask.array.all` for full documentation.

        See Also
        --------
        dask.array.all : equivalent function
        """
        from dask.array.reductions import all

        return all(self, axis=axis, keepdims=keepdims, split_every=split_every, out=out)

    def min(self, axis=None, keepdims=False, split_every=None, out=None):
        """Return the minimum along a given axis.

        Refer to :func:`dask.array.min` for full documentation.

        See Also
        --------
        dask.array.min : equivalent function
        """
        from dask.array.reductions import min

        return min(self, axis=axis, keepdims=keepdims, split_every=split_every, out=out)

    def max(self, axis=None, keepdims=False, split_every=None, out=None):
        """Return the maximum along a given axis.

        Refer to :func:`dask.array.max` for full documentation.

        See Also
        --------
        dask.array.max : equivalent function
        """
        from dask.array.reductions import max

        return max(self, axis=axis, keepdims=keepdims, split_every=split_every, out=out)

    def argmin(self, axis=None, *, keepdims=False, split_every=None, out=None):
        """Return indices of the minimum values along the given axis.

        Refer to :func:`dask.array.argmin` for full documentation.

        See Also
        --------
        dask.array.argmin : equivalent function
        """
        from dask.array.reductions import argmin

        return argmin(
            self, axis=axis, keepdims=keepdims, split_every=split_every, out=out
        )

    def argmax(self, axis=None, *, keepdims=False, split_every=None, out=None):
        """Return indices of the maximum values along the given axis.

        Refer to :func:`dask.array.argmax` for full documentation.

        See Also
        --------
        dask.array.argmax : equivalent function
        """
        from dask.array.reductions import argmax

        return argmax(
            self, axis=axis, keepdims=keepdims, split_every=split_every, out=out
        )

    def sum(self, axis=None, dtype=None, keepdims=False, split_every=None, out=None):
        """
        Return the sum of the array elements over the given axis.

        Refer to :func:`dask.array.sum` for full documentation.

        See Also
        --------
        dask.array.sum : equivalent function
        """
        from dask.array.reductions import sum

        return sum(
            self,
            axis=axis,
            dtype=dtype,
            keepdims=keepdims,
            split_every=split_every,
            out=out,
        )

    def trace(self, offset=0, axis1=0, axis2=1, dtype=None):
        """Return the sum along diagonals of the array.

        Refer to :func:`dask.array.trace` for full documentation.

        See Also
        --------
        dask.array.trace : equivalent function
        """
        from dask.array.reductions import trace

        return trace(self, offset=offset, axis1=axis1, axis2=axis2, dtype=dtype)

    def prod(self, axis=None, dtype=None, keepdims=False, split_every=None, out=None):
        """Return the product of the array elements over the given axis

        Refer to :func:`dask.array.prod` for full documentation.

        See Also
        --------
        dask.array.prod : equivalent function
        """
        from dask.array.reductions import prod

        return prod(
            self,
            axis=axis,
            dtype=dtype,
            keepdims=keepdims,
            split_every=split_every,
            out=out,
        )

    def mean(self, axis=None, dtype=None, keepdims=False, split_every=None, out=None):
        """Returns the average of the array elements along given axis.

        Refer to :func:`dask.array.mean` for full documentation.

        See Also
        --------
        dask.array.mean : equivalent function
        """
        from dask.array.reductions import mean

        return mean(
            self,
            axis=axis,
            dtype=dtype,
            keepdims=keepdims,
            split_every=split_every,
            out=out,
        )

    def std(
        self, axis=None, dtype=None, keepdims=False, ddof=0, split_every=None, out=None
    ):
        """Returns the standard deviation of the array elements along given axis.

        Refer to :func:`dask.array.std` for full documentation.

        See Also
        --------
        dask.array.std : equivalent function
        """
        from dask.array.reductions import std

        return std(
            self,
            axis=axis,
            dtype=dtype,
            keepdims=keepdims,
            ddof=ddof,
            split_every=split_every,
            out=out,
        )

    def var(
        self, axis=None, dtype=None, keepdims=False, ddof=0, split_every=None, out=None
    ):
        """Returns the variance of the array elements, along given axis.

        Refer to :func:`dask.array.var` for full documentation.

        See Also
        --------
        dask.array.var : equivalent function
        """
        from dask.array.reductions import var

        return var(
            self,
            axis=axis,
            dtype=dtype,
            keepdims=keepdims,
            ddof=ddof,
            split_every=split_every,
            out=out,
        )

    def moment(
        self,
        order,
        axis=None,
        dtype=None,
        keepdims=False,
        ddof=0,
        split_every=None,
        out=None,
    ):
        """Calculate the nth centralized moment.

        Refer to :func:`dask.array.moment` for the full documentation.

        See Also
        --------
        dask.array.moment : equivalent function
        """
        from dask.array.reductions import moment

        return moment(
            self,
            order,
            axis=axis,
            dtype=dtype,
            keepdims=keepdims,
            ddof=ddof,
            split_every=split_every,
            out=out,
        )

    @wraps(map_blocks)
    def map_blocks(self, func, *args, **kwargs):
        return map_blocks(func, self, *args, **kwargs)

    def map_overlap(self, func, depth, boundary=None, trim=True, **kwargs):
        """Map a function over blocks of the array with some overlap

        Refer to :func:`dask.array.map_overlap` for full documentation.

        See Also
        --------
        dask.array.map_overlap : equivalent function
        """
        from dask.array.overlap import map_overlap

        return map_overlap(
            func, self, depth=depth, boundary=boundary, trim=trim, **kwargs
        )

    def cumsum(self, axis, dtype=None, out=None, *, method="sequential"):
        """Return the cumulative sum of the elements along the given axis.

        Refer to :func:`dask.array.cumsum` for full documentation.

        See Also
        --------
        dask.array.cumsum : equivalent function
        """
        from dask.array.reductions import cumsum

        return cumsum(self, axis, dtype, out=out, method=method)

    def cumprod(self, axis, dtype=None, out=None, *, method="sequential"):
        """Return the cumulative product of the elements along the given axis.

        Refer to :func:`dask.array.cumprod` for full documentation.

        See Also
        --------
        dask.array.cumprod : equivalent function
        """
        from dask.array.reductions import cumprod

        return cumprod(self, axis, dtype, out=out, method=method)

    def squeeze(self, axis=None):
        """Remove axes of length one from array.

        Refer to :func:`dask.array.squeeze` for full documentation.

        See Also
        --------
        dask.array.squeeze : equivalent function
        """
        from dask.array.routines import squeeze

        return squeeze(self, axis)

    def rechunk(
        self,
        chunks="auto",
        threshold=None,
        block_size_limit=None,
        balance=False,
        method=None,
    ):
        """Convert blocks in dask array x for new chunks.

        Refer to :func:`dask.array.rechunk` for full documentation.

        See Also
        --------
        dask.array.rechunk : equivalent function
        """
        from dask.array.rechunk import rechunk  # avoid circular import

        return rechunk(self, chunks, threshold, block_size_limit, balance, method)

    def shuffle(
        self,
        indexer: list[list[int]],
        axis: int,
        chunks: Literal["auto"] = "auto",
    ):
        """Reorders one dimensions of a Dask Array based on an indexer.

        Refer to :func:`dask.array.shuffle` for full documentation.

        See Also
        --------
        dask.array.shuffle : equivalent function
        """
        from dask.array._shuffle import shuffle

        return shuffle(self, indexer, axis, chunks)

    @property
    def real(self):
        from dask.array.ufunc import real

        return real(self)

    @property
    def imag(self):
        from dask.array.ufunc import imag

        return imag(self)

    def conj(self):
        """Complex-conjugate all elements.

        Refer to :func:`dask.array.conj` for full documentation.

        See Also
        --------
        dask.array.conj : equivalent function
        """
        from dask.array.ufunc import conj

        return conj(self)

    def clip(self, min=None, max=None):
        """Return an array whose values are limited to ``[min, max]``.
        One of max or min must be given.

        Refer to :func:`dask.array.clip` for full documentation.

        See Also
        --------
        dask.array.clip : equivalent function
        """
        from dask.array.ufunc import clip

        return clip(self, min, max)

    def view(self, dtype=None, order="C"):
        """Get a view of the array as a new data type

        Parameters
        ----------
        dtype:
            The dtype by which to view the array.
            The default, None, results in the view having the same data-type
            as the original array.
        order: string
            'C' or 'F' (Fortran) ordering

        This reinterprets the bytes of the array under a new dtype.  If that
        dtype does not have the same size as the original array then the shape
        will change.

        Beware that both numpy and dask.array can behave oddly when taking
        shape-changing views of arrays under Fortran ordering.  Under some
        versions of NumPy this function will fail when taking shape-changing
        views of Fortran ordered arrays if the first dimension has chunks of
        size one.
        """
        if dtype is None:
            dtype = self.dtype
        else:
            dtype = np.dtype(dtype)
        mult = self.dtype.itemsize / dtype.itemsize

        if order == "C":
            chunks = self.chunks[:-1] + (
                tuple(ensure_int(c * mult) for c in self.chunks[-1]),
            )
        elif order == "F":
            chunks = (
                tuple(ensure_int(c * mult) for c in self.chunks[0]),
            ) + self.chunks[1:]
        else:
            raise ValueError("Order must be one of 'C' or 'F'")

        return self.map_blocks(
            chunk.view, dtype, order=order, dtype=dtype, chunks=chunks
        )

    def swapaxes(self, axis1, axis2):
        """Return a view of the array with ``axis1`` and ``axis2`` interchanged.

        Refer to :func:`dask.array.swapaxes` for full documentation.

        See Also
        --------
        dask.array.swapaxes : equivalent function
        """
        from dask.array.routines import swapaxes

        return swapaxes(self, axis1, axis2)

    def round(self, decimals=0):
        """Return array with each element rounded to the given number of decimals.

        Refer to :func:`dask.array.round` for full documentation.

        See Also
        --------
        dask.array.round : equivalent function
        """
        from dask.array.routines import round

        return round(self, decimals=decimals)

    def copy(self):
        """
        Copy array.  This is a no-op for dask.arrays, which are immutable
        """
        return Array(self.dask, self.name, self.chunks, meta=self)

    def __deepcopy__(self, memo):
        c = self.copy()
        memo[id(self)] = c
        return c

    def to_delayed(self, optimize_graph=True):
        """Convert into an array of :class:`dask.delayed.Delayed` objects, one per chunk.

        Parameters
        ----------
        optimize_graph : bool, optional
            If True [default], the graph is optimized before converting into
            :class:`dask.delayed.Delayed` objects.

        See Also
        --------
        dask.array.from_delayed
        """
        keys = self.__dask_keys__()
        graph = self.__dask_graph__()
        layer = self.__dask_layers__()[0]
        if optimize_graph:
            graph = self.__dask_optimize__(graph, keys)  # TODO, don't collape graph
            layer = "delayed-" + self.name
            graph = HighLevelGraph.from_collections(layer, graph, dependencies=())
        L = ndeepmap(self.ndim, lambda k: Delayed(k, graph, layer=layer), keys)
        return np.array(L, dtype=object)

    def repeat(self, repeats, axis=None):
        """Repeat elements of an array.

        Refer to :func:`dask.array.repeat` for full documentation.

        See Also
        --------
        dask.array.repeat : equivalent function
        """
        from dask.array.creation import repeat

        return repeat(self, repeats, axis=axis)

    def nonzero(self):
        """Return the indices of the elements that are non-zero.

        Refer to :func:`dask.array.nonzero` for full documentation.

        See Also
        --------
        dask.array.nonzero : equivalent function
        """
        from dask.array.routines import nonzero

        return nonzero(self)

    def to_zarr(self, *args, **kwargs):
        """Save array to the zarr storage format

        See https://zarr.readthedocs.io for details about the format.

        Refer to :func:`dask.array.to_zarr` for full documentation.

        See also
        --------
        dask.array.to_zarr : equivalent function
        """
        return to_zarr(self, *args, **kwargs)

    def to_tiledb(self, uri, *args, **kwargs):
        """Save array to the TileDB storage manager

        See https://docs.tiledb.io for details about the format and engine.

        See function :func:`dask.array.to_tiledb` for argument documentation.

        See also
        --------
        dask.array.to_tiledb : equivalent function
        """
        from dask.array.tiledb_io import to_tiledb

        return to_tiledb(self, uri, *args, **kwargs)


def ensure_int(f):
    i = int(f)
    if i != f:
        raise ValueError("Could not coerce %f to integer" % f)
    return i


@functools.lru_cache
def normalize_chunks_cached(
    chunks, shape=None, limit=None, dtype=None, previous_chunks=None
):
    """Cached version of normalize_chunks.

    .. note::

        chunks and previous_chunks are expected to be hashable. Dicts and lists aren't
        allowed for this function.

    See :func:`normalize_chunks` for further documentation.
    """
    return normalize_chunks(
        chunks, shape=shape, limit=limit, dtype=dtype, previous_chunks=previous_chunks
    )


def normalize_chunks(chunks, shape=None, limit=None, dtype=None, previous_chunks=None):
    """Normalize chunks to tuple of tuples

    This takes in a variety of input types and information and produces a full
    tuple-of-tuples result for chunks, suitable to be passed to Array or
    rechunk or any other operation that creates a Dask array.

    Parameters
    ----------
    chunks: tuple, int, dict, or string
        The chunks to be normalized.  See examples below for more details
    shape: Tuple[int]
        The shape of the array
    limit: int (optional)
        The maximum block size to target in bytes,
        if freedom is given to choose
    dtype: np.dtype
    previous_chunks: Tuple[Tuple[int]] optional
        Chunks from a previous array that we should use for inspiration when
        rechunking auto dimensions.  If not provided but auto-chunking exists
        then auto-dimensions will prefer square-like chunk shapes.

    Examples
    --------
    Fully explicit tuple-of-tuples

    >>> from dask.array.core import normalize_chunks
    >>> normalize_chunks(((2, 2, 1), (2, 2, 2)), shape=(5, 6))
    ((2, 2, 1), (2, 2, 2))

    Specify uniform chunk sizes

    >>> normalize_chunks((2, 2), shape=(5, 6))
    ((2, 2, 1), (2, 2, 2))

    Cleans up missing outer tuple

    >>> normalize_chunks((3, 2), (5,))
    ((3, 2),)

    Cleans up lists to tuples

    >>> normalize_chunks([[2, 2], [3, 3]])
    ((2, 2), (3, 3))

    Expands integer inputs 10 -> (10, 10)

    >>> normalize_chunks(10, shape=(30, 5))
    ((10, 10, 10), (5,))

    Expands dict inputs

    >>> normalize_chunks({0: 2, 1: 3}, shape=(6, 6))
    ((2, 2, 2), (3, 3))

    The values -1 and None get mapped to full size

    >>> normalize_chunks((5, -1), shape=(10, 10))
    ((5, 5), (10,))
    >>> normalize_chunks((5, None), shape=(10, 10))
    ((5, 5), (10,))

    Use the value "auto" to automatically determine chunk sizes along certain
    dimensions.  This uses the ``limit=`` and ``dtype=`` keywords to
    determine how large to make the chunks.  The term "auto" can be used
    anywhere an integer can be used.  See array chunking documentation for more
    information.

    >>> normalize_chunks(("auto",), shape=(20,), limit=5, dtype='uint8')
    ((5, 5, 5, 5),)
    >>> normalize_chunks("auto", (2, 3), dtype=np.int32)
    ((2,), (3,))

    You can also use byte sizes (see :func:`dask.utils.parse_bytes`) in place of
    "auto" to ask for a particular size

    >>> normalize_chunks("1kiB", shape=(2000,), dtype='float32')
    ((256, 256, 256, 256, 256, 256, 256, 208),)

    Respects null dimensions

    >>> normalize_chunks(())
    ()
    >>> normalize_chunks((), ())
    ()
    >>> normalize_chunks((1,), ())
    ()
    >>> normalize_chunks((), shape=(0, 0))
    ((0,), (0,))

    Handles NaNs

    >>> normalize_chunks((1, (np.nan,)), (1, np.nan))
    ((1,), (nan,))
    """
    if dtype and not isinstance(dtype, np.dtype):
        dtype = np.dtype(dtype)
    if chunks is None:
        raise ValueError(CHUNKS_NONE_ERROR_MESSAGE)
    if isinstance(chunks, list):
        chunks = tuple(chunks)
    if isinstance(chunks, (Number, str)):
        chunks = (chunks,) * len(shape)
    if isinstance(chunks, dict):
        chunks = tuple(chunks.get(i, None) for i in range(len(shape)))
    if isinstance(chunks, np.ndarray):
        chunks = chunks.tolist()
    if not chunks and shape and all(s == 0 for s in shape):
        chunks = ((0,),) * len(shape)

    if (
        shape
        and len(shape) == 1
        and len(chunks) > 1
        and all(isinstance(c, (Number, str)) for c in chunks)
    ):
        chunks = (chunks,)

    if shape and len(chunks) != len(shape):
        raise ValueError(
            "Chunks and shape must be of the same length/dimension. "
            "Got chunks=%s, shape=%s" % (chunks, shape)
        )
    if -1 in chunks or None in chunks:
        chunks = tuple(s if c == -1 or c is None else c for c, s in zip(chunks, shape))

    # If specifying chunk size in bytes, use that value to set the limit.
    # Verify there is only one consistent value of limit or chunk-bytes used.
    for c in chunks:
        if isinstance(c, str) and c != "auto":
            parsed = parse_bytes(c)
            if limit is None:
                limit = parsed
            elif parsed != limit:
                raise ValueError(
                    "Only one consistent value of limit or chunk is allowed."
                    "Used %s != %s" % (parsed, limit)
                )
    # Substitute byte limits with 'auto' now that limit is set.
    chunks = tuple("auto" if isinstance(c, str) and c != "auto" else c for c in chunks)

    if any(c == "auto" for c in chunks):
        chunks = auto_chunks(chunks, shape, limit, dtype, previous_chunks)

    allints = None
    if chunks and shape is not None:
        # allints: did we start with chunks as a simple tuple of ints?
        allints = all(isinstance(c, int) for c in chunks)
        chunks = _convert_int_chunk_to_tuple(shape, chunks)
    for c in chunks:
        if not c:
            raise ValueError(
                "Empty tuples are not allowed in chunks. Express "
                "zero length dimensions with 0(s) in chunks"
            )

    if not allints and shape is not None:
        if not all(
            c == s or (math.isnan(c) or math.isnan(s))
            for c, s in zip(map(sum, chunks), shape)
        ):
            raise ValueError(
                "Chunks do not add up to shape. "
                "Got chunks=%s, shape=%s" % (chunks, shape)
            )
    if allints or isinstance(sum(sum(_) for _ in chunks), int):
        # Fastpath for when we already know chunks contains only integers
        return tuple(tuple(ch) for ch in chunks)
    return tuple(
        tuple(int(x) if not math.isnan(x) else np.nan for x in c) for c in chunks
    )


def _convert_int_chunk_to_tuple(shape, chunks):
    return sum(
        (
            (
                blockdims_from_blockshape((s,), (c,))
                if not isinstance(c, (tuple, list))
                else (c,)
            )
            for s, c in zip(shape, chunks)
        ),
        (),
    )


def _compute_multiplier(limit: int, dtype, largest_block: int, result):
    """
    Utility function for auto_chunk, to fin how much larger or smaller the ideal
    chunk size is relative to what we have now.
    """
    return (
        limit
        / dtype.itemsize
        / largest_block
        / math.prod(max(r) if isinstance(r, tuple) else r for r in result.values() if r)
    )


def auto_chunks(chunks, shape, limit, dtype, previous_chunks=None):
    """Determine automatic chunks

    This takes in a chunks value that contains ``"auto"`` values in certain
    dimensions and replaces those values with concrete dimension sizes that try
    to get chunks to be of a certain size in bytes, provided by the ``limit=``
    keyword.  If multiple dimensions are marked as ``"auto"`` then they will
    all respond to meet the desired byte limit, trying to respect the aspect
    ratio of their dimensions in ``previous_chunks=``, if given.

    Parameters
    ----------
    chunks: Tuple
        A tuple of either dimensions or tuples of explicit chunk dimensions
        Some entries should be "auto"
    shape: Tuple[int]
    limit: int, str
        The maximum allowable size of a chunk in bytes
    previous_chunks: Tuple[Tuple[int]]

    See also
    --------
    normalize_chunks: for full docstring and parameters
    """
    if previous_chunks is not None:
        # rioxarray is passing ((1, ), (x,)) for shapes like (100, 5x),
        # so add this compat code for now
        # https://github.com/corteva/rioxarray/pull/820
        previous_chunks = (
            c[0] if isinstance(c, tuple) and len(c) == 1 else c for c in previous_chunks
        )
        previous_chunks = _convert_int_chunk_to_tuple(shape, previous_chunks)
    chunks = list(chunks)

    autos = {i for i, c in enumerate(chunks) if c == "auto"}
    if not autos:
        return tuple(chunks)

    if limit is None:
        limit = config.get("array.chunk-size")
    if isinstance(limit, str):
        limit = parse_bytes(limit)

    if dtype is None:
        raise TypeError("dtype must be known for auto-chunking")

    if dtype.hasobject:
        raise NotImplementedError(
            "Can not use auto rechunking with object dtype. "
            "We are unable to estimate the size in bytes of object data"
        )

    for x in tuple(chunks) + tuple(shape):
        if (
            isinstance(x, Number)
            and np.isnan(x)
            or isinstance(x, tuple)
            and np.isnan(x).any()
        ):
            raise ValueError(
                "Can not perform automatic rechunking with unknown "
                "(nan) chunk sizes.%s" % unknown_chunk_message
            )

    limit = max(1, limit)
    chunksize_tolerance = config.get("array.chunk-size-tolerance")

    largest_block = math.prod(
        cs if isinstance(cs, Number) else max(cs) for cs in chunks if cs != "auto"
    )

    if previous_chunks:
        # Base ideal ratio on the median chunk size of the previous chunks
        median_chunks = {a: np.median(previous_chunks[a]) for a in autos}
        result = {}

        # How much larger or smaller the ideal chunk size is relative to what we have now
        multiplier = _compute_multiplier(limit, dtype, largest_block, median_chunks)
        if multiplier < 1:
            # we want to update inplace, algorithm relies on it in this case
            result = median_chunks

        ideal_shape = []
        for i, s in enumerate(shape):
            chunk_frequencies = frequencies(previous_chunks[i])
            mode, count = max(chunk_frequencies.items(), key=lambda kv: kv[1])
            if mode > 1 and count >= len(previous_chunks[i]) / 2:
                ideal_shape.append(mode)
            else:
                ideal_shape.append(s)

        def _trivial_aggregate(a):
            autos.remove(a)
            del median_chunks[a]
            return True

        multiplier_remaining = True
        reduce_case = multiplier < 1
        while multiplier_remaining:  # while things change
            last_autos = set(autos)  # record previous values
            multiplier_remaining = False

            # Expand or contract each of the dimensions appropriately
            for a in sorted(autos):
                this_multiplier = multiplier ** (1 / len(last_autos))

                proposed = median_chunks[a] * this_multiplier
                this_chunksize_tolerance = chunksize_tolerance ** (1 / len(last_autos))
                max_chunk_size = proposed * this_chunksize_tolerance

                if proposed > shape[a]:  # we've hit the shape boundary
                    chunks[a] = shape[a]
                    multiplier_remaining = _trivial_aggregate(a)
                    largest_block *= shape[a]
                    result[a] = (shape[a],)
                    continue
                elif reduce_case or max(previous_chunks[a]) > max_chunk_size:
                    result[a] = round_to(proposed, ideal_shape[a])
                    if proposed < 1:
                        multiplier_remaining = True
                        autos.discard(a)
                    continue
                else:
                    dimension_result, new_chunk = [], 0
                    for c in previous_chunks[a]:
                        if c + new_chunk <= proposed:
                            # keep increasing the chunk
                            new_chunk += c
                        else:
                            # We reach the boundary so start a new chunk
                            if new_chunk > 0:
                                dimension_result.append(new_chunk)
                            new_chunk = c
                    if new_chunk > 0:
                        dimension_result.append(new_chunk)

                result[a] = tuple(dimension_result)

            # recompute how much multiplier we have left, repeat
            if multiplier_remaining or reduce_case:
                last_multiplier = multiplier
                multiplier = _compute_multiplier(
                    limit, dtype, largest_block, median_chunks
                )
                if multiplier != last_multiplier:
                    multiplier_remaining = True

        for k, v in result.items():
            chunks[k] = v if v else 0
        return tuple(chunks)

    else:
        # Check if dtype.itemsize is greater than 0
        if dtype.itemsize == 0:
            raise ValueError(
                "auto-chunking with dtype.itemsize == 0 is not supported, please pass in `chunks` explicitly"
            )
        size = (limit / dtype.itemsize / largest_block) ** (1 / len(autos))
        small = [i for i in autos if shape[i] < size]
        if small:
            for i in small:
                chunks[i] = (shape[i],)
            return auto_chunks(chunks, shape, limit, dtype)

        for i in autos:
            chunks[i] = round_to(size, shape[i])

        return tuple(chunks)


def _split_up_single_chunk(
    c: int, this_chunksize_tolerance: float, max_chunk_size: int, proposed: int
) -> list[int]:
    # Calculate by what factor we have to split this chunk
    m = c / proposed
    if math.ceil(m) / m > this_chunksize_tolerance:
        # We want to smooth things potentially if rounding up would change the result
        # by a lot
        m = math.ceil(c / max_chunk_size)
    else:
        m = math.ceil(m)
    # split the chunk
    new_c, remainder = divmod(c, min(m, c))
    x = [new_c] * min(m, c)
    for i in range(remainder):
        x[i] += 1
    return x


def round_to(c, s):
    """Return a chunk dimension that is close to an even multiple or factor

    We want values for c that are nicely aligned with s.

    If c is smaller than s we use the original chunk size and accept an
    uneven chunk at the end.

    If c is larger than s then we want the largest multiple of s that is still
    smaller than c.
    """
    if c <= s:
        return max(1, int(c))
    else:
        return c // s * s


def _get_chunk_shape(a):
    s = np.asarray(a.shape, dtype=int)
    return s[len(s) * (None,) + (slice(None),)]


def from_array(
    x,
    chunks="auto",
    name=None,
    lock=False,
    asarray=None,
    fancy=True,
    getitem=None,
    meta=None,
    inline_array=False,
):
    """Create dask array from something that looks like an array.

    Input must have a ``.shape``, ``.ndim``, ``.dtype`` and support numpy-style slicing.

    Parameters
    ----------
    x : array_like
    chunks : int, tuple
        How to chunk the array. Must be one of the following forms:

        - A blocksize like 1000.
        - A blockshape like (1000, 1000).
        - Explicit sizes of all blocks along all dimensions like
          ((1000, 1000, 500), (400, 400)).
        - A size in bytes, like "100 MiB" which will choose a uniform
          block-like shape
        - The word "auto" which acts like the above, but uses a configuration
          value ``array.chunk-size`` for the chunk size

        -1 or None as a blocksize indicate the size of the corresponding
        dimension.
    name : str or bool, optional
        The key name to use for the array. Defaults to a hash of ``x``.

        Hashing is useful if the same value of ``x`` is used to create multiple
        arrays, as Dask can then recognise that they're the same and
        avoid duplicate computations. However, it can also be slow, and if the
        array is not contiguous it is copied for hashing. If the array uses
        stride tricks (such as :func:`numpy.broadcast_to` or
        :func:`skimage.util.view_as_windows`) to have a larger logical
        than physical size, this copy can cause excessive memory usage.

        If you don't need the deduplication provided by hashing, use
        ``name=False`` to generate a random name instead of hashing, which
        avoids the pitfalls described above. Using ``name=True`` is
        equivalent to the default.

        By default, hashing uses python's standard sha1. This behaviour can be
        changed by installing cityhash, xxhash or murmurhash. If installed,
        a large-factor speedup can be obtained in the tokenisation step.

        .. note::

           Because this ``name`` is used as the key in task graphs, you should
           ensure that it uniquely identifies the data contained within. If
           you'd like to provide a descriptive name that is still unique, combine
           the descriptive name with :func:`dask.base.tokenize` of the
           ``array_like``. See :ref:`graphs` for more.

    lock : bool or Lock, optional
        If ``x`` doesn't support concurrent reads then provide a lock here, or
        pass in True to have dask.array create one for you.
    asarray : bool, optional
        If True then call np.asarray on chunks to convert them to numpy arrays.
        If False then chunks are passed through unchanged.
        If None (default) then we use True if the ``__array_function__`` method
        is undefined.

        .. note::

            Dask does not preserve the memory layout of the original array when
            the array is created using Fortran rather than C ordering.

    fancy : bool, optional
        If ``x`` doesn't support fancy indexing (e.g. indexing with lists or
        arrays) then set to False. Default is True.
    meta : Array-like, optional
        The metadata for the resulting dask array.  This is the kind of array
        that will result from slicing the input array.
        Defaults to the input array.
    inline_array : bool, default False
        How to include the array in the task graph. By default
        (``inline_array=False``) the array is included in a task by itself,
        and each chunk refers to that task by its key.

        .. code-block:: python

           >>> x = h5py.File("data.h5")["/x"]  # doctest: +SKIP
           >>> a = da.from_array(x, chunks=500)  # doctest: +SKIP
           >>> dict(a.dask)  # doctest: +SKIP
           {
              'array-original-<name>': <HDF5 dataset ...>,
              ('array-<name>', 0): (getitem, "array-original-<name>", ...),
              ('array-<name>', 1): (getitem, "array-original-<name>", ...)
           }

        With ``inline_array=True``, Dask will instead inline the array directly
        in the values of the task graph.

        .. code-block:: python

           >>> a = da.from_array(x, chunks=500, inline_array=True)  # doctest: +SKIP
           >>> dict(a.dask)  # doctest: +SKIP
           {
              ('array-<name>', 0): (getitem, <HDF5 dataset ...>, ...),
              ('array-<name>', 1): (getitem, <HDF5 dataset ...>, ...)
           }

        Note that there's no key in the task graph with just the array `x`
        anymore. Instead it's placed directly in the values.

        The right choice for ``inline_array`` depends on several factors,
        including the size of ``x``, how expensive it is to create, which
        scheduler you're using, and the pattern of downstream computations.
        As a heuristic, ``inline_array=True`` may be the right choice when
        the array ``x`` is cheap to serialize and deserialize (since it's
        included in the graph many times) and if you're experiencing ordering
        issues (see :ref:`order` for more).

        This has no effect when ``x`` is a NumPy array.

    Examples
    --------

    >>> x = h5py.File('...')['/data/path']  # doctest: +SKIP
    >>> a = da.from_array(x, chunks=(1000, 1000))  # doctest: +SKIP

    If your underlying datastore does not support concurrent reads then include
    the ``lock=True`` keyword argument or ``lock=mylock`` if you want multiple
    arrays to coordinate around the same lock.

    >>> a = da.from_array(x, chunks=(1000, 1000), lock=True)  # doctest: +SKIP

    If your underlying datastore has a ``.chunks`` attribute (as h5py and zarr
    datasets do) then a multiple of that chunk shape will be used if you
    do not provide a chunk shape.

    >>> a = da.from_array(x, chunks='auto')  # doctest: +SKIP
    >>> a = da.from_array(x, chunks='100 MiB')  # doctest: +SKIP
    >>> a = da.from_array(x)  # doctest: +SKIP

    If providing a name, ensure that it is unique

    >>> import dask.base
    >>> token = dask.base.tokenize(x)  # doctest: +SKIP
    >>> a = da.from_array('myarray-' + token)  # doctest: +SKIP

    NumPy ndarrays are eagerly sliced and then embedded in the graph.

    >>> import dask.array
    >>> a = dask.array.from_array(np.array([[1, 2], [3, 4]]), chunks=(1,1))
    >>> a.dask[a.name, 0, 0][0]
    array([1])

    Chunks with exactly-specified, different sizes can be created.

    >>> import numpy as np
    >>> import dask.array as da
    >>> rng = np.random.default_rng()
    >>> x = rng.random((100, 6))
    >>> a = da.from_array(x, chunks=((67, 33), (6,)))
    """
    if isinstance(x, Array):
        raise ValueError(
            "Array is already a dask array. Use 'asarray' or 'rechunk' instead."
        )

    if xr is not None and isinstance(x, xr.DataArray) and x.chunks is not None:
        if isinstance(x.data, Array):
            return x.data

    elif is_dask_collection(x):
        warnings.warn(
            "Passing an object to dask.array.from_array which is already a "
            "Dask collection. This can lead to unexpected behavior."
        )

    if isinstance(x, (list, tuple, memoryview) + np.ScalarType):
        x = np.array(x)

    if is_arraylike(x) and hasattr(x, "copy"):
        x = x.copy()

    if asarray is None:
        asarray = not hasattr(x, "__array_function__")

    previous_chunks = getattr(x, "chunks", None)

    chunks = normalize_chunks(
        chunks, x.shape, dtype=x.dtype, previous_chunks=previous_chunks
    )

    if name in (None, True):
        token = tokenize(x, chunks, lock, asarray, fancy, getitem, inline_array)
        name = name or "array-" + token
    elif name is False:
        name = "array-" + str(uuid.uuid1())

    if lock is True:
        lock = SerializableLock()

    is_ndarray = type(x) in (np.ndarray, np.ma.core.MaskedArray)
    is_single_block = all(len(c) == 1 for c in chunks)
    # Always use the getter for h5py etc. Not using isinstance(x, np.ndarray)
    # because np.matrix is a subclass of np.ndarray.
    if is_ndarray and not is_single_block and not lock:
        # eagerly slice numpy arrays to prevent memory blowup
        # GH5367, GH5601
        slices = slices_from_chunks(chunks)
        keys = product([name], *(range(len(bds)) for bds in chunks))
        values = [x[slc] for slc in slices]
        dsk = dict(zip(keys, values))

    elif is_ndarray and is_single_block:
        # No slicing needed
        dsk = {(name,) + (0,) * x.ndim: x}
    else:
        if getitem is None:
            if fancy:
                getitem = getter
            else:
                getitem = getter_nofancy

        dsk = graph_from_arraylike(
            x,
            chunks,
            x.shape,
            name,
            getitem=getitem,
            lock=lock,
            asarray=asarray,
            dtype=x.dtype,
            inline_array=inline_array,
        )

    # Workaround for TileDB, its indexing is 1-based,
    # and doesn't seems to support 0-length slicing
    if x.__class__.__module__.split(".")[0] == "tiledb" and hasattr(x, "_ctx_"):
        return Array(dsk, name, chunks, dtype=x.dtype)

    if meta is None:
        meta = x

    return Array(dsk, name, chunks, meta=meta, dtype=getattr(x, "dtype", None))


@lru_cache
def _zarr_v3() -> bool:
    try:
        import zarr
    except ImportError:
        return False
    else:
        return Version(zarr.__version__).major >= 3


def from_zarr(
    url,
    component=None,
    storage_options=None,
    chunks=None,
    name=None,
    inline_array=False,
    **kwargs,
):
    """Load array from the zarr storage format

    See https://zarr.readthedocs.io for details about the format.

    Parameters
    ----------
    url: Zarr Array or str or MutableMapping
        Location of the data. A URL can include a protocol specifier like s3://
        for remote data. Can also be any MutableMapping instance, which should
        be serializable if used in multiple processes.
    component: str or None
        If the location is a zarr group rather than an array, this is the
        subcomponent that should be loaded, something like ``'foo/bar'``.
    storage_options: dict
        Any additional parameters for the storage backend (ignored for local
        paths)
    chunks: tuple of ints or tuples of ints
        Passed to :func:`dask.array.from_array`, allows setting the chunks on
        initialisation, if the chunking scheme in the on-disc dataset is not
        optimal for the calculations to follow.
    name : str, optional
         An optional keyname for the array.  Defaults to hashing the input
    kwargs:
        Passed to :class:`zarr.core.Array`.
    inline_array : bool, default False
        Whether to inline the zarr Array in the values of the task graph.
        See :meth:`dask.array.from_array` for an explanation.

    See Also
    --------
    from_array
    """
    import zarr

    storage_options = storage_options or {}
    if isinstance(url, zarr.Array):
        z = url
    elif isinstance(url, (str, os.PathLike)):
        if isinstance(url, os.PathLike):
            url = os.fspath(url)
        if storage_options:
            if _zarr_v3():
                store = zarr.storage.FsspecStore.from_url(
                    url, storage_options=storage_options
                )
            else:
                store = zarr.storage.FSStore(url, **storage_options)
        else:
            store = url
        z = zarr.open_array(store=store, path=component, **kwargs)
    else:
        z = zarr.open_array(store=url, path=component, **kwargs)
    chunks = chunks if chunks is not None else z.chunks
    if name is None:
        name = "from-zarr-" + tokenize(z, component, storage_options, chunks, **kwargs)
    return from_array(z, chunks, name=name, inline_array=inline_array)


def to_zarr(
    arr,
    url,
    component=None,
    storage_options=None,
    overwrite=False,
    region=None,
    compute=True,
    return_stored=False,
    **kwargs,
):
    """Save array to the zarr storage format

    See https://zarr.readthedocs.io for details about the format.

    Parameters
    ----------
    arr: dask.array
        Data to store
    url: Zarr Array or str or MutableMapping
        Location of the data. A URL can include a protocol specifier like s3://
        for remote data. Can also be any MutableMapping instance, which should
        be serializable if used in multiple processes.
    component: str or None
        If the location is a zarr group rather than an array, this is the
        subcomponent that should be created/over-written.
    storage_options: dict
        Any additional parameters for the storage backend (ignored for local
        paths)
    overwrite: bool
        If given array already exists, overwrite=False will cause an error,
        where overwrite=True will replace the existing data.
    region: tuple of slices or None
        The region of data that should be written if ``url`` is a zarr.Array.
        Not to be used with other types of ``url``.
    compute: bool
        See :func:`~dask.array.store` for more details.
    return_stored: bool
        See :func:`~dask.array.store` for more details.
    **kwargs:
        Passed to the :func:`zarr.creation.create` function, e.g., compression options.

    Raises
    ------
    ValueError
        If ``arr`` has unknown chunk sizes, which is not supported by Zarr.
        If ``region`` is specified and ``url`` is not a zarr.Array

    See Also
    --------
    dask.array.store
    dask.array.Array.compute_chunk_sizes

    """
    import zarr

    if np.isnan(arr.shape).any():
        raise ValueError(
            "Saving a dask array with unknown chunk sizes is not "
            "currently supported by Zarr.%s" % unknown_chunk_message
        )

    if _zarr_v3():
        zarr_mem_store_types = (zarr.storage.MemoryStore,)
    else:
        zarr_mem_store_types = (dict, zarr.storage.MemoryStore, zarr.storage.KVStore)

    if isinstance(url, zarr.Array):
        z = url
        if isinstance(z.store, zarr_mem_store_types):
            try:
                from distributed import default_client

                default_client()
            except (ImportError, ValueError):
                pass
            else:
                raise RuntimeError(
                    "Cannot store into in memory Zarr Array using "
                    "the distributed scheduler."
                )

        if region is None:
            arr = arr.rechunk(z.chunks)
            regions = None
        else:
            from dask.array.slicing import new_blockdim, normalize_index

            old_chunks = normalize_chunks(z.chunks, z.shape)
            index = normalize_index(region, z.shape)
            chunks = tuple(
                tuple(new_blockdim(s, c, r))
                for s, c, r in zip(z.shape, old_chunks, index)
            )
            arr = arr.rechunk(chunks)
            regions = [region]
        return arr.store(
            z, lock=False, regions=regions, compute=compute, return_stored=return_stored
        )
    else:
        if not _check_regular_chunks(arr.chunks):
            # We almost certainly get here because auto chunking has been used
            # on irregular chunks. The max will then be smaller than auto, so using
            # max is a safe choice
            arr = arr.rechunk(tuple(map(max, arr.chunks)))

    if region is not None:
        raise ValueError("Cannot use `region` keyword when url is not a `zarr.Array`.")

    if not _check_regular_chunks(arr.chunks):
        raise ValueError(
            "Attempt to save array to zarr with irregular "
            "chunking, please call `arr.rechunk(...)` first."
        )

    storage_options = storage_options or {}

    if storage_options:
        if _zarr_v3():
            read_only = (
                kwargs["read_only"]
                if "read_only" in kwargs
                else kwargs.pop("mode", "a") == "r"
            )
            store = zarr.storage.FsspecStore.from_url(
                url, read_only=read_only, storage_options=storage_options
            )
        else:
            store = zarr.storage.FSStore(url, **storage_options)
    else:
        store = url

    chunks = [c[0] for c in arr.chunks]

    z = zarr.create(
        shape=arr.shape,
        chunks=chunks,
        dtype=arr.dtype,
        store=store,
        path=component,
        overwrite=overwrite,
        **kwargs,
    )
    return arr.store(z, lock=False, compute=compute, return_stored=return_stored)


def _check_regular_chunks(chunkset):
    """Check if the chunks are regular

    "Regular" in this context means that along every axis, the chunks all
    have the same size, except the last one, which may be smaller

    Parameters
    ----------
    chunkset: tuple of tuples of ints
        From the ``.chunks`` attribute of an ``Array``

    Returns
    -------
    True if chunkset passes, else False

    Examples
    --------
    >>> import dask.array as da
    >>> arr = da.zeros(10, chunks=(5, ))
    >>> _check_regular_chunks(arr.chunks)
    True

    >>> arr = da.zeros(10, chunks=((3, 3, 3, 1), ))
    >>> _check_regular_chunks(arr.chunks)
    True

    >>> arr = da.zeros(10, chunks=((3, 1, 3, 3), ))
    >>> _check_regular_chunks(arr.chunks)
    False
    """
    for chunks in chunkset:
        if len(chunks) == 1:
            continue
        if len(set(chunks[:-1])) > 1:
            return False
        if chunks[-1] > chunks[0]:
            return False
    return True


def from_delayed(value, shape, dtype=None, meta=None, name=None):
    """Create a dask array from a dask delayed value

    This routine is useful for constructing dask arrays in an ad-hoc fashion
    using dask delayed, particularly when combined with stack and concatenate.

    The dask array will consist of a single chunk.

    Examples
    --------
    >>> import dask
    >>> import dask.array as da
    >>> import numpy as np
    >>> value = dask.delayed(np.ones)(5)
    >>> array = da.from_delayed(value, (5,), dtype=float)
    >>> array
    dask.array<from-value, shape=(5,), dtype=float64, chunksize=(5,), chunktype=numpy.ndarray>
    >>> array.compute()
    array([1., 1., 1., 1., 1.])
    """
    from dask.delayed import Delayed, delayed

    is_future = False
    name = name or "from-value-" + tokenize(value, shape, dtype, meta)
    if isinstance(value, TaskRef):
        is_future = True
    elif not isinstance(value, Delayed) and hasattr(value, "key"):
        value = delayed(value)
    task = Alias(
        key=(name,) + (0,) * len(shape),
        target=value.key,
    )

    dsk = {task.key: task}

    if is_future:
        dsk[value.key] = value
        dependencies = []
    else:
        dependencies = [value]
    chunks = tuple((d,) for d in shape)
    # TODO: value._key may not be the name of the layer in value.dask
    # This should be fixed after we build full expression graphs
    graph = HighLevelGraph.from_collections(name, dsk, dependencies=dependencies)
    return Array(graph, name, chunks, dtype=dtype, meta=meta)


def from_func(func, shape, dtype=None, name=None, args=(), kwargs=None):
    """Create dask array in a single block by calling a function

    Calling the provided function with func(*args, **kwargs) should return a
    NumPy array of the indicated shape and dtype.

    Examples
    --------

    >>> a = from_func(np.arange, (3,), dtype='i8', args=(3,))
    >>> a.compute()
    array([0, 1, 2])

    This works particularly well when coupled with dask.array functions like
    concatenate and stack:

    >>> arrays = [from_func(np.array, (), dtype='i8', args=(n,)) for n in range(5)]
    >>> stack(arrays).compute()
    array([0, 1, 2, 3, 4])
    """
    if kwargs is None:
        kwargs = {}

    name = name or "from_func-" + tokenize(func, shape, dtype, args, kwargs)
    if args or kwargs:
        func = partial(func, *args, **kwargs)
    dsk = {(name,) + (0,) * len(shape): (func,)}
    chunks = tuple((i,) for i in shape)
    return Array(dsk, name, chunks, dtype)


def common_blockdim(blockdims):
    """Find the common block dimensions from the list of block dimensions

    Currently only implements the simplest possible heuristic: the common
    block-dimension is the only one that does not span fully span a dimension.
    This is a conservative choice that allows us to avoid potentially very
    expensive rechunking.

    Assumes that each element of the input block dimensions has all the same
    sum (i.e., that they correspond to dimensions of the same size).

    Examples
    --------
    >>> common_blockdim([(3,), (2, 1)])
    (2, 1)
    >>> common_blockdim([(1, 2), (2, 1)])
    (1, 1, 1)
    >>> common_blockdim([(2, 2), (3, 1)])  # doctest: +SKIP
    Traceback (most recent call last):
        ...
    ValueError: Chunks do not align
    """
    if not any(blockdims):
        return ()
    non_trivial_dims = {d for d in blockdims if len(d) > 1}
    if len(non_trivial_dims) == 1:
        return first(non_trivial_dims)
    if len(non_trivial_dims) == 0:
        return max(blockdims, key=first)

    if np.isnan(sum(map(sum, blockdims))):
        raise ValueError(
            "Arrays' chunk sizes (%s) are unknown.\n\n"
            "A possible solution:\n"
            "  x.compute_chunk_sizes()" % blockdims
        )

    if len(set(map(sum, non_trivial_dims))) > 1:
        raise ValueError("Chunks do not add up to same value", blockdims)

    # We have multiple non-trivial chunks on this axis
    # e.g. (5, 2) and (4, 3)

    # We create a single chunk tuple with the same total length
    # that evenly divides both, e.g. (4, 1, 2)

    # To accomplish this we walk down all chunk tuples together, finding the
    # smallest element, adding it to the output, and subtracting it from all
    # other elements and remove the element itself.  We stop once we have
    # burned through all of the chunk tuples.
    # For efficiency's sake we reverse the lists so that we can pop off the end
    rchunks = [list(ntd)[::-1] for ntd in non_trivial_dims]
    total = sum(first(non_trivial_dims))
    i = 0

    out = []
    while i < total:
        m = min(c[-1] for c in rchunks)
        out.append(m)
        for c in rchunks:
            c[-1] -= m
            if c[-1] == 0:
                c.pop()
        i += m

    return tuple(out)


def unify_chunks(*args, **kwargs):
    """
    Unify chunks across a sequence of arrays

    This utility function is used within other common operations like
    :func:`dask.array.core.map_blocks` and :func:`dask.array.core.blockwise`.
    It is not commonly used by end-users directly.

    Parameters
    ----------
    *args: sequence of Array, index pairs
        Sequence like (x, 'ij', y, 'jk', z, 'i')

    Examples
    --------
    >>> import dask.array as da
    >>> x = da.ones(10, chunks=((5, 2, 3),))
    >>> y = da.ones(10, chunks=((2, 3, 5),))
    >>> chunkss, arrays = unify_chunks(x, 'i', y, 'i')
    >>> chunkss
    {'i': (2, 3, 2, 3)}

    >>> x = da.ones((100, 10), chunks=(20, 5))
    >>> y = da.ones((10, 100), chunks=(4, 50))
    >>> chunkss, arrays = unify_chunks(x, 'ij', y, 'jk', 'constant', None)
    >>> chunkss  # doctest: +SKIP
    {'k': (50, 50), 'i': (20, 20, 20, 20, 20), 'j': (4, 1, 3, 2)}

    >>> unify_chunks(0, None)
    ({}, [0])

    Returns
    -------
    chunkss : dict
        Map like {index: chunks}.
    arrays : list
        List of rechunked arrays.

    See Also
    --------
    common_blockdim
    """
    if not args:
        return {}, []

    arginds = [
        (asanyarray(a) if ind is not None else a, ind) for a, ind in partition(2, args)
    ]  # [x, ij, y, jk]
    warn = kwargs.get("warn", True)

    arrays, inds = zip(*arginds)
    if all(ind is None for ind in inds):
        return {}, list(arrays)
    if all(ind == inds[0] for ind in inds) and all(
        a.chunks == arrays[0].chunks for a in arrays
    ):
        return dict(zip(inds[0], arrays[0].chunks)), arrays

    nameinds = []
    blockdim_dict = dict()
    max_parts = 0
    for a, ind in arginds:
        if ind is not None:
            nameinds.append((a.name, ind))
            blockdim_dict[a.name] = a.chunks
            max_parts = max(max_parts, a.npartitions)
        else:
            nameinds.append((a, ind))

    chunkss = broadcast_dimensions(nameinds, blockdim_dict, consolidate=common_blockdim)
    nparts = math.prod(map(len, chunkss.values()))

    if warn and nparts and nparts >= max_parts * 10:
        warnings.warn(
            "Increasing number of chunks by factor of %d" % (nparts / max_parts),
            PerformanceWarning,
            stacklevel=3,
        )

    arrays = []
    for a, i in arginds:
        if i is None:
            arrays.append(a)
        else:
            chunks = tuple(
                (
                    chunkss[j]
                    if a.shape[n] > 1
                    else a.shape[n] if not np.isnan(sum(chunkss[j])) else None
                )
                for n, j in enumerate(i)
            )
            if chunks != a.chunks and all(a.chunks):
                arrays.append(a.rechunk(chunks))
            else:
                arrays.append(a)
    return chunkss, arrays


def unpack_singleton(x):
    """

    >>> unpack_singleton([[[[1]]]])
    1
    >>> unpack_singleton(np.array(np.datetime64('2000-01-01')))
    array('2000-01-01', dtype='datetime64[D]')
    """
    while isinstance(x, (list, tuple)):
        try:
            x = x[0]
        except (IndexError, TypeError, KeyError):
            break
    return x


def block(arrays, allow_unknown_chunksizes=False):
    """
    Assemble an nd-array from nested lists of blocks.

    Blocks in the innermost lists are concatenated along the last
    dimension (-1), then these are concatenated along the second-last
    dimension (-2), and so on until the outermost list is reached

    Blocks can be of any dimension, but will not be broadcasted using the normal
    rules. Instead, leading axes of size 1 are inserted, to make ``block.ndim``
    the same for all blocks. This is primarily useful for working with scalars,
    and means that code like ``block([v, 1])`` is valid, where
    ``v.ndim == 1``.

    When the nested list is two levels deep, this allows block matrices to be
    constructed from their components.

    Parameters
    ----------
    arrays : nested list of array_like or scalars (but not tuples)
        If passed a single ndarray or scalar (a nested list of depth 0), this
        is returned unmodified (and not copied).

        Elements shapes must match along the appropriate axes (without
        broadcasting), but leading 1s will be prepended to the shape as
        necessary to make the dimensions match.

    allow_unknown_chunksizes: bool
        Allow unknown chunksizes, such as come from converting from dask
        dataframes.  Dask.array is unable to verify that chunks line up.  If
        data comes from differently aligned sources then this can cause
        unexpected results.

    Returns
    -------
    block_array : ndarray
        The array assembled from the given blocks.

        The dimensionality of the output is equal to the greatest of:
        * the dimensionality of all the inputs
        * the depth to which the input list is nested

    Raises
    ------
    ValueError
        * If list depths are mismatched - for instance, ``[[a, b], c]`` is
          illegal, and should be spelt ``[[a, b], [c]]``
        * If lists are empty - for instance, ``[[a, b], []]``

    See Also
    --------
    concatenate : Join a sequence of arrays together.
    stack : Stack arrays in sequence along a new dimension.
    hstack : Stack arrays in sequence horizontally (column wise).
    vstack : Stack arrays in sequence vertically (row wise).
    dstack : Stack arrays in sequence depth wise (along third dimension).
    vsplit : Split array into a list of multiple sub-arrays vertically.

    Notes
    -----

    When called with only scalars, ``block`` is equivalent to an ndarray
    call. So ``block([[1, 2], [3, 4]])`` is equivalent to
    ``array([[1, 2], [3, 4]])``.

    This function does not enforce that the blocks lie on a fixed grid.
    ``block([[a, b], [c, d]])`` is not restricted to arrays of the form::

        AAAbb
        AAAbb
        cccDD

    But is also allowed to produce, for some ``a, b, c, d``::

        AAAbb
        AAAbb
        cDDDD

    Since concatenation happens along the last axis first, `block` is _not_
    capable of producing the following directly::

        AAAbb
        cccbb
        cccDD

    Matlab's "square bracket stacking", ``[A, B, ...; p, q, ...]``, is
    equivalent to ``block([[A, B, ...], [p, q, ...]])``.
    """

    # This was copied almost verbatim from numpy.core.shape_base.block
    # See numpy license at https://github.com/numpy/numpy/blob/master/LICENSE.txt
    # or NUMPY_LICENSE.txt within this directory

    def atleast_nd(x, ndim):
        x = asanyarray(x)
        diff = max(ndim - x.ndim, 0)
        if diff == 0:
            return x
        else:
            return x[(None,) * diff + (Ellipsis,)]

    def format_index(index):
        return "arrays" + "".join(f"[{i}]" for i in index)

    rec = _Recurser(recurse_if=lambda x: type(x) is list)

    # ensure that the lists are all matched in depth
    list_ndim = None
    any_empty = False
    for index, value, entering in rec.walk(arrays):
        if type(value) is tuple:
            # not strictly necessary, but saves us from:
            #  - more than one way to do things - no point treating tuples like
            #    lists
            #  - horribly confusing behaviour that results when tuples are
            #    treated like ndarray
            raise TypeError(
                "{} is a tuple. "
                "Only lists can be used to arrange blocks, and np.block does "
                "not allow implicit conversion from tuple to ndarray.".format(
                    format_index(index)
                )
            )
        if not entering:
            curr_depth = len(index)
        elif len(value) == 0:
            curr_depth = len(index) + 1
            any_empty = True
        else:
            continue

        if list_ndim is not None and list_ndim != curr_depth:
            raise ValueError(
                "List depths are mismatched. First element was at depth {}, "
                "but there is an element at depth {} ({})".format(
                    list_ndim, curr_depth, format_index(index)
                )
            )
        list_ndim = curr_depth

    # do this here so we catch depth mismatches first
    if any_empty:
        raise ValueError("Lists cannot be empty")

    # convert all the arrays to ndarrays
    arrays = rec.map_reduce(arrays, f_map=asanyarray, f_reduce=list)

    # determine the maximum dimension of the elements
    elem_ndim = rec.map_reduce(arrays, f_map=lambda xi: xi.ndim, f_reduce=max)
    ndim = max(list_ndim, elem_ndim)

    # first axis to concatenate along
    first_axis = ndim - list_ndim

    # Make all the elements the same dimension
    arrays = rec.map_reduce(
        arrays, f_map=lambda xi: atleast_nd(xi, ndim), f_reduce=list
    )

    # concatenate innermost lists on the right, outermost on the left
    return rec.map_reduce(
        arrays,
        f_reduce=lambda xs, axis: concatenate(
            list(xs), axis=axis, allow_unknown_chunksizes=allow_unknown_chunksizes
        ),
        f_kwargs=lambda axis: dict(axis=(axis + 1)),
        axis=first_axis,
    )


def concatenate(seq, axis=0, allow_unknown_chunksizes=False):
    """
    Concatenate arrays along an existing axis

    Given a sequence of dask Arrays form a new dask Array by stacking them
    along an existing dimension (axis=0 by default)

    Parameters
    ----------
    seq: list of dask.arrays
    axis: int
        Dimension along which to align all of the arrays. If axis is None,
        arrays are flattened before use.
    allow_unknown_chunksizes: bool
        Allow unknown chunksizes, such as come from converting from dask
        dataframes.  Dask.array is unable to verify that chunks line up.  If
        data comes from differently aligned sources then this can cause
        unexpected results.

    Examples
    --------

    Create slices

    >>> import dask.array as da
    >>> import numpy as np

    >>> data = [da.from_array(np.ones((4, 4)), chunks=(2, 2))
    ...          for i in range(3)]

    >>> x = da.concatenate(data, axis=0)
    >>> x.shape
    (12, 4)

    >>> da.concatenate(data, axis=1).shape
    (4, 12)

    Result is a new dask Array

    See Also
    --------
    stack
    """
    from dask.array import wrap

    seq = [asarray(a, allow_unknown_chunksizes=allow_unknown_chunksizes) for a in seq]

    if not seq:
        raise ValueError("Need array(s) to concatenate")

    if axis is None:
        seq = [a.flatten() for a in seq]
        axis = 0

    seq_metas = [meta_from_array(s) for s in seq]
    _concatenate = concatenate_lookup.dispatch(
        type(max(seq_metas, key=lambda x: getattr(x, "__array_priority__", 0)))
    )
    meta = _concatenate(seq_metas, axis=axis)

    # Promote types to match meta
    seq = [a.astype(meta.dtype) for a in seq]

    # Find output array shape
    ndim = len(seq[0].shape)
    shape = tuple(
        sum(a.shape[i] for a in seq) if i == axis else seq[0].shape[i]
        for i in range(ndim)
    )

    # Drop empty arrays
    seq2 = [a for a in seq if a.size]
    if not seq2:
        seq2 = seq

    if axis < 0:
        axis = ndim + axis
    if axis >= ndim:
        msg = (
            "Axis must be less than than number of dimensions"
            "\nData has %d dimensions, but got axis=%d"
        )
        raise ValueError(msg % (ndim, axis))

    n = len(seq2)
    if n == 0:
        try:
            return wrap.empty_like(meta, shape=shape, chunks=shape, dtype=meta.dtype)
        except TypeError:
            return wrap.empty(shape, chunks=shape, dtype=meta.dtype)
    elif n == 1:
        return seq2[0]

    if not allow_unknown_chunksizes and not all(
        i == axis or all(x.shape[i] == seq2[0].shape[i] for x in seq2)
        for i in range(ndim)
    ):
        if any(map(np.isnan, seq2[0].shape)):
            raise ValueError(
                "Tried to concatenate arrays with unknown"
                " shape %s.\n\nTwo solutions:\n"
                "  1. Force concatenation pass"
                " allow_unknown_chunksizes=True.\n"
                "  2. Compute shapes with "
                "[x.compute_chunk_sizes() for x in seq]" % str(seq2[0].shape)
            )
        raise ValueError("Shapes do not align: %s", [x.shape for x in seq2])

    inds = [list(range(ndim)) for i in range(n)]
    for i, ind in enumerate(inds):
        ind[axis] = -(i + 1)

    uc_args = list(concat(zip(seq2, inds)))
    _, seq2 = unify_chunks(*uc_args, warn=False)

    bds = [a.chunks for a in seq2]

    chunks = (
        seq2[0].chunks[:axis]
        + (sum((bd[axis] for bd in bds), ()),)
        + seq2[0].chunks[axis + 1 :]
    )

    cum_dims = [0] + list(accumulate(add, [len(a.chunks[axis]) for a in seq2]))

    names = [a.name for a in seq2]

    name = "concatenate-" + tokenize(names, axis)
    keys = list(product([name], *[range(len(bd)) for bd in chunks]))

    values = [
        (names[bisect(cum_dims, key[axis + 1]) - 1],)
        + key[1 : axis + 1]
        + (key[axis + 1] - cum_dims[bisect(cum_dims, key[axis + 1]) - 1],)
        + key[axis + 2 :]
        for key in keys
    ]

    dsk = dict(zip(keys, values))
    graph = HighLevelGraph.from_collections(name, dsk, dependencies=seq2)

    return Array(graph, name, chunks, meta=meta)


def load_store_chunk(
    x: Any,
    out: Any,
    index: slice | None,
    region: slice | None,
    lock: Any,
    return_stored: bool,
    load_stored: bool,
) -> Any:
    """
    A function inserted in a Dask graph for storing a chunk.

    Parameters
    ----------
    x: array-like
        An array (potentially a NumPy one)
    out: array-like
        Where to store results.
    index: slice-like
        Where to store result from ``x`` in ``out``.
    lock: Lock-like or False
        Lock to use before writing to ``out``.
    return_stored: bool
        Whether to return ``out``.
    load_stored: bool
        Whether to return the array stored in ``out``.
        Ignored if ``return_stored`` is not ``True``.

    Returns
    -------

    If return_stored=True and load_stored=False
        out
    If return_stored=True and load_stored=True
        out[index]
    If return_stored=False and compute=False
        None

    Examples
    --------

    >>> a = np.ones((5, 6))
    >>> b = np.empty(a.shape)
    >>> load_store_chunk(a, b, (slice(None), slice(None)), None, False, False, False)
    """
    if region:
        # Equivalent to `out[region][index]`
        if index:
            index = fuse_slice(region, index)
        else:
            index = region
    if lock:
        lock.acquire()
    try:
        if x is not None and x.size != 0:
            if is_arraylike(x):
                out[index] = x
            else:
                out[index] = np.asanyarray(x)

        if return_stored and load_stored:
            return out[index]
        elif return_stored and not load_stored:
            return out
        else:
            return None
    finally:
        if lock:
            lock.release()


A = TypeVar("A", bound=ArrayLike)


def load_chunk(out: A, index: slice, lock: Any, region: slice | None) -> A:
    return load_store_chunk(
        None,
        out=out,
        region=region,
        index=index,
        lock=lock,
        return_stored=True,
        load_stored=True,
    )


def _as_dtype(a, dtype):
    if dtype is None:
        return a
    else:
        return a.astype(dtype)


def asarray(
    a, allow_unknown_chunksizes=False, dtype=None, order=None, *, like=None, **kwargs
):
    """Convert the input to a dask array.

    Parameters
    ----------
    a : array-like
        Input data, in any form that can be converted to a dask array. This
        includes lists, lists of tuples, tuples, tuples of tuples, tuples of
        lists and ndarrays.
    allow_unknown_chunksizes: bool
        Allow unknown chunksizes, such as come from converting from dask
        dataframes.  Dask.array is unable to verify that chunks line up.  If
        data comes from differently aligned sources then this can cause
        unexpected results.
    dtype : data-type, optional
        By default, the data-type is inferred from the input data.
    order : {‘C’, ‘F’, ‘A’, ‘K’}, optional
        Memory layout. ‘A’ and ‘K’ depend on the order of input array a.
        ‘C’ row-major (C-style), ‘F’ column-major (Fortran-style) memory
        representation. ‘A’ (any) means ‘F’ if a is Fortran contiguous, ‘C’
        otherwise ‘K’ (keep) preserve input order. Defaults to ‘C’.
    like: array-like
        Reference object to allow the creation of Dask arrays with chunks
        that are not NumPy arrays. If an array-like passed in as ``like``
        supports the ``__array_function__`` protocol, the chunk type of the
        resulting array will be defined by it. In this case, it ensures the
        creation of a Dask array compatible with that passed in via this
        argument. If ``like`` is a Dask array, the chunk type of the
        resulting array will be defined by the chunk type of ``like``.
        Requires NumPy 1.20.0 or higher.

    Returns
    -------
    out : dask array
        Dask array interpretation of a.

    Examples
    --------
    >>> import dask.array as da
    >>> import numpy as np
    >>> x = np.arange(3)
    >>> da.asarray(x)
    dask.array<array, shape=(3,), dtype=int64, chunksize=(3,), chunktype=numpy.ndarray>

    >>> y = [[1, 2, 3], [4, 5, 6]]
    >>> da.asarray(y)
    dask.array<array, shape=(2, 3), dtype=int64, chunksize=(2, 3), chunktype=numpy.ndarray>

    .. warning::
        `order` is ignored if `a` is an `Array`, has the attribute ``to_dask_array``,
        or is a list or tuple of `Array`'s.
    """
    if like is None:
        if isinstance(a, Array):
            return _as_dtype(a, dtype)
        elif hasattr(a, "to_dask_array"):
            return _as_dtype(a.to_dask_array(), dtype)
        elif type(a).__module__.split(".")[0] == "xarray" and hasattr(a, "data"):
            return _as_dtype(asarray(a.data, order=order), dtype)
        elif isinstance(a, (list, tuple)) and any(isinstance(i, Array) for i in a):
            return _as_dtype(
                stack(a, allow_unknown_chunksizes=allow_unknown_chunksizes), dtype
            )
        elif not isinstance(getattr(a, "shape", None), Iterable):
            a = np.asarray(a, dtype=dtype, order=order)
    else:
        like_meta = meta_from_array(like)
        if isinstance(a, Array):
            return a.map_blocks(
                # Pass the dtype parameter to np.asarray, not to map_blocks
                partial(asarray_safe, like=like_meta, dtype=dtype, order=order)
            )
        else:
            a = asarray_safe(a, like=like_meta, dtype=dtype, order=order)

    a = from_array(a, getitem=getter_inline, **kwargs)
    return _as_dtype(a, dtype)


def asanyarray(a, dtype=None, order=None, *, like=None, inline_array=False):
    """Convert the input to a dask array.

    Subclasses of ``np.ndarray`` will be passed through as chunks unchanged.

    Parameters
    ----------
    a : array-like
        Input data, in any form that can be converted to a dask array. This
        includes lists, lists of tuples, tuples, tuples of tuples, tuples of
        lists and ndarrays.
    dtype : data-type, optional
        By default, the data-type is inferred from the input data.
    order : {‘C’, ‘F’, ‘A’, ‘K’}, optional
        Memory layout. ‘A’ and ‘K’ depend on the order of input array a.
        ‘C’ row-major (C-style), ‘F’ column-major (Fortran-style) memory
        representation. ‘A’ (any) means ‘F’ if a is Fortran contiguous, ‘C’
        otherwise ‘K’ (keep) preserve input order. Defaults to ‘C’.
    like: array-like
        Reference object to allow the creation of Dask arrays with chunks
        that are not NumPy arrays. If an array-like passed in as ``like``
        supports the ``__array_function__`` protocol, the chunk type of the
        resulting array will be defined by it. In this case, it ensures the
        creation of a Dask array compatible with that passed in via this
        argument. If ``like`` is a Dask array, the chunk type of the
        resulting array will be defined by the chunk type of ``like``.
        Requires NumPy 1.20.0 or higher.
    inline_array:
        Whether to inline the array in the resulting dask graph. For more information,
        see the documentation for ``dask.array.from_array()``.

    Returns
    -------
    out : dask array
        Dask array interpretation of a.

    Examples
    --------
    >>> import dask.array as da
    >>> import numpy as np
    >>> x = np.arange(3)
    >>> da.asanyarray(x)
    dask.array<array, shape=(3,), dtype=int64, chunksize=(3,), chunktype=numpy.ndarray>

    >>> y = [[1, 2, 3], [4, 5, 6]]
    >>> da.asanyarray(y)
    dask.array<array, shape=(2, 3), dtype=int64, chunksize=(2, 3), chunktype=numpy.ndarray>

    .. warning::
        `order` is ignored if `a` is an `Array`, has the attribute ``to_dask_array``,
        or is a list or tuple of `Array`'s.
    """
    if like is None:
        if isinstance(a, Array):
            return _as_dtype(a, dtype)
        elif hasattr(a, "to_dask_array"):
            return _as_dtype(a.to_dask_array(), dtype)
        elif type(a).__module__.split(".")[0] == "xarray" and hasattr(a, "data"):
            return _as_dtype(asarray(a.data, order=order), dtype)
        elif isinstance(a, (list, tuple)) and any(isinstance(i, Array) for i in a):
            return _as_dtype(stack(a), dtype)
        elif not isinstance(getattr(a, "shape", None), Iterable):
            a = np.asanyarray(a, dtype=dtype, order=order)
    else:
        like_meta = meta_from_array(like)
        if isinstance(a, Array):
            return a.map_blocks(
                # Pass the dtype parameter to np.asanyarray, not to map_blocks
                partial(asanyarray_safe, like=like_meta, dtype=dtype, order=order)
            )
        else:
            a = asanyarray_safe(a, like=like_meta, dtype=dtype, order=order)

    a = from_array(
        a,
        chunks=a.shape,
        getitem=getter_inline,
        asarray=False,
        inline_array=inline_array,
    )
    return _as_dtype(a, dtype)


def is_scalar_for_elemwise(arg):
    """

    >>> is_scalar_for_elemwise(42)
    True
    >>> is_scalar_for_elemwise('foo')
    True
    >>> is_scalar_for_elemwise(True)
    True
    >>> is_scalar_for_elemwise(np.array(42))
    True
    >>> is_scalar_for_elemwise([1, 2, 3])
    True
    >>> is_scalar_for_elemwise(np.array([1, 2, 3]))
    False
    >>> is_scalar_for_elemwise(from_array(np.array(0), chunks=()))
    False
    >>> is_scalar_for_elemwise(np.dtype('i4'))
    True
    """
    # the second half of shape_condition is essentially just to ensure that
    # dask series / frame are treated as scalars in elemwise.
    maybe_shape = getattr(arg, "shape", None)
    shape_condition = not isinstance(maybe_shape, Iterable) or any(
        is_dask_collection(x) for x in maybe_shape
    )

    return (
        np.isscalar(arg)
        or shape_condition
        or isinstance(arg, np.dtype)
        or (isinstance(arg, np.ndarray) and arg.ndim == 0)
    )


def broadcast_shapes(*shapes):
    """
    Determines output shape from broadcasting arrays.

    Parameters
    ----------
    shapes : tuples
        The shapes of the arguments.

    Returns
    -------
    output_shape : tuple

    Raises
    ------
    ValueError
        If the input shapes cannot be successfully broadcast together.
    """
    if len(shapes) == 1:
        return shapes[0]
    out = []
    for sizes in zip_longest(*map(reversed, shapes), fillvalue=-1):
        if np.isnan(sizes).any():
            dim = np.nan
        else:
            dim = 0 if 0 in sizes else np.max(sizes).item()
        if any(i not in [-1, 0, 1, dim] and not np.isnan(i) for i in sizes):
            raise ValueError(
                "operands could not be broadcast together with "
                "shapes {}".format(" ".join(map(str, shapes)))
            )
        out.append(dim)
    return tuple(reversed(out))


def elemwise(op, *args, out=None, where=True, dtype=None, name=None, **kwargs):
    """Apply an elementwise ufunc-like function blockwise across arguments.

    Like numpy ufuncs, broadcasting rules are respected.

    Parameters
    ----------
    op : callable
        The function to apply. Should be numpy ufunc-like in the parameters
        that it accepts.
    *args : Any
        Arguments to pass to `op`. Non-dask array-like objects are first
        converted to dask arrays, then all arrays are broadcast together before
        applying the function blockwise across all arguments. Any scalar
        arguments are passed as-is following normal numpy ufunc behavior.
    out : dask array, optional
        If out is a dask.array then this overwrites the contents of that array
        with the result.
    where : array_like, optional
        An optional boolean mask marking locations where the ufunc should be
        applied. Can be a scalar, dask array, or any other array-like object.
        Mirrors the ``where`` argument to numpy ufuncs, see e.g. ``numpy.add``
        for more information.
    dtype : dtype, optional
        If provided, overrides the output array dtype.
    name : str, optional
        A unique key name to use when building the backing dask graph. If not
        provided, one will be automatically generated based on the input
        arguments.

    Examples
    --------
    >>> elemwise(add, x, y)  # doctest: +SKIP
    >>> elemwise(sin, x)  # doctest: +SKIP
    >>> elemwise(sin, x, out=dask_array)  # doctest: +SKIP

    See Also
    --------
    blockwise
    """
    if kwargs:
        raise TypeError(
            f"{op.__name__} does not take the following keyword arguments "
            f"{sorted(kwargs)}"
        )

    out = _elemwise_normalize_out(out)
    where = _elemwise_normalize_where(where)
    args = [np.asarray(a) if isinstance(a, (list, tuple)) else a for a in args]

    shapes = []
    for arg in args:
        shape = getattr(arg, "shape", ())
        if any(is_dask_collection(x) for x in shape):
            # Want to exclude Delayed shapes and dd.Scalar
            shape = ()
        shapes.append(shape)
    if isinstance(where, Array):
        shapes.append(where.shape)
    if isinstance(out, Array):
        shapes.append(out.shape)

    shapes = [s if isinstance(s, Iterable) else () for s in shapes]
    out_ndim = len(
        broadcast_shapes(*shapes)
    )  # Raises ValueError if dimensions mismatch
    expr_inds = tuple(range(out_ndim))[::-1]

    if dtype is not None:
        need_enforce_dtype = True
    else:
        # We follow NumPy's rules for dtype promotion, which special cases
        # scalars and 0d ndarrays (which it considers equivalent) by using
        # their values to compute the result dtype:
        # https://github.com/numpy/numpy/issues/6240
        # We don't inspect the values of 0d dask arrays, because these could
        # hold potentially very expensive calculations. Instead, we treat
        # them just like other arrays, and if necessary cast the result of op
        # to match.
        vals = [
            (
                np.empty((1,) * max(1, a.ndim), dtype=a.dtype)
                if not is_scalar_for_elemwise(a)
                else a
            )
            for a in args
        ]
        try:
            dtype = apply_infer_dtype(op, vals, {}, "elemwise", suggest_dtype=False)
        except Exception:
            return NotImplemented
        need_enforce_dtype = any(
            not is_scalar_for_elemwise(a) and a.ndim == 0 for a in args
        )

    if not name:
        name = f"{funcname(op)}-{tokenize(op, dtype, *args, where)}"

    blockwise_kwargs = dict(dtype=dtype, name=name, token=funcname(op).strip("_"))

    if where is not True:
        blockwise_kwargs["elemwise_where_function"] = op
        op = _elemwise_handle_where
        args.extend([where, out])

    if need_enforce_dtype:
        blockwise_kwargs["enforce_dtype"] = dtype
        blockwise_kwargs["enforce_dtype_function"] = op
        op = _enforce_dtype

    result = blockwise(
        op,
        expr_inds,
        *concat(
            (a, tuple(range(a.ndim)[::-1]) if not is_scalar_for_elemwise(a) else None)
            for a in args
        ),
        **blockwise_kwargs,
    )

    return handle_out(out, result)


def _elemwise_normalize_where(where):
    if where is True:
        return True
    elif where is False or where is None:
        return False
    return asarray(where)


def _elemwise_handle_where(*args, **kwargs):
    function = kwargs.pop("elemwise_where_function")
    *args, where, out = args
    if hasattr(out, "copy"):
        out = out.copy()
    return function(*args, where=where, out=out, **kwargs)


def _elemwise_normalize_out(out):
    if isinstance(out, tuple):
        if len(out) == 1:
            out = out[0]
        elif len(out) > 1:
            raise NotImplementedError("The out parameter is not fully supported")
        else:
            out = None
    if not (out is None or isinstance(out, Array)):
        raise NotImplementedError(
            f"The out parameter is not fully supported."
            f" Received type {type(out).__name__}, expected Dask Array"
        )
    return out


def handle_out(out, result):
    """Handle out parameters

    If out is a dask.array then this overwrites the contents of that array with
    the result
    """
    out = _elemwise_normalize_out(out)
    if isinstance(out, Array):
        if out.shape != result.shape:
            raise ValueError(
                "Mismatched shapes between result and out parameter. "
                "out=%s, result=%s" % (str(out.shape), str(result.shape))
            )
        out._chunks = result.chunks
        out.dask = result.dask
        out._meta = result._meta
        out._name = result.name
        return out
    else:
        return result


def _enforce_dtype(*args, **kwargs):
    """Calls a function and converts its result to the given dtype.

    The parameters have deliberately been given unwieldy names to avoid
    clashes with keyword arguments consumed by blockwise

    A dtype of `object` is treated as a special case and not enforced,
    because it is used as a dummy value in some places when the result will
    not be a block in an Array.

    Parameters
    ----------
    enforce_dtype : dtype
        Result dtype
    enforce_dtype_function : callable
        The wrapped function, which will be passed the remaining arguments
    """
    dtype = kwargs.pop("enforce_dtype")
    function = kwargs.pop("enforce_dtype_function")

    result = function(*args, **kwargs)
    if hasattr(result, "dtype") and dtype != result.dtype and dtype != object:
        if not np.can_cast(result, dtype, casting="same_kind"):
            raise ValueError(
                "Inferred dtype from function %r was %r "
                "but got %r, which can't be cast using "
                "casting='same_kind'"
                % (funcname(function), str(dtype), str(result.dtype))
            )
        if np.isscalar(result):
            # scalar astype method doesn't take the keyword arguments, so
            # have to convert via 0-dimensional array and back.
            result = result.astype(dtype)
        else:
            try:
                result = result.astype(dtype, copy=False)
            except TypeError:
                # Missing copy kwarg
                result = result.astype(dtype)
    return result


def broadcast_to(x, shape, chunks=None, meta=None):
    """Broadcast an array to a new shape.

    Parameters
    ----------
    x : array_like
        The array to broadcast.
    shape : tuple
        The shape of the desired array.
    chunks : tuple, optional
        If provided, then the result will use these chunks instead of the same
        chunks as the source array. Setting chunks explicitly as part of
        broadcast_to is more efficient than rechunking afterwards. Chunks are
        only allowed to differ from the original shape along dimensions that
        are new on the result or have size 1 the input array.
    meta : empty ndarray
        empty ndarray created with same NumPy backend, ndim and dtype as the
        Dask Array being created (overrides dtype)

    Returns
    -------
    broadcast : dask array

    See Also
    --------
    :func:`numpy.broadcast_to`
    """
    x = asarray(x)
    shape = tuple(shape)

    if meta is None:
        meta = meta_from_array(x)

    if x.shape == shape and (chunks is None or chunks == x.chunks):
        return x

    ndim_new = len(shape) - x.ndim
    if ndim_new < 0 or any(
        new != old for new, old in zip(shape[ndim_new:], x.shape) if old != 1
    ):
        raise ValueError(f"cannot broadcast shape {x.shape} to shape {shape}")

    if chunks is None:
        chunks = tuple((s,) for s in shape[:ndim_new]) + tuple(
            bd if old > 1 else (new,)
            for bd, old, new in zip(x.chunks, x.shape, shape[ndim_new:])
        )
    else:
        chunks = normalize_chunks(
            chunks, shape, dtype=x.dtype, previous_chunks=x.chunks
        )
        for old_bd, new_bd in zip(x.chunks, chunks[ndim_new:]):
            if old_bd != new_bd and old_bd != (1,):
                raise ValueError(
                    "cannot broadcast chunks %s to chunks %s: "
                    "new chunks must either be along a new "
                    "dimension or a dimension of size 1" % (x.chunks, chunks)
                )

    name = "broadcast_to-" + tokenize(x, shape, chunks)
    dsk = {}

    enumerated_chunks = product(*(enumerate(bds) for bds in chunks))
    for new_index, chunk_shape in (zip(*ec) for ec in enumerated_chunks):
        old_index = tuple(
            0 if bd == (1,) else i for bd, i in zip(x.chunks, new_index[ndim_new:])
        )
        old_key = (x.name,) + old_index
        new_key = (name,) + new_index
        dsk[new_key] = (np.broadcast_to, old_key, quote(chunk_shape))

    graph = HighLevelGraph.from_collections(name, dsk, dependencies=[x])
    return Array(graph, name, chunks, dtype=x.dtype, meta=meta)


@derived_from(np)
def broadcast_arrays(*args, subok=False):
    subok = bool(subok)

    to_array = asanyarray if subok else asarray
    args = tuple(to_array(e) for e in args)

    # Unify uneven chunking
    inds = [list(reversed(range(x.ndim))) for x in args]
    uc_args = concat(zip(args, inds))
    _, args = unify_chunks(*uc_args, warn=False)

    shape = broadcast_shapes(*(e.shape for e in args))
    chunks = broadcast_chunks(*(e.chunks for e in args))

    if NUMPY_GE_200:
        result = tuple(broadcast_to(e, shape=shape, chunks=chunks) for e in args)
    else:
        result = [broadcast_to(e, shape=shape, chunks=chunks) for e in args]

    return result


def offset_func(func, offset, *args):
    """Offsets inputs by offset

    >>> double = lambda x: x * 2
    >>> f = offset_func(double, (10,))
    >>> f(1)
    22
    >>> f(300)
    620
    """

    def _offset(*args):
        args2 = list(map(add, args, offset))
        return func(*args2)

    with contextlib.suppress(Exception):
        _offset.__name__ = "offset_" + func.__name__

    return _offset


def chunks_from_arrays(arrays):
    """Chunks tuple from nested list of arrays

    >>> x = np.array([1, 2])
    >>> chunks_from_arrays([x, x])
    ((2, 2),)

    >>> x = np.array([[1, 2]])
    >>> chunks_from_arrays([[x], [x]])
    ((1, 1), (2,))

    >>> x = np.array([[1, 2]])
    >>> chunks_from_arrays([[x, x]])
    ((1,), (2, 2))

    >>> chunks_from_arrays([1, 1])
    ((1, 1),)
    """
    if not arrays:
        return ()
    result = []
    dim = 0

    def shape(x):
        try:
            return x.shape if x.shape else (1,)
        except AttributeError:
            return (1,)

    while isinstance(arrays, (list, tuple)):
        result.append(tuple(shape(deepfirst(a))[dim] for a in arrays))
        arrays = arrays[0]
        dim += 1
    return tuple(result)


def deepfirst(seq):
    """First element in a nested list

    >>> deepfirst([[[1, 2], [3, 4]], [5, 6], [7, 8]])
    1
    """
    if not isinstance(seq, (list, tuple)):
        return seq
    else:
        return deepfirst(seq[0])


def shapelist(a):
    """Get the shape of nested list"""
    if type(a) is list:
        return tuple([len(a)] + list(shapelist(a[0])))
    else:
        return ()


def transposelist(arrays, axes, extradims=0):
    """Permute axes of nested list

    >>> transposelist([[1,1,1],[1,1,1]], [2,1])
    [[[1, 1], [1, 1], [1, 1]]]

    >>> transposelist([[1,1,1],[1,1,1]], [2,1], extradims=1)
    [[[[1], [1]], [[1], [1]], [[1], [1]]]]
    """
    if len(axes) != ndimlist(arrays):
        raise ValueError("Length of axes should equal depth of nested arrays")
    if extradims < 0:
        raise ValueError("`newdims` should be positive")
    if len(axes) > len(set(axes)):
        raise ValueError("`axes` should be unique")

    ndim = max(axes) + 1
    shape = shapelist(arrays)
    newshape = [
        shape[axes.index(i)] if i in axes else 1 for i in range(ndim + extradims)
    ]

    result = list(core.flatten(arrays))
    return reshapelist(newshape, result)


def stack(seq, axis=0, allow_unknown_chunksizes=False):
    """
    Stack arrays along a new axis

    Given a sequence of dask arrays, form a new dask array by stacking them
    along a new dimension (axis=0 by default)

    Parameters
    ----------
    seq: list of dask.arrays
    axis: int
        Dimension along which to align all of the arrays
    allow_unknown_chunksizes: bool
        Allow unknown chunksizes, such as come from converting from dask
        dataframes.  Dask.array is unable to verify that chunks line up.  If
        data comes from differently aligned sources then this can cause
        unexpected results.

    Examples
    --------

    Create slices

    >>> import dask.array as da
    >>> import numpy as np

    >>> data = [da.from_array(np.ones((4, 4)), chunks=(2, 2))
    ...         for i in range(3)]

    >>> x = da.stack(data, axis=0)
    >>> x.shape
    (3, 4, 4)

    >>> da.stack(data, axis=1).shape
    (4, 3, 4)

    >>> da.stack(data, axis=-1).shape
    (4, 4, 3)

    Result is a new dask Array

    See Also
    --------
    concatenate
    """
    from dask.array import wrap

    seq = [asarray(a, allow_unknown_chunksizes=allow_unknown_chunksizes) for a in seq]

    if not seq:
        raise ValueError("Need array(s) to stack")
    if not allow_unknown_chunksizes and not all(x.shape == seq[0].shape for x in seq):
        idx = first(i for i in enumerate(seq) if i[1].shape != seq[0].shape)
        raise ValueError(
            "Stacked arrays must have the same shape. The first array had shape "
            f"{seq[0].shape}, while array {idx[0] + 1} has shape {idx[1].shape}."
        )

    meta = np.stack([meta_from_array(a) for a in seq], axis=axis)
    seq = [x.astype(meta.dtype) for x in seq]

    ndim = meta.ndim - 1
    if axis < 0:
        axis = ndim + axis + 1
    shape = tuple(
        (
            len(seq)
            if i == axis
            else (seq[0].shape[i] if i < axis else seq[0].shape[i - 1])
        )
        for i in range(meta.ndim)
    )

    seq2 = [a for a in seq if a.size]
    if not seq2:
        seq2 = seq

    n = len(seq2)
    if n == 0:
        try:
            return wrap.empty_like(meta, shape=shape, chunks=shape, dtype=meta.dtype)
        except TypeError:
            return wrap.empty(shape, chunks=shape, dtype=meta.dtype)

    ind = list(range(ndim))
    uc_args = list(concat((x, ind) for x in seq2))
    _, seq2 = unify_chunks(*uc_args)

    assert len({a.chunks for a in seq2}) == 1  # same chunks
    chunks = seq2[0].chunks[:axis] + ((1,) * n,) + seq2[0].chunks[axis:]

    names = [a.name for a in seq2]
    name = "stack-" + tokenize(names, axis)
    keys = list(product([name], *[range(len(bd)) for bd in chunks]))

    inputs = [
        (names[key[axis + 1]],) + key[1 : axis + 1] + key[axis + 2 :] for key in keys
    ]
    values = [
        (
            getitem,
            inp,
            (slice(None, None, None),) * axis
            + (None,)
            + (slice(None, None, None),) * (ndim - axis),
        )
        for inp in inputs
    ]

    layer = dict(zip(keys, values))
    graph = HighLevelGraph.from_collections(name, layer, dependencies=seq2)

    return Array(graph, name, chunks, meta=meta)


def concatenate_shaped(arrays, shape):
    shaped = reshapelist(shape, arrays)
    return concatenate3(shaped)


def concatenate3(arrays):
    """Recursive np.concatenate

    Input should be a nested list of numpy arrays arranged in the order they
    should appear in the array itself.  Each array should have the same number
    of dimensions as the desired output and the nesting of the lists.

    >>> x = np.array([[1, 2]])
    >>> concatenate3([[x, x, x], [x, x, x]])
    array([[1, 2, 1, 2, 1, 2],
           [1, 2, 1, 2, 1, 2]])

    >>> concatenate3([[x, x], [x, x], [x, x]])
    array([[1, 2, 1, 2],
           [1, 2, 1, 2],
           [1, 2, 1, 2]])
    """
    # We need this as __array_function__ may not exist on older NumPy versions.
    # And to reduce verbosity.
    NDARRAY_ARRAY_FUNCTION = getattr(np.ndarray, "__array_function__", None)

    arrays = concrete(arrays)
    if not arrays or all(el is None for el in flatten(arrays)):
        return np.empty(0)

    advanced = max(
        core.flatten(arrays, container=(list, tuple)),
        key=lambda x: getattr(x, "__array_priority__", 0),
    )

    if not all(
        NDARRAY_ARRAY_FUNCTION
        is getattr(type(arr), "__array_function__", NDARRAY_ARRAY_FUNCTION)
        for arr in core.flatten(arrays, container=(list, tuple))
    ):
        try:
            x = unpack_singleton(arrays)
            return _concatenate2(arrays, axes=tuple(range(x.ndim)))
        except TypeError:
            pass

    if concatenate_lookup.dispatch(type(advanced)) is not np.concatenate:
        x = unpack_singleton(arrays)
        return _concatenate2(arrays, axes=list(range(x.ndim)))

    ndim = ndimlist(arrays)
    if not ndim:
        return arrays
    chunks = chunks_from_arrays(arrays)
    shape = tuple(map(sum, chunks))

    def dtype(x):
        try:
            return x.dtype
        except AttributeError:
            return type(x)

    result = np.empty(shape=shape, dtype=dtype(deepfirst(arrays)))

    for idx, arr in zip(
        slices_from_chunks(chunks), core.flatten(arrays, container=(list, tuple))
    ):
        if hasattr(arr, "ndim"):
            while arr.ndim < ndim:
                arr = arr[None, ...]
        result[idx] = arr

    return result


def concatenate_axes(arrays, axes):
    """Recursively call np.concatenate along axes"""
    if len(axes) != ndimlist(arrays):
        raise ValueError("Length of axes should equal depth of nested arrays")

    extradims = max(0, deepfirst(arrays).ndim - (max(axes) + 1))
    return concatenate3(transposelist(arrays, axes, extradims=extradims))


def to_hdf5(filename, *args, chunks=True, **kwargs):
    """Store arrays in HDF5 file

    This saves several dask arrays into several datapaths in an HDF5 file.
    It creates the necessary datasets and handles clean file opening/closing.

    Parameters
    ----------
    chunks: tuple or ``True``
        Chunk shape, or ``True`` to pass the chunks from the dask array.
        Defaults to ``True``.

    Examples
    --------

    >>> da.to_hdf5('myfile.hdf5', '/x', x)  # doctest: +SKIP

    or

    >>> da.to_hdf5('myfile.hdf5', {'/x': x, '/y': y})  # doctest: +SKIP

    Optionally provide arguments as though to ``h5py.File.create_dataset``

    >>> da.to_hdf5('myfile.hdf5', '/x', x, compression='lzf', shuffle=True)  # doctest: +SKIP

    >>> da.to_hdf5('myfile.hdf5', '/x', x, chunks=(10,20,30))  # doctest: +SKIP

    This can also be used as a method on a single Array

    >>> x.to_hdf5('myfile.hdf5', '/x')  # doctest: +SKIP

    See Also
    --------
    da.store
    h5py.File.create_dataset
    """
    if len(args) == 1 and isinstance(args[0], dict):
        data = args[0]
    elif len(args) == 2 and isinstance(args[0], str) and isinstance(args[1], Array):
        data = {args[0]: args[1]}
    else:
        raise ValueError("Please provide {'/data/path': array} dictionary")

    import h5py

    with h5py.File(filename, mode="a") as f:
        dsets = [
            f.require_dataset(
                dp,
                shape=x.shape,
                dtype=x.dtype,
                chunks=tuple(c[0] for c in x.chunks) if chunks is True else chunks,
                **kwargs,
            )
            for dp, x in data.items()
        ]
        store(list(data.values()), dsets)


def interleave_none(a, b):
    """

    >>> interleave_none([0, None, 2, None], [1, 3])
    (0, 1, 2, 3)
    """
    result = []
    i = j = 0
    n = len(a) + len(b)
    while i + j < n:
        if a[i] is not None:
            result.append(a[i])
            i += 1
        else:
            result.append(b[j])
            i += 1
            j += 1
    return tuple(result)


def keyname(name, i, okey):
    """

    >>> keyname('x', 3, [None, None, 0, 2])
    ('x', 3, 0, 2)
    """
    return (name, i) + tuple(k for k in okey if k is not None)


def _vindex(x, *indexes):
    """Point wise indexing with broadcasting.

    >>> x = np.arange(56).reshape((7, 8))
    >>> x
    array([[ 0,  1,  2,  3,  4,  5,  6,  7],
           [ 8,  9, 10, 11, 12, 13, 14, 15],
           [16, 17, 18, 19, 20, 21, 22, 23],
           [24, 25, 26, 27, 28, 29, 30, 31],
           [32, 33, 34, 35, 36, 37, 38, 39],
           [40, 41, 42, 43, 44, 45, 46, 47],
           [48, 49, 50, 51, 52, 53, 54, 55]])

    >>> d = from_array(x, chunks=(3, 4))
    >>> result = _vindex(d, [0, 1, 6, 0], [0, 1, 0, 7])
    >>> result.compute()
    array([ 0,  9, 48,  7])
    """
    indexes = replace_ellipsis(x.ndim, indexes)

    nonfancy_indexes = []
    reduced_indexes = []
    for ind in indexes:
        if isinstance(ind, Number):
            nonfancy_indexes.append(ind)
        elif isinstance(ind, slice):
            nonfancy_indexes.append(ind)
            reduced_indexes.append(slice(None))
        else:
            nonfancy_indexes.append(slice(None))
            reduced_indexes.append(ind)

    nonfancy_indexes = tuple(nonfancy_indexes)
    reduced_indexes = tuple(reduced_indexes)

    x = x[nonfancy_indexes]

    array_indexes = {}
    for i, (ind, size) in enumerate(zip(reduced_indexes, x.shape)):
        if not isinstance(ind, slice):
            ind = np.array(ind, copy=True)
            if ind.dtype.kind == "b":
                raise IndexError("vindex does not support indexing with boolean arrays")
            if ((ind >= size) | (ind < -size)).any():
                raise IndexError(
                    "vindex key has entries out of bounds for "
                    "indexing along axis %s of size %s: %r" % (i, size, ind)
                )
            ind %= size
            array_indexes[i] = ind

    if array_indexes:
        x = _vindex_array(x, array_indexes)

    return x


def _vindex_array(x, dict_indexes):
    """Point wise indexing with only NumPy Arrays."""

    token = tokenize(x, dict_indexes)
    try:
        broadcast_shape = np.broadcast_shapes(
            *(arr.shape for arr in dict_indexes.values())
        )
    except ValueError as e:
        # note: error message exactly matches numpy
        shapes_str = " ".join(str(a.shape) for a in dict_indexes.values())
        raise IndexError(
            "shape mismatch: indexing arrays could not be "
            "broadcast together with shapes " + shapes_str
        ) from e
    npoints = math.prod(broadcast_shape)
    axes = [i for i in range(x.ndim) if i in dict_indexes]

    def _subset_to_indexed_axes(iterable):
        for i, elem in enumerate(iterable):
            if i in axes:
                yield elem

    bounds2 = tuple(
        np.array(cached_cumsum(c, initial_zero=True))
        for c in _subset_to_indexed_axes(x.chunks)
    )
    axis = _get_axis(tuple(i if i in axes else None for i in range(x.ndim)))
    out_name = "vindex-merge-" + token

    # Now compute indices of each output element within each input block
    # The index is relative to the block, not the array.
    block_idxs = tuple(
        np.searchsorted(b, ind, side="right") - 1
        for b, ind in zip(bounds2, dict_indexes.values())
    )
    starts = (b[i] for i, b in zip(block_idxs, bounds2))
    inblock_idxs = []
    for idx, start in zip(dict_indexes.values(), starts):
        a = idx - start
        if len(a) > 0:
            dtype = np.min_scalar_type(np.max(a, axis=None))
            inblock_idxs.append(a.astype(dtype, copy=False))
        else:
            inblock_idxs.append(a)

    inblock_idxs = np.broadcast_arrays(*inblock_idxs)

    chunks = [c for i, c in enumerate(x.chunks) if i not in axes]
    # determine number of points in one single output block.
    # Use the input chunk size to determine this.
    max_chunk_point_dimensions = reduce(
        mul, map(cached_max, _subset_to_indexed_axes(x.chunks))
    )

    n_chunks, remainder = divmod(npoints, max_chunk_point_dimensions)
    chunks.insert(
        0,
        (
            (max_chunk_point_dimensions,) * n_chunks
            + ((remainder,) if remainder > 0 else ())
            if npoints > 0
            else (0,)
        ),
    )
    chunks = tuple(chunks)

    if npoints > 0:
        other_blocks = product(
            *[
                range(len(c)) if i not in axes else [None]
                for i, c in enumerate(x.chunks)
            ]
        )

        full_slices = [
            slice(None, None) if i not in axes else None for i in range(x.ndim)
        ]

        # The output is constructed as a new dimension and then reshaped
        # So the index of the output point is simply an `arange`
        outinds = np.arange(npoints).reshape(broadcast_shape)
        # Which output block is the point going to, and what is the index within that block?
        outblocks, outblock_idx = np.divmod(outinds, max_chunk_point_dimensions)

        # Now we try to be clever. We need to construct a graph where
        # if input chunk (0,0) contributes to output chunks 0 and 2 then
        # (0,0) => (0, 0, 0) and (2, 0, 0).
        # The following is a groupby over output key by using ravel_multi_index to convert
        # the (output_block, *input_block) tuple to a unique code.
        ravel_shape = (n_chunks + 1, *_subset_to_indexed_axes(x.numblocks))
        keys = np.ravel_multi_index([outblocks, *block_idxs], ravel_shape)
        # by sorting the data that needs to be inserted in the graph here,
        # we can slice in the hot loop below instead of using fancy indexing which will
        # always copy inside the hot loop.
        sortidx = np.argsort(keys, axis=None)
        sorted_keys = keys.flat[sortidx]  # flattens
        sorted_inblock_idxs = [_.flat[sortidx] for _ in inblock_idxs]
        sorted_outblock_idx = outblock_idx.flat[sortidx]
        dtype = np.min_scalar_type(max_chunk_point_dimensions)
        sorted_outblock_idx = sorted_outblock_idx.astype(dtype, copy=False)
        # Determine the start and end of each unique key. We will loop over this
        flag = np.concatenate([[True], sorted_keys[1:] != sorted_keys[:-1], [True]])
        (key_bounds,) = flag.nonzero()

        name = "vindex-slice-" + token
        vindex_merge_name = "vindex-merge-" + token
        dsk = {}
        for okey in other_blocks:
            merge_inputs = defaultdict(list)
            merge_indexer = defaultdict(list)
            for i, (start, stop) in enumerate(
                zip(key_bounds[:-1], key_bounds[1:], strict=True)
            ):
                slicer = slice(start, stop)
                key = sorted_keys[start]
                outblock, *input_blocks = np.unravel_index(key, ravel_shape)
                inblock = [_[slicer] for _ in sorted_inblock_idxs]
                k = keyname(name, i, okey)
                dsk[k] = Task(
                    k,
                    _vindex_slice_and_transpose,
                    TaskRef((x.name,) + interleave_none(okey, input_blocks)),
                    interleave_none(full_slices, inblock),
                    axis,
                )
                merge_inputs[outblock].append(TaskRef(keyname(name, i, okey)))
                merge_indexer[outblock].append(sorted_outblock_idx[slicer])

            for i in merge_inputs.keys():
                k = keyname(vindex_merge_name, i, okey)
                dsk[k] = Task(
                    k,
                    _vindex_merge,
                    merge_indexer[i],
                    List(merge_inputs[i]),
                )

        result_1d = Array(
            HighLevelGraph.from_collections(out_name, dsk, dependencies=[x]),
            out_name,
            chunks,
            x.dtype,
            meta=x._meta,
        )
        return result_1d.reshape(broadcast_shape + result_1d.shape[1:])

    # output has a zero dimension, just create a new zero-shape array with the
    # same dtype
    from dask.array.wrap import empty

    result_1d = empty(
        tuple(map(sum, chunks)), chunks=chunks, dtype=x.dtype, name=out_name
    )
    return result_1d.reshape(broadcast_shape + result_1d.shape[1:])


def _get_axis(indexes):
    """Get axis along which point-wise slicing results lie

    This is mostly a hack because I can't figure out NumPy's rule on this and
    can't be bothered to go reading.

    >>> _get_axis([[1, 2], None, [1, 2], None])
    0
    >>> _get_axis([None, [1, 2], [1, 2], None])
    1
    >>> _get_axis([None, None, [1, 2], [1, 2]])
    2
    """
    ndim = len(indexes)
    indexes = [slice(None, None) if i is None else [0] for i in indexes]
    x = np.empty((2,) * ndim)
    x2 = x[tuple(indexes)]
    return x2.shape.index(1)


def _vindex_slice_and_transpose(block, points, axis):
    """Pull out point-wise slices from block and rotate block so that
    points are on the first dimension"""
    points = [p if isinstance(p, slice) else list(p) for p in points]
    block = block[tuple(points)]
    axes = [axis] + list(range(axis)) + list(range(axis + 1, block.ndim))
    return block.transpose(axes)


def _vindex_merge(locations, values):
    """

    >>> locations = [0], [2, 1]
    >>> values = [np.array([[1, 2, 3]]),
    ...           np.array([[10, 20, 30], [40, 50, 60]])]

    >>> _vindex_merge(locations, values)
    array([[ 1,  2,  3],
           [40, 50, 60],
           [10, 20, 30]])
    """
    locations = list(map(list, locations))
    values = list(values)

    n = sum(map(len, locations))

    shape = list(values[0].shape)
    shape[0] = n
    shape = tuple(shape)

    dtype = values[0].dtype

    x = np.empty_like(values[0], dtype=dtype, shape=shape)

    ind = [slice(None, None) for i in range(x.ndim)]
    for loc, val in zip(locations, values):
        ind[0] = loc
        x[tuple(ind)] = val

    return x


def to_npy_stack(dirname, x, axis=0):
    """Write dask array to a stack of .npy files

    This partitions the dask.array along one axis and stores each block along
    that axis as a single .npy file in the specified directory

    Examples
    --------
    >>> x = da.ones((5, 10, 10), chunks=(2, 4, 4))  # doctest: +SKIP
    >>> da.to_npy_stack('data/', x, axis=0)  # doctest: +SKIP

    The ``.npy`` files store numpy arrays for ``x[0:2], x[2:4], and x[4:5]``
    respectively, as is specified by the chunk size along the zeroth axis::

        $ tree data/
        data/
        |-- 0.npy
        |-- 1.npy
        |-- 2.npy
        |-- info

    The ``info`` file stores the dtype, chunks, and axis information of the array.
    You can load these stacks with the :func:`dask.array.from_npy_stack` function.

    >>> y = da.from_npy_stack('data/')  # doctest: +SKIP

    See Also
    --------
    from_npy_stack
    """

    chunks = tuple((c if i == axis else (sum(c),)) for i, c in enumerate(x.chunks))
    xx = x.rechunk(chunks)

    if not os.path.exists(dirname):
        os.mkdir(dirname)

    meta = {"chunks": chunks, "dtype": x.dtype, "axis": axis}

    with open(os.path.join(dirname, "info"), "wb") as f:
        pickle.dump(meta, f)

    name = "to-npy-stack-" + str(uuid.uuid1())
    dsk = {
        (name, i): (np.save, os.path.join(dirname, "%d.npy" % i), key)
        for i, key in enumerate(core.flatten(xx.__dask_keys__()))
    }

    graph = HighLevelGraph.from_collections(name, dsk, dependencies=[xx])
    compute_as_if_collection(Array, graph, list(dsk))


def from_npy_stack(dirname, mmap_mode="r"):
    """Load dask array from stack of npy files

    Parameters
    ----------
    dirname: string
        Directory of .npy files
    mmap_mode: (None or 'r')
        Read data in memory map mode

    See Also
    --------
    to_npy_stack
    """
    with open(os.path.join(dirname, "info"), "rb") as f:
        info = pickle.load(f)

    dtype = info["dtype"]
    chunks = info["chunks"]
    axis = info["axis"]

    name = "from-npy-stack-%s" % dirname
    keys = list(product([name], *[range(len(c)) for c in chunks]))
    values = [
        (np.load, os.path.join(dirname, "%d.npy" % i), mmap_mode)
        for i in range(len(chunks[axis]))
    ]
    dsk = dict(zip(keys, values))

    return Array(dsk, name, chunks, dtype)


def new_da_object(dsk, name, chunks, meta=None, dtype=None):
    """Generic constructor for dask.array or dask.dataframe objects.

    Decides the appropriate output class based on the type of `meta` provided.
    """
    if is_dataframe_like(meta) or is_series_like(meta) or is_index_like(meta):
        from dask.dataframe import from_graph

        assert all(len(c) == 1 for c in chunks[1:])
        divisions = [None] * (len(chunks[0]) + 1)
        return from_graph(dict(dsk), meta, divisions, dsk.layers[name].keys(), name)
    else:
        return Array(dsk, name=name, chunks=chunks, meta=meta, dtype=dtype)


class BlockView:
    """An array-like interface to the blocks of an array.

    ``BlockView`` provides an array-like interface
    to the blocks of a dask array.  Numpy-style indexing of a
     ``BlockView`` returns a selection of blocks as a new dask array.

    You can index ``BlockView`` like a numpy array of shape
    equal to the number of blocks in each dimension, (available as
    array.blocks.size).  The dimensionality of the output array matches
    the dimension of this array, even if integer indices are passed.
    Slicing with ``np.newaxis`` or multiple lists is not supported.

    Examples
    --------
    >>> import dask.array as da
    >>> from dask.array.core import BlockView
    >>> x = da.arange(8, chunks=2)
    >>> bv = BlockView(x)
    >>> bv.shape # aliases x.numblocks
    (4,)
    >>> bv.size
    4
    >>> bv[0].compute()
    array([0, 1])
    >>> bv[:3].compute()
    array([0, 1, 2, 3, 4, 5])
    >>> bv[::2].compute()
    array([0, 1, 4, 5])
    >>> bv[[-1, 0]].compute()
    array([6, 7, 0, 1])
    >>> bv.ravel()  # doctest: +NORMALIZE_WHITESPACE
    [dask.array<blocks, shape=(2,), dtype=int64, chunksize=(2,), chunktype=numpy.ndarray>,
     dask.array<blocks, shape=(2,), dtype=int64, chunksize=(2,), chunktype=numpy.ndarray>,
     dask.array<blocks, shape=(2,), dtype=int64, chunksize=(2,), chunktype=numpy.ndarray>,
     dask.array<blocks, shape=(2,), dtype=int64, chunksize=(2,), chunktype=numpy.ndarray>]

    Returns
    -------
    An instance of ``da.array.Blockview``
    """

    def __init__(self, array: Array):
        self._array = array

    def __getitem__(self, index: Any) -> Array:
        from dask.array.slicing import normalize_index

        if not isinstance(index, tuple):
            index = (index,)
        if sum(isinstance(ind, (np.ndarray, list)) for ind in index) > 1:
            raise ValueError("Can only slice with a single list")
        if any(ind is None for ind in index):
            raise ValueError("Slicing with np.newaxis or None is not supported")
        index = normalize_index(index, self._array.numblocks)
        index = tuple(
            slice(k, k + 1) if isinstance(k, Number) else k  # type: ignore
            for k in index
        )

        name = "blocks-" + tokenize(self._array, index)

        new_keys = self._array._key_array[index]

        chunks = tuple(
            tuple(np.array(c)[i].tolist()) for c, i in zip(self._array.chunks, index)
        )

        keys = product(*(range(len(c)) for c in chunks))

        graph: Graph = {(name,) + key: tuple(new_keys[key].tolist()) for key in keys}

        hlg = HighLevelGraph.from_collections(name, graph, dependencies=[self._array])
        return Array(hlg, name, chunks, meta=self._array)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, BlockView):
            return self._array is other._array
        else:
            return NotImplemented

    @property
    def size(self) -> int:
        """
        The total number of blocks in the array.
        """
        return math.prod(self.shape)

    @property
    def shape(self) -> tuple[int, ...]:
        """
        The number of blocks per axis. Alias of ``dask.array.numblocks``.
        """
        return self._array.numblocks

    def ravel(self) -> list[Array]:
        """
        Return a flattened list of all the blocks in the array in C order.
        """
        return [self[idx] for idx in np.ndindex(self.shape)]


def _numpy_vindex(indexer, arr):
    return arr[indexer]


from dask.array.blockwise import blockwise
