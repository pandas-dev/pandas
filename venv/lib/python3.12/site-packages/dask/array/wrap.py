from __future__ import annotations

from functools import partial
from itertools import product

import numpy as np
from tlz import curry

from dask.array.backends import array_creation_dispatch
from dask.array.core import Array, normalize_chunks
from dask.array.utils import meta_from_array
from dask.base import tokenize
from dask.blockwise import blockwise as core_blockwise
from dask.layers import ArrayChunkShapeDep
from dask.utils import funcname


def _parse_wrap_args(func, args, kwargs, shape):
    if isinstance(shape, np.ndarray):
        shape = shape.tolist()

    if not isinstance(shape, (tuple, list)):
        shape = (shape,)

    name = kwargs.pop("name", None)
    chunks = kwargs.pop("chunks", "auto")

    dtype = kwargs.pop("dtype", None)
    if dtype is None:
        dtype = func(shape, *args, **kwargs).dtype
    dtype = np.dtype(dtype)

    chunks = normalize_chunks(chunks, shape, dtype=dtype)

    name = name or funcname(func) + "-" + tokenize(
        func, shape, chunks, dtype, args, kwargs
    )

    return {
        "shape": shape,
        "dtype": dtype,
        "kwargs": kwargs,
        "chunks": chunks,
        "name": name,
    }


def wrap_func_shape_as_first_arg(func, *args, **kwargs):
    """
    Transform np creation function into blocked version
    """
    if "shape" not in kwargs:
        shape, args = args[0], args[1:]
    else:
        shape = kwargs.pop("shape")

    if isinstance(shape, Array):
        raise TypeError(
            "Dask array input not supported. "
            "Please use tuple, list, or a 1D numpy array instead."
        )

    parsed = _parse_wrap_args(func, args, kwargs, shape)
    shape = parsed["shape"]
    dtype = parsed["dtype"]
    chunks = parsed["chunks"]
    name = parsed["name"]
    kwargs = parsed["kwargs"]
    func = partial(func, dtype=dtype, **kwargs)

    out_ind = dep_ind = tuple(range(len(shape)))
    graph = core_blockwise(
        func,
        name,
        out_ind,
        ArrayChunkShapeDep(chunks),
        dep_ind,
        numblocks={},
    )

    return Array(graph, name, chunks, dtype=dtype, meta=kwargs.get("meta", None))


def wrap_func_like(func, *args, **kwargs):
    """
    Transform np creation function into blocked version
    """
    x = args[0]
    meta = meta_from_array(x)
    shape = kwargs.get("shape", x.shape)

    parsed = _parse_wrap_args(func, args, kwargs, shape)
    shape = parsed["shape"]
    dtype = parsed["dtype"]
    chunks = parsed["chunks"]
    name = parsed["name"]
    kwargs = parsed["kwargs"]

    keys = product([name], *[range(len(bd)) for bd in chunks])
    shapes = product(*chunks)
    shapes = list(shapes)
    kw = [kwargs for _ in shapes]
    for i, s in enumerate(list(shapes)):
        kw[i]["shape"] = s
    vals = ((partial(func, dtype=dtype, **k),) + args for (k, s) in zip(kw, shapes))

    dsk = dict(zip(keys, vals))

    return Array(dsk, name, chunks, meta=meta.astype(dtype))


@curry
def wrap(wrap_func, func, func_like=None, **kwargs):
    if func_like is None:
        f = partial(wrap_func, func, **kwargs)
    else:
        f = partial(wrap_func, func_like, **kwargs)
    template = """
    Blocked variant of %(name)s

    Follows the signature of %(name)s exactly except that it also features
    optional keyword arguments ``chunks: int, tuple, or dict`` and ``name: str``.

    Original signature follows below.
    """
    if func.__doc__ is not None:
        f.__doc__ = template % {"name": func.__name__} + func.__doc__
        f.__name__ = "blocked_" + func.__name__
    return f


w = wrap(wrap_func_shape_as_first_arg)


@curry
def _broadcast_trick_inner(func, shape, meta=(), *args, **kwargs):
    # cupy-specific hack. numpy is happy with hardcoded shape=().
    null_shape = () if shape == () else 1

    return np.broadcast_to(func(meta, *args, shape=null_shape, **kwargs), shape)


def broadcast_trick(func):
    """
    Provide a decorator to wrap common numpy function with a broadcast trick.

    Dask arrays are currently immutable; thus when we know an array is uniform,
    we can replace the actual data by a single value and have all elements point
    to it, thus reducing the size.

    >>> x = np.broadcast_to(1, (100,100,100))
    >>> x.base.nbytes
    8

    Those array are not only more efficient locally, but dask serialisation is
    aware of the _real_ size of those array and thus can send them around
    efficiently and schedule accordingly.

    Note that those array are read-only and numpy will refuse to assign to them,
    so should be safe.
    """
    inner = _broadcast_trick_inner(func)
    inner.__doc__ = func.__doc__
    inner.__name__ = func.__name__
    return inner


ones = array_creation_dispatch.register_inplace(
    backend="numpy",
    name="ones",
)(w(broadcast_trick(np.ones_like), dtype="f8"))


zeros = array_creation_dispatch.register_inplace(
    backend="numpy",
    name="zeros",
)(w(broadcast_trick(np.zeros_like), dtype="f8"))


empty = array_creation_dispatch.register_inplace(
    backend="numpy",
    name="empty",
)(w(broadcast_trick(np.empty_like), dtype="f8"))


w_like = wrap(wrap_func_like)


empty_like = w_like(np.empty, func_like=np.empty_like)


# full and full_like require special casing due to argument check on fill_value
# Generate wrapped functions only once
_full = array_creation_dispatch.register_inplace(
    backend="numpy",
    name="full",
)(w(broadcast_trick(np.full_like)))
_full_like = w_like(np.full, func_like=np.full_like)


# workaround for numpy doctest failure: https://github.com/numpy/numpy/pull/17472
if _full.__doc__ is not None:
    _full.__doc__ = _full.__doc__.replace(
        "array([0.1,  0.1,  0.1,  0.1,  0.1,  0.1])",
        "array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1])",
    )
    _full.__doc__ = _full.__doc__.replace(
        ">>> np.full_like(y, [0, 0, 255])",
        ">>> np.full_like(y, [0, 0, 255])  # doctest: +NORMALIZE_WHITESPACE",
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


def full_like(a, fill_value, *args, **kwargs):
    if np.ndim(fill_value) != 0:
        raise ValueError(
            f"fill_value must be scalar. Received {type(fill_value).__name__} instead."
        )
    return _full_like(
        *args,
        a=a,
        fill_value=fill_value,
        **kwargs,
    )


full.__doc__ = _full.__doc__
full_like.__doc__ = _full_like.__doc__
