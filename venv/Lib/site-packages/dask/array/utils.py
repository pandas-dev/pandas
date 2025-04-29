from __future__ import annotations

import contextlib
import functools
import itertools
import math
import numbers
import warnings

import numpy as np
from tlz import concat, frequencies

from dask.array.numpy_compat import AxisError
from dask.base import is_dask_collection, tokenize
from dask.highlevelgraph import HighLevelGraph
from dask.utils import has_keyword, is_arraylike, is_cupy_type, typename


def normalize_to_array(x):
    if is_cupy_type(x):
        return x.get()
    else:
        return x


def meta_from_array(x, ndim=None, dtype=None):
    """Normalize an array to appropriate meta object

    Parameters
    ----------
    x: array-like, callable
        Either an object that looks sufficiently like a Numpy array,
        or a callable that accepts shape and dtype keywords
    ndim: int
        Number of dimensions of the array
    dtype: Numpy dtype
        A valid input for ``np.dtype``

    Returns
    -------
    array-like with zero elements of the correct dtype
    """
    # If using x._meta, x must be a Dask Array, some libraries (e.g. zarr)
    # implement a _meta attribute that are incompatible with Dask Array._meta
    if hasattr(x, "_meta") and is_dask_collection(x) and is_arraylike(x):
        x = x._meta

    if dtype is None and x is None:
        raise ValueError("You must specify the meta or dtype of the array")

    if np.isscalar(x):
        x = np.array(x)

    if x is None:
        x = np.ndarray
    elif dtype is None and hasattr(x, "dtype"):
        dtype = x.dtype

    if isinstance(x, type):
        x = x(shape=(0,) * (ndim or 0), dtype=dtype)

    if isinstance(x, list) or isinstance(x, tuple):
        ndims = [
            (
                0
                if isinstance(a, numbers.Number)
                else a.ndim if hasattr(a, "ndim") else len(a)
            )
            for a in x
        ]
        a = [a if nd == 0 else meta_from_array(a, nd) for a, nd in zip(x, ndims)]
        return a if isinstance(x, list) else tuple(x)

    if (
        not hasattr(x, "shape")
        or not hasattr(x, "dtype")
        or not isinstance(x.shape, tuple)
    ):
        return x

    if ndim is None:
        ndim = x.ndim

    try:
        meta = x[tuple(slice(0, 0, None) for _ in range(x.ndim))]
        if meta.ndim != ndim:
            if ndim > x.ndim:
                meta = meta[(Ellipsis,) + tuple(None for _ in range(ndim - meta.ndim))]
                meta = meta[tuple(slice(0, 0, None) for _ in range(meta.ndim))]
            elif ndim == 0:
                meta = meta.sum()
            else:
                meta = meta.reshape((0,) * ndim)
        if meta is np.ma.masked:
            meta = np.ma.array(np.empty((0,) * ndim, dtype=dtype or x.dtype), mask=True)
    except Exception:
        meta = np.empty((0,) * ndim, dtype=dtype or x.dtype)

    if np.isscalar(meta):
        meta = np.array(meta)

    if dtype and meta.dtype != dtype:
        try:
            meta = meta.astype(dtype)
        except ValueError as e:
            if (
                any(
                    s in str(e)
                    for s in [
                        "invalid literal",
                        "could not convert string to float",
                    ]
                )
                and meta.dtype.kind in "SU"
            ):
                meta = np.array([]).astype(dtype)
            else:
                raise e

    return meta


def compute_meta(func, _dtype, *args, **kwargs):
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)

        args_meta = [meta_from_array(x) if is_arraylike(x) else x for x in args]
        kwargs_meta = {
            k: meta_from_array(v) if is_arraylike(v) else v for k, v in kwargs.items()
        }

        # todo: look for alternative to this, causes issues when using map_blocks()
        # with np.vectorize, such as dask.array.routines._isnonzero_vec().
        if isinstance(func, np.vectorize):
            meta = func(*args_meta)
        else:
            try:
                # some reduction functions need to know they are computing meta
                if has_keyword(func, "computing_meta"):
                    kwargs_meta["computing_meta"] = True
                meta = func(*args_meta, **kwargs_meta)
            except TypeError as e:
                if any(
                    s in str(e)
                    for s in [
                        "unexpected keyword argument",
                        "is an invalid keyword for",
                        "Did not understand the following kwargs",
                    ]
                ):
                    raise
                else:
                    return None
            except ValueError as e:
                # min/max functions have no identity, just use the same input type when there's only one
                if len(
                    args_meta
                ) == 1 and "zero-size array to reduction operation" in str(e):
                    meta = args_meta[0]
                else:
                    return None
            except Exception:
                return None

        if _dtype and getattr(meta, "dtype", None) != _dtype:
            with contextlib.suppress(AttributeError):
                meta = meta.astype(_dtype)

        if np.isscalar(meta):
            meta = np.array(meta)

        return meta


def allclose(a, b, equal_nan=False, **kwargs):
    a = normalize_to_array(a)
    b = normalize_to_array(b)
    if getattr(a, "dtype", None) != "O":
        if hasattr(a, "mask") or hasattr(b, "mask"):
            return np.ma.allclose(a, b, masked_equal=True, **kwargs)
        else:
            return np.allclose(a, b, equal_nan=equal_nan, **kwargs)
    if equal_nan:
        return a.shape == b.shape and all(
            np.isnan(b) if np.isnan(a) else a == b for (a, b) in zip(a.flat, b.flat)
        )
    return (a == b).all()


def same_keys(a, b):
    def key(k):
        if isinstance(k, str):
            return (k, -1, -1, -1)
        else:
            return k

    return sorted(a.dask, key=key) == sorted(b.dask, key=key)


def _not_empty(x):
    return x.shape and 0 not in x.shape


def _check_dsk(dsk):
    """Check that graph is well named and non-overlapping"""
    if not isinstance(dsk, HighLevelGraph):
        return

    dsk.validate()
    assert all(isinstance(k, (tuple, str)) for k in dsk.layers)
    freqs = frequencies(concat(dsk.layers.values()))
    non_one = {k: v for k, v in freqs.items() if v != 1}
    key_collisions = set()
    # Allow redundant keys if the values are equivalent
    collisions = set()
    for k in non_one.keys():
        for layer in dsk.layers.values():
            try:
                key_collisions.add(tokenize(layer[k]))
            except KeyError:
                pass
        if len(key_collisions) >= 2:
            collisions.add(k)
        key_collisions.clear()

    assert len(collisions) == 0, {k: non_one[k] for k in collisions}


def assert_eq_shape(a, b, check_ndim=True, check_nan=True):
    if check_ndim:
        assert len(a) == len(b)

    for aa, bb in zip(a, b):
        if math.isnan(aa) or math.isnan(bb):
            if check_nan:
                assert math.isnan(aa) == math.isnan(bb)
        else:
            assert aa == bb


def _check_chunks(x, check_ndim=True, scheduler=None):
    x = x.persist(scheduler=scheduler)
    for idx in itertools.product(*(range(len(c)) for c in x.chunks)):
        chunk = x.dask[(x.name,) + idx]
        if hasattr(chunk, "result"):  # it's a future
            chunk = chunk.result()
        if not hasattr(chunk, "dtype"):
            chunk = np.array(chunk, dtype="O")
        expected_shape = tuple(c[i] for c, i in zip(x.chunks, idx))
        assert_eq_shape(
            expected_shape, chunk.shape, check_ndim=check_ndim, check_nan=False
        )
        assert (
            chunk.dtype == x.dtype
        ), "maybe you forgot to pass the scheduler to `assert_eq`?"
    return x


def _get_dt_meta_computed(
    x,
    check_shape=True,
    check_graph=True,
    check_chunks=True,
    check_ndim=True,
    scheduler=None,
):
    x_original = x
    x_meta = None
    x_computed = None

    if is_dask_collection(x) and is_arraylike(x):
        assert x.dtype is not None
        adt = x.dtype
        if check_graph:
            _check_dsk(x.dask)
        x_meta = getattr(x, "_meta", None)
        if check_chunks:
            # Replace x with persisted version to avoid computing it twice.
            x = _check_chunks(x, check_ndim=check_ndim, scheduler=scheduler)
        x = x.compute(scheduler=scheduler)
        x_computed = x
        if hasattr(x, "todense"):
            x = x.todense()
        if not hasattr(x, "dtype"):
            x = np.array(x, dtype="O")
        if _not_empty(x):
            assert x.dtype == x_original.dtype
        if check_shape:
            assert_eq_shape(x_original.shape, x.shape, check_nan=False)
    else:
        if not hasattr(x, "dtype"):
            x = np.array(x, dtype="O")
        adt = getattr(x, "dtype", None)

    return x, adt, x_meta, x_computed


def assert_eq(
    a,
    b,
    check_shape=True,
    check_graph=True,
    check_meta=True,
    check_chunks=True,
    check_ndim=True,
    check_type=True,
    check_dtype=True,
    equal_nan=True,
    scheduler="sync",
    **kwargs,
):
    a_original = a
    b_original = b

    if isinstance(a, (list, int, float)):
        a = np.array(a)
    if isinstance(b, (list, int, float)):
        b = np.array(b)

    a, adt, a_meta, a_computed = _get_dt_meta_computed(
        a,
        check_shape=check_shape,
        check_graph=check_graph,
        check_chunks=check_chunks,
        check_ndim=check_ndim,
        scheduler=scheduler,
    )
    b, bdt, b_meta, b_computed = _get_dt_meta_computed(
        b,
        check_shape=check_shape,
        check_graph=check_graph,
        check_chunks=check_chunks,
        check_ndim=check_ndim,
        scheduler=scheduler,
    )

    if check_dtype and str(adt) != str(bdt):
        raise AssertionError(f"a and b have different dtypes: (a: {adt}, b: {bdt})")

    try:
        assert (
            a.shape == b.shape
        ), f"a and b have different shapes (a: {a.shape}, b: {b.shape})"
        if check_type:
            _a = a if a.shape else a.item()
            _b = b if b.shape else b.item()
            assert type(_a) == type(
                _b
            ), f"a and b have different types (a: {type(_a)}, b: {type(_b)})"
        if check_meta:
            if hasattr(a, "_meta") and hasattr(b, "_meta"):
                assert_eq(a._meta, b._meta)
            if hasattr(a_original, "_meta"):
                msg = (
                    f"compute()-ing 'a' changes its number of dimensions "
                    f"(before: {a_original._meta.ndim}, after: {a.ndim})"
                )
                assert a_original._meta.ndim == a.ndim, msg
                if a_meta is not None:
                    msg = (
                        f"compute()-ing 'a' changes its type "
                        f"(before: {type(a_original._meta)}, after: {type(a_meta)})"
                    )
                    assert type(a_original._meta) == type(a_meta), msg
                    if not (np.isscalar(a_meta) or np.isscalar(a_computed)):
                        msg = (
                            f"compute()-ing 'a' results in a different type than implied by its metadata "
                            f"(meta: {type(a_meta)}, computed: {type(a_computed)})"
                        )
                        assert type(a_meta) == type(a_computed), msg
            if hasattr(b_original, "_meta"):
                msg = (
                    f"compute()-ing 'b' changes its number of dimensions "
                    f"(before: {b_original._meta.ndim}, after: {b.ndim})"
                )
                assert b_original._meta.ndim == b.ndim, msg
                if b_meta is not None:
                    msg = (
                        f"compute()-ing 'b' changes its type "
                        f"(before: {type(b_original._meta)}, after: {type(b_meta)})"
                    )
                    assert type(b_original._meta) == type(b_meta), msg
                    if not (np.isscalar(b_meta) or np.isscalar(b_computed)):
                        msg = (
                            f"compute()-ing 'b' results in a different type than implied by its metadata "
                            f"(meta: {type(b_meta)}, computed: {type(b_computed)})"
                        )
                        assert type(b_meta) == type(b_computed), msg
        msg = "found values in 'a' and 'b' which differ by more than the allowed amount"
        assert allclose(a, b, equal_nan=equal_nan, **kwargs), msg
        return True
    except TypeError:
        pass

    c = a == b

    if isinstance(c, np.ndarray):
        assert c.all()
    else:
        assert c

    return True


def safe_wraps(wrapped, assigned=functools.WRAPPER_ASSIGNMENTS):
    """Like functools.wraps, but safe to use even if wrapped is not a function.

    Only needed on Python 2.
    """
    if all(hasattr(wrapped, attr) for attr in assigned):
        return functools.wraps(wrapped, assigned=assigned)
    else:
        return lambda x: x


def _dtype_of(a):
    """Determine dtype of an array-like."""
    try:
        # Check for the attribute before using asanyarray, because some types
        # (notably sparse arrays) don't work with it.
        return a.dtype
    except AttributeError:
        return np.asanyarray(a).dtype


def arange_safe(*args, like, **kwargs):
    """
    Use the `like=` from `np.arange` to create a new array dispatching
    to the downstream library. If that fails, falls back to the
    default NumPy behavior, resulting in a `numpy.ndarray`.
    """
    if like is None:
        return np.arange(*args, **kwargs)
    else:
        try:
            return np.arange(*args, like=meta_from_array(like), **kwargs)
        except TypeError:
            return np.arange(*args, **kwargs)


def _array_like_safe(np_func, da_func, a, like, **kwargs):
    from dask.array import Array

    if like is a and hasattr(a, "__array_function__"):
        return a

    if isinstance(like, Array):
        return da_func(a, **kwargs)

    if isinstance(a, Array) and is_cupy_type(a._meta):
        a = a.compute(scheduler="sync")

    if hasattr(like, "__array_function__"):
        return np_func(a, like=like, **kwargs)

    if type(like).__module__.startswith("scipy.sparse"):
        # e.g. scipy.sparse.csr_matrix
        kwargs.pop("order", None)
        return type(like)(a, **kwargs)

    # Unknown namespace with no __array_function__ support.
    # Quietly disregard like= parameter.
    # FIXME is there a use case for this?
    return np_func(a, **kwargs)


def array_safe(a, like, **kwargs):
    """
    If `a` is `dask.array`, return `dask.array.asarray(a, **kwargs)`,
    otherwise return `np.asarray(a, like=like, **kwargs)`, dispatching
    the call to the library that implements the like array. Note that
    when `a` is a `dask.Array` backed by `cupy.ndarray` but `like`
    isn't, this function will call `a.compute(scheduler="sync")`
    before `np.array`, as downstream libraries are unlikely to know how
    to convert a `dask.Array` and CuPy doesn't implement `__array__` to
    prevent implicit copies to host.
    """
    from dask.array.routines import array

    return _array_like_safe(np.array, array, a, like, **kwargs)


def asarray_safe(a, like, **kwargs):
    """
    If a is dask.array, return dask.array.asarray(a, **kwargs),
    otherwise return np.asarray(a, like=like, **kwargs), dispatching
    the call to the library that implements the like array. Note that
    when a is a dask.Array but like isn't, this function will call
    a.compute(scheduler="sync") before np.asarray, as downstream
    libraries are unlikely to know how to convert a dask.Array.
    """
    from dask.array.core import asarray

    return _array_like_safe(np.asarray, asarray, a, like, **kwargs)


def asanyarray_safe(a, like, **kwargs):
    """
    If a is dask.array, return dask.array.asanyarray(a, **kwargs),
    otherwise return np.asanyarray(a, like=like, **kwargs), dispatching
    the call to the library that implements the like array. Note that
    when a is a dask.Array but like isn't, this function will call
    a.compute(scheduler="sync") before np.asanyarray, as downstream
    libraries are unlikely to know how to convert a dask.Array.
    """
    from dask.array.core import asanyarray

    return _array_like_safe(np.asanyarray, asanyarray, a, like, **kwargs)


def validate_axis(axis, ndim):
    """Validate an input to axis= keywords"""
    if isinstance(axis, (tuple, list)):
        return tuple(validate_axis(ax, ndim) for ax in axis)
    if not isinstance(axis, numbers.Integral):
        raise TypeError("Axis value must be an integer, got %s" % axis)
    if axis < -ndim or axis >= ndim:
        raise AxisError(
            "Axis %d is out of bounds for array of dimension %d" % (axis, ndim)
        )
    if axis < 0:
        axis += ndim
    return axis


def svd_flip(u, v, u_based_decision=False):
    """Sign correction to ensure deterministic output from SVD.

    This function is useful for orienting eigenvectors such that
    they all lie in a shared but arbitrary half-space. This makes
    it possible to ensure that results are equivalent across SVD
    implementations and random number generator states.

    Parameters
    ----------

    u : (M, K) array_like
        Left singular vectors (in columns)
    v : (K, N) array_like
        Right singular vectors (in rows)
    u_based_decision: bool
        Whether or not to choose signs based
        on `u` rather than `v`, by default False

    Returns
    -------

    u : (M, K) array_like
        Left singular vectors with corrected sign
    v:  (K, N) array_like
        Right singular vectors with corrected sign
    """
    # Determine half-space in which all singular vectors
    # lie relative to an arbitrary vector; summation
    # equivalent to dot product with row vector of ones
    if u_based_decision:
        dtype = u.dtype
        signs = np.sum(u, axis=0, keepdims=True)
    else:
        dtype = v.dtype
        signs = np.sum(v, axis=1, keepdims=True).T
    signs = 2.0 * ((signs >= 0) - 0.5).astype(dtype)
    # Force all singular vectors into same half-space
    u, v = u * signs, v * signs.T
    return u, v


def scipy_linalg_safe(func_name, *args, **kwargs):
    # need to evaluate at least the first input array
    # for gpu/cpu checking
    a = args[0]
    if is_cupy_type(a):
        import cupyx.scipy.linalg

        func = getattr(cupyx.scipy.linalg, func_name)
    else:
        import scipy.linalg

        func = getattr(scipy.linalg, func_name)

    return func(*args, **kwargs)


def solve_triangular_safe(a, b, lower=False):
    return scipy_linalg_safe("solve_triangular", a, b, lower=lower)


def __getattr__(name):
    # Can't use the @_deprecated decorator as it would not work on `except AxisError`
    if name == "AxisError":
        warnings.warn(
            "AxisError was deprecated after version 2021.10.0 and will be removed in a "
            f"future release. Please use {typename(AxisError)} instead.",
            category=FutureWarning,
            stacklevel=2,
        )
        return AxisError
    else:
        raise AttributeError(f"module {__name__} has no attribute {name}")
