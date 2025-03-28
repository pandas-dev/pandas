from __future__ import annotations

import operator

import numpy as np
import toolz

from dask._collections import new_collection
from dask.array._array_expr._blockwise import Blockwise, Elemwise
from dask.array._array_expr._expr import ArrayExpr
from dask.array._array_expr._io import FromGraph
from dask.array.core import finalize
from dask.base import DaskMethodsMixin, named_schedulers
from dask.core import flatten
from dask.utils import key_split


class Array(DaskMethodsMixin):
    __dask_scheduler__ = staticmethod(
        named_schedulers.get("threads", named_schedulers["sync"])
    )
    __dask_optimize__ = staticmethod(lambda dsk, keys, **kwargs: dsk)

    def __init__(self, expr):
        self._expr = expr

    @property
    def expr(self) -> ArrayExpr:
        return self._expr

    def __dask_postcompute__(self):
        return finalize, ()

    def __dask_postpersist__(self):
        state = self.expr.lower_completely()
        return from_graph, (
            state._meta,
            state.chunks,
            list(flatten(state.__dask_keys__())),
            key_split(state._name),
        )

    @property
    def dask(self):
        return self.__dask_graph__()

    def __dask_graph__(self):
        out = self.expr.lower_completely()
        return out.__dask_graph__()

    def __dask_keys__(self):
        out = self.expr.lower_completely()
        return out.__dask_keys__()

    def compute(self, **kwargs):
        return DaskMethodsMixin.compute(self.optimize(), **kwargs)

    def persist(self, **kwargs):
        return DaskMethodsMixin.persist(self.optimize(), **kwargs)

    def optimize(self):
        return new_collection(self.expr.optimize())

    def simplify(self):
        return new_collection(self.expr.simplify())

    @property
    def _meta(self):
        return self.expr._meta

    @property
    def dtype(self):
        return self.expr.dtype

    @property
    def shape(self):
        return self.expr.shape

    @property
    def chunks(self):
        return self.expr.chunks

    @property
    def ndim(self):
        return self.expr.ndim

    @property
    def name(self):
        return self.expr.name

    def __add__(self, other):
        return elemwise(operator.add, self, other)

    def __radd__(self, other):
        return elemwise(operator.add, other, self)

    def __mul__(self, other):
        return elemwise(operator.mul, self, other)

    def __rmul__(self, other):
        return elemwise(operator.mul, other, self)

    def __sub__(self, other):
        return elemwise(operator.sub, self, other)

    def __rsub__(self, other):
        return elemwise(operator.sub, other, self)

    def __pow__(self, other):
        return elemwise(operator.pow, self, other)

    def __rpow__(self, other):
        return elemwise(operator.pow, other, self)

    def __truediv__(self, other):
        return elemwise(operator.truediv, self, other)

    def __rtruediv__(self, other):
        return elemwise(operator.truediv, other, self)

    def __floordiv__(self, other):
        return elemwise(operator.floordiv, self, other)

    def __rfloordiv__(self, other):
        return elemwise(operator.floordiv, other, self)

    def __array_function__(self, func, types, args, kwargs):
        # TODO(expr-soon): Not done yet, but needed for assert_eq to identify us as an Array
        raise NotImplementedError


def from_graph(layer, _meta, chunks, keys, name_prefix):
    return new_collection(
        FromGraph(
            layer=layer,
            _meta=_meta,
            chunks=chunks,
            keys=keys,
            name_prefix=name_prefix,
        )
    )


def blockwise(
    func,
    out_ind,
    *args,
    name=None,
    token=None,
    dtype=None,
    adjust_chunks=None,
    new_axes=None,
    align_arrays=False,  # TODO(expr-soon): Change to True
    concatenate=None,
    meta=None,
    **kwargs,
):
    """Tensor operation: Generalized inner and outer products

    A broad class of blocked algorithms and patterns can be specified with a
    concise multi-index notation.  The ``blockwise`` function applies an in-memory
    function across multiple blocks of multiple inputs in a variety of ways.
    Many dask.array operations are special cases of blockwise including
    elementwise, broadcasting, reductions, tensordot, and transpose.

    Parameters
    ----------
    func : callable
        Function to apply to individual tuples of blocks
    out_ind : iterable
        Block pattern of the output, something like 'ijk' or (1, 2, 3)
    *args : sequence of Array, index pairs
        You may also pass literal arguments, accompanied by None index
        e.g. (x, 'ij', y, 'jk', z, 'i', some_literal, None)
    **kwargs : dict
        Extra keyword arguments to pass to function
    dtype : np.dtype
        Datatype of resulting array.
    concatenate : bool, keyword only
        If true concatenate arrays along dummy indices, else provide lists
    adjust_chunks : dict
        Dictionary mapping index to function to be applied to chunk sizes
    new_axes : dict, keyword only
        New indexes and their dimension lengths
    align_arrays: bool
        Whether or not to align chunks along equally sized dimensions when
        multiple arrays are provided.  This allows for larger chunks in some
        arrays to be broken into smaller ones that match chunk sizes in other
        arrays such that they are compatible for block function mapping. If
        this is false, then an error will be thrown if arrays do not already
        have the same number of blocks in each dimension.

    Examples
    --------
    2D embarrassingly parallel operation from two arrays, x, and y.

    >>> import operator, numpy as np, dask.array as da
    >>> x = da.from_array([[1, 2],
    ...                    [3, 4]], chunks=(1, 2))
    >>> y = da.from_array([[10, 20],
    ...                    [0, 0]])
    >>> z = blockwise(operator.add, 'ij', x, 'ij', y, 'ij', dtype='f8')
    >>> z.compute()
    array([[11, 22],
           [ 3,  4]])

    Outer product multiplying a by b, two 1-d vectors

    >>> a = da.from_array([0, 1, 2], chunks=1)
    >>> b = da.from_array([10, 50, 100], chunks=1)
    >>> z = blockwise(np.outer, 'ij', a, 'i', b, 'j', dtype='f8')
    >>> z.compute()
    array([[  0,   0,   0],
           [ 10,  50, 100],
           [ 20, 100, 200]])

    z = x.T

    >>> z = blockwise(np.transpose, 'ji', x, 'ij', dtype=x.dtype)
    >>> z.compute()
    array([[1, 3],
           [2, 4]])

    The transpose case above is illustrative because it does transposition
    both on each in-memory block by calling ``np.transpose`` and on the order
    of the blocks themselves, by switching the order of the index ``ij -> ji``.

    We can compose these same patterns with more variables and more complex
    in-memory functions

    z = X + Y.T

    >>> z = blockwise(lambda x, y: x + y.T, 'ij', x, 'ij', y, 'ji', dtype='f8')
    >>> z.compute()
    array([[11,  2],
           [23,  4]])

    Any index, like ``i`` missing from the output index is interpreted as a
    contraction (note that this differs from Einstein convention; repeated
    indices do not imply contraction.)  In the case of a contraction the passed
    function should expect an iterable of blocks on any array that holds that
    index.  To receive arrays concatenated along contracted dimensions instead
    pass ``concatenate=True``.

    Inner product multiplying a by b, two 1-d vectors

    >>> def sequence_dot(a_blocks, b_blocks):
    ...     result = 0
    ...     for a, b in zip(a_blocks, b_blocks):
    ...         result += a.dot(b)
    ...     return result

    >>> z = blockwise(sequence_dot, '', a, 'i', b, 'i', dtype='f8')
    >>> z.compute()
    np.int64(250)

    Add new single-chunk dimensions with the ``new_axes=`` keyword, including
    the length of the new dimension.  New dimensions will always be in a single
    chunk.

    >>> def f(a):
    ...     return a[:, None] * np.ones((1, 5))

    >>> z = blockwise(f, 'az', a, 'a', new_axes={'z': 5}, dtype=a.dtype)

    New dimensions can also be multi-chunk by specifying a tuple of chunk
    sizes.  This has limited utility as is (because the chunks are all the
    same), but the resulting graph can be modified to achieve more useful
    results (see ``da.map_blocks``).

    >>> z = blockwise(f, 'az', a, 'a', new_axes={'z': (5, 5)}, dtype=x.dtype)
    >>> z.chunks
    ((1, 1, 1), (5, 5))

    If the applied function changes the size of each chunk you can specify this
    with a ``adjust_chunks={...}`` dictionary holding a function for each index
    that modifies the dimension size in that index.

    >>> def double(x):
    ...     return np.concatenate([x, x])

    >>> y = blockwise(double, 'ij', x, 'ij',
    ...               adjust_chunks={'i': lambda n: 2 * n}, dtype=x.dtype)
    >>> y.chunks
    ((2, 2), (2,))

    Include literals by indexing with None

    >>> z = blockwise(operator.add, 'ij', x, 'ij', 1234, None, dtype=x.dtype)
    >>> z.compute()
    array([[1235, 1236],
           [1237, 1238]])
    """
    assert not align_arrays  # TODO(expr-soon): Remove
    new_axes = new_axes or {}

    # Input Validation
    if len(set(out_ind)) != len(out_ind):
        raise ValueError(
            "Repeated elements not allowed in output index",
            [k for k, v in toolz.frequencies(out_ind).items() if v > 1],
        )
    new = (
        set(out_ind)
        - {a for arg in args[1::2] if arg is not None for a in arg}
        - set(new_axes or ())
    )
    if new:
        raise ValueError("Unknown dimension", new)

    return new_collection(
        Blockwise(
            func,
            out_ind,
            name,
            token,
            dtype,
            adjust_chunks,
            new_axes,
            align_arrays,
            concatenate,
            meta,
            kwargs,
            *args,
        )
    )


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
    if where is not True:
        # TODO(expr-soon): Need asarray for this
        where = True

    if out is not None:
        raise NotImplementedError("elemwise does not support out=")

    args = [np.asarray(a) if isinstance(a, (list, tuple)) else a for a in args]

    return new_collection(Elemwise(op, dtype, name, where, *args))
