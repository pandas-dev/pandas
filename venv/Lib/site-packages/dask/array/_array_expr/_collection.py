from __future__ import annotations

import operator
import warnings
from collections.abc import Iterable

import numpy as np
import toolz
from toolz import concat, first

from dask._collections import new_collection
from dask.array import chunk
from dask.array._array_expr._blockwise import Blockwise, Elemwise, Transpose
from dask.array._array_expr._expr import (
    ArrayExpr,
    Concatenate,
    Stack,
    unify_chunks_expr,
)
from dask.array._array_expr._io import FromArray, FromGraph
from dask.array.core import T_IntOrNaN, finalize, getter_inline
from dask.array.dispatch import concatenate_lookup
from dask.array.utils import meta_from_array
from dask.base import DaskMethodsMixin, is_dask_collection, named_schedulers
from dask.core import flatten
from dask.utils import derived_from, is_arraylike, key_split


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

    @property
    def _name(self):
        return self.expr._name

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

    def __dask_tokenize__(self):
        return "Array", self.expr._name

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
    def numblocks(self):
        return self.expr.numblocks

    @property
    def size(self) -> T_IntOrNaN:
        return self.expr.size

    @property
    def name(self):
        return self.expr.name

    def __len__(self):
        return self.expr.__len__()

    def __getitem__(self, index):
        # Field access, e.g. x['a'] or x[['a', 'b']]
        if isinstance(index, str) or (
            isinstance(index, list) and index and all(isinstance(i, str) for i in index)
        ):
            # TODO(expr-soon): needs map_blocks that we don't support yet,
            #  but implementation is trivial after we have that
            raise NotImplementedError()

        if not isinstance(index, tuple):
            index = (index,)

        from dask.array._array_expr._slicing import (
            slice_array,
            slice_with_int_dask_array,
        )
        from dask.array.slicing import normalize_index

        index2 = normalize_index(index, self.shape)
        dependencies = {self.name}
        for i in index2:
            if isinstance(i, Array):
                dependencies.add(i.name)

        if any(isinstance(i, Array) and i.dtype.kind in "iu" for i in index2):
            self, index2 = slice_with_int_dask_array(self, index2)
        if any(isinstance(i, Array) and i.dtype == bool for i in index2):
            # TODO(expr-soon): This is simple but needs ravel which needs reshape,
            # which is not simple. Trivial to add after we have reshape
            raise NotImplementedError

        if all(isinstance(i, slice) and i == slice(None) for i in index2):
            return self

        result = slice_array(self.expr, index2)
        return new_collection(result)

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

    def transpose(self, axes=None):
        if axes:
            if len(axes) != self.ndim:
                raise ValueError("axes don't match array")
            axes = tuple(d + self.ndim if d < 0 else d for d in axes)
        else:
            axes = tuple(range(self.ndim))[::-1]

        return new_collection(Transpose(self, axes))

    @property
    def T(self):
        return self.transpose()

    def rechunk(
        self,
        chunks="auto",
        threshold=None,
        block_size_limit=None,
        balance=False,
        method=None,
    ):
        return rechunk(self, chunks, threshold, block_size_limit, balance, method)

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

    def map_blocks(self, func, *args, **kwargs):
        from dask.array._array_expr._map_blocks import map_blocks

        return map_blocks(func, self, *args, **kwargs)


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
    align_arrays=True,
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


def rechunk(
    x,
    chunks="auto",
    threshold=None,
    block_size_limit=None,
    balance=False,
    method=None,
):
    """
    Convert blocks in dask array x for new chunks.

    Parameters
    ----------
    x: dask array
        Array to be rechunked.
    chunks:  int, tuple, dict or str, optional
        The new block dimensions to create. -1 indicates the full size of the
        corresponding dimension. Default is "auto" which automatically
        determines chunk sizes.
    threshold: int, optional
        The graph growth factor under which we don't bother introducing an
        intermediate step.
    block_size_limit: int, optional
        The maximum block size (in bytes) we want to produce
        Defaults to the configuration value ``array.chunk-size``
    balance : bool, default False
        If True, try to make each chunk to be the same size.

        This means ``balance=True`` will remove any small leftover chunks, so
        using ``x.rechunk(chunks=len(x) // N, balance=True)``
        will almost certainly result in ``N`` chunks.
    method: {'tasks', 'p2p'}, optional.
        Rechunking method to use.


    Examples
    --------
    >>> import dask.array as da
    >>> x = da.ones((1000, 1000), chunks=(100, 100))

    Specify uniform chunk sizes with a tuple

    >>> y = x.rechunk((1000, 10))

    Or chunk only specific dimensions with a dictionary

    >>> y = x.rechunk({0: 1000})

    Use the value ``-1`` to specify that you want a single chunk along a
    dimension or the value ``"auto"`` to specify that dask can freely rechunk a
    dimension to attain blocks of a uniform block size

    >>> y = x.rechunk({0: -1, 1: 'auto'}, block_size_limit=1e8)

    If a chunk size does not divide the dimension then rechunk will leave any
    unevenness to the last chunk.

    >>> x.rechunk(chunks=(400, -1)).chunks
    ((400, 400, 200), (1000,))

    However if you want more balanced chunks, and don't mind Dask choosing a
    different chunksize for you then you can use the ``balance=True`` option.

    >>> x.rechunk(chunks=(400, -1), balance=True).chunks
    ((500, 500), (1000,))
    """

    return new_collection(
        x.expr.rechunk(chunks, threshold, block_size_limit, balance, method)
    )


def from_array(
    x,
    chunks="auto",
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
    elif is_dask_collection(x):
        warnings.warn(
            "Passing an object to dask.array.from_array which is already a "
            "Dask collection. This can lead to unexpected behavior."
        )

    if isinstance(x, (list, tuple, memoryview) + np.ScalarType):
        x = np.array(x)

    if is_arraylike(x) and hasattr(x, "copy"):
        x = x.copy()

    return new_collection(
        FromArray(
            x,
            chunks,
            lock=lock,
            asarray=asarray,
            fancy=fancy,
            getitem=getitem,
            meta=meta,
            inline_array=inline_array,
        )
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
            return a.map_blocks(np.asarray, like=like_meta, dtype=dtype, order=order)
        else:
            a = np.asarray(a, like=like_meta, dtype=dtype, order=order)

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
            return a.map_blocks(np.asanyarray, like=like_meta, dtype=dtype, order=order)
        else:
            a = np.asanyarray(a, like=like_meta, dtype=dtype, order=order)

    a = from_array(
        a,
        chunks=a.shape,
        getitem=getter_inline,
        asarray=False,
        inline_array=inline_array,
    )
    return _as_dtype(a, dtype)


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
    uc_args = list(concat((x.expr, ind) for x in seq2))
    _, seq2, _ = unify_chunks_expr(*uc_args)

    assert len({a.chunks for a in seq2}) == 1  # same chunks

    return new_collection(Stack(seq2[0], axis, meta, *seq2[1:]))


@derived_from(np)
def array(x, dtype=None, ndmin=None, *, like=None):
    x = asarray(x, like=like)
    while ndmin is not None and x.ndim < ndmin:
        x = x[None, :]
    if dtype is not None and x.dtype != dtype:
        x = x.astype(dtype)
    return x


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

    seq_tmp = [s.expr for s in seq2]
    uc_args = list(concat((s, i) for s, i in zip(seq_tmp, inds)))
    _, seq2, _ = unify_chunks_expr(*uc_args)
    return new_collection(Concatenate(seq2[0], axis, meta, *seq2[1:]))
