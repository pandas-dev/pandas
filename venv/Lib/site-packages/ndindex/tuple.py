import itertools

from .ndindex import NDIndexCommon, ndindex
from .subindex_helpers import subindex_slice
from .shapetools import asshape, broadcast_shapes

from ._tuple import _Tuple

class Tuple(_Tuple, NDIndexCommon):
    """
    Represents a tuple of single-axis indices.

    Valid single axis indices are

    - :class:`Integer`
    - :class:`Slice`
    - :class:`ellipsis`
    - :class:`Newaxis`
    - :class:`IntegerArray`
    - :class:`BooleanArray`

    See :doc:`../indexing-guide/multidimensional-indices/tuples` for a
    description of the semantics of tuple indices.

    `Tuple(x1, x2, …, xn)` represents the index `a[x1, x2, …, xn]` or,
    equivalently, `a[(x1, x2, …, xn)]`. `Tuple()` with no arguments is the
    empty tuple index, `a[()]`, which returns `a` unchanged.

    >>> from ndindex import Tuple, Slice
    >>> import numpy as np
    >>> idx = Tuple(0, Slice(2, 4))
    >>> a = np.arange(10).reshape((2, 5))
    >>> a
    array([[0, 1, 2, 3, 4],
           [5, 6, 7, 8, 9]])
    >>> a[0, 2:4]
    array([2, 3])
    >>> a[idx.raw]
    array([2, 3])

    .. note::

       `Tuple` does *not* represent a tuple, but rather an *tuple index*. It
       does not have most methods that `tuple` has, and should not be used in
       non-indexing contexts. See the document on :any:`type-confusion` for
       more details.

    """
    __slots__ = ()

    def __hash__(self):
        # Since self.args is itself a tuple, it will match the hash of
        # self.raw when it is hashable.
        return hash(self.args)

    def __repr__(self):
        # Since tuples are nested, we can print the raw form of the args to
        # make them a little more readable.
        def _repr(s):
            if s == ...:
                return '...'
            if isinstance(s, ArrayIndex):
                if s.shape and 0 not in s.shape:
                    return repr(s.array.tolist())
                return repr(s)
            return repr(s.raw)
        return f"{self.__class__.__name__}({', '.join(map(_repr, self.args))})"

    def __str__(self):
        # Since tuples are nested, we can print the raw form of the args to
        # make them a little more readable.
        def _str(s):
            if s == ...:
                return '...'
            if isinstance(s, ArrayIndex):
                return str(s)
            return str(s.raw)
        return f"{self.__class__.__name__}({', '.join(map(_str, self.args))})"

    @property
    def has_ellipsis(self):
        """
        Returns True if self has an ellipsis
        """
        return ... in self.args

    @property
    def ellipsis_index(self):
        """
        Give the index i of `self.args` where the ellipsis is.

        If `self` doesn't have an ellipsis, it gives `len(self.args)`, since
        tuple indices without an ellipsis always implicitly end in an
        ellipsis.

        The resulting value `i` is such that `self.args[:i]` indexes the
        beginning axes of an array and `self.args[i+1:]` indexes the end axes
        of an array.

        >>> from ndindex import Tuple
        >>> idx = Tuple(0, 1, ..., 2, 3)
        >>> i = idx.ellipsis_index
        >>> i
        2
        >>> idx.args[:i]
        (Integer(0), Integer(1))
        >>> idx.args[i+1:]
        (Integer(2), Integer(3))

        >>> Tuple(0, 1).ellipsis_index
        2

        """
        if self.has_ellipsis:
            return self.args.index(...)
        return len(self.args)

    def reduce(self, shape=None, *, negative_int=False):
        r"""
        Reduce a Tuple index on an array of shape `shape`

        A `Tuple` with a single argument is always reduced to that single
        argument (because `a[idx,]` is the same as `a[idx]`).

        >>> from ndindex import Tuple

        >>> Tuple(slice(2, 4)).reduce()
        Slice(2, 4, 1)

        If an explicit array shape is given, the result will either be
        `IndexError` if the index is invalid for the given shape, or an index
        that is as simple as possible:

        - All the elements of the :any:`Tuple` are recursively :any:`reduced
          <NDIndex.reduce>`.

        - Any axes that can be merged into an :any:`ellipsis` are removed.
          This includes the implicit ellipsis at the end of a Tuple that
          doesn't contain any explicit ellipses.

        - :any:`Ellipses <ellipsis>` that don't match any axes are removed.

        - An :any:`ellipsis` at the end of the :any:`Tuple` is removed.

        - Scalar :any:`BooleanArray` arguments (`True` or `False`) are
          combined into a single term (the first boolean scalar is replaced
          with the AND of all the boolean scalars).

        - If the resulting :any:`Tuple` would have a single argument, that
          argument is returned.

        >>> idx = Tuple(0, ..., slice(0, 3))
        >>> idx.reduce((5, 4))
        Tuple(0, slice(0, 3, 1))
        >>> idx.reduce((5, 3))
        Integer(0)

        >>> idx = Tuple(slice(0, 10), -3)
        >>> idx.reduce((5,))
        Traceback (most recent call last):
        ...
        IndexError: too many indices for array: array is 1-dimensional, but 2 were indexed
        >>> idx.reduce((5, 2))
        Traceback (most recent call last):
        ...
        IndexError: index -3 is out of bounds for axis 1 with size 2

        Note
        ====

        ndindex presently does not distinguish between scalar objects and
        0-D arrays. It is possible for the original index to produce one
        and the reduced index to produce the other. In particular, the
        presence of a redundant ellipsis forces NumPy to return a 0-D array
        instead of a scalar.

        >>> import numpy as np
        >>> a = np.array([0, 1])
        >>> Tuple(..., 1).reduce(a.shape)
        Integer(1)
        >>> a[..., 1]
        array(1)
        >>> a[1] # doctest: +SKIPNP1
        np.int64(1)

        See https://github.com/Quansight-Labs/ndindex/issues/22.

        See Also
        ========

        .Tuple.expand
        .NDIndex.reduce
        .Slice.reduce
        .Integer.reduce
        .ellipsis.reduce
        .Newaxis.reduce
        .IntegerArray.reduce
        .BooleanArray.reduce

        """
        args = list(self.args)

        boolean_scalars = [i for i in args if _is_boolean_scalar(i)]
        if len(boolean_scalars) > 1:
            _args = []
            seen_boolean_scalar = False
            for s in args:
                if _is_boolean_scalar(s):
                    if seen_boolean_scalar:
                        continue
                    _args.append(BooleanArray(all(i == True for i in boolean_scalars)))
                    seen_boolean_scalar = True
                else:
                    _args.append(s)
            return type(self)(*_args).reduce(shape, negative_int=negative_int)

        arrays = []
        for i in args:
            if _is_boolean_scalar(i):
                continue
            elif isinstance(i, IntegerArray):
                arrays.append(i.raw)
            elif isinstance(i, BooleanArray):
                # TODO: Avoid explicitly calling nonzero
                arrays.extend(i.raw.nonzero())
        if arrays:
            from numpy import broadcast_to

            broadcast_shape = broadcast_shapes(*[a.shape for a in arrays])
        else:
            broadcast_shape = ()

        # If the broadcast shape is empty, out of bounds indices in
        # non-empty arrays are ignored, e.g., ([], [10]) would broadcast to
        # ([], []), so the bounds for 10 are not checked. Thus, we must do
        # this before calling reduce() on the arguments. This rule, however,
        # is *not* followed for scalar integer indices.
        if 0 in broadcast_shape:
            for i in range(len(args)):
                s = args[i]
                if isinstance(s, IntegerArray):
                    if s.ndim == 0:
                        args[i] = Integer(s.raw)
                    else:
                        # broadcast_to(x) gives a readonly view on x, which is also
                        # readonly, so set _copy=False to avoid representing the full
                        # broadcasted array in memory.
                        args[i] = type(s)(broadcast_to(s.raw, broadcast_shape),
                                          _copy=False)

        if shape is not None:
            # assert self.args.count(...) == 1
            # assert self.args.count(False) <= 1
            # assert self.args.count(True) <= 1
            n_newaxis = self.args.count(None)
            n_boolean = sum(i.ndim - 1 for i in args if
                            isinstance(i, BooleanArray) and not _is_boolean_scalar(i))
            if True in args or False in args:
                n_boolean -= 1
            indexed_args = len(args) + n_boolean - n_newaxis - 1 # -1 for the

            shape = asshape(shape, axis=indexed_args - 1)

        ellipsis_i = self.ellipsis_index

        preargs = []
        removable = shape is not None
        begin_offset = args[:ellipsis_i].count(None)
        begin_offset -= sum(j.ndim - 1 for j in args[:ellipsis_i] if
                            isinstance(j, BooleanArray))
        for i, s in enumerate(reversed(args[:ellipsis_i]), start=1):
            if s == None:
                begin_offset -= 1
            elif isinstance(s, BooleanArray):
                begin_offset += s.ndim - 1
            axis = ellipsis_i - i - begin_offset
            reduced = s.reduce(shape, axis=axis, negative_int=negative_int)
            if (removable
                and isinstance(reduced, Slice)
                and reduced == Slice(0, shape[axis], 1)):
                continue
            else:
                removable = False
                preargs.insert(0, reduced)

        if shape is None:
            endargs = [s.reduce(negative_int=negative_int) for s in args[ellipsis_i+1:]]
        else:
            endargs = []
            end_offset = 0
            for i, s in enumerate(reversed(args[ellipsis_i+1:]), start=1):
                if isinstance(s, BooleanArray):
                    end_offset -= s.ndim - 1
                elif s == None:
                    end_offset += 1
                axis = len(shape) - i + end_offset
                if not (isinstance(s, IntegerArray) and (0 in broadcast_shape or
                                                         False in args)):
                    # Array bounds are not checked when the broadcast shape is empty
                    s = s.reduce(shape, axis=axis, negative_int=negative_int)
                endargs.insert(0, s)

        if shape is not None:
            # Remove redundant slices
            axis = len(shape) - len(endargs) + end_offset
            for i, s in enumerate(endargs):
                axis += i
                if (isinstance(s, Slice)
                    and s == Slice(0, shape[axis], 1)):
                    i += 1
                    continue
                else:
                    break
            if endargs:
                endargs = endargs[i:]

        if shape is None or (endargs and len(preargs) + len(endargs)
                             < len(shape) + args.count(None) - n_boolean):
            preargs = preargs + [...]

        newargs = preargs + endargs

        if newargs and newargs[-1] == ...:
            newargs = newargs[:-1]

        if len(newargs) == 1:
            return newargs[0]

        return type(self)(*newargs)

    def broadcast_arrays(self):
        args = self.args
        boolean_scalars = [i for i in args if _is_boolean_scalar(i)]
        if len(boolean_scalars) > 1:
            _args = []
            seen_boolean_scalar = False
            for s in args:
                if _is_boolean_scalar(s):
                    if seen_boolean_scalar:
                        continue
                    _args.append(BooleanArray(all(i == True for i in boolean_scalars)))
                    seen_boolean_scalar = True
                else:
                    _args.append(s)
            return type(self)(*_args).broadcast_arrays()

        # Broadcast all array indices. Note that broadcastability is checked
        # in the Tuple constructor, so this should not fail.
        boolean_nonzero = {}
        arrays = []
        for s in args:
            if _is_boolean_scalar(s):
                continue
            elif isinstance(s, IntegerArray):
                arrays.append(s.raw)
            elif isinstance(s, BooleanArray):
                nz = s.raw.nonzero()
                arrays.extend(nz)
                boolean_nonzero[s] = nz
        if not arrays:
            return self

        from numpy import array, broadcast_to, intp

        broadcast_shape = broadcast_shapes(*[a.shape for a in arrays])

        newargs = []
        for s in args:
            if isinstance(s, BooleanArray):
                if not _is_boolean_scalar(s):
                    newargs.extend([IntegerArray(broadcast_to(i, broadcast_shape))
                                    for i in boolean_nonzero[s]])
            elif isinstance(s, Integer):
                # broadcast_to(x) gives a readonly view on x, which is also
                # readonly, so set _copy=False to avoid representing the full
                # broadcasted array in memory.
                newargs.append(IntegerArray(broadcast_to(array(s.raw, dtype=intp),
                                              broadcast_shape), _copy=False))
            elif isinstance(s, IntegerArray):
                newargs.append(IntegerArray(broadcast_to(s.raw, broadcast_shape),
                                            _copy=False))
            else:
                newargs.append(s)
        return Tuple(*newargs)

    def expand(self, shape):
        # The expand() docstring is on NDIndex.expand()
        args = list(self.args)
        if ... not in args:
            return type(self)(*args, ...).expand(shape)

        # TODO: Use broadcast_arrays here. The challenge is that we still need
        # to do bounds checks on nonscalar integer arrays that get broadcast
        # away.
        boolean_scalars = [i for i in args if _is_boolean_scalar(i)]
        if len(boolean_scalars) > 1:
            _args = []
            seen_boolean_scalar = False
            for s in args:
                if _is_boolean_scalar(s):
                    if seen_boolean_scalar:
                        continue
                    _args.append(BooleanArray(all(i == True for i in boolean_scalars)))
                    seen_boolean_scalar = True
                else:
                    _args.append(s)
            return type(self)(*_args).expand(shape)

        # Broadcast all array indices. Note that broadcastability is checked
        # in the Tuple constructor, so this should not fail.
        arrays = []
        for i in args:
            if _is_boolean_scalar(i):
                continue
            elif isinstance(i, IntegerArray):
                arrays.append(i.raw)
            elif isinstance(i, BooleanArray):
                # TODO: Avoid calling nonzero twice
                arrays.extend(i.raw.nonzero())

        if arrays:
            from numpy import broadcast_to, array, intp

            broadcast_shape = broadcast_shapes(*[a.shape for a in arrays])
            # If the broadcast shape is empty, out of bounds indices in
            # non-empty arrays are ignored, e.g., ([], [10]) would broadcast to
            # ([], []), so the bounds for 10 are not checked. Thus, we must do
            # this before calling reduce() on the arguments. This rule, however,
            # is *not* followed for scalar integer indices.

            for i in range(len(args)):
                s = args[i]
                if isinstance(s, IntegerArray):
                    if s.ndim == 0:
                        args[i] = Integer(s.raw)
                    else:
                        # broadcast_to(x) gives a readonly view on x, which is also
                        # readonly, so set _copy=False to avoid representing the full
                        # broadcasted array in memory.
                        args[i] = type(s)(broadcast_to(s.raw, broadcast_shape),
                                          _copy=False)

        # assert args.count(...) == 1
        # assert args.count(False) <= 1
        # assert args.count(True) <= 1
        n_newaxis = args.count(None)
        n_boolean = sum(i.ndim - 1 for i in args if
                        isinstance(i, BooleanArray) and not _is_boolean_scalar(i))
        if True in args or False in args:
            n_boolean -= 1
        indexed_args = len(args) + n_boolean - n_newaxis - 1 # -1 for the ellipsis
        shape = asshape(shape, axis=indexed_args - 1)

        ellipsis_i = self.ellipsis_index

        startargs = []
        begin_offset = 0
        for i, s in enumerate(args[:ellipsis_i]):
            axis = i + begin_offset
            if not (isinstance(s, IntegerArray) and (0 in broadcast_shape or
                                                     False in args)):
                s = s.reduce(shape, axis=axis)
            if isinstance(s, ArrayIndex):
                if isinstance(s, BooleanArray):
                    begin_offset += s.ndim - 1
                    if not _is_boolean_scalar(s):
                        s = s.reduce(shape, axis=axis)
                        startargs.extend([IntegerArray(broadcast_to(i,
                                                                  broadcast_shape))
                                        for i in s.array.nonzero()])
                        continue
            elif arrays and isinstance(s, Integer):
                s = IntegerArray(broadcast_to(array(s.raw, dtype=intp),
                                              broadcast_shape), _copy=False)
            elif s == None:
                begin_offset -= 1
            startargs.append(s)

        # TODO: Merge this with the above loop
        endargs = []
        end_offset = 0
        for i, s in enumerate(reversed(args[ellipsis_i+1:]), start=1):
            if isinstance(s, ArrayIndex):
                if isinstance(s, BooleanArray):
                    end_offset -= s.ndim - 1
                    if not _is_boolean_scalar(s):
                        s = s.reduce(shape, axis=len(shape) - i + end_offset)
                        endargs.extend([IntegerArray(broadcast_to(i,
                                                                  broadcast_shape))
                                        for i in reversed(s.array.nonzero())])
                        continue
            elif arrays and isinstance(s, Integer):
                if (0 in broadcast_shape or False in args):
                    s = s.reduce(shape, axis=len(shape)-i+end_offset)
                s = IntegerArray(broadcast_to(array(s.raw, dtype=intp),
                                              broadcast_shape), _copy=False)
            elif s == None:
                end_offset += 1
            axis = len(shape) - i + end_offset
            assert axis >= 0
            if not (isinstance(s, IntegerArray) and (0 in broadcast_shape or
                                                     False in args)):
                # Array bounds are not checked when the broadcast shape is empty
                s = s.reduce(shape, axis=axis)
            endargs.append(s)

        idx_offset = begin_offset - end_offset

        midargs = [Slice(None).reduce(shape, axis=i + ellipsis_i + begin_offset) for
                        i in range(len(shape) - len(args) + 1 - idx_offset)]


        newargs = startargs + midargs + endargs[::-1]

        return type(self)(*newargs)

    def newshape(self, shape):
        # The docstring for this method is on the NDIndex base class
        shape = asshape(shape)

        if self == Tuple():
            return shape

        # This will raise any IndexErrors
        self = self.expand(shape)

        newshape = []
        axis = 0
        arrays = False
        for i, s in enumerate(self.args):
            if s == None:
                newshape.append(1)
                axis -= 1
            # After expand(), there will be at most one boolean scalar
            elif s == True:
                newshape.append(1)
                axis -= 1
            elif s == False:
                newshape.append(0)
                axis -= 1
            elif isinstance(s, ArrayIndex):
                if not arrays:
                    # Multiple arrays are all broadcast together (in expand())
                    # and iterated as one, so we only need to get the shape
                    # for the first array we see. Note that arrays separated
                    # by ellipses, slices, or newaxes affect the shape
                    # differently, but these are currently unsupported (see
                    # the comments in the Tuple constructor).

                    # expand() should remove all non scalar boolean arrays
                    assert not isinstance(s, BooleanArray)

                    newshape.extend(list(s.newshape(shape[axis])))
                    arrays = True
            else:
                newshape.extend(list(s.newshape(shape[axis])))
            axis += 1
        return tuple(newshape)

    def as_subindex(self, index):
        index = ndindex(index).reduce().broadcast_arrays()

        self = self.broadcast_arrays()

        if ... in self.args:
            raise NotImplementedError("Tuple.as_subindex() is not yet implemented for tuples with ellipses")

        if isinstance(index, (Integer, ArrayIndex, Slice)):
            index = Tuple(index)
        if isinstance(index, Tuple):
            new_args = []
            boolean_arrays = []
            integer_arrays = []
            if any(isinstance(i, Slice) and i.step < 0 for i in index.args):
                    raise NotImplementedError("Tuple.as_subindex() is only implemented on slices with positive steps")
            if ... in index.args:
                raise NotImplementedError("Tuple.as_subindex() is not yet implemented for tuples with ellipses")
            for self_arg, index_arg in zip(self.args, index.args):
                if (isinstance(self_arg, IntegerArray) and
                    isinstance(index_arg, Slice)):
                    if (self_arg.array < 0).any():
                        raise NotImplementedError("IntegerArray.as_subindex() is only implemented for arrays with all nonnegative entries. Try calling reduce() with a shape first.")
                    if index_arg.step < 0:
                        raise NotImplementedError("IntegerArray.as_subindex(Slice) is only implemented for slices with positive steps")

                    # After reducing, start is not None when step > 0
                    if index_arg.stop is None or index_arg.start < 0 or index_arg.stop < 0:
                        raise NotImplementedError("IntegerArray.as_subindex(Slice) is only implemented for slices with nonnegative start and stop. Try calling reduce() with a shape first.")

                    s = self_arg.array
                    start, stop, step = subindex_slice(
                        s, s+1, 1, index_arg.start, index_arg.stop, index_arg.step)
                    if (stop <= 0).all():
                        raise ValueError("Indices do not intersect")
                    if start.shape == ():
                        if start >= stop:
                            raise ValueError("Indices do not intersect")

                    integer_arrays.append((start, stop))
                    # Placeholder. We need to mask out the stops below.
                    new_args.append(IntegerArray(start))
                else:
                    subindex = self_arg.as_subindex(index_arg)
                    if isinstance(subindex, Tuple):
                        assert subindex == ()
                        subindex # Workaround https://github.com/nedbat/coveragepy/issues/1029
                        continue
                    if isinstance(subindex, BooleanArray):
                        boolean_arrays.append(subindex)
                    new_args.append(subindex)
            args_remainder = self.args[min(len(self.args), len(index.args)):]
            index_remainder = index.args[min(len(self.args), len(index.args)):]
            if any(isinstance(i, ArrayIndex) and i.isempty() for i in
                   index_remainder):
                raise ValueError("Indices do not intersect")
            for arg in args_remainder:
                if isinstance(arg, BooleanArray):
                    boolean_arrays.append(arg)
                if isinstance(arg, IntegerArray):
                    integer_arrays.append((arg.array, arg.array+1))
                new_args.append(arg)
            # Replace all boolean arrays with the logical AND of them.
            if any(i.isempty() for i in boolean_arrays):
                raise ValueError("Indices do not intersect")
            if boolean_arrays:
                if len(boolean_arrays) > 1:
                    from numpy import logical_and
                    new_array = BooleanArray(logical_and.reduce([i.array for i in boolean_arrays]))
                else:
                    new_array = boolean_arrays[0]
                new_args2 = []
                first = True
                for arg in new_args:
                    if arg in boolean_arrays:
                        if first:
                            new_args2.append(new_array)
                            first = False
                    else:
                        new_args2.append(arg)
                new_args = new_args2

            # Mask out integer arrays to only where the start is less than the
            # stop for all arrays.
            if integer_arrays:
                from numpy import array, broadcast_arrays, logical_and
                starts, stops = zip(*integer_arrays)
                starts = array(broadcast_arrays(*starts))
                stops = array(broadcast_arrays(*stops))
                mask = logical_and.reduce(starts < stops, axis=0)
                new_args2 = []
                i = 0
                for arg in new_args:
                    if isinstance(arg, IntegerArray):
                        if mask.ndim == 0:
                            # Integer arrays always result in a 1 dimensional
                            # result, except when we have a scalar, we want to
                            # have a 0 dimensional result to match Integer().
                            new_args2.append(IntegerArray(starts[i]))
                        elif mask.all():
                            new_args2.append(IntegerArray(starts[i]))
                        else:
                            new_args2.append(IntegerArray(starts[i, mask]))
                        if new_args2[-1].isempty():
                            raise ValueError("Indices do not intersect")
                        i += 1
                    else:
                        new_args2.append(arg)
                new_args = new_args2
            return Tuple(*new_args)
        raise NotImplementedError(f"Tuple.as_subindex() is not implemented for type '{type(index).__name__}")

    def isempty(self, shape=None):
        if shape is not None:
            return 0 in self.newshape(shape)

        return any(i.isempty() for i in self.args)

    def selected_indices(self, shape):
        shape = asshape(shape)
        idx = self.expand(shape)

        def _zipped_array_indices(array_indices, shape, axis=0):
            return zip(*[i.selected_indices(shape, axis=axis+j)
                            for j, i in enumerate(array_indices)])

        def _flatten(l):
            for element in l:
                if isinstance(element, tuple):
                    yield from element
                else:
                    yield element

        # We need to zip all array indices into a single iterator.
        iterators = []
        array_indices = []
        axis = 0
        for i in idx.args:
            if i == False:
                return
            elif i == True:
                pass
            elif isinstance(i, IntegerArray):
                array_indices.append(i)
            else:
                # Tuples do not support array indices separated by slices,
                # newaxes, or ellipses. Furthermore, if there are (non-scalar
                # boolean) array indices, any Integer and BooleanArray indices
                # are converted to IntegerArray. So we can assume all array
                # indices are together in a single block, and this is the end
                # of it.
                if array_indices:
                    iterators.append(_zipped_array_indices(array_indices,
                                                           shape, axis=axis))
                    axis += len(array_indices)
                    array_indices.clear()
                if i != None:
                    iterators.append(i.selected_indices(shape, axis=axis))
                    axis += 1
        if idx.args and isinstance(idx.args[-1], IntegerArray):
            iterators.append(_zipped_array_indices(array_indices,
                                                   shape, axis=axis))

        for i in itertools.product(*iterators):
            yield Tuple(*_flatten(i)).reduce()

# Imports at the bottom to avoid circular import issues
from .array import ArrayIndex
from .slice import Slice
from .integer import Integer
from .booleanarray import BooleanArray, _is_boolean_scalar
from .integerarray import IntegerArray
