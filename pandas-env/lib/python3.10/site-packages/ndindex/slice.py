from .ndindex import NDIndexCommon
from .subindex_helpers import subindex_slice
from .shapetools import asshape
from ._slice import _Slice

class default:
    """
    A default keyword argument value.

    Used as the default value for keyword arguments where `None` is also a
    meaningful value but not the default.

    """
    pass

class Slice(_Slice, NDIndexCommon):
    """
    Represents a slice on an axis of an nd-array.

    `Slice(x)` with one argument is equivalent to `Slice(None, x)`. `Slice(x,
    y)` with two arguments is equivalent to `Slice(x, y, None)`.

    `start` and `stop` can be any integer, or `None`. `step` can be any
    nonzero integer or `None`.

    `Slice(a, b)` is the same as the syntax `a:b` in an index and `Slice(a, b,
    c)` is the same as `a:b:c`. An argument being `None` is equivalent to the
    syntax where the item is omitted, for example, `Slice(None, None, k)` is
    the same as the syntax `::k`.

    `Slice.args` always has three arguments, and does not make any distinction
    between, for instance, `Slice(x, y)` and `Slice(x, y, None)`. This is
    because Python itself does not make the distinction between `x:y` and
    `x:y:` syntactically.

    See :doc:`../indexing-guide/slices` for a complete description of the
    semantics of slices.

    Slice has attributes `start`, `stop`, and `step` to access the
    corresponding attributes.

    >>> from ndindex import Slice
    >>> s = Slice(10)
    >>> s
    Slice(None, 10, None)
    >>> print(s.start)
    None
    >>> s.args
    (None, 10, None)
    >>> s.raw
    slice(None, 10, None)

    For most use cases, it's more convenient to create Slice objects using
    `ndindex[slice]`, which allows using `a:b` slicing syntax:

    >>> from ndindex import ndindex
    >>> ndindex[0:10]
    Slice(0, 10, None)

    """
    __slots__ = ()

    def __repr__(self):
        return f"{self.__class__.__name__}({', '.join(map(repr, self.args))})"

    def __hash__(self):
        # Slices are only hashable in Python 3.12+
        try:
            return hash(self.raw)
        except TypeError: # pragma: no cover
            return hash(self.args)

    def __len__(self):
        """
        `len()` gives the maximum size of an axis sliced with `self`.

        An actual array may produce a smaller size if it is smaller than the
        bounds of the slice. For instance, `[0, 1, 2][2:4]` only has 1 element
        but the maximum length of the slice `2:4` is 2.

        >>> from ndindex import Slice
        >>> [0, 1, 2][2:4]
        [2]
        >>> len(Slice(2, 4))
        2
        >>> [0, 1, 2, 3][2:4]
        [2, 3]

        If there is no such maximum, it raises `ValueError`.

        >>> # From the second element to the end, which could have any size
        >>> len(Slice(1, None))
        Traceback (most recent call last):
        ...
        ValueError: Cannot determine max length of slice

        The :meth:`Slice.reduce` method with a `shape` argument returns a
        `Slice` that always has a correct `len` which doesn't raise
        `ValueError`.

        >>> Slice(2, 4).reduce(3)
        Slice(2, 3, 1)
        >>> len(_)
        1

        Be aware that `len(Slice)` only gives the size of the axis being
        sliced. It does not say anything about the total shape of the array.
        In particular, the array may be empty after slicing if one of its
        dimensions is 0, but the other dimensions may be nonzero. To check if
        an array will empty after indexing, use :meth:`isempty`.

        See Also
        ========
        isempty

        """
        s = self
        if None in self.args or self.start < 0 or self.stop < 0:
            s = s.reduce()
        start, stop, step = s.args
        error = ValueError("Cannot determine max length of slice")
        # We reuse the logic in range.__len__. However, it is only correct if
        # start and stop are nonnegative.
        if step > 0:
            # start cannot be None
            if stop is None:
                if start >= 0:
                    # a[n:]. Extends to the end of the array.
                    raise error
                else:
                    # a[-n:]. From n from the end to the end. Same as
                    # range(-n, 0).
                    stop = 0
            elif start < 0 and stop >= 0:
                # a[-n:m] indexes from nth element from the end to the
                # m-1th element from the beginning.
                start, stop = 0, min(-start, stop)
            elif start >=0 and stop < 0:
                # a[n:-m]. The max length depends on the size of the array.
                raise error
        else:
            if stop is None:
                if start >= 0:
                    # a[n::-1] (start != None by above). Same as range(n, -1, -1)
                    stop = -1
                else:
                    # a[-n::-1]. From n from the end to the beginning of the
                    # array backwards. The max length depends on the size of
                    # the array.
                    raise error
            elif start < 0 and stop >= 0:
                # a[-n:m:-1]. The max length depends on the size of the array
                raise error
            elif start >=0 and stop < 0:
                # a[n:-m:-1] indexes from the nth element backwards to the mth
                # element from the end.
                start, stop = 0, min(start+1, -stop - 1)
                step = -step

        return len(range(start, stop, step))

    def reduce(self, shape=None, *, axis=0, negative_int=False):
        """
        `Slice.reduce` returns a slice that is canonicalized for an array of the
        given shape, or for any shape if `shape` is `None` (the default).

        `Slice.reduce` is a perfect canonicalization, meaning that two slices
        are equal---for all array shapes if `shape=None` or for arrays of
        shape `shape` otherwise---if and only if they `reduce` to the same
        `Slice` object. Note that ndindex objects do not simplify
        automatically, and `==` only does exact equality comparison, so to
        test that two slices are equal, use `slice1.reduce(shape) ==
        slice2.reduce(shape)`.

        - If `shape` is `None`, the following properties hold after calling
          `reduce()`:

          - `start` is not `None`.

          - `stop` is not `None`, when possible. The reduced `stop` can only
            be `None` if the original `stop` is.

          - `step` is not `None`.

          - `step` is as close to 0 as possible.

          - If the slice is always empty, the resulting slice will be
            `Slice(0, 0, 1)`. However, one should prefer the :any:`isempty`
            method to test if a slice is always empty.

          In particular, `stop` may be `None`, even after canonicalization
          with `reduce()` with no `shape`. This is because some slices are
          impossible to represent without `None` without making assumptions
          about the array shape. For example, `Slice(0, None)` cannot be
          equivalent to a slice with `stop != None` for all array shapes. To
          get a slice where the `start`, `stop`, and `step` are always
          integers, use `reduce(shape)` with an explicit array shape.

          Note that `Slice` objects that index a single element are not
          canonicalized to `Integer`, because integer indices always remove an
          axis whereas slices keep the axis. Furthermore, slices cannot raise
          `IndexError` except on arrays with shape equal to `()`.

          >>> from ndindex import Slice
          >>> Slice(10).reduce()
          Slice(0, 10, 1)
          >>> Slice(1, 3, 3).reduce()
          Slice(1, 2, 1)

        - If an explicit shape is given, the following properties are true
          after calling `Slice.reduce(shape)`:

          - `start`, `stop`, and `step` are not `None`,

          - `start` is nonnegative.

          - `stop` is nonnegative whenever possible. In particular, `stop` is
            only negative when it has to be to represent the given slice,
            i.e., a slice with negative `step` that indexes more than 1
            element and indexes the first (index `0`) element (in this case,
            it will be `-n - 1` where `n` is the size of the axis being
            sliced).

          - `stop` is as small as possible for positive `step` or large as
            possible for negative `step`.

          - `step` is as close to 0 as possible.

          - If the slice is empty for the given shape, the resulting slice
            will be `Slice(0, 0, 1)`. However, one should prefer the
            :any:`isempty` method to test if a slice is always empty.

          - If the slice indexes a single element, the resulting slice will be
            of the form `Slice(i, i+1, 1)`. However, one should prefer using
            `len(s.reduce(shape)) == 1` to test if a slice indexes exactly 1
            element.

          - :any:`len() <Slice.__len__>` gives the true size of the axis for a
            sliced array of the given shape, and never raises `ValueError`.

          The `axis` argument can be used to specify an axis of the shape (by
          default, `axis=0`). For convenience, `shape` can be passed as an integer
          for a single dimension.


          >>> from ndindex import Slice
          >>> Slice(1, 10).reduce(3)
          Slice(1, 3, 1)
          >>> Slice(-1, 1, -2).reduce(4)
          Slice(3, 4, 1)
          >>> Slice(1, 10, 3).reduce((4, 5), axis=0)
          Slice(1, 2, 1)
          >>> Slice(1, 10, 3).reduce((4, 5), axis=1)
          Slice(1, 5, 3)

          >>> s = Slice(2, None)
          >>> len(s)
          Traceback (most recent call last):
          ...
          ValueError: Cannot determine max length of slice
          >>> s.reduce((5,))
          Slice(2, 5, 1)
          >>> len(_)
          3

        See Also
        ========

        .NDIndex.reduce
        .Tuple.reduce
        .Integer.reduce
        .ellipsis.reduce
        .Newaxis.reduce
        .IntegerArray.reduce
        .BooleanArray.reduce

        """
        if self._reduced and shape is None:
            return self

        start, stop, step = self.args

        # Canonicalize with no shape

        if step is None:
            step = 1
        if start is None:
            if step > 0:
                start = 0
            else: # step < 0
                start = -1

        if start is not None and stop is not None:
            if start >= 0 and stop >= 0 or start < 0 and stop < 0:
                if step > 0:
                    if stop <= start:
                        start, stop, step = 0, 0, 1
                    elif start >= 0 and start + step >= stop:
                        # Indexes 1 element. Start has to be >= 0 because a
                        # negative start could be less than the size of the
                        # axis, in which case it will clip and the single
                        # element will be element 0. We can only do that
                        # reduction if we know the shape.

                        # Note that returning Integer here is wrong, because
                        # slices keep the axis and integers remove it.
                        stop, step = start + 1, 1
                    elif start < 0 and start + step > stop:
                        # The exception is this case where stop is already
                        # start + 1.
                        step = stop - start
                    if start >= 0:
                        stop -= (stop - start - 1) % step
                else: # step < 0
                    if stop >= start:
                        start, stop, step = 0, 0, 1
                    elif start < 0 and start + step <= stop:
                        if start < -1:
                            stop, step = start + 1, 1
                        else: # start == -1
                            stop, step = start - 1, -1
                    elif stop == start - 1:
                        stop, step = start + 1, 1
                    elif start >= 0 and start + step <= stop:
                        # Indexes 0 or 1 elements. We can't change stop
                        # because start might clip to a smaller true start if
                        # the axis is smaller than it, and increasing stop
                        # would prevent it from indexing an element in that
                        # case. The exception is the case right before this
                        # one (stop == start - 1). In that case start cannot
                        # clip past the stop (it always indexes the same one
                        # element in the cases where it indexes anything at
                        # all).
                        step = stop - start
                    if start < 0:
                        stop -= (stop - start + 1) % step
            elif start >= 0 and stop < 0 and step < 0 and (start < -step or
                                                           -stop - 1 < -step):
                if stop == -1:
                    start, stop, step = 0, 0, 1
                else:
                    step = max(-start - 1, stop + 1)
            elif start < 0 and stop == 0 and step > 0:
                start, stop, step = 0, 0, 1
            elif start < 0 and stop >= 0 and step >= min(-start, stop):
                step = min(-start, stop)
                if start == -1 or stop == 1:
                    # Can only index 0 or 1 elements. We can either pick a
                    # version with positive start and negative step, or
                    # negative start and positive step. We prefer the former
                    # as it matches what is done for reduce() with a shape
                    # (start is always nonnegative).
                    assert step == 1
                    start, stop, step = stop - 1, start - 1, -1
        elif start is not None and stop is None:
            if start == -1 and step > 0:
                start, stop, step = (-1, -2, -1)
            elif start < 0 and step >= -start:
                step = -start
            elif step < 0:
                if start == 0:
                    start, stop, step = 0, 1, 1
                elif 0 <= start < -step:
                    step = -start - 1
        if shape is None:
            return type(self)(start, stop, step, _reduced=True)

        # Further canonicalize with an explicit array shape

        shape = asshape(shape, axis=axis)
        size = shape[axis]

        if stop is None:
            if step > 0:
                stop = size
            else:
                stop = -size - 1

        if stop < -size:
            stop = -size - 1

        if size == 0:
            start, stop, step = 0, 0, 1
        elif step > 0:
            # start cannot be None
            if start < 0:
                start = size + start
            if start < 0:
                start = 0
            if start >= size:
                start, stop, step = 0, 0, 1

            if stop < 0:
                stop = size + stop
                if stop < 0:
                    stop = 0
            else:
                stop = min(stop, size)
            stop -= (stop - start - 1) % step

            if stop - start == 1:
                # Indexes 1 element.
                step = 1
            elif stop - start <= 0:
                start, stop, step = 0, 0, 1
        else:
            if start < 0:
                if start >= -size:
                    start = size + start
                else:
                    start, stop = 0, 0
            if start >= 0:
                start = min(size - 1, start)

            if -size <= stop < 0:
                stop += size

            if stop >= 0:
                if start - stop == 1:
                    stop, step = start + 1, 1
                elif start - stop <= 0:
                    start, stop, step = 0, 0, 1
                else:
                    stop += (start - stop - 1) % -step

            # start >= 0
            if (stop < 0 and start - size - stop <= -step
                or stop >= 0 and start - stop <= -step):
                stop, step = start + 1, 1
            if stop < 0 and start % step != 0:
                # At this point, negative stop is only necessary to index the
                # first element. If that element isn't actually indexed, we
                # prefer a nonnegative stop. Otherwise, stop will be -size - 1.
                stop = start % -step - 1
        return self.__class__(start, stop, step, _reduced=True)

    def isvalid(self, shape):
        # The docstring for this method is on the NDIndex base class
        shape = asshape(shape)

        # All slices are valid as long as there is at least one dimension
        return bool(shape)

    def newshape(self, shape):
        # The docstring for this method is on the NDIndex base class
        shape = asshape(shape)

        idx = self.reduce(shape)

        # len() won't raise an error after reducing with a shape
        return (len(idx),) + shape[1:]

    # TODO: Better name?
    def as_subindex(self, index):
        # The docstring of this method is currently on NDindex.as_subindex, as
        # this is the only method that is actually implemented so far.
        index = ndindex(index)
        index_orig = index

        s = self.reduce()
        index = index.reduce()

        if isinstance(index, Tuple):
            return Tuple(self).as_subindex(index)

        if isinstance(index, Integer):
            if index == -1:
                s = self.as_subindex(Slice(index.args[0], None))
            else:
                s = self.as_subindex(Slice(index.args[0], index.args[0] + 1))
            if s == Slice(0, 0, 1):
                # There is no index that we can return here. The intersection
                # of `self` and `index` is empty. Ideally we want to give an
                # index that gives an empty array, but we cannot make the
                # shape match. If a is dimension 1, then a[index] is dimension
                # 0, so a[index][slice(0, 0)] will not work. A possibility
                # would be to return False, which would add a length-0
                # dimension to the array. But
                #
                # 1. this isn't implemented yet, and
                # 2. a False can only add a length-0 dimension once, so it
                #    still wouldn't work in every case. For example,
                #    Tuple(slice(0), slice(0)).as_subindex((0, 0)) would need
                #    to return an index that replaces the first two
                #    dimensions with length-0 dimensions.
                raise ValueError(f"{self} and {index_orig} do not intersect")
            assert len(s) == 1
            return Tuple()

        if s.step < 0:
            raise NotImplementedError("Slice.as_subindex() is only implemented for slices with positive steps")

        # After reducing, start is not None when step > 0
        if s.stop is None or s.start < 0 or s.stop < 0:
            raise NotImplementedError("Slice.as_subindex() is only implemented for slices with nonnegative start and stop. Try calling reduce() with a shape first.")

        if isinstance(index, IntegerArray):
            idx = index.array
            if (idx < 0).any():
                raise NotImplementedError("Slice.as_subindex(IntegerArray) is not yet implemented for arrays with negative values. Try calling reduce with a shape first.")
            start, stop, step = subindex_slice(s.start, s.stop, s.step,
                                                     idx, idx+1, 1)
            res = BooleanArray(start < stop)

            if not res.count_nonzero:
                raise ValueError("Indices do not intersect")

            return res

        if not isinstance(index, Slice):
            raise NotImplementedError("Slice.as_subindex() is only implemented for tuples, integers, arrays and slices")

        if index.step < 0:
            raise NotImplementedError("Slice.as_subindex() is only implemented for slices with positive steps")

        # After reducing, start is not None when step > 0
        if index.stop is None or index.start < 0 or index.stop < 0:
            raise NotImplementedError("Slice.as_subindex() is only implemented for slices with nonnegative start and stop. Try calling reduce() with a shape first.")

        return Slice(*subindex_slice(s.start, s.stop, s.step, index.start,
                                     index.stop, index.step)).reduce()

    def isempty(self, shape=None):
        if shape is not None:
            return 0 in self.newshape(shape)

        try:
            l = len(self)
        except (TypeError, ValueError):
            return False
        return l == 0

    def selected_indices(self, shape, axis=None):
        if axis is None:
            yield from self.expand(shape).selected_indices(shape)
        else:
            shape = asshape(shape, axis=axis)
            for i in range(shape[axis])[self.raw]:
                yield Integer(i)

# Imports at the bottom to avoid circular import issues
from .ndindex import ndindex
from .tuple import Tuple
from .integer import Integer
from .integerarray import IntegerArray
from .booleanarray import BooleanArray
