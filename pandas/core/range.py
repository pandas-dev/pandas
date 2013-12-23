import pandas as pd
import numpy as np
import pandas.algos as _algos
from pandas.core.index import Int64Index, Index
import pandas.core.index as pdindex

EMPTY_RANGE = lambda: pd.Index([], dtype='int64')

def _delegate_to_int64(func_name, doc=None):
    def wrapper(self, *args, **kwargs):
        return getattr(self.as_int64index, func_name)(*args, **kwargs)
    wrapper.__name__ = func_name
    if hasattr(Int64Index, func_name):
        doc = doc or getattr(Int64Index, func_name).__doc__
    wrapper.__doc__ = doc or ''
    return wrapper

def _not_implemented(self, *args, **kwargs):
    raise NotImplementedError('method not implemented')

class RangeIndex(Index):
    """Represents a range with left-open interval.

    Parameters
    ----------
    left, right : int
        start and end of range (i.e., ``s[left:right]``). If left > right,
        assumes reversed (i.e., ``s[left:right:-1]``)
    """
    _groupby = _algos.groupby_range
    # My thinking on this structure:
    # - from a 'set theoretic' standpoint, order doesn't matter, so ascending
    #   vs. descending is generally an implementation detail
    # - left, right are necessary for comparison, but making start and end
    #   separate makes it much easier to work with
    # - prohibiting empty RangeIndex helps simplify edge cases considerably.
    # - From pandas' perspective, RangeIndex should behave *exactly* the same
    #   as an Int64Index. (except for ops where it can yield more RangeIndexes)
    # - supporting steps might be possible, but less simple (and less clear
    #   that it would be useful to pandas proper, given that you'd have to know
    #   ahead of time that you could convert to a stepped range). Plus, only
    #   helpful when you specifically have a consistent step
    # - Certain operations with RangeIndex could benefit from allowing nested
    #   ranges, i.e. the union of RangeIndex(10, 5) and RangeIndex(3, 7) could be
    #   something like: [RangeIndex(10,5), RangeIndex(3, 7)] and then could just
    #   iterate over that. But that's for after all of this starts working.
    # - It would be nice if groupby() accepted an Index
    # - It might be valuable to cache the values property of a RangeIndex, but
    #   I'm not totally convinced that's the best strategy (yet!), unless RangeIndex
    #   could *replace* itself on its parent. (tradeoff between instantiation time and
    #   ability to gc values when they aren't needed anymore)
    # TODO: Block setting of start and end
    def __new__(cls, left, right=None, step=1, name=None):
        if step not in (1, -1):
            raise ValueError("Invalid step %s" % step)
        if right is None:
            left, right = 0, left
        # RangeIndex can never be empty (makes checking things simpler)
        if left == right:
            return EMPTY_RANGE()
        self = np.array([], dtype='int64').view(RangeIndex)
        # shouldn't happen
        if left == right:
            raise ValueError("Can't have empty range")

        # TODO: builtin range function only accepts integers, could make more
        # sense to just do a single isinstance check there. Not sure...
        # want to coerce where possible, but not change value
        l, r = left, right
        left, right = int(left), int(right)

        if left != l or right != r:
            raise TypeError("Need to pass integral values")

        self.left = left
        self.right = right
        self.ascending = left < right
        self.name = name
        if self.ascending:
            self.start, self.stop = left, right
            self.step = 1
        else:
            # because non-inclusive. want start and stop to be the actual
            # bounds of the range if the range were ascending.
            # e.g., range(10, 5, -1) == range(6, 11)[::-1]
            self.start, self.stop = right + 1, left + 1
            self.step = -1
        return self

    @property
    def _constructor(self):
        # functions depend on being able to pass a list to this function,
        # so makes more sense to use Index [could potentially use Int64Index
        # instead]
        return Index

    # Stuff that also needs to / could be overriden:
        # _has_valid_type funcs (_convert_scalar_indexer,
        #                        _convert_slice_indexer, etc)

    def map(self, mapper):
        return _algos.arrmap_range(self, mapper)

    def groupby(self, to_groupby):
        return self._groupby(self, to_groupby)

    def __array_finalize__(self, obj):
        if not isinstance(obj, type(self)):
            return
        for attr in ('ascending', 'left', 'right', 'start',
                     'stop', 'step'):
            if hasattr(obj, attr):
                setattr(self, attr, getattr(obj, attr))
        self.name = getattr(obj, 'name', None)

    join = _delegate_to_int64('join')
    to_series = _delegate_to_int64('to_series')
    astype = _delegate_to_int64('astype')
    to_datetime = _delegate_to_int64('to_datetime')
    _format_native_types = _delegate_to_int64('_format_native_types')
    argsort = _delegate_to_int64('argsort')
    asof = _not_implemented
    asof_locs = _not_implemented
    nlevels = 1
    is_integer = lambda self: True
    is_floating = lambda self: False
    is_numeric = lambda self: True
    is_mixed = lambda self: False
    holds_integer = lambda self: True
    is_all_dates = lambda self: False
    is_unique = True
    get_duplicates = lambda self: []
    inferred_type = 'integer'

    def order(self, return_indexers=False, ascending=True):
        result = self
        left, right = self.left, self.right
        if ascending != self.ascending:
            result = result[::-1]
        if return_indexers:
            if ascending != self.ascending:
                indexer = np.arange(len(self) - 1, -1, -1)
            else:
                indexer = np.arange(len(self))
            return result, indexer
        else:
            return result

    @property
    def values(self):
        if self.ascending:
            vals = np.arange(self.left, self.right, 1, dtype='int64')
        else:
            vals = np.arange(self.left, self.right, -1, dtype='int64')
        return vals

    @property
    def as_int64index(self):
        # TODO: Maybe fastpath this!
        return Int64Index(self.values, name=self.name)

    def union(self, other):
        """Union this with another RangeIndex. Always returns ascending RangeIndex."""
        if not isinstance(other, RangeIndex):
            return self.as_int64index.union(other)

        if not self._overlaps(other):
            return self.values | other.values

        start = min(self.start, other.start)
        stop = max(self.stop, other.stop)
        return RangeIndex(start, stop)

    def intersection(self, other):
        if not isinstance(other, RangeIndex):
            return self.as_int64index.intersection(other)
        # not overlapping or start touches end or vice versa
        if not self._overlaps(other) or (self.start == other.stop) or (self.stop == other.start):
            return EMPTY_RANGE()
        else:
            return RangeIndex(max(self.start, other.start), min(self.stop, other.stop))

    def _shallow_copy(self):
        # recursion issue: index view() calls _shallow_copy(), probably need to
        # decide if _shallow_copy() is necessary.
        return RangeIndex(self.left, self.right, self.step)

    def view(self, *args, **kwargs):
        if not args and not kwargs:
            return self._shallow_copy()
        else:
            return self.as_int64index.view(*args, **kwargs)

    def difference(self, other):
        if not isinstance(other, RangeIndex):
            return self.as_int64index.difference(other)

        if not self._overlaps(other) or self.start == other.stop or self.stop == other.start:
            return self.view()
        # completely contained
        elif self.start >= other.start and self.stop <= other.stop:
            return EMPTY_RANGE()
        elif self.start < other.start:
            return RangeIndex(self.start, other.start)
        # starts within other [because must overlap]
        elif self.start > other.start:
            assert other.stop > self.stop, (self, other)
            return RangeIndex(other.stop, self.stop)
        assert False, "Shouldn't get to here"

    @property
    def empty(self):
        return False

    @property
    def is_monotonic(self):
        # monotonically increasing
        return self.ascending

    def all(self):
        # False only if it *spans* zero
        return not (self.start <= 0 and self.stop > 0)

    def any(self):
        # if includes any number other than zero than any is True
        return len(self) > 1 or self.start != 0

    def __array__(self):
        return self.values

    # TODO: Probably remove these functions when Index is no longer a subclass of
    # ndarray [need to override them for now to make them work with np.asarray
    # and buddies].
    @property
    def __array_interface__(self):
        raise AttributeError("No attribute __array_interface__")

    @property
    def __array_struct__(self):
        raise AttributeError("No attribute __array_struct__ [disabled]")

    def __or__(self, other):
        return self.union(other)

    def __and__(self, other):
        return self.intersection(other)

    def __sub__(self, other):
        return self.difference(other)

    def equals(self, other):
        if not isinstance(other, RangeIndex):
            return self.as_int64index.equals(other)
        return self.ascending == other.ascending and self.start == other.start and self.stop == other.stop

    def identical(self, other):
        other = pdindex._ensure_index(other)
        return self.equals(other) and self.name == other.name

    def _overlaps(self, other):
        # cheers to Ned Batchelder on this
        # only overlaps if each ranges' end is beyond or *at* the other ranges' start.
        # touching does not count as overlapping
        return other.stop > self.start and self.stop > other.start

    def nonzero(self):
        if self.start > 0:
            return (np.arange(len(self)),)
        else:
            # need to skip when self is zero
            res = range(len(self))
            res.pop(0 - self.start * self.step)
            return (np.array(res),)

    def __contains__(self, val):
        # can only hold integers
        try:
            v = val
            val = int(val)
        except (TypeError, ValueError):
            return False

        # pd.isnull(val)?
        if v != val or val != val:
            return False

        return self.start <= val < self.stop

    def __iter__(self):
        return iter(xrange(self.left, self.right, self.step))

    def __len__(self):
        return self.stop - self.start

    def __str__(self):
        return str(self.as_int64index)

    def __repr__(self):
        # TODO: Either change to Int64Repr OR to RangeIndex(left, right)
        return "RangeIndex(%s, %s, %s)" % (self.left, self.right, self.step)

    def get_indexer(self, arr, method=None):
        """Returns indexer (i.e., matching indices between the index and
        arr)."""
        # bfill : will fill everything < start as 0; > stop filled as -1
        # ffill : will fill everything < start as -1; > stop filled as
        #         len(self) - 1
        if method not in (None, 'ffill', 'pad', 'bfill', 'backfill'):
            raise ValueError("Unknown method: %r" % method)
        if method and not self.is_monotonic:
            kind = 'forward' if method in ('ffill', 'pad') else 'backward'
            raise ValueError("Must be monotonic for %s fill." % kind)

        arr = np.asarray(arr, dtype='int64')
        indexer = arr - self.start
        if not self.ascending:
            indexer = (len(self) - 1) - indexer

        if method in ('ffill', 'pad'):
            # next valid observation always 0
            min_fill, max_fill = 0, -1
        elif method in ('bfill', 'backfill'):
            # last valid observation is len(self) - 1
            min_fill, max_fill = -1, len(self) - 1
        else:
            min_fill = max_fill = -1

        indexer[indexer < 0] = min_fill
        indexer[indexer >= len(self)] = max_fill
        return indexer

    def get_loc(self, val):
        if val not in self:
            return -1
        if self.ascending:
            return val - self.start
        else:
            return self.stop - val

    def __getitem__(self, val):
        # only want to handle the simple stuff here, otherwise let Int64Index
        # handle it
        if isinstance(val, slice) and val.step in (1, -1):
            # Step 1 - convert the slice to be forward index (i.e. step == -1
            #          --> step == 1)
            v_start, v_stop = _get_forward_indices(len(self), val)
            left, right = self.left, self.right
            step = 1 if self.ascending else -1
            if v_start is None:
                # empty range
                return EMPTY_RANGE()

            # Step 3 - increment left by start of slice
            left += v_start * step

            # Step 4 - set right to min(right, stop)

            # Step 5 - flip bounds if they were reversed
            return RangeIndex(start, stop, step)
        elif np.isscalar(val):
            if -len(self) <= val < 0:
                val = len(self) + val
            if val > len(self):
                raise IndexError("%d out of range" % val)
            step = 1 if self.ascending else -1
            return self.left + val * step
        else:
            return self.as_int64index[val]

    def __bool__(self):
        raise ValueError("The truth value of an array is ambiguous...") # blah blah blah

    __nonzero__ = __bool__

    # don't need to override slice_indexer
    def slice_locs(self, start=None, end=None):
        start = self.get_loc(start) if start is not None else 0
        end = self.get_loc(end) + 1 if end is not None else len(self)
        return start, end

    def get_indexer_non_unique(self, arr):
        return self.get_indexer(self, arr), np.array([], dtype='int64')

def _flip_bounds(start, stop, step):
    """Returns bounds and step for reversed range (where end is non-inclusive):
        >>> range(3, 6)
        [3, 4, 5]
        >>> _flip_bounds(3, 6, 1)
        (5, 2, -1)
        >>> range(*_flip_bounds(3, 6, 1))
        [5, 4, 3]
    """
    return stop - step, start - step, step * -1


def _get_forward_indices(length, slc):
    """Converts given slice to positive, forward step indices.
    Returns (None, None) if not possible to convert.

    >>> _get_forward_indices(10, slice(5, 1, -2))
    (2, 6)
    >>> _get_forward_indices(10, slice(-100, -90, 5))
    (None, None)
    >>> _get_forward_indices(5, slice(3, 4, 1))
    (3, 4)
    """
    if slc.step == 0 or length == 0:
        return None, None
    start, stop = slc.start, slc.stop
    if slc.step < 0:
        # when you flip direction, need to increment edges
        # e.g., [6:2:-1] --> [3:7][::-1]
        start = start + 1 if start is not None else length - 1
        stop = stop + 1 if stop is not None else 0
        start, stop = stop, start
    else:
        start = start if start is not None else 0
        stop = stop if stop is not None else length - 1

    if start >= stop or start > length or stop == 0 or stop < -length:
        return (None, None)
    if start < 0:
        start = length + start if start > -length else 0
    if stop < 0:
        stop = length + stop

    return start, min(length, stop)
