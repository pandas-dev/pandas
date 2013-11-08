import pandas as pd
import numpy as np


class RangeIndex(object):
    """Represents a range with left-open interval.

    Parameters
    ----------
    left, right : int
        start and end of range (i.e., ``s[left:right]``). If left > right,
        assumes reversed (i.e., ``s[left:right:-1]``)
    """
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
    def __new__(cls, left, right):
        if left == right:
            return pd.Index([], dtype='int64')
        else:
            return object.__new__(cls)
    def __init__(self, left, right):
        # shouldn't happen
        if left == right:
            raise ValueError("Can't have empty range")

        l, r = left, right
        left, right = int(left), int(right)

        if left != l or right != r:
            raise ValueError("Need to pass integral values")

        self.left = left
        self.right = right
        self.ascending = left < right
        if self.ascending:
            self.start, self.end = left, right
            self.step = 1
        else:
            self.start, self.end = right, left
            self.step = -1

    @property
    def values(self):
        return np.arange(self.left, self.right, (1 if self.ascending else -1), dtype='int64')

    def union(self, other):
        """Union this with another RangeIndex. Always returns ascending RangeIndex."""
        if not isinstance(other, RangeIndex):
            raise NotImplementedError("Other not range index")

        if not self._overlaps(other):
            return self.values | other.values

        start = min(self.start, other.start)
        end = max(self.end, other.end)
        return RangeIndex(start, end)

    def intersection(self, other):
        if not isinstance(other, RangeIndex):
            raise NotImplementedError("Other not range index")
        # not overlapping or start touches end or vice versa
        if not self._overlaps(other) or (self.start == other.end) or (self.end == other.start):
            return pd.Index([], dtype='int64')
        else:
            return RangeIndex(max(self.start, other.start), min(self.end, other.end))

    def view(self, other):
        return self

    def difference(self, other):
        if not isinstance(other, RangeIndex):
            raise NotImplementedError("Other not range index")

        if not self._overlaps(other) or self.start == other.end or self.end == other.start:
            return self.view()
        # completely contained
        elif self.start >= other.start and self.end <= other.end:
            return pd.Index([], dtype='int64')
        elif self.start < other.start:
            return RangeIndex(self.start, other.start)
        # starts within other [because must overlap]
        elif self.start > other.start:
            assert other.end > self.end, (self, other)
            return RangeIndex(other.end, self.end)
        assert False, "Shouldn't get to here"

    @property
    def empty(self):
        return False

    def all(self):
        return True

    def any(self):
        return True

    def __array__(self):
        return self.values

    def __or__(self, other):
        return self.union(other)


    __add__ = __or__

    def __and__(self, other):
        return self.intersection(other)

    def __sub__(self, other):
        return self.difference(other)

    def equals(self, other):
        return self.left == other.left and self.right == other.right

    def _overlaps(self, other):
        # cheers to Ned Batchelder on this
        # only overlaps if each ranges' end is beyond or *at* the other ranges' start.
        # touching does not count as overlapping
        return other.end > self.start and self.end > other.start

        # # starts before or on other's start and ends after or on other's start
        # return ((self.start <= other.start and self.end >= other.start) or
        # # starts within other
        # (self.start > other.start and self.start <= other.end))
    def nonzero(self):
        if self.start > 0:
            return np.arange(len(self))
        else:
            # need to skip when self is zero
            res = range(len(self))
            res.pop(0 - self.start * self.step)
            return np.array(res)

    def __contains__(self, val):
        # can only hold integers
        try:
            val = int(val)
        except (TypeError, ValueError):
            return False

        if val != val:
            return False

        return self.start <= val < self.end

    def __iter__(self):
        return iter(xrange(self.left, self.right, self.step))

    def __len__(self):
        return self.end - self.start

    def __repr__(self):
        # TODO: Either change to Int64Repr OR to RangeIndex(left, right)
        return "RangeIndex(%r)" % (dict(left=self.left, right=self.right, start=self.start, end=self.end, ascending=self.ascending))

    def get_indexer(self, arr, method=None):
        arr = np.asarray(arr, dtype='int64')
        indexer = arr - self.start
        if not self.ascending:
            indexer = (len(self) - 1) - indexer
        indexer[(indexer < 0) | (indexer >= len(self) )] = -1
        return indexer

    def get_loc(self, val):
        if val in self:
            return val - (self.start if self.ascending else self.end)
        else:
            return -1

    def __getitem__(self, val):
        if isinstance(val, slice):
            if slice.step not in (1, -1):
                return self.values[val]
            if slice.start >= 0 and slice.end >= 0:
                start = slice.start if slice.start is None or slice.start > self.start else self.start
                end = slice.end if slice.end is None or slice.end < self.end else self.end

            if self.step != slice.step:
                start, end = end, start

            return RangeIndex(start, end)
        else:
            if 0 <= val < len(self):
                return self.left + val * self.step
            elif -len(self) <= val < 0:
                return self.right + val * self.step
            else:
                raise IndexError("%d out of range" % val)

    def __bool__(self):
        raise ValueError("The truth value of an array is ambiguous...") # blah blah blah

    __nonzero__ = __bool__

    def slice_locs(self, start=None, end=None):
        pass

    def get_indexer_non_unique(self, arr):
        return self.get_indexer(self, arr), np.array([], dtype='int64')
