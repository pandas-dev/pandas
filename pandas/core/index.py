# pylint: disable=E1101,E1103,W0232

from datetime import time
from itertools import izip

import numpy as np

from pandas.core.common import _format, adjoin, _ensure_index, _is_bool_indexer
import pandas.core.common as common
import pandas._tseries as _tseries

__all__ = ['Index']

def _indexOp(opname):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    def wrapper(self, other):
        func = getattr(self.view(np.ndarray), opname)
        return func(other)
    return wrapper

class Index(np.ndarray):
    """
    Immutable ndarray implementing an ordered, sliceable set

    Parameters
    ----------
    data : array-like (1-dimensional)
    dtype : NumPy dtype (default: object)
    copy : bool
        Make a copy of input ndarray

    Note
    ----
    An Index instance can **only** contain immutable objects for
    reasons of hashability.
    """
    def __new__(cls, data, dtype=object, copy=False):
        if isinstance(data, np.ndarray):
            subarr = np.array(data, dtype=dtype, copy=copy)
        elif np.isscalar(data):
            raise ValueError('Index(...) must be called with a collection '
                             'of some kind, %s was passed' % repr(data))
        else:
            subarr = np.empty(len(data), dtype=dtype)
            subarr[:] = data

        assert(subarr.ndim == 1)
        return subarr.view(cls)

    def summary(self):
        if len(self) > 0:
            index_summary = ', %s to %s' % (str(self[0]), str(self[-1]))
        else:
            index_summary = ''
        return 'Index: %s entries%s' % (len(self), index_summary)

    @property
    def indexMap(self):
        if not hasattr(self, '_cache_indexMap'):
            self._cache_indexMap = _tseries.map_indices_buf(self)
            self._verify_integrity()

        return self._cache_indexMap

    def is_all_dates(self):
        if not hasattr(self, '_cache_allDates'):
            self._cache_allDates = _tseries.isAllDates(self)

        return self._cache_allDates

    def _verify_integrity(self):
        if len(self.indexMap) < len(self):
            raise Exception('Index cannot contain duplicate values!')

    def __iter__(self):
        return iter(self.view(np.ndarray))

    def __setstate__(self, state):
        """Necessary for making this object picklable"""
        np.ndarray.__setstate__(self, state)

    def __deepcopy__(self, memo={}):
        """
        Index is not mutable, so disabling deepcopy
        """
        return self

    def __contains__(self, key):
        return key in self.indexMap

    def __hash__(self):
        return hash(self.view(np.ndarray))

    def __setitem__(self, key, value):
        """Disable the setting of values."""
        raise Exception(str(self.__class__) + ' object is immutable' )

    def __getitem__(self, key):
        """Override numpy.ndarray's __getitem__ method to work as desired"""
        arr_idx = self.view(np.ndarray)
        if np.isscalar(key):
            return arr_idx[key]
        else:
            if _is_bool_indexer(key):
                key = np.asarray(key)

            return Index(arr_idx[key])

    def format(self):
        """
        Render a string representation of the Index
        """
        if self.is_all_dates():
            to_join = []
            zero_time = time(0, 0)
            for dt in self:
                if dt.time() != zero_time or dt.tzinfo is not None:
                    return '\n'.join(str(x) for x in self)
                to_join.append(dt.strftime("%Y-%m-%d"))
            return '\n'.join(to_join)

        return '\n'.join(str(x) for x in self)

    def equals(self, other):
        """
        Determines if two Index objects contain the same elements.
        """
        if self is other:
            return True

        if not isinstance(other, Index):
            return False

        return np.array_equal(self, other)

    def asOfDate(self, date):
        if date not in self.indexMap:
            loc = self.searchsorted(date, side='left')
            if loc > 0:
                return self[loc-1]
            else:
                return None

        return date

    def sort(self, *args, **kwargs):
        raise Exception('Cannot sort an Index object')

    def shift(self, periods, offset):
        if periods == 0:
            # OK because immutable
            return self

        offset = periods * offset
        return Index([idx + offset for idx in self])

    def argsort(self, *args, **kwargs):
        return self.view(np.ndarray).argsort(*args, **kwargs)

    def __add__(self, other):
        if isinstance(other, Index):
            return self.union(other)
        else:
            return Index(self.view(np.ndarray) + other)

    __eq__ = _indexOp('__eq__')
    __ne__ = _indexOp('__ne__')
    __lt__ = _indexOp('__lt__')
    __gt__ = _indexOp('__gt__')
    __le__ = _indexOp('__le__')
    __ge__ = _indexOp('__ge__')

    def union(self, other):
        """
        Form the union of two Index objects and sorts if possible

        Parameters
        ----------
        other : Index or array-like

        Returns
        -------
        Index
        """
        if not hasattr(other, '__iter__'):
            raise Exception('Input must be iterable!')

        if len(other) == 0 or self.equals(other):
            return self

        new_seq = np.concatenate((self, other))
        try:
            new_seq = np.unique(new_seq)
        except Exception:
            # Not sortable / multiple types
            pass
        return Index(new_seq)

    def intersection(self, other):
        """
        Form the intersection of two Index objects and sorts if possible

        Parameters
        ----------
        other : Index or array-like

        Returns
        -------
        Index
        """
        if not hasattr(other, '__iter__'):
            raise Exception('Input must be iterable!')

        if self.equals(other):
            return self

        theIntersection = sorted(set(self) & set(other))
        return Index(theIntersection)

    def diff(self, other):
        if not hasattr(other, '__iter__'):
            raise Exception('Input must be iterable!')

        if self.equals(other):
            return Index([])

        otherArr = np.asarray(other)
        theDiff = sorted(set(self) - set(otherArr))
        return Index(theDiff)

    __sub__ = lambda self, other: self.diff(other)

    def take(self, *args, **kwargs):
        taken = self.view(np.ndarray).take(*args, **kwargs)
        return Index(taken)

    def get_loc(self, key):
        return self.indexMap[key]

    def get_indexer(self, target, method=None):
        """

        Parameters
        ----------
        target : Index
        method :

        Returns
        -------
        (indexer, mask)
        """
        if method:
            method = method.upper()

        aliases = {
            'FFILL' : 'PAD',
            'BFILL' : 'BACKFILL'
        }

        target = _ensure_index(target)

        method = aliases.get(method, method)
        indexer, mask = _tseries.getFillVec(self, target, self.indexMap,
                                            target.indexMap, method)
        return indexer, mask

    def slice_locs(self, start=None, end=None):
        """


        Returns
        -------
        (begin, end) : tuple

        Notes
        -----
        This function assumes that the data is sorted, so use at your own peril
        """
        if start is None:
            beg_slice = 0
        elif start in self:
            beg_slice = self.indexMap[start]
        else:
            beg_slice = self.searchsorted(start, side='left')

        if end is None:
            end_slice = len(self)
        elif end in self.indexMap:
            end_slice = self.indexMap[end] + 1
        else:
            end_slice = self.searchsorted(end, side='right')

        return beg_slice, end_slice


class DateIndex(Index):
    pass


class Factor(object):
    """
    Represents a categorical variable in classic R / S-plus fashion
    """
    def __init__(self, labels, levels):
        self.labels = labels
        self.levels = levels

    @classmethod
    def fromarray(cls, values):
        values = np.asarray(values, dtype=object)
        levels = Index(sorted(set(values)))
        labels, _ = _tseries.getMergeVec(values, levels.indexMap)
        return Factor(labels, levels=levels)

    def asarray(self):
        return self.levels[self.labels]

    def __len__(self):
        return len(self.labels)

    def __repr__(self):
        temp = 'Factor:\n%s\nLevels (%d): %s'
        values = self.asarray()
        return temp % (repr(values), len(self.levels), self.levels)

    def __getitem__(self, key):
        if isinstance(key, (int, np.integer)):
            i = self.labels[key]
            return self.levels[i]
        else:
            new_labels = self.labels[key]
            return Factor(new_labels, self.levels)


"""
something like this?

            ----------------------------------------------------
            | common                     | uncommon            |
            ----------------------------------------------------
            | foo     | bar     | baz    | qux     | wibble    |
            ----------------------------------------------------
A       1     ...       ...       ...      ...       ...
        2
        3
B       1
        2
C       1
        2
        3     ...       ...       ...      ...       ...


 common                uncommon
 ------                --------
|foo    bar    baz    |qux    wibble

"""

class MultiIndex(Index):
    """
    Implements multi-level, a.k.a. hierarchical, index object for pandas objects

    Parameters
    ----------
    levels : list or tuple of arrays
    labels :

    """

    def __new__(cls, levels=None, labels=None, sortorder=None):
        return np.arange(len(labels[0]), dtype=object).view(cls)

    def __init__(self, levels, labels, sortorder=None):
        self.levels = [_ensure_index(lev) for lev in levels]
        self.labels = [np.asarray(labs, dtype=np.int32) for labs in labels]

        if sortorder is not None:
            self.sortorder = int(sortorder)
        else:
            self.sortorder = sortorder

    def __iter__(self):
        values = [np.asarray(lev).take(lab)
                  for lev, lab in zip(self.levels, self.labels)]
        return izip(*values)

    def __contains__(self, key):
        try:
            label_key = self._get_label_key(key)
            return label_key in self.indexMap
        except Exception:
            return False

    def format(self, space=2):
        stringified_levels = [lev.format().split('\n') for lev in self.levels]

        result_levels = []
        for lab, lev in zip(self.labels, stringified_levels):
            taken = np.array(lev, dtype=object).take(lab)
            result_levels.append(taken)

        return adjoin(space, *result_levels)

    def is_all_dates(self):
        return False

    @classmethod
    def from_arrays(cls, arrays, sortorder=None):
        """
        Convert arrays to MultiIndex

        Parameters
        ----------
        arrays : list / sequence
        sortorder : int or None
            Level of sortedness (must be lexicographically sorted by that level)

        Returns
        -------
        index : MultiIndex
        """
        levels = []
        labels = []
        for arr in arrays:
            factor = Factor.fromarray(arr)
            levels.append(factor.levels)
            labels.append(factor.labels)

        return MultiIndex(levels=levels, labels=labels, sortorder=sortorder)

    @classmethod
    def from_tuples(cls, tuples, sortorder=None):
        arrays = zip(*tuples)
        return MultiIndex.from_arrays(arrays, sortorder=sortorder)

    @property
    def indexMap(self):
        if not hasattr(self, '_cache_indexMap'):
            zipped = zip(*self.labels)
            self._cache_indexMap = _tseries.map_indices_list(zipped)
            self._verify_integrity()

        return self._cache_indexMap

    @property
    def nlevels(self):
        return len(self.levels)

    @property
    def levshape(self):
        return tuple(len(x) for x in self.levels)

    def __reduce__(self):
        """Necessary for making this object picklable"""
        object_state = list(np.ndarray.__reduce__(self))
        subclass_state = (self.levels, self.labels)
        object_state[2] = (object_state[2], subclass_state)
        return tuple(object_state)

    def __setstate__(self, state):
        """Necessary for making this object picklable"""
        nd_state, own_state = state
        np.ndarray.__setstate__(self, nd_state)
        levels, labels, = own_state

        self.levels = [Index(x) for x in levels]
        self.labels = labels

    def __getitem__(self, key):
        arr_idx = self.view(np.ndarray)
        if np.isscalar(key):
            return tuple(lev[lab[key]]
                         for lev, lab in zip(self.levels, self.labels))
        else:
            if _is_bool_indexer(key):
                key = np.asarray(key)
                sortorder = self.sortorder
            else:
                # cannot be sure whether the result will be sorted
                sortorder = None

            new_tuples = arr_idx[key]
            new_labels = [lab[key] for lab in self.labels]

            # an optimization
            result = new_tuples.view(MultiIndex)
            result.levels = self.levels
            result.labels = new_labels
            result.sortorder = sortorder

            return result

    def __getslice__(self, i, j):
        return self.__getitem__(slice(i, j))

    def sortlevel(self, level=0, ascending=True):
        """

        """
        labels = list(self.labels)
        primary = labels.pop(level)

        # Lexsort starts from END
        indexer = np.lexsort(tuple(labels[::-1]) + (primary,))

        if not ascending:
            indexer = indexer[::-1]

        new_labels = [lab.take(indexer) for lab in self.labels]
        new_index = MultiIndex(levels=self.levels, labels=new_labels,
                               sortorder=level)

        return new_index, indexer

    def get_indexer(self, target, method=None):
        """

        Parameters
        ----------
        target : Index
        method :

        Returns
        -------
        (indexer, mask)
        """
        if method:
            method = method.upper()

        aliases = {
            'FFILL' : 'PAD',
            'BFILL' : 'BACKFILL'
        }
        method = aliases.get(method, method)

        if not isinstance(target, MultiIndex):
            raise TypeError('Can only align with other MultiIndex objects')

        self_index = self.get_tuple_index()
        target_index = target.get_tuple_index()

        indexer, mask = _tseries.getFillVec(self_index, target_index,
                                            self_index.indexMap,
                                            target_index.indexMap, method)
        return indexer, mask

    def get_tuple_index(self):
        to_join = []
        for lev, lab in zip(self.levels, self.labels):
            to_join.append(np.asarray(lev).take(lab))

        return Index(zip(*to_join))

    def slice_locs(self, start=None, end=None):
        """

        Returns
        -------

        Notes
        -----
        This function assumes that the data is sorted by the first level
        """
        assert(self.sortorder == 0)

        level0 = self.levels[0]

        if start is None:
            start_slice = 0
        else:
            if not isinstance(start, tuple):
                start = start,
            start_slice = self._partial_tup_index(start, side='left')

        if end is None:
            end_slice = len(self)
        else:
            if not isinstance(end, tuple):
                end = end,
            end_slice = self._partial_tup_index(end, side='right')

        return start_slice, end_slice

    def _partial_tup_index(self, tup, side='left'):
        assert(self.sortorder == 0)
        tup_labels = self._get_label_key_approx(tup, side=side)

        n = len(tup_labels)
        start, end = 0, len(self)
        for k, (idx, labs) in enumerate(zip(tup_labels, self.labels)):
            section = labs[start:end]
            if k < n - 1:
                start = start + section.searchsorted(idx, side='left')
                end = start + section.searchsorted(idx, side='right')
            else:
                return start + labs.searchsorted(idx, side=side)

    def _get_label_key_approx(self, tup, side='left'):
        result = []
        for lev, v in zip(self.levels, tup):
            try:
                label = lev.get_loc(v)
            except KeyError:
                label = lev.searchsorted(v, side=side)

                if side == 'right' and label > 0:
                    label -= 1

            result.append(label)
        return tuple(result)

    def get_loc(self, key):
        if isinstance(key, tuple):
            return self._get_tuple_loc(key)
        else:
            assert(self.sortorder == 0)
            # slice level 0
            level = self.levels[0]
            labels = self.labels[0]

            loc = level.get_loc(key)
            i = labels.searchsorted(loc, side='left')
            j = labels.searchsorted(loc, side='right')
            return slice(i, j)

    def _get_tuple_loc(self, tup):
        indexer = self._get_label_key(tup)
        try:
            return self.indexMap[indexer]
        except KeyError:
            raise KeyError(str(tup))

    def _get_label_key(self, tup):
        return tuple(lev.get_loc(v) for lev, v in zip(self.levels, tup))

    def truncate(self, before=None, after=None):
        """
        Slice index between two major axis values, return new MultiIndex

        Parameters
        ----------
        before : type of major_axis values or None, default None
            None defaults to start of panel

        after : type of major_axis values or None, default None
            None defaults to after of panel

        Returns
        -------
        MultiIndex
        """
        i, j = self._get_axis_bounds(before, after)
        left, right = self._get_label_bounds(i, j)

        new_levels = list(self.levels)
        new_levels[0] = new_levels[0][i:j]

        new_labels = [lab[left:right] for lab in self.labels]
        new_labels[0] = new_labels[0] - i

        return MultiIndex(levels=new_levels, labels=new_labels)

    def equals(self, other):
        """
        Determines if two MultiIndex objects are the same
        """
        if self is other:
            return True

        if not isinstance(other, MultiIndex):
            return False

        if self.nlevels != other.nlevels:
            return False

        if len(self) != len(other):
            return False

        # if not self.equal_levels(other):
        #     return False

        for i in xrange(self.nlevels):
            if not self.levels[i].equals(other.levels[i]):
                return False
            if not np.array_equal(self.labels[i], other.labels[i]):
                return False
        return True

    def equal_levels(self, other):
        if self.nlevels != other.nlevels:
            return False

        for i in xrange(self.nlevels):
            if not self.levels[i].equals(other.levels[i]):
                return False
        return True

    def union(self, other):
        """
        Form the union of two MultiIndex objects

        Parameters
        ----------
        other : Index or array-like

        Returns
        -------
        Index
        """
        self._assert_can_do_setop(other)

        if len(other) == 0 or self.equals(other):
            return self

        # TODO: optimize / make less wasteful
        self_tuples = self.get_tuple_index()
        other_tuples = other.get_tuple_index()
        uniq_tuples = np.unique(np.concatenate((self_tuples, other_tuples)))
        return MultiIndex.from_arrays(zip(*uniq_tuples), sortorder=0)

    def intersection(self, other):
        """
        Form the intersection of two Index objects and sorts if possible

        Parameters
        ----------
        other : Index or array-like

        Returns
        -------
        Index
        """
        self._assert_can_do_setop(other)

        if self.equals(other):
            return self

        # TODO: optimize / make less wasteful
        self_tuples = self.get_tuple_index()
        other_tuples = other.get_tuple_index()
        uniq_tuples = sorted(set(self_tuples) & set(other_tuples))
        return MultiIndex.from_arrays(zip(*uniq_tuples), sortorder=0)

    def _assert_can_do_setop(self, other):
        if not isinstance(other, MultiIndex):
            raise TypeError('can only call with other hierarchical '
                            'index objects')

        assert(self.nlevels == other.nlevels)

    def get_major_bounds(self, begin=None, end=None):
        """
        Return index bounds for slicing LongPanel labels and / or
        values

        Parameters
        ----------
        begin : axis value or None
        end : axis value or None

        Returns
        -------
        y : tuple
            (left, right) absolute bounds on LongPanel values
        """
        i, j = self._get_axis_bounds(begin, end)
        left, right = self._get_label_bounds(i, j)

        return left, right

    def _get_axis_bounds(self, begin, end):
        """
        Return major axis locations corresponding to interval values
        """
        if begin is not None:
            i = self.levels[0].indexMap.get(begin)
            if i is None:
                i = self.levels[0].searchsorted(begin, side='right')
        else:
            i = 0

        if end is not None:
            j = self.levels[0].indexMap.get(end)
            if j is None:
                j = self.levels[0].searchsorted(end)
            else:
                j = j + 1
        else:
            j = len(self.levels[0])

        if i > j:
            raise ValueError('Must have begin <= end!')

        return i, j

    def _get_label_bounds(self, i, j):
        "Return slice points between two major axis locations"

        left = self._bounds[i]

        if j >= len(self.levels[0]):
            right = len(self.labels[0])
        else:
            right = self._bounds[j]

        return left, right

    __bounds = None
    @property
    def _bounds(self):
        "Return or compute and return slice points for major axis"
        if self.__bounds is None:
            inds = np.arange(len(self.levels[0]))
            self.__bounds = self.labels[0].searchsorted(inds)

        return self.__bounds

# For utility purposes

NULL_INDEX = Index([])

