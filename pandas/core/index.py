# pylint: disable=E1101,E1103,W0232

import numpy as np

from pandas.core.common import _ensure_index, _is_bool_indexer
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
        else:
            subarr = np.empty(len(data), dtype=dtype)
            subarr[:] = data

        if subarr.ndim == 0:
            raise Exception('Index(...) must be called with a collection '
                            'of some kind, %s was passed' % repr(data))

        return subarr.view(cls)

    def __array_finalize__(self, obj):
        if self.ndim == 0:
            # tolist will cause a bus error if this is not here, hmm
            return self.item()
            # raise Exception('Cannot create 0-dimensional Index!')

        return

        # New instance creation
        if obj is None:
            pass
        # New from template / slicing
        elif isinstance(obj, type(self)) and len(self) != len(obj.indexMap):
            pass
        # View casting
        else:
            pass

    def summary(self):
        if len(self) > 0:
            index_summary = ', %s to %s' % (self[0], self[-1])
        else:
            index_summary = ''
        return 'Index: %s entries%s' % (len(self), index_summary)

    @property
    def indexMap(self):
        if not hasattr(self, '_cache_indexMap'):
            self._cache_indexMap = _tseries.map_indices(self)
            self._verify_integrity()

        return self._cache_indexMap

    @property
    def _allDates(self):
        if not hasattr(self, '_cache_allDates'):
            self._cache_allDates = _tseries.isAllDates(self)

        return self._cache_allDates

    def is_all_dates(self):
        return self._allDates

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

    def __contains__(self, date):
        return date in self.indexMap

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

            # easier to ask forgiveness than permission
            try:
                return Index(arr_idx[key])
            except Exception, e1:
                try:
                    return Index(arr_idx[np.asarray(key)])
                except Exception, e2: # pragma: no cover
                    raise e1

    def format(self):
        from datetime import time

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

    __sub__ = diff

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



class Factor(object):
    """
    Represents a categorical variable in classic R / S-plus fashion
    """
    def __init__(self, labels, levels):
        self.labels = labels
        self.levels = levels

    @classmethod
    def fromarray(cls, values):
        levels = np.array(sorted(set(values)), dtype=object)
        labels = levels.searchsorted(values)

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

class MultiLevelIndex(Index):
    """
    Implements multi-level, a.k.a. hierarchical, index object for pandas objects


    Parameters
    ----------
    levels : list or tuple of arrays
    labels :

    """

    def __new__(cls, levels=None, labels=None):
        arr = np.empty(len(labels[0]), dtype=object)
        arr[:] = zip(*labels)
        arr = arr.view(cls)
        return arr

    def __init__(self, levels=None, labels=None):
        self.levels = [_ensure_index(lev) for lev in levels]
        self.labels = [np.asarray(labs, dtype=np.int32) for labs in labels]
        self._verify_integrity()

    def __array_finalize__(self, obj):
        self.labels = getattr(obj, 'labels', None)
        self.levels = getattr(obj, 'levels', None)

    @property
    def nlevels(self):
        return len(levels)

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

    def format(self, space=2):
        from pandas.core.common import _format, adjoin

        stringified_levels = [lev.format().split('\n') for lev in self.levels]

        padded_levels = []
        for lab, lev in zip(self.labels, stringified_levels):
            maxlen = max(len(x) for x in lev)
            padded = [x.ljust(maxlen) for x in lev]
            padded = np.array(padded, dtype=object).take(lab)
            padded_levels.append(padded)

        return adjoin(2, *padded_levels)

    def sort(self, bylevel=0):
        pass

    @classmethod
    def from_arrays(cls, *arrays):
        levels = []
        labels = []
        return cls(levels, labels)

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

class LongPanelIndex(MultiLevelIndex):
    """
    Holds axis indexing information for a LongPanel instance

    Parameters
    ----------
    major_axis : Index-like
    minor_axis : Index-like
    major_labels : ndarray
    minor_labels : ndarray
    mask : ndarray (bool), optional
        observation selection vector using major and minor labels, for
        converting to wide format.
    """

    @property
    def major_axis(self):
        return self.levels[0]

    @property
    def minor_axis(self):
        return self.levels[1]


    @property
    def major_labels(self):
        return self.labels[0]

    @property
    def minor_labels(self):
        return self.labels[1]

    def truncate(self, before=None, after=None):
        """
        Slice index between two major axis values, return new
        LongPanelIndex

        Parameters
        ----------
        before : type of major_axis values or None, default None
            None defaults to start of panel

        after : type of major_axis values or None, default None
            None defaults to after of panel

        Returns
        -------
        LongPanelIndex
        """
        i, j = self._get_axis_bounds(before, after)
        left, right = self._get_label_bounds(i, j)

        return LongPanelIndex([self.major_axis[i : j],
                               self.minor_axis],
                              [self.major_labels[left : right] - i,
                               self.minor_labels[left : right]])

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
            i = self.major_axis.indexMap.get(begin)
            if i is None:
                i = self.major_axis.searchsorted(begin, side='right')
        else:
            i = 0

        if end is not None:
            j = self.major_axis.indexMap.get(end)
            if j is None:
                j = self.major_axis.searchsorted(end)
            else:
                j = j + 1
        else:
            j = len(self.major_axis)

        if i > j:
            raise ValueError('Must have begin <= end!')

        return i, j

    def _get_label_bounds(self, i, j):
        "Return slice points between two major axis locations"

        left = self._bounds[i]

        if j >= len(self.major_axis):
            right = len(self.major_labels)
        else:
            right = self._bounds[j]

        return left, right

    __bounds = None
    @property
    def _bounds(self):
        "Return or compute and return slice points for major axis"
        if self.__bounds is None:
            inds = np.arange(len(self.major_axis))
            self.__bounds = self.major_labels.searchsorted(inds)

        return self.__bounds

    @property
    def mask(self):
        return self._make_mask()
        # if self._mask is None:
        #     self._mask = self._make_mask()

        # return self._mask

    def _make_mask(self):
        """
        Create observation selection vector using major and minor
        labels, for converting to wide format.
        """
        N, K = self.levshape
        selector = self.minor_labels + K * self.major_labels

        mask = np.zeros(N * K, dtype=bool)
        mask[selector] = True

        return mask

    # @property
    # def shape(self):
    #     return len(self.major_axis), len(self.minor_axis)


# For utility purposes

NULL_INDEX = Index([])

