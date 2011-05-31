# pylint: disable=E1101,E1103,W0232

import numpy as np
import pandas.lib.tseries as _tseries

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
        subarr = np.array(data, dtype=dtype, copy=copy)

        if subarr.ndim == 0:
            raise Exception('Index(...) must be called with a collection '
                            'of some kind, %s was passed' % repr(data))

        return subarr.view(cls)

    def __array_finalize__(self, obj):
        if self.ndim == 0:
            # tolist will cause a bus error if this is not here, hmm
            return self.item()
            # raise Exception('Cannot create 0-dimensional Index!')

        # New instance creation
        if obj is None:
            pass

        # New from template / slicing
        elif isinstance(obj, type(self)) and len(self) != len(obj.indexMap):
            pass

        # View casting
        else:
            if hasattr(obj, '_cache_indexMap'):
                self._cache_indexMap = obj._cache_indexMap
                self._cache_allDates = getattr(obj, '_cache_allDates', None)

        self._checkForDuplicates()

    @property
    def indexMap(self):
        if not hasattr(self, '_cache_indexMap'):
            self._cache_indexMap = _tseries.map_indices(self)

        return self._cache_indexMap

    @property
    def _allDates(self):
        if not hasattr(self, '_cache_allDates'):
            self._cache_allDates = _tseries.isAllDates(self)

        return self._cache_allDates

    def is_all_dates(self):
        return self._allDates

    def _checkForDuplicates(self):
        if len(self.indexMap) < len(self):
            raise Exception('Index cannot contain duplicate values!')

    def __iter__(self):
        return iter(self.view(np.ndarray))

    def __setstate__(self, state):
        """Necessary for making this object picklable"""
        np.ndarray.__setstate__(self, state)
        self._cache_indexMap = _tseries.map_indices(self)
        self._cache_allDates = _tseries.isAllDates(self)

    def __deepcopy__(self, memo={}):
        """
        Index is not mutable, so disabling deepcopy
        """
        return self

    def __contains__(self, date):
        return date in self.indexMap

    def __setitem__(self, key, value):
        """Disable the setting of values."""
        raise Exception(str(self.__class__) + ' object is immutable' )

    def __getitem__(self, key):
        """Override numpy.ndarray's __getitem__ method to work as desired"""
        arr_idx = self.view(np.ndarray)
        if np.isscalar(key):
            return arr_idx[key]
        else:
            # easier to ask forgiveness than permission
            try:
                return Index(arr_idx[key])
            except Exception, e1:
                try:
                    return Index(arr_idx[np.asarray(key)])
                except Exception, e2: # pragma: no cover
                    raise e1

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
        import pandas.lib.tseries as tseries

        if method:
            method = method.upper()

        aliases = {
            'FFILL' : 'PAD',
            'BFILL' : 'BACKFILL'
        }

        method = aliases.get(method, method)
        indexer, mask = tseries.getFillVec(self, target, self.indexMap,
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

# For utility purposes

NULL_INDEX = Index([])

