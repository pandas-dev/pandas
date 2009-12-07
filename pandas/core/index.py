# pylint: disable-msg=E1101
# pylint: disable-msg=E1103
# pylint: disable-msg=W0232

import numpy as np
from pandas.lib.tdates import isAllDates
from pandas.lib.tseries import map_indices

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
    """Extension of numpy-array to represent a series index,
    dates or otherwise.

    Index is immutable always (don't even try to change elements!).

    Note that the Index can ONLY contain immutable objects. Mutable objects are not
    hashable, and that's bad!
    """
    __md5 = None
    def __new__(cls, data, dtype=object, copy=False):
        subarr = np.array(data, dtype=dtype, copy=copy)

        if subarr.ndim == 0:
            raise Exception('Index(...) must be called with a collection '
                            'of some kind, %s was passed' % repr(data))

        subarr = subarr.view(cls)
        return subarr

    def __array_finalize__(self, obj):
        if self.ndim == 0:
            # tolist will cause a bus error if this is not here, hmm
            
            return self.item()

        
            # raise Exception('Cannot create 0-dimensional Index!')

        # New instance creation
        if obj is None:
            self.indexMap = map_indices(self)
            self._allDates = isAllDates(self)

        # New from template / slicing
        elif isinstance(obj, type(self)) and len(self) != len(obj.indexMap):
            self.indexMap = map_indices(self)
            self._allDates = isAllDates(self)

        # View casting
        else:
            self.indexMap = getattr(obj, 'indexMap', map_indices(self))
            self._allDates = getattr(obj, '_allDates', isAllDates(self))

        self._checkForDuplicates()

    def _checkForDuplicates(self):
        if len(self.indexMap) < len(self):
            raise Exception('Index cannot contain duplicate values!')

    def __iter__(self):
        return iter(self.view(np.ndarray))

    def __setstate__(self,state):
        """Necessary for making this object picklable"""
        np.ndarray.__setstate__(self,state)
        self.indexMap = map_indices(self)
        self._allDates = isAllDates(self)

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
        if np.isscalar(key):
            return np.ndarray.__getitem__(self, key)
        else:
            return Index(self.view(np.ndarray)[key])

    def equals(self, other):
        """
        Determines if two Index objects contain the same elements.

        If the compared object is of the right type and length, the MD5
        checksum is compared
        """
        if self is other:
            return True

        if not isinstance(other, Index):
            return False

        return np.array_equal(self, other)

    def _computeMD5(self):
        import hashlib
        return hashlib.md5(self.data).hexdigest()

    @property
    def _md5(self):
        """
        Return MD5 hex-digested hash for the Index elements. Note that
        this quantity is only computed once.
        """
        if self.__md5 is None:
            self.__md5 = self._computeMD5()

        return self.__md5

    def asOfDate(self, date):
        if date not in self.indexMap:
            loc = self.searchsorted(date, side='left')
            if loc > 0:
                return self[loc-1]
            else:
                return None

        return date

    def sort(self, *args, **kwargs):
        raise Exception('Tried to sort an Index object, too dangerous to be OK!')

    def argsort(self, *args, **kwargs):
        return self.view(np.ndarray).argsort(*args, **kwargs)

    def __add__(self, other):
        if isinstance(other, Index):
            return self.union(other)
        else:
            return np.ndarray.__add__(self, other)

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
        other: Index or array-like

        Returns
        -------
        Index
        """
        if not hasattr(other, '__iter__'):
            raise Exception('Input must be iterable!')

        if other is self:
            return self

        if isinstance(other, Index):
            if self.equals(other):
                return self

        f = self.indexMap.__contains__
        newElts = [x for x in other if not f(x)]
        if len(newElts) > 0:
            newSeq = np.concatenate((self, newElts))
            try:
                newSeq = np.unique(newSeq)
            except Exception:
                # Not sortable / multiple types
                pass
            return Index(newSeq)
        else:
            return self

    def intersection(self, other):
        """
        Form the intersection of two Index objects and sorts if possible

        Parameters
        ----------
        other: Index or array-like

        Returns
        -------
        Index
        """
        if not hasattr(other, '__iter__'):
            raise Exception('Input must be iterable!')

        if other is self:
            return self

        theIntersection = sorted(set(self) & set(other))
        return Index(theIntersection)

    def diff(self, other):
        if not hasattr(other, '__iter__'):
            raise Exception('Input must be iterable!')

        if other is self:
            return Index([])

        otherArr = np.asarray(other)
        theDiff = sorted(set(self) - set(otherArr))
        return Index(theDiff)

    __sub__ = diff

# For utility purposes

NULL_INDEX = Index([])

