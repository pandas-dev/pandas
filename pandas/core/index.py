import numpy as np
from pandas.lib.tdates import isAllDates
from pandas.lib.tseries import map_indices


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
            return self.item()

        if len(self) > 0:
            self.indexMap = map_indices(self)

            if hasattr(obj, '_allDates'):
                self._allDates = obj._allDates
            else:
                self._allDates = isAllDates(self)
        else:
            self.indexMap = {}
            self._allDates = False

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
        result = self.view(np.ndarray)[key]
        if isinstance(result, np.ndarray):
            return Index(result)
        else:
            return result

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

        if len(self) != len(other):
            return False

        return self._md5 == other._md5

    def _computeMD5(self):
        import hashlib
        m = hashlib.md5(self.tostring())
        return m.hexdigest()

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
        import bisect

        if date not in self.indexMap:
            loc = bisect.bisect(self, date)
            if loc > 0:
                return self[loc-1]
            else:
                return None
        return date

    def sort(self, *args, **kwargs):
        raise Exception('Tried to sort an Index object, too dangerous to be OK!')

    def __add__(self, other):
        if isinstance(other, Index):
            return self.union(other)
        else:
            return np.ndarray.__add__(self, other)

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
        newElts = filter(lambda x: x not in self.indexMap, other)
        if len(newElts) > 0:
            newSeq = np.concatenate((self, newElts))
            try:
                newSeq = np.unique(newSeq)
            except:
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
        otherArr = np.asarray(other)
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

