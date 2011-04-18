import numpy as np

import operator

from pandas.core.series import Series
from pandas.core.frame import DataFrame

import pandas.lib.sparse as splib

def make_sparse(series, kind='block', sparse_value=np.NaN):
    """

    kind : {'block', 'integer'}
    """

class SparseSeries(Series):
    """
    Data structure for labeled, sparse floating point data
    """

    def __new__(cls, data, index=None, copy=False):
        if isinstance(data, Series):
            if index is None:
                index = data.index
        elif isinstance(data, dict):
            if index is None:
                index = Index(sorted(data.keys()))
            data = [data[idx] for idx in index]

        # Create array, do *not* copy data by default, infer type
        try:
            subarr = np.array(data, dtype=dtype, copy=copy)
        except ValueError:
            if dtype:
                raise

            subarr = np.array(data, dtype=object)

        if subarr.ndim == 0:
            if isinstance(data, list): # pragma: no cover
                subarr = np.array(data, dtype=object)
            elif index is not None:
                value = data

                # If we create an empty array using a string to infer
                # the dtype, NumPy will only allocate one character per entry
                # so this is kind of bad. Alternately we could use np.repeat
                # instead of np.empty (but then you still don't want things
                # coming out as np.str_!
                if isinstance(value, basestring) and dtype is None:
                    dtype = np.object_

                if dtype is None:
                    subarr = np.empty(len(index), dtype=type(value))
                else:
                    subarr = np.empty(len(index), dtype=dtype)
                subarr.fill(value)
            else:
                return subarr.item()

        elif subarr.ndim > 1:
            raise Exception('Data must be 1-dimensional')

        if index is None:
            raise Exception('Index cannot be None!')

        # This is to prevent mixed-type Series getting all casted to
        # NumPy string type, e.g. NaN --> '-1#IND'.
        if issubclass(subarr.dtype.type, basestring):
            subarr = np.array(data, dtype=object, copy=copy)

        # Change the class of the array to be the subclass type.
        subarr = subarr.view(cls)
        subarr.index = index

        if subarr.index._allDates:
            subarr = subarr.view(TimeSeries)

        return subarr

    def __hash__(self):
        raise TypeError('unhashable type')

    _index = None
    def _get_index(self):
        return self._index

    def _set_index(self, index):
        indexTypes = ndarray, Index, list, tuple
        if not isinstance(index, indexTypes):
            raise TypeError("Expected index to be in %s; was %s."
                            % (indexTypes, type(index)))

        if len(self) != len(index):
            raise AssertionError('Lengths of index and values did not match!')

        if not isinstance(index, Index):
            index = Index(index)

        self._index = index

    index = property(fget=_get_index, fset=_set_index)

    def __array_finalize__(self, obj):
        """
        Gets called after any ufunc or other array operations, necessary
        to pass on the index.
        """
        self._index = getattr(obj, '_index', None)

    def toDict(self):
        return dict(self.iteritems())

    @classmethod
    def fromValue(cls, value=np.NaN, index=None, dtype=None): # pragma: no cover
        warnings.warn("'fromValue', can call Series(value, index=index) now",
                      FutureWarning)

        return Series(value, index=index, dtype=dtype)

    def __contains__(self, key):
        return key in self.index

    def __reduce__(self):
        """Necessary for making this object picklable"""
        object_state = list(ndarray.__reduce__(self))
        subclass_state = (self.index, )
        object_state[2] = (object_state[2], subclass_state)
        return tuple(object_state)

    def __setstate__(self, state):
        """Necessary for making this object picklable"""
        nd_state, own_state = state
        ndarray.__setstate__(self, nd_state)
        index, = own_state
        self.index = index

    def __getitem__(self, key):
        """
        Returns item(s) for requested index/sequence, overrides default behavior
        for series[key].

        Logic is as follows:
            - If key is in the index, return the value corresponding
              to that index
            - Otherwise, use key (presumably one integer or a sequence
              of integers) to obtain values from the series. In the case
              of a sequence, a 'slice' of the series (with corresponding dates)
              will be returned, otherwise a single value.
        """
        values = self.values

        try:
            # Check that we can even look for this in the index
            return values[self.index.indexMap[key]]
        except KeyError:
            if isinstance(key, (int, np.integer)):
                return values[key]
            raise Exception('Requested index not in this series!')
        except TypeError:
            # Could not hash item
            pass

        # is there a case where this would NOT be an ndarray?
        # need to find an example, I took out the case for now

        dataSlice = values[key]
        indices = Index(self.index.view(ndarray)[key])
        return Series(dataSlice, index=indices)

    def get(self, key, default=None):
        """
        Returns value occupying requested index, default to specified
        missing value if not present

        Parameters
        ----------
        key : object
            Index value looking for
        default : object, optional
            Value to return if key not in index

        Returns
        -------
        y : scalar
        """
        if key in self.index:
            return ndarray.__getitem__(self, self.index.indexMap[key])
        else:
            return default

    def __getslice__(self, i, j):
        """
        Returns a slice of the Series.

        Note that the underlying values are COPIES.

        The reason that the getslice returns copies is that otherwise you
        will have a reference to the original series which could be
        inadvertently changed if the slice were altered (made mutable).
        """
        newArr = self.values[i:j].copy()
        newIndex = self.index[i:j]

        return Series(newArr, index=newIndex)

    def __setitem__(self, key, value):
        """
        If this series is mutable, set specified indices equal to given values.
        """
        try:
            loc = self.index.indexMap[key]
            ndarray.__setitem__(self, loc, value)
        except Exception:
            values = self.values
            values[key] = value

    def __setslice__(self, i, j, value):
        """Set slice equal to given value(s)"""
        ndarray.__setslice__(self, i, j, value)

    def __repr__(self):
        """Clean string representation of a Series"""
        vals = self.values
        index = self.index

        if len(index) > 500:
            head = _seriesRepr(index[:50], vals[:50])
            tail = _seriesRepr(index[-50:], vals[-50:])
            return head + '\n...\n' + tail + '\nlength: %d' % len(vals)
        elif len(index) > 0:
            return _seriesRepr(index, vals)
        else:
            return '%s' % ndarray.__repr__(self)

    def toString(self, buffer=sys.stdout, nanRep='NaN'):
        print >> buffer, _seriesRepr(self.index, self.values,
                                     nanRep=nanRep)

    def __str__(self):
        return repr(self)

    def __iter__(self):
        return iter(self.values)

    def to_dense(self):
        pass

    def copy(self):
        vec_copy = self._vector.copy()
        return SparseSeries(vec_copy, index=self.index)

class SparseDataFrame(DataFrame):
    pass
