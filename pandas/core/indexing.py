from pandas.core.common import _asarray_tuplesafe
from pandas.core.index import MultiIndex

import numpy as np

class _SeriesIndexer(object):
    """
    Class to support fancy indexing, potentially using labels

    Notes
    -----
    Indexing based on labels is INCLUSIVE
    Slicing uses PYTHON SEMANTICS (endpoint is excluded)

    If Index contains int labels, these will be used rather than the locations,
    so be very careful (ambiguous).

    Examples
    --------
    >>> ts.ix[5:10] # equivalent to ts[5:10]
    >>> ts.ix[[date1, date2, date3]]
    >>> ts.ix[date1:date2] = 0
    """
    def __init__(self, series):
        self.series = series

    def __getitem__(self, key):
        op = self._fancy_index(key, operation='get')
        return op()

    def __setitem__(self, key, value):
        op = self._fancy_index(key, value, operation='set')
        op()

    def _fancy_index(self, key, value=None, operation='get'):
        from pandas.core.series import Series

        # going to great lengths to avoid code dup
        series = self.series

        if operation == 'get':
            def do_default():
                return series[key]

            def do_list_like():
                if isinstance(series.index, MultiIndex):
                    try:
                        return series[key]
                    except (KeyError, TypeError, IndexError):
                        pass
                return series.reindex(key)
        else:
            def do_default():
                series[key] = value

            def do_list_like():
                inds, mask = series.index.get_indexer(key)
                if not mask.all():
                    raise Exception('Indices %s not found' % key[-mask])
                series.put(inds, value)
        op = do_default
        if _isboolarr(key):
            if isinstance(key, Series):
                if not key.index.equals(series.index):
                    raise Exception('Cannot use boolean index with misaligned '
                                    'or unequal labels')
        elif isinstance(key, slice):
            if _is_label_slice(series.index, key):
                i, j = series.index.slice_locs(key.start, key.stop)
                key = slice(i, j)
        elif _is_list_like(key):
            op = do_list_like
        return op


class AmbiguousIndexError(Exception):
    pass

class _DataFrameIndexer(object):
    """
    Class to support fancy indexing, potentially using labels of DataFrame

    Notes
    -----
    Indexing based on labels is INCLUSIVE
    Slicing uses PYTHON SEMANTICS (endpoint is excluded)

    If Index contains int labels, these will be used rather than the locations,
    so be very careful (ambiguous).

    Examples
    --------
    >>> frame.ix[5:10, ['A', 'B']]
    >>> frame.ix[date1:date2, 'A']
    """

    def __init__(self, frame):
        self.frame = frame

    def __getitem__(self, key):
        if isinstance(key, slice):
            return self._fancy_getitem_axis(key, axis=0)
        elif isinstance(key, tuple):
            return self._getitem_tuple(key)
        elif _is_list_like(key):
            return self._fancy_getitem_iterable(key, axis=0)
        else:
            return self._fancy_getitem_axis(key, axis=0)

    def _getitem_tuple(self, key):
        if isinstance(self.frame.index, MultiIndex):
            try:
                return self.frame.xs(key)
            except (KeyError, TypeError):
                # could do something more intelligent here? like raising the
                # exception if each tuple value are in the levels?
                pass

        if len(key) != 2:
            raise Exception('only length 2 tuple supported')
        return self._fancy_getitem_tuple(*key)

    def __setitem__(self, key, value):
        # also has the side effect of consolidating in-place
        if self.frame._is_mixed_type:
            raise Exception('setting on mixed-type frames not yet supported')

        if isinstance(key, tuple):
            if len(key) != 2:
                raise Exception('only length 2 tuple supported')

            x, y = key
            xidx = self._convert_to_indexer(x, axis=0)
            yidx = self._convert_to_indexer(y, axis=1)
            indexer = _maybe_convert_ix(xidx, yidx)
        else:
            indexer = self._convert_to_indexer(key)

        self.frame.values[indexer] = value

    def _convert_to_indexer(self, obj, axis=0):
        """
        Convert indexing key into something we can use to do actual fancy
        indexing on an ndarray

        Examples
        ix[:5] -> slice(0, 5)
        ix[[1,2,3]] -> [1,2,3]
        ix[['foo', 'bar', 'baz']] -> [i, j, k] (indices of foo, bar, baz)

        Going by Zen of Python?
        "In the face of ambiguity, refuse the temptation to guess."
        raise AmbiguousIndexError with integer labels?
        - No, prefer label-based indexing
        """
        index = self.frame._get_axis(axis)
        is_int_index = _is_integer_index(index)
        if isinstance(obj, slice):
            if _is_label_slice(index, obj):
                i, j = index.slice_locs(obj.start, obj.stop)
                return slice(i, j)
            else:
                return obj
        elif _is_list_like(obj):
            objarr = _asarray_tuplesafe(obj)

            if objarr.dtype == np.bool_:
                if not obj.index.equals(index):
                    raise Exception('Cannot use boolean index with misaligned '
                                    'or unequal labels')
                return objarr
            else:
                # If have integer labels, defer to label-based indexing
                if _is_integer_dtype(objarr) and not is_int_index:
                    return objarr

                indexer, mask = index.get_indexer(objarr)
                if not mask.all():
                    raise KeyError('%s not in index' % objarr[-mask])

                return indexer
        else:
            if _is_int_like(obj) and not is_int_index:
                return obj
            return index.get_loc(obj)

    def _fancy_getitem_tuple(self, rowkey, colkey):
        # to avoid wasted computation
        # df.ix[d1:d2, 0] -> columns first (True)
        # df.ix[0, ['C', 'B', A']] -> rows first (False)
        if _is_label_like(colkey):
            return self._fancy_getitem_axis(colkey, axis=1).ix[rowkey]
        elif _is_label_like(rowkey):
            return self._fancy_getitem_axis(rowkey, axis=0).ix[colkey]

        result = self._fancy_getitem_axis(colkey, axis=1)
        return result.ix[rowkey]

    def _fancy_getitem_axis(self, key, axis=0):
        if isinstance(key, slice):
            return self._get_slice_axis(key, axis=axis)
        elif _is_list_like(key):
            return self._fancy_getitem_iterable(key, axis=axis)
        elif axis == 0:
            is_int_index = _is_integer_index(self.frame.index)

            idx = key
            if _is_int_like(key) and not is_int_index:
                idx = self.frame.index[key]

            if self.frame._is_mixed_type:
                return self.frame.xs(idx)
            else:
                # get a view if possible
                return self.frame.xs(idx, copy=False)
        else:
            col = key
            if _is_int_like(key) and not _is_integer_index(self.frame.columns):
                col = self.frame.columns[key]
            return self.frame[col]

    def _fancy_getitem_iterable(self, key, axis=0):
        from pandas.core.series import Series

        labels = self.frame._get_axis(axis)
        axis_name = self.frame._get_axis_name(axis)

        # asarray can be unsafe, NumPy strings are weird
        keyarr = _asarray_tuplesafe(key)

        if keyarr.dtype == np.bool_:
            if isinstance(key, Series):
                if not key.index.equals(labels):
                    raise Exception('Cannot use boolean index with misaligned '
                                    'or unequal labels')
            return self.frame.reindex(**{axis_name : labels[np.asarray(key)]})
        else:
            if _is_integer_dtype(keyarr) and not _is_integer_index(labels):
                keyarr = labels.take(keyarr)

            return self.frame.reindex(**{axis_name : keyarr})

    def _get_slice_axis(self, slice_obj, axis=0):
        frame = self.frame

        axis_name = frame._get_axis_name(axis)
        labels = getattr(frame, axis_name)
        if _is_label_slice(labels, slice_obj):
            i, j = labels.slice_locs(slice_obj.start, slice_obj.stop)
            slicer = slice(i, j)
        else:
            slicer = slice_obj

        if not _need_slice(slice_obj):
            return frame

        return frame._slice(slicer, axis=axis)

def _maybe_convert_ix(*args):
    """
    We likely want to take the cross-product
    """
    ixify = True
    for arg in args:
        if not isinstance(arg, (np.ndarray, list)):
            ixify = False

    if ixify:
        return np.ix_(*args)
    else:
        return args

def _is_integer_dtype(arr):
    return issubclass(arr.dtype.type, np.integer)

def _is_integer_index(index):
    # make an educated and not too intelligent guess
    if len(index) == 0: # pragma: no cover
        return False
    else:
        return _is_int_like(index[0])

def _is_int_like(val):
    return isinstance(val, (int, np.integer))

def _is_label_like(key):
    # select a label or row
    return not isinstance(key, slice) and not _is_list_like(key)

def _is_list_like(obj):
    return np.iterable(obj) and not isinstance(obj, basestring)

def _is_label_slice(labels, obj):
    def crit(x):
        if x in labels:
            return False
        else:
            return isinstance(x, int) or x is None
    return not crit(obj.start) or not crit(obj.stop)

def _need_slice(obj):
    return (obj.start is not None or
            obj.stop is not None or
            (obj.step is not None and obj.step != 1))

def _maybe_droplevels(index, key):
    # drop levels
    if isinstance(key, tuple):
        for _ in key:
            index = index.droplevel(0)
    else:
        index = index.droplevel(0)

    return index

_isboolarr = lambda x: np.asarray(x).dtype == np.bool_

