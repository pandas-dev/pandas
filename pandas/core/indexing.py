# pylint: disable=W0223

from pandas.core.common import _asarray_tuplesafe
from pandas.core.index import Index, MultiIndex

import numpy as np

# "null slice"
_NS = slice(None, None)


class IndexingError(Exception):
    pass


class AmbiguousIndexError(Exception):
    pass


class _NDFrameIndexer(object):

    def __init__(self, obj):
        self.obj = obj
        self.ndim = obj.ndim

    def __getitem__(self, key):
        if isinstance(key, tuple):
            return self._getitem_tuple(key)
        else:
            return self._getitem_axis(key, axis=0)

    def _getitem_xs(self, idx, axis=0):
        try:
            return self.obj.xs(idx, axis=axis, copy=False)
        except Exception:
            return self.obj.xs(idx, axis=axis, copy=True)

    def __setitem__(self, key, value):
        if isinstance(key, tuple):
            if len(key) > self.ndim:
                raise IndexingError('only tuples of length <= %d supported',
                                    self.ndim)

            keyidx = []
            for i, k in enumerate(key):
                idx = self._convert_to_indexer(k, axis=i)
                keyidx.append(idx)
            indexer = _maybe_convert_ix(*keyidx)
        else:
            indexer = self._convert_to_indexer(key)

        self._setitem_with_indexer(indexer, value)

    def _setitem_with_indexer(self, indexer, value):
        # also has the side effect of consolidating in-place
        if self.obj._is_mixed_type:
            if not isinstance(indexer, tuple):
                indexer = self._tuplify(indexer)

            het_axis = self.obj._het_axis
            het_idx = indexer[het_axis]

            if isinstance(het_idx, (int, long)):
                het_idx = [het_idx]

            if not np.isscalar(value):
                raise IndexingError('setting on mixed-type frames only '
                                    'allowed with scalar values')

            plane_indexer = indexer[:het_axis] + indexer[het_axis+1:]
            item_labels = self.obj._get_axis(het_axis)
            for item in item_labels[het_idx]:
                data = self.obj[item]
                data.values[plane_indexer] = value
        else:
            self.obj.values[indexer] = value

    def _getitem_tuple(self, tup):
        # a bit kludgy
        if isinstance(self.obj._get_axis(0), MultiIndex):
            try:
                return self._getitem_xs(tup, axis=0)
            except (KeyError, TypeError):
                pass

        try:
            return self._getitem_lowerdim(tup)
        except IndexingError:
            pass

        # no shortcut needed
        retval = self.obj
        for i, key in enumerate(tup):
            # hack?
            retval = retval.ix._getitem_axis(key, axis=i)

        return retval

    def _getitem_lowerdim(self, tup):
        from pandas.core.frame import DataFrame

        # to avoid wasted computation
        # df.ix[d1:d2, 0] -> columns first (True)
        # df.ix[0, ['C', 'B', A']] -> rows first (False)
        for i, key in enumerate(tup):
            if _is_label_like(key):
                section = self._getitem_axis(key, axis=i)

                # might have been a MultiIndex
                if section.ndim == self.ndim:
                    new_key = tup[:i] + (_NS,) + tup[i+1:]
                else:
                    new_key = tup[:i] + tup[i+1:]

                    # unfortunately need an odious kludge here because of
                    # DataFrame transposing convention
                    if (isinstance(section, DataFrame) and i > 0
                        and len(new_key) == 2):
                        a, b = new_key
                        new_key = b, a

                    if len(new_key) == 1:
                        new_key, = new_key

                return section.ix[new_key]

        raise IndexingError('not applicable')

    def _getitem_axis(self, key, axis=0):
        if isinstance(key, slice):
            return self._get_slice_axis(key, axis=axis)
        elif _is_list_like(key):
            return self._getitem_iterable(key, axis=axis)
        elif axis == 0:
            labels = self.obj._get_axis(0)
            is_int_index = _is_integer_index(labels)

            idx = key
            if _is_int_like(key):
                if isinstance(labels, MultiIndex):
                    try:
                        return self._getitem_xs(key, axis=0)
                    except (KeyError, TypeError):
                        if _is_integer_index(self.obj.index.levels[0]):
                            raise

                if not is_int_index:
                    idx = labels[key]

            return self._getitem_xs(idx, axis=0)
        else:
            labels = self.obj._get_axis(axis)
            lab = key
            if _is_int_like(key) and not _is_integer_index(labels):
                lab = labels[key]
            return self._getitem_xs(lab, axis=axis)

    def _getitem_iterable(self, key, axis=0):
        labels = self.obj._get_axis(axis)
        axis_name = self.obj._get_axis_name(axis)

        # asarray can be unsafe, NumPy strings are weird
        if isinstance(key, Index):
            # want Index objects to pass through untouched
            keyarr = key
        else:
            keyarr = _asarray_tuplesafe(key)

        if keyarr.dtype == np.bool_:
            if _is_series(key):
                if not key.index.equals(labels):
                    raise IndexingError('Cannot use boolean index with '
                                        'misaligned or unequal labels')
            return self.obj.reindex(**{axis_name : labels[np.asarray(key)]})
        else:
            if _is_integer_dtype(keyarr) and not _is_integer_index(labels):
                keyarr = labels.take(keyarr)

            return self.obj.reindex(**{axis_name : keyarr})

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
        index = self.obj._get_axis(axis)
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
                    raise IndexingError('Cannot use boolean index with '
                                        'misaligned or unequal labels')
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

    def _tuplify(self, loc):
        tup = [slice(None, None) for _ in range(self.ndim)]
        tup[0] = loc
        return tuple(tup)

    def _get_slice_axis(self, slice_obj, axis=0):
        obj = self.obj

        axis_name = obj._get_axis_name(axis)
        labels = getattr(obj, axis_name)
        if _is_label_slice(labels, slice_obj):
            i, j = labels.slice_locs(slice_obj.start, slice_obj.stop)
            slicer = slice(i, j)
        else:
            slicer = slice_obj

        if not _need_slice(slice_obj):
            return obj

        return obj._slice(slicer, axis=axis)

class _SeriesIndexer(_NDFrameIndexer):
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

    def __getitem__(self, key):
        op = self._fancy_index(key, operation='get')
        return op()

    def __setitem__(self, key, value):
        op = self._fancy_index(key, value, operation='set')
        op()

    def _fancy_index(self, key, value=None, operation='get'):
        # going to great lengths to avoid code dup
        obj = self.obj

        if operation == 'get':
            def do_default():
                return obj[key]

            def do_list_like():
                if isinstance(obj.index, MultiIndex):
                    try:
                        return obj[key]
                    except (KeyError, TypeError, IndexError):
                        pass
                return obj.reindex(key)
        else:
            def do_default():
                obj[key] = value

            def do_list_like():
                inds, mask = obj.index.get_indexer(key)
                if not mask.all():
                    raise IndexingError('Indices %s not found' % key[-mask])
                obj.put(inds, value)
        op = do_default
        if _isboolarr(key):
            if _is_series(key):
                if not key.index.equals(obj.index):
                    raise IndexingError('Cannot use boolean index with '
                                        'misaligned or unequal labels')
        elif isinstance(key, slice):
            if _is_label_slice(obj.index, key):
                i, j = obj.index.slice_locs(key.start, key.stop)
                key = slice(i, j)
        elif _is_list_like(key):
            op = do_list_like
        return op

def _is_series(obj):
    from pandas.core.series import Series
    return isinstance(obj, Series)

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
