# pylint: disable=W0223

from pandas.core.common import _asarray_tuplesafe
from pandas.core.index import Index, MultiIndex
import pandas.core.common as com

import numpy as np

# "null slice"
_NS = slice(None, None)


class IndexingError(Exception):
    pass


class _NDFrameIndexer(object):

    def __init__(self, obj):
        self.obj = obj
        self.ndim = obj.ndim

    def __getitem__(self, key):
        if isinstance(key, tuple):
            try:
                return self.obj.get_value(*key)
            except Exception:
                pass

            return self._getitem_tuple(key)
        else:
            return self._getitem_axis(key, axis=0)

    def _get_label(self, label, axis=0):
        try:
            return self.obj.xs(label, axis=axis, copy=False)
        except Exception:
            return self.obj.xs(label, axis=axis, copy=True)

    def _slice(self, obj, axis=0):
        return self.obj._slice(obj, axis=axis)

    def __setitem__(self, key, value):
        # kludgetastic
        ax = self.obj._get_axis(0)
        if isinstance(ax, MultiIndex):
            try:
                indexer = ax.get_loc(key)
                self._setitem_with_indexer(indexer, value)
                return
            except Exception:
                pass

        if isinstance(key, tuple):
            if len(key) > self.ndim:
                raise IndexingError('only tuples of length <= %d supported',
                                    self.ndim)
            indexer = self._convert_tuple(key)
        else:
            indexer = self._convert_to_indexer(key)

        self._setitem_with_indexer(indexer, value)

    def _convert_tuple(self, key):
        keyidx = []
        for i, k in enumerate(key):
            idx = self._convert_to_indexer(k, axis=i)
            keyidx.append(idx)
        return tuple(keyidx)

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

            plane_indexer = indexer[:het_axis] + indexer[het_axis + 1:]
            item_labels = self.obj._get_axis(het_axis)
            for item in item_labels[het_idx]:
                data = self.obj[item]
                data.values[plane_indexer] = value
        else:
            if isinstance(indexer, tuple):
                indexer = _maybe_convert_ix(*indexer)
            self.obj.values[indexer] = value

    def _getitem_tuple(self, tup):
        try:
            return self._getitem_lowerdim(tup)
        except IndexingError:
            pass

        # no shortcut needed
        retval = self.obj
        for i, key in enumerate(tup):
            if i >= self.obj.ndim:
                raise IndexingError('Too many indexers')

            if _is_null_slice(key):
                continue

            retval = retval.ix._getitem_axis(key, axis=i)

        return retval

    def _getitem_lowerdim(self, tup):
        from pandas.core.frame import DataFrame

        ax0 = self.obj._get_axis(0)
        # a bit kludgy
        if isinstance(ax0, MultiIndex):
            try:
                return self._get_label(tup, axis=0)
            except TypeError:
                # slices are unhashable
                pass
            except Exception:
                if isinstance(tup[0], slice):
                    raise IndexingError
                if tup[0] not in ax0:
                    raise

        # to avoid wasted computation
        # df.ix[d1:d2, 0] -> columns first (True)
        # df.ix[0, ['C', 'B', A']] -> rows first (False)
        for i, key in enumerate(tup):
            if _is_label_like(key):
                section = self._getitem_axis(key, axis=i)

                # might have been a MultiIndex
                if section.ndim == self.ndim:
                    new_key = tup[:i] + (_NS,) + tup[i + 1:]
                    # new_key = tup[:i] + tup[i+1:]
                else:
                    new_key = tup[:i] + tup[i + 1:]

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
        labels = self.obj._get_axis(axis)
        if isinstance(key, slice):
            return self._get_slice_axis(key, axis=axis)
        elif _is_list_like(key) and not (isinstance(key, tuple) and
                                         isinstance(labels, MultiIndex)):

            return self._getitem_iterable(key, axis=axis)
        elif axis == 0:
            is_int_index = _is_integer_index(labels)

            idx = key
            if com.is_integer(key):
                if isinstance(labels, MultiIndex):
                    try:
                        return self._get_label(key, axis=0)
                    except (KeyError, TypeError):
                        if _is_integer_index(self.obj.index.levels[0]):
                            raise

                if not is_int_index:
                    idx = labels[key]

            return self._get_label(idx, axis=0)
        else:
            labels = self.obj._get_axis(axis)
            lab = key
            if com.is_integer(key) and not _is_integer_index(labels):
                lab = labels[key]
            return self._get_label(lab, axis=axis)

    def _getitem_iterable(self, key, axis=0):
        labels = self.obj._get_axis(axis)
        axis_name = self.obj._get_axis_name(axis)

        if com._is_bool_indexer(key):
            key = _check_bool_indexer(labels, key)
            return self.obj.reindex(**{axis_name: labels[np.asarray(key)]})
        else:
            if isinstance(key, Index):
                # want Index objects to pass through untouched
                keyarr = key
            else:
                # asarray can be unsafe, NumPy strings are weird
                keyarr = _asarray_tuplesafe(key)
            if _is_integer_dtype(keyarr) and not _is_integer_index(labels):
                keyarr = labels.take(keyarr)

            return self.obj.reindex(**{axis_name: keyarr})

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
        labels = self.obj._get_axis(axis)
        is_int_index = _is_integer_index(labels)

        if com.is_integer(obj) and not is_int_index:
            return obj

        try:
            return labels.get_loc(obj)
        except (KeyError, TypeError):
            pass

        if isinstance(obj, slice):

            int_slice = _is_index_slice(obj)
            null_slice = obj.start is None and obj.stop is None
            # could have integers in the first level of the MultiIndex
            position_slice = (int_slice
                              and not labels.inferred_type == 'integer'
                              and not isinstance(labels, MultiIndex))

            if null_slice or position_slice:
                slicer = obj
            else:
                try:
                    i, j = labels.slice_locs(obj.start, obj.stop)
                    slicer = slice(i, j, obj.step)
                except Exception:
                    if _is_index_slice(obj):
                        if labels.inferred_type == 'integer':
                            raise
                        slicer = obj
                    else:
                        raise

            return slicer

        elif _is_list_like(obj):
            if com._is_bool_indexer(obj):
                objarr = _check_bool_indexer(labels, obj)
                return objarr
            else:
                objarr = _asarray_tuplesafe(obj)

                # If have integer labels, defer to label-based indexing
                if _is_integer_dtype(objarr) and not is_int_index:
                    return objarr

                indexer = labels.get_indexer(objarr)
                mask = indexer == -1
                if mask.any():
                    raise KeyError('%s not in index' % objarr[mask])

                return indexer
        else:
            return labels.get_loc(obj)

    def _tuplify(self, loc):
        tup = [slice(None, None) for _ in range(self.ndim)]
        tup[0] = loc
        return tuple(tup)

    def _get_slice_axis(self, slice_obj, axis=0):
        obj = self.obj

        axis_name = obj._get_axis_name(axis)
        labels = getattr(obj, axis_name)

        int_slice = _is_index_slice(slice_obj)

        null_slice = slice_obj.start is None and slice_obj.stop is None
        # could have integers in the first level of the MultiIndex
        position_slice = (int_slice and not labels.inferred_type == 'integer'
                          and not isinstance(labels, MultiIndex))
        if null_slice or position_slice:
            slicer = slice_obj
        else:
            try:
                i, j = labels.slice_locs(slice_obj.start, slice_obj.stop)
                slicer = slice(i, j, slice_obj.step)
            except Exception:
                if _is_index_slice(slice_obj):
                    if labels.inferred_type == 'integer':
                        raise
                    slicer = slice_obj
                else:
                    raise

        if not _need_slice(slice_obj):
            return obj

        return self._slice(slicer, axis=axis)


def _is_index_slice(obj):
    def _is_valid_index(x):
        return com.is_integer(x) or com.is_float(x) and np.allclose(x, int(x))

    def _crit(v):
        return v is None or _is_valid_index(v)

    both_none = obj.start is None and obj.stop is None

    return not both_none and (_crit(obj.start) and _crit(obj.stop))


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

    def _get_label(self, key, axis=0):
        return self.obj[key]

    def _slice(self, indexer, axis=0):
        return self.obj._get_values(indexer)

    def _setitem_with_indexer(self, indexer, value):
        self.obj._set_values(indexer, value)


def _check_bool_indexer(ax, key):
    # boolean indexing, need to check that the data are aligned, otherwise
    # disallowed
    result = key
    if _is_series(key) and key.dtype == np.bool_:
        if not key.index.equals(ax):
            result = key.reindex(ax)

    if isinstance(result, np.ndarray) and result.dtype == np.object_:
        mask = com.isnull(result)
        if mask.any():
            raise IndexingError('cannot index with vector containing '
                                'NA / NaN values')

    return result


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


def _is_null_slice(obj):
    return (isinstance(obj, slice) and obj.start is None and
            obj.stop is None and obj.step is None)


def _is_integer_dtype(arr):
    return issubclass(arr.dtype.type, np.integer)


def _is_integer_index(index):
    return index.inferred_type == 'integer'


def _is_label_like(key):
    # select a label or row
    return not isinstance(key, slice) and not _is_list_like(key)


def _is_list_like(obj):
    return np.iterable(obj) and not isinstance(obj, basestring)


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
