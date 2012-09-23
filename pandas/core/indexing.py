# pylint: disable=W0223

from pandas.core.common import _asarray_tuplesafe
from pandas.core.index import Index, MultiIndex
import pandas.core.common as com
import pandas.lib as lib

import numpy as np

# "null slice"
_NS = slice(None, None)

def _is_sequence(x):
    try:
        iter(x)
        assert(not isinstance(x, basestring))
        return True
    except Exception:
        return False

class IndexingError(Exception):
    pass


class _NDFrameIndexer(object):

    def __init__(self, obj):
        self.obj = obj
        self.ndim = obj.ndim

    def __iter__(self):
        raise NotImplementedError('ix is not iterable')

    def __getitem__(self, key):
        if type(key) is tuple:
            try:
                return self.obj.get_value(*key)
            except Exception:
                pass

            return self._getitem_tuple(key)
        else:
            return self._getitem_axis(key, axis=0)

    def _get_label(self, label, axis=0):
        # ueber-hack
        if (isinstance(label, tuple) and
            isinstance(label[axis], slice)):

            raise IndexingError('no slices here')

        try:
            return self.obj.xs(label, axis=axis, copy=False)
        except Exception:
            return self.obj.xs(label, axis=axis, copy=True)

    def _get_loc(self, key, axis=0):
        return self.obj._ixs(key, axis=axis)

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
        from pandas.core.frame import DataFrame, Series

        # also has the side effect of consolidating in-place

        # mmm, spaghetti

        if self.obj._is_mixed_type:
            if not isinstance(indexer, tuple):
                indexer = self._tuplify(indexer)

            if isinstance(value, Series):
                value = self._align_series(indexer, value)

            het_axis = self.obj._het_axis
            het_idx = indexer[het_axis]

            if isinstance(het_idx, (int, long)):
                het_idx = [het_idx]

            plane_indexer = indexer[:het_axis] + indexer[het_axis + 1:]
            item_labels = self.obj._get_axis(het_axis)

            if isinstance(value, (np.ndarray, DataFrame)) and value.ndim > 1:
                raise ValueError('Setting mixed-type DataFrames with '
                                 'array/DataFrame pieces not yet supported')

            try:
                for item in item_labels[het_idx]:
                    data = self.obj[item]
                    data.values[plane_indexer] = value
            except ValueError:
                for item, v in zip(item_labels[het_idx], value):
                    data = self.obj[item]
                    data.values[plane_indexer] = v
        else:
            if isinstance(indexer, tuple):
                indexer = _maybe_convert_ix(*indexer)

            if isinstance(value, Series):
                value = self._align_series(indexer, value)

            if isinstance(value, DataFrame):
                value = self._align_frame(indexer, value)

            self.obj.values[indexer] = value

    def _align_series(self, indexer, ser):
        # indexer to assign Series can be tuple or scalar
        if isinstance(indexer, tuple):
            for i, idx in enumerate(indexer):
                ax = self.obj.axes[i]
                if _is_sequence(idx) or isinstance(idx, slice):
                    new_ix = ax[idx]
                    if ser.index.equals(new_ix):
                        return ser.values.copy()
                    return ser.reindex(new_ix).values

        elif np.isscalar(indexer):
            ax = self.obj._get_axis(1)

            if ser.index.equals(ax):
                return ser.values.copy()

            return ser.reindex(ax).values

        raise ValueError('Incompatible indexer with Series')

    def _align_frame(self, indexer, df):
        from pandas import DataFrame
        is_frame = isinstance(self.obj, DataFrame)
        if not is_frame:
            df = df.T
        if isinstance(indexer, tuple):
            idx, cols = None, None
            for i, ix in enumerate(indexer):
                ax = self.obj.axes[i]
                if _is_sequence(ix) or isinstance(ix, slice):
                    if idx is None:
                        idx = ax[ix]
                    elif cols is None:
                        cols = ax[ix]
                    else:
                        break

            if idx is not None and cols is not None:
                if df.index.equals(idx) and df.columns.equals(cols):
                    val = df.copy().values
                else:
                    val = df.reindex(idx, columns=cols).values
                return val

        elif ((isinstance(indexer, slice) or com.is_list_like(indexer))
              and is_frame):
            ax = self.obj.index[indexer]
            if df.index.equals(ax):
                val = df.copy().values
            else:
                val = df.reindex(ax).values
            return val

        elif np.isscalar(indexer) and not is_frame:
            idx = self.obj.axes[1]
            cols = self.obj.axes[2]

            if idx.equals(df.index) and cols.equals(df.columns):
                return df.copy().values
            return df.reindex(idx, columns=cols).values

        raise ValueError('Incompatible indexer with DataFrame')

    def _getitem_tuple(self, tup):
        try:
            return self._getitem_lowerdim(tup)
        except IndexingError:
            pass

        # ugly hack for GH #836
        if self._multi_take_opportunity(tup):
            return self._multi_take(tup)

        # no shortcut needed
        retval = self.obj
        for i, key in enumerate(tup):
            if i >= self.obj.ndim:
                raise IndexingError('Too many indexers')

            if _is_null_slice(key):
                continue

            retval = retval.ix._getitem_axis(key, axis=i)

        return retval

    def _multi_take_opportunity(self, tup):
        from pandas.core.generic import NDFrame

        # ugly hack for GH #836
        if not isinstance(self.obj, NDFrame):
            return False

        if not all(_is_list_like(x) for x in tup):
            return False

        # just too complicated
        for ax in self.obj._data.axes:
            if isinstance(ax, MultiIndex):
                return False

        return True

    def _multi_take(self, tup):
        from pandas.core.frame import DataFrame
        from pandas.core.panel import Panel

        if isinstance(self.obj, DataFrame):
            index = self._convert_for_reindex(tup[0], axis=0)
            columns = self._convert_for_reindex(tup[1], axis=1)
            return self.obj.reindex(index=index, columns=columns)
        elif isinstance(self.obj, Panel):
            conv = [self._convert_for_reindex(x, axis=i)
                    for i, x in enumerate(tup)]
            return self.obj.reindex(items=tup[0], major=tup[1], minor=tup[2])

    def _convert_for_reindex(self, key, axis=0):
        labels = self.obj._get_axis(axis)

        if com._is_bool_indexer(key):
            key = _check_bool_indexer(labels, key)
            return labels[np.asarray(key)]
        else:
            if isinstance(key, Index):
                # want Index objects to pass through untouched
                keyarr = key
            else:
                # asarray can be unsafe, NumPy strings are weird
                keyarr = _asarray_tuplesafe(key)

            if _is_integer_dtype(keyarr) and not _is_integer_index(labels):
                return labels.take(keyarr)

            return keyarr

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
            except Exception, e1:
                if isinstance(tup[0], slice):
                    raise IndexingError
                try:
                    loc = ax0.get_loc(tup[0])
                except KeyError:
                    raise e1

        # to avoid wasted computation
        # df.ix[d1:d2, 0] -> columns first (True)
        # df.ix[0, ['C', 'B', A']] -> rows first (False)
        for i, key in enumerate(tup):
            if _is_label_like(key) or isinstance(key, tuple):
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

            if hasattr(key, 'ndim') and key.ndim > 1:
                raise ValueError('Cannot index with multidimensional key')

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
                    return self._get_loc(key, axis=0)

            return self._get_label(idx, axis=0)
        else:
            labels = self.obj._get_axis(axis)
            lab = key
            if com.is_integer(key) and not _is_integer_index(labels):
                return self._get_loc(key, axis=axis)
            return self._get_label(lab, axis=axis)

    def _getitem_iterable(self, key, axis=0):
        labels = self.obj._get_axis(axis)

        def _reindex(keys, level=None):
            try:
                return self.obj.reindex_axis(keys, axis=axis, level=level)
            except AttributeError:
                # Series
                assert(axis == 0)
                return self.obj.reindex(keys, level=level)

        if com._is_bool_indexer(key):
            key = _check_bool_indexer(labels, key)
            inds, = np.asarray(key, dtype=bool).nonzero()
            return self.obj.take(inds, axis=axis)
        else:
            was_index = isinstance(key, Index)
            if was_index:
                # want Index objects to pass through untouched
                keyarr = key
            else:
                # asarray can be unsafe, NumPy strings are weird
                keyarr = _asarray_tuplesafe(key)

            if _is_integer_dtype(keyarr):
                if labels.inferred_type != 'integer':
                    keyarr = np.where(keyarr < 0,
                                      len(labels) + keyarr, keyarr)

                if labels.inferred_type == 'mixed-integer':
                    indexer = labels.get_indexer(keyarr)
                    if (indexer >= 0).all():
                        self.obj.take(indexer, axis=axis)
                    else:
                        return self.obj.take(keyarr, axis=axis)
                elif not labels.inferred_type == 'integer':

                    return self.obj.take(keyarr, axis=axis)

            # this is not the most robust, but...
            if (isinstance(labels, MultiIndex) and
                not isinstance(keyarr[0], tuple)):
                level = 0
            else:
                level = None

            if labels.is_unique:
                return _reindex(keyarr, level=level)
            else:
                mask = labels.isin(keyarr)
                return self.obj.take(mask.nonzero()[0], axis=axis)

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
            ltype = labels.inferred_type

            if ltype == 'floating':
                int_slice = _is_int_slice(obj)
            else:
                # floats that are within tolerance of int used
                int_slice = _is_index_slice(obj)

            null_slice = obj.start is None and obj.stop is None
            # could have integers in the first level of the MultiIndex
            position_slice = (int_slice
                              and not ltype == 'integer'
                              and not isinstance(labels, MultiIndex))

            start, stop = obj.start, obj.stop

            # last ditch effort: if we are mixed and have integers
            try:
                if 'mixed' in ltype and int_slice:
                    if start is not None:
                        i = labels.get_loc(start)
                    if stop is not None:
                        j = labels.get_loc(stop)
                    position_slice = False
            except KeyError:
                if ltype == 'mixed-integer-float':
                    raise

            if null_slice or position_slice:
                slicer = obj
            else:
                try:
                    i, j = labels.slice_locs(start, stop)
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
                if isinstance(obj, Index):
                    objarr = obj.values
                else:
                    objarr = _asarray_tuplesafe(obj)

                # If have integer labels, defer to label-based indexing
                if _is_integer_dtype(objarr) and not is_int_index:
                    if labels.inferred_type != 'integer':
                        objarr = np.where(objarr < 0,
                                          len(labels) + objarr, objarr)
                    return objarr

                # this is not the most robust, but...
                if (isinstance(labels, MultiIndex) and
                    not isinstance(objarr[0], tuple)):
                    level = 0
                    _, indexer = labels.reindex(objarr, level=level)

                    check = labels.levels[0].get_indexer(objarr)
                else:
                    level = None
                    # XXX
                    if labels.is_unique:
                        indexer = check = labels.get_indexer(objarr)
                    else:
                        mask = np.zeros(len(labels), dtype=bool)
                        lvalues = labels.values
                        for x in objarr:
                            # ugh
                            to_or = lib.map_infer(lvalues, x.__eq__)
                            if not to_or.any():
                                raise KeyError('%s not in index' % str(x))
                            mask |= to_or

                        indexer = check = mask.nonzero()[0]

                mask = check == -1
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

        start = slice_obj.start
        stop = slice_obj.stop

        # in case of providing all floats, use label-based indexing
        float_slice = (labels.inferred_type == 'floating'
                       and _is_float_slice(slice_obj))

        null_slice = slice_obj.start is None and slice_obj.stop is None

        # could have integers in the first level of the MultiIndex, in which
        # case we wouldn't want to do position-based slicing
        position_slice = (int_slice
                          and labels.inferred_type != 'integer'
                          and not isinstance(labels, MultiIndex)
                          and not float_slice)

        # last ditch effort: if we are mixed and have integers
        try:
            if 'mixed' in labels.inferred_type and int_slice:
                if start is not None:
                    i = labels.get_loc(start)
                if stop is not None:
                    j = labels.get_loc(stop)
                position_slice = False
        except KeyError:
            if labels.inferred_type == 'mixed-integer-float':
                raise

        if null_slice or position_slice:
            slicer = slice_obj
        else:
            try:
                i, j = labels.slice_locs(start, stop)
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

# 32-bit floating point machine epsilon
_eps = np.finfo('f4').eps

def _is_index_slice(obj):
    def _is_valid_index(x):
        return (com.is_integer(x) or com.is_float(x)
                and np.allclose(x, int(x), rtol=_eps, atol=0))

    def _crit(v):
        return v is None or _is_valid_index(v)

    both_none = obj.start is None and obj.stop is None

    return not both_none and (_crit(obj.start) and _crit(obj.stop))

def _is_int_slice(obj):
    def _is_valid_index(x):
        return com.is_integer(x)

    def _crit(v):
        return v is None or _is_valid_index(v)

    both_none = obj.start is None and obj.stop is None

    return not both_none and (_crit(obj.start) and _crit(obj.stop))

def _is_float_slice(obj):
    def _is_valid_index(x):
        return com.is_float(x)

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

    def _get_loc(self, key, axis=0):
        return self.obj.values[key]

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
    return (issubclass(arr.dtype.type, np.integer) and
            not arr.dtype.type == np.datetime64)


def _is_integer_index(index):
    return index.inferred_type == 'integer'


def _is_label_like(key):
    # select a label or row
    return not isinstance(key, slice) and not _is_list_like(key)


def _is_list_like(obj):
    # Consider namedtuples to be not list like as they are useful as indices
    return (np.iterable(obj)
            and not isinstance(obj, basestring)
            and not (isinstance(obj, tuple) and type(obj) is not tuple))


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
