# pylint: disable=W0223

from datetime import datetime
from pandas.core.common import _asarray_tuplesafe
from pandas.core.index import Index, MultiIndex, _ensure_index
import pandas.core.common as com
import pandas.lib as lib

import numpy as np

# the supported indexers
def get_indexers_list():

    return [
        ('ix'  ,_NDFrameIndexer),
        ('iloc',_iLocIndexer   ),
        ('loc' ,_LocIndexer    ),
        ('at'  ,_AtIndexer     ),
        ('iat' ,_iAtIndexer    ),
        ]

# "null slice"
_NS = slice(None, None)


class IndexingError(Exception):
    pass


class _NDFrameIndexer(object):

    def __init__(self, obj, name):
        self.obj = obj
        self.ndim = obj.ndim
        self.name = name

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
            return self.obj._xs(label, axis=axis, copy=False)
        except Exception:
            return self.obj._xs(label, axis=axis, copy=True)

    def _get_loc(self, key, axis=0):
        return self.obj._ixs(key, axis=axis)

    def _slice(self, obj, axis=0, raise_on_error=False):
        return self.obj._slice(obj, axis=axis, raise_on_error=raise_on_error)

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

    def _has_valid_tuple(self, key):
        pass

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

            if com.is_integer(het_idx):
                het_idx = [het_idx]

            plane_indexer = indexer[:het_axis] + indexer[het_axis + 1:]
            item_labels = self.obj._get_axis(het_axis)

            def setter(item, v):
                data = self.obj[item]
                values = data.values
                if np.prod(values.shape):
                    result, changed = com._maybe_upcast_indexer(values,plane_indexer,v,dtype=getattr(data,'dtype',None))
                    self.obj[item] = result

            labels = item_labels[het_idx]

            if _is_list_like(value):

                # we have an equal len Frame
                if isinstance(value, DataFrame) and value.ndim > 1:

                    for item in labels:

                        # align to
                        if item in value:
                            v = value[item]
                            v = v.reindex(self.obj[item].index & v.index)
                            setter(item, v.values)
                        else:
                            setter(item, np.nan)

                # we have an equal len ndarray to our labels
                elif isinstance(value, np.ndarray) and value.ndim == 2:
                    if len(labels) != value.shape[1]:
                        raise ValueError('Must have equal len keys and value when'
                                         ' setting with an ndarray')

                    for i, item in enumerate(labels):
                        setter(item, value[:,i])

                # we have an equal len list/ndarray
                elif len(labels) == 1 and (
                    len(self.obj[labels[0]]) == len(value) or len(plane_indexer[0]) == len(value)):
                    setter(labels[0], value)

                # per label values
                else:

                    for item, v in zip(labels, value):
                        setter(item, v)
            else:

                # scalar
                for item in labels:
                    setter(item, value)

        else:
            if isinstance(indexer, tuple):
                indexer = _maybe_convert_ix(*indexer)

            if isinstance(value, Series):
                value = self._align_series(indexer, value)

            if isinstance(value, DataFrame):
                value = self._align_frame(indexer, value)

            # 2096
            values = self.obj.values
            if np.prod(values.shape):
                values[indexer] = value

    def _align_series(self, indexer, ser):
        # indexer to assign Series can be tuple or scalar
        if isinstance(indexer, tuple):
            for i, idx in enumerate(indexer):
                ax = self.obj.axes[i]
                if com._is_sequence(idx) or isinstance(idx, slice):
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
                if com._is_sequence(ix) or isinstance(ix, slice):
                    if idx is None:
                        idx = ax[ix].ravel()
                    elif cols is None:
                        cols = ax[ix].ravel()
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

        # no multi-index, so validate all of the indexers
        self._has_valid_tuple(tup)

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

            retval = getattr(retval,self.name)._getitem_axis(key, axis=i)

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
        from pandas.core.panel4d import Panel4D

        if isinstance(self.obj, DataFrame):
            index = self._convert_for_reindex(tup[0], axis=0)
            columns = self._convert_for_reindex(tup[1], axis=1)
            return self.obj.reindex(index=index, columns=columns)
        elif isinstance(self.obj, Panel4D):
            conv = [self._convert_for_reindex(x, axis=i)
                    for i, x in enumerate(tup)]
            return self.obj.reindex(labels=tup[0], items=tup[1], major=tup[2], minor=tup[3])
        elif isinstance(self.obj, Panel):
            conv = [self._convert_for_reindex(x, axis=i)
                    for i, x in enumerate(tup)]
            return self.obj.reindex(items=tup[0], major=tup[1], minor=tup[2])

    def _convert_for_reindex(self, key, axis=0):
        labels = self.obj._get_axis(axis)

        if com._is_bool_indexer(key):
            key = _check_bool_indexer(labels, key)
            return labels[key]
        else:
            if isinstance(key, Index):
                # want Index objects to pass through untouched
                keyarr = key
            else:
                # asarray can be unsafe, NumPy strings are weird
                keyarr = _asarray_tuplesafe(key)

            if _is_integer_dtype(keyarr) and not _is_integer_index(labels):
                keyarr = com._ensure_platform_int(keyarr)
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
                if isinstance(tup[0], (slice, Index)):
                    raise IndexingError

                # raise the error if we are not sorted
                if not ax0.is_lexsorted_for_tuple(tup):
                    raise e1
                try:
                    loc = ax0.get_loc(tup[0])
                except KeyError:
                    raise e1

        if len(tup) > self.obj.ndim:
            raise IndexingError

        # to avoid wasted computation
        # df.ix[d1:d2, 0] -> columns first (True)
        # df.ix[0, ['C', 'B', A']] -> rows first (False)
        for i, key in enumerate(tup):
            if _is_label_like(key) or isinstance(key, tuple):
                section = self._getitem_axis(key, axis=i)

                # we have yielded a scalar ?
                if not _is_list_like(section):
                    return section

                # might have been a MultiIndex
                elif section.ndim == self.ndim:
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

                return getattr(section,self.name)[new_key]

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
        else:
            if com.is_integer(key):
                if axis == 0 and isinstance(labels, MultiIndex):
                    try:
                        return self._get_label(key, axis=axis)
                    except (KeyError, TypeError):
                        if _is_integer_index(self.obj.index.levels[0]):
                            raise

                if not _is_integer_index(labels):
                    return self._get_loc(key, axis=axis)

            return self._get_label(key, axis=axis)

    def _getitem_iterable(self, key, axis=0):
        labels = self.obj._get_axis(axis)

        def _reindex(keys, level=None):
            try:
                return self.obj.reindex_axis(keys, axis=axis, level=level)
            except AttributeError:
                # Series
                if axis != 0:
                    raise AssertionError('axis must be 0')
                return self.obj.reindex(keys, level=level)

        if com._is_bool_indexer(key):
            key = _check_bool_indexer(labels, key)
            inds, = key.nonzero()
            return self.obj.take(inds, axis=axis, convert=False)
        else:
            if isinstance(key, Index):
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
                        self.obj.take(indexer, axis=axis, convert=True)
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

            if labels.is_unique and Index(keyarr).is_unique:
                return _reindex(keyarr, level=level)
            else:
                indexer, missing = labels.get_indexer_non_unique(keyarr)
                check = indexer != -1
                result = self.obj.take(indexer[check], axis=axis, convert=False)

                # need to merge the result labels and the missing labels
                if len(missing):
                    l = np.arange(len(indexer))

                    missing = com._ensure_platform_int(missing)
                    missing_labels = keyarr.take(missing)
                    missing_indexer = com._ensure_int64(l[~check])
                    cur_labels = result._get_axis(axis).values
                    cur_indexer = com._ensure_int64(l[check])

                    new_labels = np.empty(tuple([len(indexer)]),dtype=object)
                    new_labels[cur_indexer]     = cur_labels
                    new_labels[missing_indexer] = missing_labels
                    new_indexer = (Index(cur_indexer) + Index(missing_indexer)).values
                    new_indexer[missing_indexer] = -1

                    # reindex with the specified axis
                    ndim = self.obj.ndim
                    if axis+1 > ndim:
                        raise AssertionError("invalid indexing error with non-unique index")

                    args = [None] * (2*ndim)
                    args[2*axis] = new_labels
                    args[2*axis+1] = new_indexer

                    result = result._reindex_with_indexers(*args, copy=False, fill_value=np.nan)

                return result

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

            # in case of providing all floats, use label-based indexing
            float_slice = (labels.inferred_type == 'floating'
                           and _is_float_slice(obj))

            # floats that are within tolerance of int used as positions
            int_slice = _is_index_slice(obj)

            null_slice = obj.start is None and obj.stop is None

            # could have integers in the first level of the MultiIndex,
            # in which case we wouldn't want to do position-based slicing
            position_slice = (int_slice
                              and not ltype == 'integer'
                              and not isinstance(labels, MultiIndex)
                              and not float_slice)

            start, stop = obj.start, obj.stop

            # last ditch effort: if we are mixed and have integers
            try:
                if position_slice and 'mixed' in ltype:
                    if start is not None:
                        i = labels.get_loc(start)
                    if stop is not None:
                        j = labels.get_loc(stop)
                    position_slice = False
            except KeyError:
                if ltype == 'mixed-integer-float':
                    raise

            if null_slice or position_slice:
                indexer = obj
            else:
                try:
                    indexer = labels.slice_indexer(start, stop, obj.step)
                except Exception:
                    if _is_index_slice(obj):
                        if ltype == 'integer':
                            raise
                        indexer = obj
                    else:
                        raise

            return indexer

        elif _is_list_like(obj):
            if com._is_bool_indexer(obj):
                obj = _check_bool_indexer(labels, obj)
                inds, = obj.nonzero()
                return inds
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

                    # unique index
                    if labels.is_unique:
                        indexer = check = labels.get_indexer(objarr)

                    # non-unique (dups)
                    else:
                        indexer, missing = labels.get_indexer_non_unique(objarr)
                        check = indexer

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

        if not _need_slice(slice_obj):
            return obj

        labels = obj._get_axis(axis)

        ltype = labels.inferred_type

        # in case of providing all floats, use label-based indexing
        float_slice = (labels.inferred_type == 'floating'
                       and _is_float_slice(slice_obj))

        # floats that are within tolerance of int used as positions
        int_slice = _is_index_slice(slice_obj)

        null_slice = slice_obj.start is None and slice_obj.stop is None

        # could have integers in the first level of the MultiIndex,
        # in which case we wouldn't want to do position-based slicing
        position_slice = (int_slice
                          and not ltype == 'integer'
                          and not isinstance(labels, MultiIndex)
                          and not float_slice)

        start, stop = slice_obj.start, slice_obj.stop

        # last ditch effort: if we are mixed and have integers
        try:
            if position_slice and 'mixed' in ltype:
                if start is not None:
                    i = labels.get_loc(start)
                if stop is not None:
                    j = labels.get_loc(stop)
                position_slice = False
        except KeyError:
            if ltype == 'mixed-integer-float':
                raise

        if null_slice or position_slice:
            indexer = slice_obj
        else:
            try:
                indexer = labels.slice_indexer(start, stop, slice_obj.step)
            except Exception:
                if _is_index_slice(slice_obj):
                    if ltype == 'integer':
                        raise
                    indexer = slice_obj
                else:
                    raise

        if isinstance(indexer, slice):
            return self._slice(indexer, axis=axis)
        else:
            return self.obj.take(indexer, axis=axis)

class _LocationIndexer(_NDFrameIndexer):
    _valid_types = None
    _exception   = Exception

    def _has_valid_type(self, k, axis):
        raise NotImplementedError()

    def _has_valid_tuple(self, key):
        """ check the key for valid keys across my indexer """
        for i, k in enumerate(key):
            if i >= self.obj.ndim:
                raise ValueError('Too many indexers')
            if not self._has_valid_type(k,i):
                raise ValueError("Location based indexing can only have [%s] types" % self._valid_types)

    def __getitem__(self, key):
        if type(key) is tuple:
            return self._getitem_tuple(key)
        else:
            return self._getitem_axis(key, axis=0)

    def _getitem_axis(self, key, axis=0):
        raise NotImplementedError()

    def _getbool_axis(self, key, axis=0):
            labels = self.obj._get_axis(axis)
            key = _check_bool_indexer(labels, key)
            inds, = key.nonzero()
            try:
                return self.obj.take(inds, axis=axis, convert=False)
            except (Exception), detail:
                raise self._exception(detail)
    def _get_slice_axis(self, slice_obj, axis=0):
        """ this is pretty simple as we just have to deal with labels """
        obj = self.obj
        if not _need_slice(slice_obj):
            return obj

        labels = obj._get_axis(axis)
        indexer = labels.slice_indexer(slice_obj.start, slice_obj.stop, slice_obj.step)

        if isinstance(indexer, slice):
            return self._slice(indexer, axis=axis)
        else:
            return self.obj.take(indexer, axis=axis)

class _LocIndexer(_LocationIndexer):
    """ purely label based location based indexing """
    _valid_types = "labels (MUST BE IN THE INDEX), slices of labels (BOTH endpoints included! Can be slices of integers if the index is integers), listlike of labels, boolean"
    _exception   = KeyError

    def _has_valid_type(self, key, axis):
        ax = self.obj._get_axis(axis)

        # valid for a label where all labels are in the index
        # slice of lables (where start-end in labels)
        # slice of integers (only if in the lables)
        # boolean

        if isinstance(key, slice):

            if key.start is not None:
                if key.start not in ax:
                    raise KeyError("start bound [%s] is not the [%s]" % (key.start,self.obj._get_axis_name(axis)))
            if key.stop is not None:
                if key.stop not in ax:
                    raise KeyError("stop bound [%s] is not in the [%s]" % (key.stop,self.obj._get_axis_name(axis)))

        elif com._is_bool_indexer(key):
                return True

        elif _is_list_like(key):

            # require all elements in the index
            idx = _ensure_index(key)
            if not idx.isin(ax).all():
                raise KeyError("[%s] are not in ALL in the [%s]" % (key,self.obj._get_axis_name(axis)))

            return True

        else:

            # if its empty we want a KeyError here
            if not len(ax):
                raise KeyError("The [%s] axis is empty" % self.obj._get_axis_name(axis))

            if not key in ax:
                raise KeyError("the label [%s] is not in the [%s]" % (key,self.obj._get_axis_name(axis)))

        return True

    def _getitem_axis(self, key, axis=0):
        labels = self.obj._get_axis(axis)

        if isinstance(key, slice):
            self._has_valid_type(key,axis)
            return self._get_slice_axis(key, axis=axis)
        elif com._is_bool_indexer(key):
            return self._getbool_axis(key, axis=axis)
        elif _is_list_like(key) and not (isinstance(key, tuple) and
                                         isinstance(labels, MultiIndex)):

            if hasattr(key, 'ndim') and key.ndim > 1:
                raise ValueError('Cannot index with multidimensional key')

            return self._getitem_iterable(key, axis=axis)
        else:
            return self._get_label(key, axis=axis)

class _iLocIndexer(_LocationIndexer):
    """ purely integer based location based indexing """
    _valid_types = "integer, integer slice (START point is INCLUDED, END point is EXCLUDED), listlike of integers, boolean array"
    _exception   = IndexError

    def _has_valid_type(self, key, axis):
        if com._is_bool_indexer(key):
            if hasattr(key,'index') and isinstance(key.index,Index):
                if key.index.inferred_type == 'integer':
                    raise NotImplementedError("iLocation based boolean indexing on an integer type is not available")
                raise ValueError("iLocation based boolean indexing cannot use an indexable as a mask")
            return True

        return isinstance(key, slice) or com.is_integer(key) or _is_list_like(key)

    def _getitem_tuple(self, tup):

        self._has_valid_tuple(tup)
        try:
            return self._getitem_lowerdim(tup)
        except:
            pass

        retval = self.obj
        for i, key in enumerate(tup):
            if i >= self.obj.ndim:
                raise IndexingError('Too many indexers')

            if _is_null_slice(key):
                continue

            retval = getattr(retval,self.name)._getitem_axis(key, axis=i)

        return retval

    def _get_slice_axis(self, slice_obj, axis=0):
        obj = self.obj

        if not _need_slice(slice_obj):
            return obj

        if isinstance(slice_obj, slice):
            return self._slice(slice_obj, axis=axis, raise_on_error=True)
        else:
            return self.obj.take(slice_obj, axis=axis)

    def _getitem_axis(self, key, axis=0):

        if isinstance(key, slice):
            self._has_valid_type(key,axis)
            return self._get_slice_axis(key, axis=axis)

        elif com._is_bool_indexer(key):
            self._has_valid_type(key,axis)
            return self._getbool_axis(key, axis=axis)

        # a single integer or a list of integers
        else:

            if not (com.is_integer(key) or _is_list_like(key)):
                raise ValueError("Cannot index by location index with a non-integer key")

            return self._get_loc(key,axis=axis)

    def _convert_to_indexer(self, obj, axis=0):
        """ much simpler as we only have to deal with our valid types """
        if self._has_valid_type(obj,axis):
            return obj

        raise ValueError("Can only index by location with a [%s]" % self._valid_types)


class _ScalarAccessIndexer(_NDFrameIndexer):
    """ access scalars quickly """

    def _convert_key(self, key):
        return list(key)

    def __getitem__(self, key):
        if not isinstance(key, tuple):

            # we could have a convertible item here (e.g. Timestamp)
            if not _is_list_like(key):
                key = tuple([ key ])
            else:
                raise ValueError('Invalid call for scalar access (getting)!')

        if len(key) != self.obj.ndim:
            raise ValueError('Not enough indexers for scalar access (getting)!')
        key = self._convert_key(key)
        return self.obj.get_value(*key)

    def __setitem__(self, key, value):
        if not isinstance(key, tuple):
            raise ValueError('Invalid call for scalar access (setting)!')
        if len(key) != self.obj.ndim:
            raise ValueError('Not enough indexers for scalar access (setting)!')
        key = self._convert_key(key)
        key.append(value)
        self.obj.set_value(*key)

class _AtIndexer(_ScalarAccessIndexer):
    """ label based scalar accessor """
    pass

class _iAtIndexer(_ScalarAccessIndexer):
    """ integer based scalar accessor """

    def _convert_key(self, key):
        """ require  integer args (and convert to label arguments) """
        ckey = []
        for a, i in zip(self.obj.axes,key):
            if not com.is_integer(i):
                raise ValueError("iAt based indexing can only have integer indexers")
            ckey.append(a[i])
        return ckey

# 32-bit floating point machine epsilon
_eps = np.finfo('f4').eps


def _convert_to_index_sliceable(obj, key):
    """ if we are index sliceable, then return my slicer, otherwise return None """
    idx = obj.index
    if isinstance(key, slice):
        idx_type = idx.inferred_type
        if idx_type == 'floating':
            indexer = obj.ix._convert_to_indexer(key, axis=0)
        elif idx_type == 'integer' or _is_index_slice(key):
            indexer = key
        else:
            indexer = obj.ix._convert_to_indexer(key, axis=0)
        return indexer

    elif isinstance(key, basestring):

        # we are an actual column
        if key in obj._data.items:
            return None

        # we need a timelike key here
        if idx.is_all_dates:
            try:
                return idx._get_string_slice(key)
            except:
                return None

    return None

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

    # this function assumes that com._is_bool_indexer(key) == True

    result = key
    if _is_series(key) and not key.index.equals(ax):
        result = result.reindex(ax)
        mask = com.isnull(result)
        if mask.any():
            raise IndexingError('Unalignable boolean Series key provided')

    # com._is_bool_indexer has already checked for nulls in the case of an
    # object array key, so no check needed here
    result = np.asarray(result, dtype=bool)
    return result

def _is_series(obj):
    from pandas.core.series import Series
    return isinstance(obj, Series)


def _maybe_convert_indices(indices, n):
    """ if we have negative indicies, translate to postive here
        if have indicies that are out-of-bounds, raise an IndexError """
    if isinstance(indices, list):
        indices = np.array(indices)

    mask = indices<0
    if mask.any():
        indices[mask] += n
    mask = (indices>=n) | (indices<0)
    if mask.any():
        raise IndexError("indices are out-of-bounds")
    return indices

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


def _check_slice_bounds(slobj, values):
    l = len(values)
    start = slobj.start
    if start is not None:
        if start < -l or start > l-1:
            raise IndexError("out-of-bounds on slice (start)")
    stop = slobj.stop
    if stop is not None:
        if stop < -l-1 or stop > l:
            raise IndexError("out-of-bounds on slice (end)")

def _maybe_droplevels(index, key):
    # drop levels
    original_index = index
    if isinstance(key, tuple):
        for _ in key:
            try:
                index = index.droplevel(0)
            except:
                # we have dropped too much, so back out
                return original_index
    else:
        index = index.droplevel(0)

    return index
