"""
Data structures for sparse float data. Life is made simpler by dealing only
with float64 data
"""

# pylint: disable=E1101,E1103,W0231

from pandas.compat import range, lrange, zip
from pandas import compat
import numpy as np

from pandas.core.index import Index, MultiIndex, _ensure_index
from pandas.core.frame import DataFrame
from pandas.core.panel import Panel
from pandas.sparse.frame import SparseDataFrame
from pandas.util.decorators import deprecate

import pandas.core.common as com
import pandas.core.ops as ops


class SparsePanelAxis(object):

    def __init__(self, cache_field, frame_attr):
        self.cache_field = cache_field
        self.frame_attr = frame_attr

    def __get__(self, obj, type=None):
        return getattr(obj, self.cache_field, None)

    def __set__(self, obj, value):
        value = _ensure_index(value)

        if isinstance(value, MultiIndex):
            raise NotImplementedError

        for v in compat.itervalues(obj._frames):
            setattr(v, self.frame_attr, value)

        setattr(obj, self.cache_field, value)


class SparsePanel(Panel):

    """
    Sparse version of Panel

    Parameters
    ----------
    frames : dict of DataFrame objects
    items : array-like
    major_axis : array-like
    minor_axis : array-like
    default_kind : {'block', 'integer'}, default 'block'
        Default sparse kind for converting Series to SparseSeries. Will not
        override SparseSeries passed into constructor
    default_fill_value : float
        Default fill_value for converting Series to SparseSeries. Will not
        override SparseSeries passed in

    Notes
    -----
    """
    ndim = 3
    _typ = 'panel'
    _subtyp = 'sparse_panel'

    def __init__(self, frames, items=None, major_axis=None, minor_axis=None,
                 default_fill_value=np.nan, default_kind='block',
                 copy=False):
        if isinstance(frames, np.ndarray):
            new_frames = {}
            for item, vals in zip(items, frames):
                new_frames[item] = \
                    SparseDataFrame(vals, index=major_axis,
                                    columns=minor_axis,
                                    default_fill_value=default_fill_value,
                                    default_kind=default_kind)
            frames = new_frames

        if not isinstance(frames, dict):
            raise TypeError('input must be a dict, a %r was passed' %
                            type(frames).__name__)

        self.default_fill_value = fill_value = default_fill_value
        self.default_kind = kind = default_kind

        # pre-filter, if necessary
        if items is None:
            items = Index(sorted(frames.keys()))
        items = _ensure_index(items)

        (clean_frames,
         major_axis,
         minor_axis) = _convert_frames(frames, major_axis,
                                       minor_axis, kind=kind,
                                       fill_value=fill_value)

        self._frames = clean_frames

        # do we want to fill missing ones?
        for item in items:
            if item not in clean_frames:
                raise ValueError('column %r not found in data' % item)

        self._items = items
        self.major_axis = major_axis
        self.minor_axis = minor_axis

    def _consolidate_inplace(self):  # pragma: no cover
        # do nothing when DataFrame calls this method
        pass

    def __array_wrap__(self, result):
        return SparsePanel(result, items=self.items,
                           major_axis=self.major_axis,
                           minor_axis=self.minor_axis,
                           default_kind=self.default_kind,
                           default_fill_value=self.default_fill_value)

    @classmethod
    def from_dict(cls, data):
        """
        Analogous to Panel.from_dict
        """
        return SparsePanel(data)

    def to_dense(self):
        """
        Convert SparsePanel to (dense) Panel

        Returns
        -------
        dense : Panel
        """
        return Panel(self.values, self.items, self.major_axis,
                     self.minor_axis)

    def as_matrix(self):
        return self.values

    @property
    def values(self):
        # return dense values
        return np.array([self._frames[item].values
                         for item in self.items])

    # need a special property for items to make the field assignable

    _items = None

    def _get_items(self):
        return self._items

    def _set_items(self, new_items):
        new_items = _ensure_index(new_items)
        if isinstance(new_items, MultiIndex):
            raise NotImplementedError

        # need to create new frames dict

        old_frame_dict = self._frames
        old_items = self._items
        self._frames = dict((new_k, old_frame_dict[old_k])
                            for new_k, old_k in zip(new_items, old_items))
        self._items = new_items
    items = property(fget=_get_items, fset=_set_items)

    # DataFrame's index
    major_axis = SparsePanelAxis('_major_axis', 'index')

    # DataFrame's columns / "items"
    minor_axis = SparsePanelAxis('_minor_axis', 'columns')

    def _ixs(self, i, axis=0):
        """
        for compat as we don't support Block Manager here
        i : int, slice, or sequence of integers
        axis : int
        """

        key = self._get_axis(axis)[i]

        # xs cannot handle a non-scalar key, so just reindex here
        if com.is_list_like(key):
            return self.reindex(**{self._get_axis_name(axis): key})

        return self.xs(key, axis=axis)

    def _slice(self, slobj, axis=0, raise_on_error=False, typ=None):
        """
        for compat as we don't support Block Manager here
        """
        axis = self._get_axis_name(axis)
        index = self._get_axis(axis)

        return self.reindex(**{axis: index[slobj]})

    def _get_item_cache(self, key):
        return self._frames[key]

    def __setitem__(self, key, value):
        if isinstance(value, DataFrame):
            value = value.reindex(index=self.major_axis,
                                  columns=self.minor_axis)
            if not isinstance(value, SparseDataFrame):
                value = value.to_sparse(fill_value=self.default_fill_value,
                                        kind=self.default_kind)
        else:
            raise ValueError('only DataFrame objects can be set currently')

        self._frames[key] = value

        if key not in self.items:
            self._items = Index(list(self.items) + [key])

    def set_value(self, item, major, minor, value):
        """
        Quickly set single value at (item, major, minor) location

        Parameters
        ----------
        item : item label (panel item)
        major : major axis label (panel item row)
        minor : minor axis label (panel item column)
        value : scalar

        Notes
        -----
        This method *always* returns a new object. It is not particularly
        efficient but is provided for API compatibility with Panel

        Returns
        -------
        panel : SparsePanel
        """
        dense = self.to_dense().set_value(item, major, minor, value)
        return dense.to_sparse(kind=self.default_kind,
                               fill_value=self.default_fill_value)

    def __delitem__(self, key):
        loc = self.items.get_loc(key)
        indices = lrange(loc) + lrange(loc + 1, len(self.items))
        del self._frames[key]
        self._items = self._items.take(indices)

    def __getstate__(self):
        # pickling
        return (self._frames, com._pickle_array(self.items),
                com._pickle_array(self.major_axis),
                com._pickle_array(self.minor_axis),
                self.default_fill_value, self.default_kind)

    def __setstate__(self, state):
        frames, items, major, minor, fv, kind = state

        self.default_fill_value = fv
        self.default_kind = kind
        self._items = _ensure_index(com._unpickle_array(items))
        self._major_axis = _ensure_index(com._unpickle_array(major))
        self._minor_axis = _ensure_index(com._unpickle_array(minor))
        self._frames = frames

    def copy(self, deep=True):
        """
        Make a copy of the sparse panel

        Returns
        -------
        copy : SparsePanel
        """

        d = self._construct_axes_dict()
        if deep:
            new_data = dict((k, v.copy(deep=True)) for k, v in compat.iteritems(self._frames))
            d = dict((k, v.copy(deep=True)) for k, v in compat.iteritems(d))
        else:
            new_data = self._frames.copy()
        d['default_fill_value']=self.default_fill_value
        d['default_kind']=self.default_kind

        return SparsePanel(new_data, **d)

    def to_frame(self, filter_observations=True):
        """
        Convert SparsePanel to (dense) DataFrame

        Returns
        -------
        frame : DataFrame
        """
        if not filter_observations:
            raise TypeError('filter_observations=False not supported for '
                            'SparsePanel.to_long')

        I, N, K = self.shape
        counts = np.zeros(N * K, dtype=int)

        d_values = {}
        d_indexer = {}

        for item in self.items:
            frame = self[item]

            values, major, minor = _stack_sparse_info(frame)

            # values are stacked column-major
            indexer = minor * N + major
            counts.put(indexer, counts.take(indexer) + 1)  # cuteness

            d_values[item] = values
            d_indexer[item] = indexer

        # have full set of observations for each item
        mask = counts == I

        # for each item, take mask values at index locations for those sparse
        # values, and use that to select values
        values = np.column_stack([d_values[item][mask.take(d_indexer[item])]
                                  for item in self.items])

        inds, = mask.nonzero()

        # still column major
        major_labels = inds % N
        minor_labels = inds // N

        index = MultiIndex(levels=[self.major_axis, self.minor_axis],
                           labels=[major_labels, minor_labels],
                           verify_integrity=False)

        df = DataFrame(values, index=index, columns=self.items)
        return df.sortlevel(level=0)

    to_long = deprecate('to_long', to_frame)
    toLong = deprecate('toLong', to_frame)

    def reindex(self, major=None, items=None, minor=None, major_axis=None,
                minor_axis=None, copy=False):
        """
        Conform / reshape panel axis labels to new input labels

        Parameters
        ----------
        major : array-like, default None
        items : array-like, default None
        minor : array-like, default None
        copy : boolean, default False
            Copy underlying SparseDataFrame objects

        Returns
        -------
        reindexed : SparsePanel
        """
        major = com._mut_exclusive(major=major, major_axis=major_axis)
        minor = com._mut_exclusive(minor=minor, minor_axis=minor_axis)

        if com._all_none(items, major, minor):
            raise ValueError('Must specify at least one axis')

        major = self.major_axis if major is None else major
        minor = self.minor_axis if minor is None else minor

        if items is not None:
            new_frames = {}
            for item in items:
                if item in self._frames:
                    new_frames[item] = self._frames[item]
                else:
                    raise NotImplementedError('Reindexing with new items not yet '
                                              'supported')
        else:
            new_frames = self._frames

        if copy:
            new_frames = dict((k, v.copy()) for k, v in compat.iteritems(new_frames))

        return SparsePanel(new_frames, items=items,
                           major_axis=major,
                           minor_axis=minor,
                           default_fill_value=self.default_fill_value,
                           default_kind=self.default_kind)

    def _combine(self, other, func, axis=0):
        if isinstance(other, DataFrame):
            return self._combineFrame(other, func, axis=axis)
        elif isinstance(other, Panel):
            return self._combinePanel(other, func)
        elif np.isscalar(other):
            new_frames = dict((k, func(v, other))
                              for k, v in compat.iteritems(self))
            return self._new_like(new_frames)

    def _combineFrame(self, other, func, axis=0):
        index, columns = self._get_plane_axes(axis)
        axis = self._get_axis_number(axis)

        other = other.reindex(index=index, columns=columns)

        if axis == 0:
            new_values = func(self.values, other.values)
        elif axis == 1:
            new_values = func(self.values.swapaxes(0, 1), other.values.T)
            new_values = new_values.swapaxes(0, 1)
        elif axis == 2:
            new_values = func(self.values.swapaxes(0, 2), other.values)
            new_values = new_values.swapaxes(0, 2)

        # TODO: make faster!
        new_frames = {}
        for item, item_slice in zip(self.items, new_values):
            old_frame = self[item]
            ofv = old_frame.default_fill_value
            ok = old_frame.default_kind
            new_frames[item] = SparseDataFrame(item_slice,
                                               index=self.major_axis,
                                               columns=self.minor_axis,
                                               default_fill_value=ofv,
                                               default_kind=ok)

        return self._new_like(new_frames)

    def _new_like(self, new_frames):
        return SparsePanel(new_frames, self.items, self.major_axis,
                           self.minor_axis,
                           default_fill_value=self.default_fill_value,
                           default_kind=self.default_kind)

    def _combinePanel(self, other, func):
        items = self.items + other.items
        major = self.major_axis + other.major_axis
        minor = self.minor_axis + other.minor_axis

        # could check that everything's the same size, but forget it

        this = self.reindex(items=items, major=major, minor=minor)
        other = other.reindex(items=items, major=major, minor=minor)

        new_frames = {}
        for item in items:
            new_frames[item] = func(this[item], other[item])

        if not isinstance(other, SparsePanel):
            new_default_fill = self.default_fill_value
        else:
            # maybe unnecessary
            new_default_fill = func(self.default_fill_value,
                                    other.default_fill_value)

        return SparsePanel(new_frames, items, major, minor,
                           default_fill_value=new_default_fill,
                           default_kind=self.default_kind)

    def major_xs(self, key):
        """
        Return slice of panel along major axis

        Parameters
        ----------
        key : object
            Major axis label

        Returns
        -------
        y : DataFrame
            index -> minor axis, columns -> items
        """
        slices = dict((k, v.xs(key)) for k, v in compat.iteritems(self))
        return DataFrame(slices, index=self.minor_axis, columns=self.items)

    def minor_xs(self, key):
        """
        Return slice of panel along minor axis

        Parameters
        ----------
        key : object
            Minor axis label

        Returns
        -------
        y : SparseDataFrame
            index -> major axis, columns -> items
        """
        slices = dict((k, v[key]) for k, v in compat.iteritems(self))
        return SparseDataFrame(slices, index=self.major_axis,
                               columns=self.items,
                               default_fill_value=self.default_fill_value,
                               default_kind=self.default_kind)

    # TODO: allow SparsePanel to work with flex arithmetic.
    # pow and mod only work for scalars for now
    def pow(self, val, *args, **kwargs):
        """wrapper around `__pow__` (only works for scalar values)"""
        return self.__pow__(val)

    def mod(self, val, *args, **kwargs):
        """wrapper around `__mod__` (only works for scalar values"""
        return self.__mod__(val)

# Sparse objects opt out of numexpr
SparsePanel._add_aggregate_operations(use_numexpr=False)
ops.add_special_arithmetic_methods(SparsePanel, use_numexpr=False, **ops.panel_special_funcs)
SparseWidePanel = SparsePanel


def _convert_frames(frames, index, columns, fill_value=np.nan, kind='block'):
    from pandas.core.panel import _get_combined_index
    output = {}
    for item, df in compat.iteritems(frames):
        if not isinstance(df, SparseDataFrame):
            df = SparseDataFrame(df, default_kind=kind,
                                 default_fill_value=fill_value)

        output[item] = df

    if index is None:
        all_indexes = [df.index for df in output.values()]
        index = _get_combined_index(all_indexes)
    if columns is None:
        all_columns = [df.columns for df in output.values()]
        columns = _get_combined_index(all_columns)

    index = _ensure_index(index)
    columns = _ensure_index(columns)

    for item, df in compat.iteritems(output):
        if not (df.index.equals(index) and df.columns.equals(columns)):
            output[item] = df.reindex(index=index, columns=columns)

    return output, index, columns


def _stack_sparse_info(frame):
    lengths = [s.sp_index.npoints for _, s in compat.iteritems(frame)]

    # this is pretty fast
    minor_labels = np.repeat(np.arange(len(frame.columns)), lengths)

    inds_to_concat = []
    vals_to_concat = []
    for col in frame.columns:
        series = frame[col]

        if not np.isnan(series.fill_value):
            raise TypeError('This routine assumes NaN fill value')

        int_index = series.sp_index.to_int_index()
        inds_to_concat.append(int_index.indices)
        vals_to_concat.append(series.sp_values)

    major_labels = np.concatenate(inds_to_concat)
    sparse_values = np.concatenate(vals_to_concat)

    return sparse_values, major_labels, minor_labels
