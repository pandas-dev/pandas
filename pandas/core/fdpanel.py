""" FDPanel: a 4-d dict like collection of panels """

import operator
import sys
import numpy as np

from pandas.core.common import (PandasError, _mut_exclusive,
                                _try_sort, _default_index, _infer_dtype)
from pandas.core.index import (Index, MultiIndex, _ensure_index,
                               _get_combined_index)
from pandas.core.indexing import _NDFrameIndexer, _maybe_droplevels, _is_null_slice
from pandas.core.internals import BlockManager, make_block, form_blocks
from pandas.core.frame import DataFrame
from pandas.core.generic import NDFrame
from pandas.core.panel import Panel
from pandas.util import py3compat
from pandas.util.decorators import deprecate, Appender, Substitution
import pandas.core.common as com
import pandas.core.nanops as nanops
import pandas.lib as lib

class FDPanel(Panel):
    _AXIS_NUMBERS = {
        'labels' : 0,
        'items'  : 1,
        'major_axis' : 2,
        'minor_axis' : 3
        }

    _AXIS_ALIASES = {
        'major' : 'major_axis',
        'minor' : 'minor_axis'
    }

    _AXIS_NAMES = {
        0 : 'labels',
        1 : 'items',
        2 : 'major_axis',
        3 : 'minor_axis'
    }

    # major
    _default_stat_axis = 2

    labels     = lib.AxisProperty(0)
    items      = lib.AxisProperty(1)
    major_axis = lib.AxisProperty(2)
    minor_axis = lib.AxisProperty(3)

    def __init__(self, data=None, labels=None, items=None, major_axis=None, minor_axis=None, copy=False, dtype=None):
        """
        Represents a 4 dimensonal structured

        Parameters
        ----------
        data : ndarray (labels x items x major x minor), or dict of Panels

        labels : Index or array-like : axis=0
        items  : Index or array-like : axis=1
        major_axis : Index or array-like: axis=2
        minor_axis : Index or array-like: axis=3

        dtype : dtype, default None
            Data type to force, otherwise infer
        copy : boolean, default False
            Copy data from inputs. Only affects DataFrame / 2d ndarray input
        """
        if data is None:
            data = {}

        passed_axes = [labels,items, major_axis, minor_axis]
        axes = None
        if isinstance(data, BlockManager):
            if any(x is not None for x in passed_axes):
                axes = [x if x is not None else y
                        for x, y in zip(passed_axes, data.axes)]
            mgr = data
        elif isinstance(data, dict):
            mgr = self._init_dict(data, passed_axes, dtype=dtype)
            copy = False
            dtype = None
        elif isinstance(data, (np.ndarray, list)):
            mgr = self._init_matrix(data, passed_axes, dtype=dtype, copy=copy)
            copy = False
            dtype = None
        else: # pragma: no cover
            raise PandasError('FDPanel constructor not properly called!')

        NDFrame.__init__(self, mgr, axes=axes, copy=copy, dtype=dtype)

    @classmethod
    def from_dict(cls, data, intersect=False, orient='items', dtype=None):
        """ not supporting intersect/orient arguments """
        return cls(data, dtype = dtype)

    def _init_dict(self, data, axes, dtype=None):
        labels, items, major, minor = axes

        # prefilter if labels passed
        if labels is not None:
            labels = _ensure_index(labels)
            data = dict((k, v) for k, v in data.iteritems() if k in labels)
        else:
            labels = Index(_try_sort(data.keys()))

        for k, v in data.iteritems():
            if isinstance(v, dict):
                data[k] = Panel(v)
                
        if items is None:
            items = _extract_axis(data, axis=0)
            
        if major is None:
            major = _extract_axis(data, axis=1)

        if minor is None:
            minor = _extract_axis(data, axis=2)

        axes = [labels, items, major, minor]
        reshaped_data = data.copy() # shallow

        label_shape = len(items), len(major), len(minor)
        for label in labels:
            v = values = data.get(label)
            if v is None:
                values = np.empty(item_shape, dtype=dtype)
                values.fill(np.nan)
            elif isinstance(v, Panel):
                v = v.reindex(items=items, major_axis=major, minor_axis=minor, copy=False)
                if dtype is not None:
                    v = v.astype(dtype)
                values = v.values
            reshaped_data[label] = values

        # segregates dtypes and forms blocks matching to columns
        blocks = form_blocks(reshaped_data, axes)
        mgr = BlockManager(blocks, axes).consolidate()
        return mgr

    def _init_matrix(self, data, axes, dtype=None, copy=False):
        values = _prep_ndarray(data, copy=copy)

        if dtype is not None:
            try:
                values = values.astype(dtype)
            except Exception:
                raise ValueError('failed to cast to %s' % dtype)

        shape = values.shape
        fixed_axes = []
        for i, ax in enumerate(axes):
            if ax is None:
                ax = _default_index(shape[i])
            else:
                ax = _ensure_index(ax)
            fixed_axes.append(ax)

        items = fixed_axes[0]
        block = make_block(values, items, items)
        return BlockManager([block], fixed_axes)

    @property
    def shape(self):
        return len(self.labels), len(self.items), len(self.major_axis), len(self.minor_axis)

    def __array_wrap__(self, result):
        return self._constructor(result, 
                                 labels    =self.labels,
                                 items     =self.items,
                                 major_axis=self.major_axis,
                                 minor_axis=self.minor_axis, copy=False)
    
    #----------------------------------------------------------------------
    # Magic methods

    def __repr__(self):
        class_name = str(self.__class__)

        L, I, N, K = len(self.labels), len(self.items), len(self.major_axis), len(self.minor_axis)

        dims = 'Dimensions: %d (labels) x %d (items) x %d (major) x %d (minor)' % (L, I, N, K)

        if len(self.major_axis) > 0:
            major = 'Major axis: %s to %s' % (self.major_axis[0],
                                              self.major_axis[-1])
        else:
            major = 'Major axis: None'

        if len(self.minor_axis) > 0:
            minor = 'Minor axis: %s to %s' % (self.minor_axis[0],
                                              self.minor_axis[-1])
        else:
            minor = 'Minor axis: None'

        if len(self.items) > 0:
            items = 'Items: %s to %s' % (self.items[0], self.items[-1])
        else:
            items = 'Items: None'

        if len(self.labels) > 0:
            labels= 'Labels: %s to %s' % (self.labels[0], self.labels[-1])
        else:
            labels = 'Labels: None'

        output = '%s\n%s\n%s\n%s\n%s\n%s' % (class_name, dims, labels, items, major, minor)

        return output

    def __iter__(self):
        return iter(self.labels)

    def iteritems(self):
        for label in self.labels:
            yield label, self[label]

    iterkv = iteritems

    #----------------------------------------------------------------------
    # Getting and setting elements

    def get_value(self, label, item, major, minor):
        """
        Quickly retrieve single value at (labe, item, major, minor) location

        Parameters
        ----------
        label : label (fdpanel item)
        item  : item label (fdpanel item)
        major : major axis label (fdpanel item row)
        minor : minor axis label (fdpanel item column)

        Returns
        -------
        value : scalar value
        """
        # hm, two layers to the onion
        p = self._get_item_cache(label)
        return p.get_value(item, major, minor)

    def set_value(self, label, item, major, minor, value):
        """
        Quickly set single value at (labe, item, major, minor) location

        Parameters
        ----------
        label : label (fdpanel item)
        item  : item label (fdpanel item)
        major : major axis label (fdpanel item row)
        minor : minor axis label (fdpanel item column)

        Returns
        -------
        label : FDPanel
            If label combo is contained, will be reference to calling Panel,
            otherwise a new object
        """
        try:
            p = self._get_item_cache(label)
            p.set_value(item, major, minor, value)
            return self
        except KeyError:
            ax1, ax2, ax3, ax4 = self._expand_axes((label,item, major, minor))
            result = self.reindex(labels = ax1, items=ax2, major=ax3, minor=ax4, copy=False)

            likely_dtype = com._infer_dtype(value)
            made_bigger = not np.array_equal(ax1, self.labels)
            # how to make this logic simpler?
            if made_bigger:
                com._possibly_cast_item(result, label, likely_dtype)

            return result.set_value(label, item, major, minor, value)

    def _box_item_values(self, key, values):
        return Panel(values, items=self.items, major_axis=self.major_axis, minor_axis=self.minor_axis)

    def __getattr__(self, name):
        """After regular attribute access, try looking up the name of an item.
        This allows simpler access to items for interactive use."""
        if name in self.labels:
            return self[name]
        raise AttributeError("'%s' object has no attribute '%s'" %
                             (type(self).__name__, name))

    def __setitem__(self, key, value):
        _, I, N, K = self.shape
        if isinstance(value, Panel):
            value = value.reindex(items     =self.items,
                                  major_axis=self.major_axis,
                                  minor_axis=self.minor_axis)
            mat = value.values
        elif isinstance(value, np.ndarray):
            assert(value.shape == (I, N, K))
            mat = np.asarray(value)
        elif np.isscalar(value):
            dtype = _infer_dtype(value)
            mat = np.empty((I, N, K), dtype=dtype)
            mat.fill(value)

        mat = mat.reshape((1, I, N, K))
        NDFrame._set_item(self, key, mat)

    def _get_plane_axes(self, axis):
        axis = self._get_axis_name(axis)

        if axis == 'major_axis':
            items = self.labels
            major = self.minor_axis
            minor = self.items
        elif axis == 'minor_axis':
            items = self.labels
            major = self.major_axis
            minor = self.items
        elif axis == 'items':
            items = self.labels
            major = self.major_axis
            minor = self.minor_axis
        elif axis == 'labels':
            items = self.items
            major = self.major_axis
            minor = self.minor_axis

        return items, major, minor

    def conform(self, panel, axis='labels'):
        """
        Conform input Panel to align with chosen axis pair.

        Parameters
        ----------
        panel : Panel
        axis : {'labels', 'items', 'major', 'minor'}

        Returns
        -------
        Panel
        """
        items, major, minor = self._get_plane_axes(axis)
        return panel.reindex(items=items,major_axis=major,minor_axis=minor)

    def reindex(self, labels=None, major=None, items=None, minor=None, method=None,
                major_axis=None, minor_axis=None, copy=True):
        """
        Conform fdpanel to new axis or axes

        Parameters
        ----------
        labels: Index or sequence, default None
        items : Index or sequence, default None
        major : Index or sequence, default None
            Can also use 'major_axis' keyword
        minor : Index or sequence, default None
            Can also use 'minor_axis' keyword
        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default None
            Method to use for filling holes in reindexed Series

            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap
        copy : boolean, default True
            Return a new object, even if the passed indexes are the same

        Returns
        -------
        FDPanel (new object)
        """
        result = self

        major = _mut_exclusive(major, major_axis)
        minor = _mut_exclusive(minor, minor_axis)

        if major is not None:
            result = result._reindex_axis(major, method, 2, copy)

        if minor is not None:
            result = result._reindex_axis(minor, method, 3, copy)

        if items is not None:
            result = result._reindex_axis(items, method, 1, copy)

        if labels is not None:
            result = result._reindex_axis(labels, method, 0, copy)

        if result is self and copy:
            raise ValueError('Must specify at least one axis')

        return result

    def reindex_like(self, other, method=None):
        """
        Reindex FDPanel to match indices of another Panel

        Parameters
        ----------
        other : FDPanel
        method : string or None

        Returns
        -------
        reindexed : FDPanel
        """
        # todo: object columns
        return self.reindex(labels=other.labels, major=other.major_axis, items=other.items, minor=other.minor_axis, method=method)

    def swapaxes(self, axis1='major', axis2='minor'):
        """
        Interchange axes and swap values axes appropriately

        Returns
        -------
        y : FDPanel (new object)
        """
        i = self._get_axis_number(axis1)
        j = self._get_axis_number(axis2)

        if i == j:
            raise ValueError('Cannot specify the same axis')

        mapping = {i : j, j : i}

        new_axes = (self._get_axis(mapping.get(k, k))
                    for k in range(4))
        new_values = self.values.swapaxes(i, j).copy()

        return self._constructor(new_values, *new_axes)

    def xs(self, key, axis=2, copy=True):
        """
        Return slice of fdpanel along selected axis

        Parameters
        ----------
        key : object
            Label
        axis : {'labels', 'items', 'major', 'minor}, default 1/'major'

        Returns
        -------
        y : Panel
        """
        if axis == 0:
            data = self[key]
            if copy:
                data = data.copy()
            return data

        self._consolidate_inplace()
        axis_number = self._get_axis_number(axis)
        new_data = self._data.xs(key, axis=axis_number, copy=copy)
        return Panel(new_data)

    ### remove operations ####
    def major_xs(self, *args, **kwargs):
        raise NotImplementedError
    def minor_xs(self, *args, **kwargs):
        raise NotImplementedError
    def to_frame(self, *args, **kwargs):
        raise NotImplementedError
    def to_excel(self, *args, **kwargs):
        raise NotImplementedError

def _prep_ndarray(values, copy=True):
    if not isinstance(values, np.ndarray):
        values = np.asarray(values)
            # NumPy strings are a pain, convert to object
        if issubclass(values.dtype.type, basestring):
            values = np.array(values, dtype=object, copy=True)
    else:
        if copy:
            values = values.copy()
    assert(values.ndim == 4)
    return values

def _extract_axis(data, axis=0, intersect=False):
    from pandas.core.index import _union_indexes

    if len(data) == 0:
        index = Index([])
    elif len(data) > 0:
        raw_lengths = []
        indexes = []

        have_raw_arrays = False
        have_panels = False

        for v in data.values():
            if isinstance(v, Panel):
                have_panels = True
                indexes.append(v._get_axis(axis))
            else:
                have_raw_arrays = True
                raw_lengths.append(v.shape[axis])

        if have_panels:
            index = _get_combined_index(indexes, intersect=intersect)

        if have_raw_arrays:
            lengths = list(set(raw_lengths))
            if len(lengths) > 1:
                raise ValueError('ndarrays must match shape on axis %d' % axis)

            if have_panels:
                assert(lengths[0] == len(index))
            else:
                index = Index(np.arange(lengths[0]))

    return _ensure_index(index)
