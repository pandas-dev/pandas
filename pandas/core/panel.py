"""
Contains data structures designed for manipulating panel (3-dimensional) data
"""
# pylint: disable=E1103,W0231,W0212,W0621

import operator
import sys
import numpy as np

from pandas.core.common import (PandasError, _mut_exclusive,
                                _try_sort, _default_index, _infer_dtype)
from pandas.core.index import (Factor, Index, MultiIndex, _ensure_index,
                               _get_combined_index, NULL_INDEX)
from pandas.core.indexing import _NDFrameIndexer
from pandas.core.internals import BlockManager, make_block, form_blocks
from pandas.core.frame import DataFrame
from pandas.core.generic import NDFrame
from pandas.util import py3compat
from pandas.util.decorators import deprecate, Appender, Substitution
import pandas.core.common as com
import pandas.core.nanops as nanops
import pandas._tseries as lib


def _ensure_like_indices(time, panels):
    """
    Makes sure that time and panels are conformable
    """
    n_time = len(time)
    n_panel = len(panels)
    u_panels = np.unique(panels) # this sorts!
    u_time = np.unique(time)
    if len(u_time) == n_time:
        time = np.tile(u_time, len(u_panels))
    if len(u_panels) == n_panel:
        panels = np.repeat(u_panels, len(u_time))
    return time, panels

def panel_index(time, panels, names=['time', 'panel']):
    """
    Returns a multi-index suitable for a panel-like DataFrame

    Parameters
    ----------
    time : array-like
        Time index, does not have to repeat
    panels : array-like
        Panel index, does not have to repeat
    names : list, optional
        List containing the names of the indices

    Returns
    -------
    multi_index : MultiIndex
        Time index is the first level, the panels are the second level.

    Examples
    --------
    >>> years = range(1960,1963)
    >>> panels = ['A', 'B', 'C']
    >>> panel_idx = panel_index(years, panels)
    >>> panel_idx
    MultiIndex([(1960, 'A'), (1961, 'A'), (1962, 'A'), (1960, 'B'), (1961, 'B'),
       (1962, 'B'), (1960, 'C'), (1961, 'C'), (1962, 'C')], dtype=object)

    or

    >>> import numpy as np
    >>> years = np.repeat(range(1960,1963), 3)
    >>> panels = np.tile(['A', 'B', 'C'], 3)
    >>> panel_idx = panel_index(years, panels)
    >>> panel_idx
    MultiIndex([(1960, 'A'), (1960, 'B'), (1960, 'C'), (1961, 'A'), (1961, 'B'),
       (1961, 'C'), (1962, 'A'), (1962, 'B'), (1962, 'C')], dtype=object)
    """
    time, panels = _ensure_like_indices(time, panels)
    time_factor = Factor(time)
    panel_factor = Factor(panels)

    labels = [time_factor.labels, panel_factor.labels]
    levels = [time_factor.levels, panel_factor.levels]
    return MultiIndex(levels, labels, sortorder=None, names=names)

class PanelError(Exception):
    pass

def _arith_method(func, name):
    # work only for scalars

    def f(self, other):
        if not np.isscalar(other):
            raise ValueError('Simple arithmetic with Panel can only be '
                            'done with scalar values')

        return self._combine(other, func)
    f.__name__ = name
    return f

def _panel_arith_method(op, name):
    @Substitution(op)
    def f(self, other, axis='items'):
        """
        Wrapper method for %s

        Parameters
        ----------
        other : DataFrame or Panel class
        axis : {'items', 'major', 'minor'}
            Axis to broadcast over

        Returns
        -------
        Panel
        """
        return self._combine(other, op, axis=axis)

    f.__name__ = name
    return f


_agg_doc = """
Return %(desc)s over requested axis

Parameters
----------
axis : {'items', 'major', 'minor'} or {0, 1, 2}
skipna : boolean, default True
    Exclude NA/null values. If an entire row/column is NA, the result
    will be NA

Returns
-------
%(outname)s : DataFrame
"""

_na_info = """

NA/null values are %s.
If all values are NA, result will be NA"""


class Panel(NDFrame):
    _AXIS_NUMBERS = {
        'items' : 0,
        'major_axis' : 1,
        'minor_axis' : 2
    }

    _AXIS_ALIASES = {
        'major' : 'major_axis',
        'minor' : 'minor_axis'
    }

    _AXIS_NAMES = {
        0 : 'items',
        1 : 'major_axis',
        2 : 'minor_axis'
    }

    # major
    _default_stat_axis = 1
    _het_axis = 0

    items = lib.AxisProperty(0)
    major_axis = lib.AxisProperty(1)
    minor_axis = lib.AxisProperty(2)

    __add__ = _arith_method(operator.add, '__add__')
    __sub__ = _arith_method(operator.sub, '__sub__')
    __truediv__ = _arith_method(operator.truediv, '__truediv__')
    __floordiv__ = _arith_method(operator.floordiv, '__floordiv__')
    __mul__ = _arith_method(operator.mul, '__mul__')
    __pow__ = _arith_method(operator.pow, '__pow__')

    __radd__ = _arith_method(operator.add, '__radd__')
    __rmul__ = _arith_method(operator.mul, '__rmul__')
    __rsub__ = _arith_method(lambda x, y: y - x, '__rsub__')
    __rtruediv__ = _arith_method(lambda x, y: y / x, '__rtruediv__')
    __rfloordiv__ = _arith_method(lambda x, y: y // x, '__rfloordiv__')
    __rpow__ = _arith_method(lambda x, y: y ** x, '__rpow__')

    if not py3compat.PY3:
        __div__ = _arith_method(operator.div, '__div__')
        __rdiv__ = _arith_method(lambda x, y: y / x, '__rdiv__')

    def __init__(self, data=None, items=None, major_axis=None, minor_axis=None,
                 copy=False, dtype=None):
        """
        Represents wide format panel data, stored as 3-dimensional array

        Parameters
        ----------
        data : ndarray (items x major x minor), or dict of DataFrames
        items : Index or array-like
            axis=1
        major_axis : Index or array-like
            axis=1
        minor_axis : Index or array-like
            axis=2
        dtype : dtype, default None
            Data type to force, otherwise infer
        copy : boolean, default False
            Copy data from inputs. Only affects DataFrame / 2d ndarray input
        """
        if data is None:
            data = {}

        passed_axes = [items, major_axis, minor_axis]
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
            raise PandasError('Panel constructor not properly called!')

        NDFrame.__init__(self, mgr, axes=axes, copy=copy, dtype=dtype)

    @classmethod
    def _from_axes(cls, data, axes):
        # for construction from BlockManager
        if isinstance(data, BlockManager):
            return cls(data)
        else:
            items, major, minor = axes
            return cls(data, items=items, major_axis=major,
                       minor_axis=minor, copy=False)

    def _init_dict(self, data, axes, dtype=None):
        items, major, minor = axes

        # prefilter if items passed
        if items is not None:
            items = _ensure_index(items)
            data = dict((k, v) for k, v in data.iteritems() if k in items)
        else:
            items = Index(_try_sort(data.keys()))

        for k, v in data.iteritems():
            if isinstance(v, dict):
                data[k] = DataFrame(v)

        if major is None:
            major = _extract_axis(data, axis=0)

        if minor is None:
            minor = _extract_axis(data, axis=1)

        axes = [items, major, minor]
        reshaped_data = data.copy() # shallow

        item_shape = len(major), len(minor)
        for item in items:
            v = values = data.get(item)
            if v is None:
                values = np.empty(item_shape, dtype=dtype)
                values.fill(np.nan)
            elif isinstance(v, DataFrame):
                v = v.reindex(index=major, columns=minor, copy=False)
                if dtype is not None:
                    v = v.astype(dtype)
                values = v.values
            reshaped_data[item] = values

        # segregates dtypes and forms blocks matching to columns
        blocks = form_blocks(reshaped_data, axes)
        mgr = BlockManager(blocks, axes).consolidate()
        return mgr

    @property
    def shape(self):
        return len(self.items), len(self.major_axis), len(self.minor_axis)

    @classmethod
    def from_dict(cls, data, intersect=False, orient='items', dtype=None):
        """
        Construct Panel from dict of DataFrame objects

        Parameters
        ----------
        data : dict
            {field : DataFrame}
        intersect : boolean
            Intersect indexes of input DataFrames
        orient : {'items', 'minor'}, default 'items'
            The "orientation" of the data. If the keys of the passed dict
            should be the items of the result panel, pass 'items'
            (default). Otherwise if the columns of the values of the passed
            DataFrame objects should be the items (which in the case of
            mixed-dtype data you should do), instead pass 'minor'


        Returns
        -------
        Panel
        """
        from collections import defaultdict

        orient = orient.lower()
        if orient == 'minor':
            new_data = defaultdict(dict)
            for col, df in data.iteritems():
                for item, s in df.iteritems():
                    new_data[item][col] = s
            data = new_data
        elif orient != 'items':  # pragma: no cover
            raise ValueError('only recognize items or minor for orientation')

        data, index, columns = _homogenize_dict(data, intersect=intersect,
                                                dtype=dtype)
        items = Index(sorted(data.keys()))
        return Panel(data, items, index, columns)

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


    #----------------------------------------------------------------------
    # Array interface

    def __array__(self, dtype=None):
        return self.values

    def __array_wrap__(self, result):
        return self._constructor(result, items=self.items,
                                 major_axis=self.major_axis,
                                 minor_axis=self.minor_axis, copy=False)

    #----------------------------------------------------------------------
    # Magic methods

    def __repr__(self):
        class_name = str(self.__class__)

        I, N, K = len(self.items), len(self.major_axis), len(self.minor_axis)

        dims = 'Dimensions: %d (items) x %d (major) x %d (minor)' % (I, N, K)

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

        output = '%s\n%s\n%s\n%s\n%s' % (class_name, dims, items, major, minor)

        return output

    def __iter__(self):
        return iter(self.items)

    def iteritems(self):
        for item in self.items:
            yield item, self[item]

    # Name that won't get automatically converted to items by 2to3. items is
    # already in use for the first axis.
    iterkv = iteritems

    def _get_plane_axes(self, axis):
        """

        """
        axis = self._get_axis_name(axis)

        if axis == 'major_axis':
            index = self.minor_axis
            columns = self.items
        if axis == 'minor_axis':
            index = self.major_axis
            columns = self.items
        elif axis == 'items':
            index = self.major_axis
            columns = self.minor_axis

        return index, columns

    @property
    def _constructor(self):
        return Panel

    # Fancy indexing
    _ix = None

    @property
    def ix(self):
        if self._ix is None:
            self._ix = _NDFrameIndexer(self)

        return self._ix

    def _wrap_array(self, arr, axes, copy=False):
        items, major, minor = axes
        return self._constructor(arr, items=items, major_axis=major,
                                 minor_axis=minor, copy=copy)

    fromDict = from_dict

    def to_sparse(self, fill_value=None, kind='block'):
        """
        Convert to SparsePanel

        Parameters
        ----------
        fill_value : float, default NaN
        kind : {'block', 'integer'}

        Returns
        -------
        y : SparseDataFrame
        """
        from pandas.core.sparse import SparsePanel
        frames = dict(self.iterkv())
        return SparsePanel(frames, items=self.items,
                           major_axis=self.major_axis,
                           minor_axis=self.minor_axis,
                           default_kind=kind,
                           default_fill_value=fill_value)

    def to_excel(self, path, na_rep=''):
        """
        Write each DataFrame in Panel to a separate excel sheet 

        Parameters
        ----------
        excel_writer : string or ExcelWriter object
            File path or existing ExcelWriter 
        na_rep : string, default ''
            Missing data rep'n
        """
        from pandas.io.parsers import ExcelWriter
        writer = ExcelWriter(path)
        for item, df in self.iteritems():
            name = str(item)
            df.to_excel(writer, name, na_rep=na_rep)
        writer.save()

    # TODO: needed?
    def keys(self):
        return list(self.items)

    def _get_values(self):
        self._consolidate_inplace()
        return self._data.as_matrix()

    values = property(fget=_get_values)

    #----------------------------------------------------------------------
    # Getting and setting elements

    def get_value(self, item, major, minor):
        """
        Quickly retrieve single value at (item, major, minor) location

        Parameters
        ----------
        item : item label (panel item)
        major : major axis label (panel item row)
        minor : minor axis label (panel item column)

        Returns
        -------
        value : scalar value
        """
        # hm, two layers to the onion
        frame = self._get_item_cache(item)
        return frame.get_value(major, minor)

    def set_value(self, item, major, minor, value):
        """
        Quickly set single value at (item, major, minor) location

        Parameters
        ----------
        item : item label (panel item)
        major : major axis label (panel item row)
        minor : minor axis label (panel item column)
        value : scalar

        Returns
        -------
        panel : Panel
            If label combo is contained, will be reference to calling Panel,
            otherwise a new object
        """
        try:
            frame = self._get_item_cache(item)
            frame.set_value(major, minor, value)
            return self
        except KeyError:
            ax1, ax2, ax3 = self._expand_axes((item, major, minor))
            result = self.reindex(items=ax1, major=ax2, minor=ax3, copy=False)

            likely_dtype = com._infer_dtype(value)
            made_bigger = not np.array_equal(ax1, self.items)
            # how to make this logic simpler?
            if made_bigger:
                com._possibly_cast_item(result, item, likely_dtype)

            return result.set_value(item, major, minor, value)

    def _box_item_values(self, key, values):
        return DataFrame(values, index=self.major_axis, columns=self.minor_axis)

    def __getattr__(self, name):
        """After regular attribute access, try looking up the name of an item.
        This allows simpler access to items for interactive use."""
        if name in self.items:
            return self[name]
        raise AttributeError("'%s' object has no attribute '%s'" %
                             (type(self).__name__, name))

    def _slice(self, slobj, axis=0):
        new_data = self._data.get_slice(slobj, axis=axis)
        return self._constructor(new_data)

    def __setitem__(self, key, value):
        _, N, K = self.shape
        if isinstance(value, DataFrame):
            value = value.reindex(index=self.major_axis,
                                  columns=self.minor_axis)
            mat = value.values
        elif isinstance(value, np.ndarray):
            assert(value.shape == (N, K))
            mat = np.asarray(value)
        elif np.isscalar(value):
            dtype = _infer_dtype(value)
            mat = np.empty((N, K), dtype=dtype)
            mat.fill(value)

        mat = mat.reshape((1, N, K))
        NDFrame._set_item(self, key, mat)

    def pop(self, item):
        """
        Return item slice from panel and delete from panel

        Parameters
        ----------
        key : object
            Must be contained in panel's items

        Returns
        -------
        y : DataFrame
        """
        return NDFrame.pop(self, item)

    def __getstate__(self):
        "Returned pickled representation of the panel"
        return self._data

    def __setstate__(self, state):
        # old Panel pickle
        if isinstance(state, BlockManager):
            self._data = state
        elif len(state) == 4: # pragma: no cover
            self._unpickle_panel_compat(state)
        else: # pragma: no cover
            raise ValueError('unrecognized pickle')
        self._item_cache = {}

    def _unpickle_panel_compat(self, state): # pragma: no cover
        "Unpickle the panel"
        _unpickle = com._unpickle_array
        vals, items, major, minor = state

        items = _unpickle(items)
        major = _unpickle(major)
        minor = _unpickle(minor)
        values = _unpickle(vals)
        wp = Panel(values, items, major, minor)
        self._data = wp._data

    def conform(self, frame, axis='items'):
        """
        Conform input DataFrame to align with chosen axis pair.

        Parameters
        ----------
        frame : DataFrame
        axis : {'items', 'major', 'minor'}

            Axis the input corresponds to. E.g., if axis='major', then
            the frame's columns would be items, and the index would be
            values of the minor axis

        Returns
        -------
        DataFrame
        """
        index, columns = self._get_plane_axes(axis)
        return frame.reindex(index=index, columns=columns)

    def reindex(self, major=None, items=None, minor=None, method=None,
                major_axis=None, minor_axis=None, copy=True):
        """
        Conform panel to new axis or axes

        Parameters
        ----------
        major : Index or sequence, default None
            Can also use 'major_axis' keyword
        items : Index or sequence, default None
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
        Panel (new object)
        """
        result = self

        major = _mut_exclusive(major, major_axis)
        minor = _mut_exclusive(minor, minor_axis)

        if major is not None:
            result = result._reindex_axis(major, method, 1, copy)

        if minor is not None:
            result = result._reindex_axis(minor, method, 2, copy)

        if items is not None:
            result = result._reindex_axis(items, method, 0, copy)

        if result is self and copy:
            raise ValueError('Must specify at least one axis')

        return result

    def reindex_like(self, other, method=None):
        """
        Reindex Panel to match indices of another Panel

        Parameters
        ----------
        other : Panel
        method : string or None

        Returns
        -------
        reindexed : Panel
        """
        # todo: object columns
        return self.reindex(major=other.major_axis, items=other.items,
                            minor=other.minor_axis, method=method)

    def _combine(self, other, func, axis=0):
        if isinstance(other, Panel):
            return self._combine_panel(other, func)
        elif isinstance(other, DataFrame):
            return self._combine_frame(other, func, axis=axis)
        elif np.isscalar(other):
            new_values = func(self.values, other)
            return Panel(new_values, self.items, self.major_axis,
                             self.minor_axis)

    def __neg__(self):
        return -1 * self

    def _combine_frame(self, other, func, axis=0):
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

        return Panel(new_values, self.items, self.major_axis,
                     self.minor_axis)

    def _combine_panel(self, other, func):
        items = self.items + other.items
        major = self.major_axis + other.major_axis
        minor = self.minor_axis + other.minor_axis

        # could check that everything's the same size, but forget it
        this = self.reindex(items=items, major=major, minor=minor)
        other = other.reindex(items=items, major=major, minor=minor)

        result_values = func(this.values, other.values)

        return Panel(result_values, items, major, minor)

    def fillna(self, value=None, method='pad'):
        """
        Fill NaN values using the specified method.

        Member Series / TimeSeries are filled separately.

        Parameters
        ----------
        value : any kind (should be same type as array)
            Value to use to fill holes (e.g. 0)

        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default 'pad'
            Method to use for filling holes in reindexed Series

            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap

        Returns
        -------
        y : DataFrame

        See also
        --------
        DataFrame.reindex, DataFrame.asfreq
        """
        if value is None:
            result = {}
            for col, s in self.iterkv():
                result[col] = s.fillna(method=method, value=value)

            return Panel.from_dict(result)
        else:
            new_data = self._data.fillna(value)
            return Panel(new_data)

    add = _panel_arith_method(operator.add, 'add')
    subtract = sub = _panel_arith_method(operator.sub, 'subtract')
    multiply = mul = _panel_arith_method(operator.mul, 'multiply')

    try:
        divide = div = _panel_arith_method(operator.div, 'divide')
    except AttributeError:  # pragma: no cover
        # Python 3
        divide = div = _panel_arith_method(operator.truediv, 'divide')

    def major_xs(self, key, copy=True):
        """
        Return slice of panel along major axis

        Parameters
        ----------
        key : object
            Major axis label
        copy : boolean, default False
            Copy data

        Returns
        -------
        y : DataFrame
            index -> minor axis, columns -> items
        """
        return self.xs(key, axis=1, copy=copy)

    def minor_xs(self, key, copy=True):
        """
        Return slice of panel along minor axis

        Parameters
        ----------
        key : object
            Minor axis label
        copy : boolean, default False
            Copy data

        Returns
        -------
        y : DataFrame
            index -> major axis, columns -> items
        """
        return self.xs(key, axis=2, copy=copy)

    def xs(self, key, axis=1, copy=True):
        """
        Return slice of panel along selected axis

        Parameters
        ----------
        key : object
            Label
        axis : {'items', 'major', 'minor}, default 1/'major'

        Returns
        -------
        y : DataFrame
        """
        if axis == 0:
            data = self[key]
            if copy:
                data = data.copy()
            return data

        self._consolidate_inplace()
        axis_number = self._get_axis_number(axis)
        new_data = self._data.xs(key, axis=axis_number, copy=copy)
        return DataFrame(new_data)

    def groupby(self, function, axis='major'):
        """
        Group data on given axis, returning GroupBy object

        Parameters
        ----------
        function : callable
            Mapping function for chosen access
        axis : {'major', 'minor', 'items'}, default 'major'

        Returns
        -------
        grouped : PanelGroupBy
        """
        from pandas.core.groupby import PanelGroupBy
        axis = self._get_axis_number(axis)
        return PanelGroupBy(self, function, axis=axis)

    def swapaxes(self, axis1='major', axis2='minor'):
        """
        Interchange axes and swap values axes appropriately

        Returns
        -------
        y : Panel (new object)
        """
        i = self._get_axis_number(axis1)
        j = self._get_axis_number(axis2)

        if i == j:
            raise ValueError('Cannot specify the same axis')

        mapping = {i : j, j : i}

        new_axes = (self._get_axis(mapping.get(k, k))
                    for k in range(3))
        new_values = self.values.swapaxes(i, j).copy()

        return Panel(new_values, *new_axes)

    def to_frame(self, filter_observations=True):
        """
        Transform wide format into long (stacked) format as DataFrame

        Parameters
        ----------
        filter_observations : boolean, default True
            Drop (major, minor) pairs without a complete set of observations
            across all the items

        Returns
        -------
        y : DataFrame
        """
        _, N, K = self.shape

        if filter_observations:
            mask = com.notnull(self.values).all(axis=0)
            # size = mask.sum()
            selector = mask.ravel()
        else:
            # size = N * K
            selector = slice(None, None)

        data = {}
        for item in self.items:
            data[item] = self[item].values.ravel()[selector]

        major_labels = np.arange(N).repeat(K)[selector]

        # Anyone think of a better way to do this? np.repeat does not
        # do what I want
        minor_labels = np.arange(K).reshape(1, K)[np.zeros(N, dtype=int)]
        minor_labels = minor_labels.ravel()[selector]

        index = MultiIndex(levels=[self.major_axis, self.minor_axis],
                           labels=[major_labels, minor_labels],
                           names=['major', 'minor'])

        return DataFrame(data, index=index, columns=self.items)

    to_long = deprecate('to_long', to_frame)
    toLong = deprecate('toLong', to_frame)

    def filter(self, items):
        """
        Restrict items in panel to input list

        Parameters
        ----------
        items : sequence

        Returns
        -------
        y : Panel
        """
        intersection = self.items.intersection(items)
        return self.reindex(items=intersection)

    def apply(self, func, axis='major'):
        """
        Apply

        Parameters
        ----------
        func : numpy function
            Signature should match numpy.{sum, mean, var, std} etc.
        axis : {'major', 'minor', 'items'}
        fill_value : boolean, default True
            Replace NaN values with specified first

        Returns
        -------
        result : DataFrame or Panel
        """
        i = self._get_axis_number(axis)
        result = np.apply_along_axis(func, i, self.values)
        return self._wrap_result(result, axis=axis)

    def _reduce(self, op, axis=0, skipna=True):
        axis_name = self._get_axis_name(axis)
        axis_number = self._get_axis_number(axis_name)
        f = lambda x: op(x, axis=axis_number, skipna=skipna)

        result = f(self.values)

        index, columns = self._get_plane_axes(axis_name)
        if axis_name != 'items':
            result = result.T

        return DataFrame(result, index=index, columns=columns)

    def _wrap_result(self, result, axis):
        axis = self._get_axis_name(axis)
        index, columns = self._get_plane_axes(axis)

        if axis != 'items':
            result = result.T

        return DataFrame(result, index=index, columns=columns)

    def count(self, axis='major'):
        """
        Return number of observations over requested axis.

        Parameters
        ----------
        axis : {'items', 'major', 'minor'} or {0, 1, 2}

        Returns
        -------
        count : DataFrame
        """
        i = self._get_axis_number(axis)

        values = self.values
        mask = np.isfinite(values)
        result = mask.sum(axis=i)

        return self._wrap_result(result, axis)

    @Substitution(desc='sum', outname='sum')
    @Appender(_agg_doc)
    def sum(self, axis='major', skipna=True):
        return self._reduce(nanops.nansum, axis=axis, skipna=skipna)

    @Substitution(desc='mean', outname='mean')
    @Appender(_agg_doc)
    def mean(self, axis='major', skipna=True):
        return self._reduce(nanops.nanmean, axis=axis, skipna=skipna)

    @Substitution(desc='unbiased variance', outname='variance')
    @Appender(_agg_doc)
    def var(self, axis='major', skipna=True):
        return self._reduce(nanops.nanvar, axis=axis, skipna=skipna)

    @Substitution(desc='unbiased standard deviation', outname='stdev')
    @Appender(_agg_doc)
    def std(self, axis='major', skipna=True):
        return self.var(axis=axis, skipna=skipna).apply(np.sqrt)

    @Substitution(desc='unbiased skewness', outname='skew')
    @Appender(_agg_doc)
    def skew(self, axis='major', skipna=True):
        return self._reduce(nanops.nanskew, axis=axis, skipna=skipna)

    @Substitution(desc='product', outname='prod')
    @Appender(_agg_doc)
    def prod(self, axis='major', skipna=True):
        return self._reduce(nanops.nanprod, axis=axis, skipna=skipna)

    @Substitution(desc='compounded percentage', outname='compounded')
    @Appender(_agg_doc)
    def compound(self, axis='major', skipna=True):
        return (1 + self).prod(axis=axis, skipna=skipna) - 1

    @Substitution(desc='median', outname='median')
    @Appender(_agg_doc)
    def median(self, axis='major', skipna=True):
        return self._reduce(nanops.nanmedian, axis=axis, skipna=skipna)

    @Substitution(desc='maximum', outname='maximum')
    @Appender(_agg_doc)
    def max(self, axis='major', skipna=True):
        return self._reduce(nanops.nanmax, axis=axis, skipna=skipna)

    @Substitution(desc='minimum', outname='minimum')
    @Appender(_agg_doc)
    def min(self, axis='major', skipna=True):
        return self._reduce(nanops.nanmin, axis=axis, skipna=skipna)

    def shift(self, lags, axis='major'):
        """
        Shift major or minor axis by specified number of lags. Drops periods

        Parameters
        ----------
        lags : int
            Needs to be a positive number currently
        axis : {'major', 'minor'}

        Returns
        -------
        shifted : Panel
        """
        values = self.values
        items = self.items
        major_axis = self.major_axis
        minor_axis = self.minor_axis

        if axis == 'major':
            values = values[:, :-lags, :]
            major_axis = major_axis[lags:]
        elif axis == 'minor':
            values = values[:, :, :-lags]
            minor_axis = minor_axis[lags:]
        else:
            raise ValueError('Invalid axis')

        return Panel(values, items=items, major_axis=major_axis,
                         minor_axis=minor_axis)

    def truncate(self, before=None, after=None, axis='major'):
        """Function truncates a sorted Panel before and/or after some
        particular values on the requested axis

        Parameters
        ----------
        before : date
            Left boundary
        after : date
            Right boundary
        axis : {'major', 'minor', 'items'}

        Returns
        -------
        Panel
        """
        axis = self._get_axis_name(axis)
        index = self._get_axis(axis)

        beg_slice, end_slice = index.slice_locs(before, after)
        new_index = index[beg_slice:end_slice]

        return self.reindex(**{axis : new_index})

    def join(self, other, how='left', lsuffix='', rsuffix=''):
        """
        Join items with other Panel either on major and minor axes column

        Parameters
        ----------
        other : Panel or list of Panels
            Index should be similar to one of the columns in this one
        how : {'left', 'right', 'outer', 'inner'}
            How to handle indexes of the two objects. Default: 'left'
            for joining on index, None otherwise
            * left: use calling frame's index
            * right: use input frame's index
            * outer: form union of indexes
            * inner: use intersection of indexes
        lsuffix : string
            Suffix to use from left frame's overlapping columns
        rsuffix : string
            Suffix to use from right frame's overlapping columns

        Returns
        -------
        joined : Panel
        """
        from pandas.tools.merge import concat

        if isinstance(other, Panel):
            join_major, join_minor = self._get_join_index(other, how)
            this = self.reindex(major=join_major, minor=join_minor)
            other = other.reindex(major=join_major, minor=join_minor)
            merged_data = this._data.merge(other._data, lsuffix, rsuffix)
            return self._constructor(merged_data)
        else:
            if lsuffix or rsuffix:
                raise ValueError('Suffixes not supported when passing multiple '
                                 'panels')

            if how == 'left':
                how = 'outer'
                join_axes = [self.major_axis, self.minor_axis]
            elif how == 'right':
                raise ValueError('Right join not supported with multiple '
                                 'panels')
            else:
                join_axes = None

            return concat([self] + list(other), axis=0, join=how,
                          join_axes=join_axes, verify_integrity=True)

    def _get_join_index(self, other, how):
        if how == 'left':
            join_major, join_minor = self.major_axis, self.minor_axis
        elif how == 'right':
            join_major, join_minor = other.major_axis, other.minor_axis
        elif how == 'inner':
            join_major = self.major_axis.intersection(other.major_axis)
            join_minor = self.minor_axis.intersection(other.minor_axis)
        elif how == 'outer':
            join_major = self.major_axis.union(other.major_axis)
            join_minor = self.minor_axis.union(other.minor_axis)
        return join_major, join_minor

WidePanel = Panel
LongPanel = DataFrame

def _prep_ndarray(values, copy=True):
    if not isinstance(values, np.ndarray):
        values = np.asarray(values)
        # NumPy strings are a pain, convert to object
        if issubclass(values.dtype.type, basestring):
            values = np.array(values, dtype=object, copy=True)
    else:
        if copy:
            values = values.copy()
    assert(values.ndim == 3)
    return values

def _homogenize_dict(frames, intersect=True, dtype=None):
    """
    Conform set of DataFrame-like objects to either an intersection
    of indices / columns or a union.

    Parameters
    ----------
    frames : dict
    intersect : boolean, default True

    Returns
    -------
    dict of aligned frames, index, columns
    """
    result = {}

    adj_frames = {}
    for k, v in frames.iteritems():
        if isinstance(v, dict):
            adj_frames[k] = DataFrame(v)
        else:
            adj_frames[k] = v

    index = _extract_axis(adj_frames, axis=0, intersect=intersect)
    columns = _extract_axis(adj_frames, axis=1, intersect=intersect)

    for key, frame in adj_frames.iteritems():
        result[key] = frame.reindex(index=index, columns=columns,
                                    copy=False)

    return result, index, columns


def _extract_axis(data, axis=0, intersect=False):
    from pandas.core.index import _union_indexes

    if len(data) == 0:
        index = NULL_INDEX
    elif len(data) > 0:
        raw_lengths = []
        indexes = []

        have_raw_arrays = False
        have_frames = False

        for v in data.values():
            if isinstance(v, DataFrame):
                have_frames = True
                indexes.append(v._get_axis(axis))
            else:
                have_raw_arrays = True
                raw_lengths.append(v.shape[axis])

        if have_frames:
            index = _get_combined_index(indexes, intersect=intersect)

        if have_raw_arrays:
            lengths = list(set(raw_lengths))
            if len(lengths) > 1:
                raise ValueError('ndarrays must match shape on axis %d' % axis)

            if have_frames:
                assert(lengths[0] == len(index))
            else:
                index = Index(np.arange(lengths[0]))

    if len(index) == 0:
        index = NULL_INDEX

    return _ensure_index(index)


def _monotonic(arr):
    return not (arr[1:] < arr[:-1]).any()

def install_ipython_completers():  # pragma: no cover
    """Register the Panel type with IPython's tab completion machinery, so
    that it knows about accessing column names as attributes."""
    from IPython.utils.generics import complete_object

    @complete_object.when_type(Panel)
    def complete_dataframe(obj, prev_completions):
        return prev_completions + [c for c in obj.items \
                    if isinstance(c, basestring) and py3compat.isidentifier(c)]

# Importing IPython brings in about 200 modules, so we want to avoid it unless
# we're in IPython (when those modules are loaded anyway).
if "IPython" in sys.modules:  # pragma: no cover
    try:
        install_ipython_completers()
    except Exception:
        pass

