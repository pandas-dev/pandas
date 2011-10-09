"""
Contains data structures designed for manipulating panel (3-dimensional) data
"""
# pylint: disable=E1103,W0231,W0212,W0621

import operator
import numpy as np

from pandas.core.common import (PandasError, _mut_exclusive,
                                _try_sort, _default_index, _infer_dtype)
from pandas.core.index import Factor, Index, MultiIndex, _ensure_index
from pandas.core.indexing import _NDFrameIndexer
from pandas.core.internals import BlockManager, make_block, form_blocks
from pandas.core.frame import DataFrame, _union_indexes
from pandas.core.generic import AxisProperty, NDFrame
from pandas.core.series import Series
from pandas.util.decorators import deprecate
from pandas.util import py3compat
import pandas.core.common as common
import pandas._tseries as _tseries


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

    return f

def _panel_arith_method(op, name):
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
        LongPanel
        """
        return self._combine(other, op, axis=axis)

    f.__name__ = name
    f.__doc__ = f.__doc__ % str(op)

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


def _add_docs(method, desc, outname):
    doc = _agg_doc % {'desc' : desc,
                      'outname' : outname}
    method.__doc__ = doc

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

    items = AxisProperty(0)
    major_axis = AxisProperty(1)
    minor_axis = AxisProperty(2)

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

    def __init__(self, data, items=None, major_axis=None, minor_axis=None,
                 copy=False, dtype=None):
        """
        Represents wide format panel data, stored as 3-dimensional array

        Parameters
        ----------
        values : ndarray (items x major x minor)
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
        passed_axes = [items, major_axis, minor_axis]
        if isinstance(data, BlockManager):
            mgr = data
            if copy and dtype is None:
                mgr = mgr.copy()
            elif dtype is not None:
                # no choice but to copy
                mgr = mgr.astype(dtype)
        elif isinstance(data, dict):
            mgr = self._init_dict(data, passed_axes, dtype=dtype)
        elif isinstance(data, (np.ndarray, list)):
            mgr = self._init_matrix(data, passed_axes, dtype=dtype, copy=copy)
        else: # pragma: no cover
            raise PandasError('Panel constructor not properly called!')

        self._data = mgr

    def _init_dict(self, data, axes, dtype=None):
        items, major, minor = axes

        # prefilter if items passed
        if items is not None:
            items = _ensure_index(items)
            data = dict((k, v) for k, v in data.iteritems() if k in items)
        else:
            items = Index(_try_sort(data.keys()))

        for k, v in data.iteritems():
            if not isinstance(v, DataFrame):
                data[k] = DataFrame(v)

        if major is None:
            indexes = [v.index for v in data.values()]
            major = _union_indexes(indexes)

        if minor is None:
            indexes = [v.columns for v in data.values()]
            minor = _union_indexes(indexes)

        axes = [items, major, minor]

        reshaped_data = data.copy() # shallow
        # homogenize
        for k, v in data.iteritems():
            v = v.reindex(index=major, columns=minor, copy=False)
            if dtype is not None:
                v = v.astype(dtype)
            values = v.values
            shape = values.shape
            reshaped_data[k] = values.reshape((1,) + shape)

        # segregates dtypes and forms blocks matching to columns
        blocks = form_blocks(reshaped_data, axes)
        mgr = BlockManager(blocks, axes).consolidate()
        return mgr

    @property
    def shape(self):
        return len(self.items), len(self.major_axis), len(self.minor_axis)

    @classmethod
    def from_dict(cls, data, intersect=False, dtype=float):
        """
        Construct Panel from dict of DataFrame objects

        Parameters
        ----------
        data : dict
            {field : DataFrame}
        intersect : boolean

        Returns
        -------
        Panel
        """
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

    def __repr__(self):
        class_name = str(self.__class__)

        I, N, K = len(self.items), len(self.major_axis), len(self.minor_axis)

        dims = 'Dimensions: %d (items) x %d (major) x %d (minor)' % (I, N, K)

        major = 'Major axis: %s to %s' % (self.major_axis[0],
                                          self.major_axis[-1])

        minor = 'Minor axis: %s to %s' % (self.minor_axis[0],
                                          self.minor_axis[-1])

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

    # TODO: needed?
    def keys(self):
        return list(self.items)

    def _get_values(self):
        self._consolidate_inplace()
        return self._data.as_matrix()

    values = property(fget=_get_values)

    def __getitem__(self, key):
        mat = self._data.get(key)
        return DataFrame(mat, index=self.major_axis, columns=self.minor_axis)

    def _slice(self, slobj, axis=0):
        new_data = self._data.get_slice(slobj, axis=axis)
        return self._constructor(new_data)

    def __setitem__(self, key, value):
        _, N, K = self.shape

        # XXX
        if isinstance(value, LongPanel):
            if len(value.items) != 1:
                raise ValueError('Input panel must have only one item!')

            value = value.to_wide()[value.items[0]]

        if isinstance(value, DataFrame):
            value = value.reindex(index=self.major_axis,
                                  columns=self.minor_axis)
            mat = value.values

        elif np.isscalar(value):
            dtype = _infer_dtype(value)
            mat = np.empty((N, K), dtype=dtype)
            mat.fill(value)

        mat = mat.reshape((1, N, K))
        self._data.set(key, mat)

    def __delitem__(self, key):
        self._data.delete(key)

    def pop(self, key):
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
        result = self[key]
        del self[key]
        return result

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

    def _unpickle_panel_compat(self, state): # pragma: no cover
        "Unpickle the panel"
        _unpickle = common._unpickle_array
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
        if isinstance(other, (Panel, LongPanel)):
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
        if isinstance(other, LongPanel):
            other = other.to_wide()

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
    except AttributeError:   # Python 3
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

    def to_long(self, filter_observations=True):
        """
        Transform wide format into long (stacked) format

        Parameters
        ----------
        filter_observations : boolean, default True
            Drop (major, minor) pairs without a complete set of observations
            across all the items

        Returns
        -------
        y : LongPanel
        """
        I, N, K = self.shape

        if filter_observations:
            mask = common.notnull(self.values).all(axis=0)
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
                           labels=[major_labels, minor_labels])

        return LongPanel(data, index=index, columns=self.items)

    toLong = to_long

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

    def _array_method(self, func, axis='major', fill_value=None, skipna=True):
        """
        Parameters
        ----------
        func : numpy function
            Signature should match numpy.{sum, mean, var, std} etc.
        axis : {'major', 'minor', 'items'}
        fill_value : boolean, default True
            Replace NaN values with specified first
        skipna : boolean, default True
            Exclude NA/null values. If an entire row/column is NA, the result
            will be NA

        Returns
        -------
        y : DataFrame
        """
        result = self._values_aggregate(func, axis, fill_value, skipna=skipna)
        return self._wrap_result(result, axis=axis)

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

    def sum(self, axis='major', skipna=True):
        return self._array_method(np.sum, axis=axis, fill_value=0,
                                  skipna=skipna)

    _add_docs(sum, 'sum', 'sum')

    def mean(self, axis='major', skipna=True):
        the_sum = self.sum(axis=axis, skipna=skipna)
        the_count = self.count(axis=axis)
        return the_sum / the_count

    _add_docs(mean, 'mean', 'mean')

    def var(self, axis='major', skipna=True):
        i = self._get_axis_number(axis)
        y = np.array(self.values)
        mask = np.isnan(y)

        count = (-mask).sum(axis=i).astype(float)

        if skipna:
            y[mask] = 0

        X = y.sum(axis=i)
        XX = (y ** 2).sum(axis=i)

        theVar = (XX - X**2 / count) / (count - 1)

        return self._wrap_result(theVar, axis)

    _add_docs(var, 'unbiased variance', 'variance')

    def std(self, axis='major', skipna=True):
        return self.var(axis=axis, skipna=skipna).apply(np.sqrt)

    _add_docs(std, 'unbiased standard deviation', 'stdev')

    def skew(self, axis='major', skipna=True):
        raise NotImplementedError

    def prod(self, axis='major', skipna=True):
        return self._array_method(np.prod, axis=axis, fill_value=1,
                                  skipna=skipna)

    _add_docs(prod, 'product', 'prod')

    def compound(self, axis='major', skipna=True):
        return (1 + self).prod(axis=axis, skipna=skipna) - 1

    _add_docs(compound, 'compounded percentage', 'compounded')

    def median(self, axis='major', skipna=True):
        def f(arr):
            mask = common.notnull(arr)
            if skipna:
                return _tseries.median(arr[mask])
            else:
                if not mask.all():
                    return np.nan
                return _tseries.median(arr)
        return self.apply(f, axis=axis)

    _add_docs(median, 'median', 'median')

    def max(self, axis='major', skipna=True):
        i = self._get_axis_number(axis)

        y = np.array(self.values)
        if skipna:
            np.putmask(y, np.isnan(y), -np.inf)

        result = y.max(axis=i)
        result = np.where(np.isneginf(result), np.nan, result)
        return self._wrap_result(result, axis)

    _add_docs(max, 'maximum', 'maximum')

    def min(self, axis='major', skipna=True):
        i = self._get_axis_number(axis)

        y = np.array(self.values)
        if skipna:
            np.putmask(y, np.isnan(y), np.inf)

        result = y.min(axis=i)
        result = np.where(np.isinf(result), np.nan, result)
        return self._wrap_result(result, axis)

    _add_docs(min, 'minimum', 'minimum')

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

    def join(self, other, how=None, lsuffix='', rsuffix=''):
        """
        Join items with other Panel either on major and minor axes column

        Parameters
        ----------
        other : Panel
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
        if how is None:
            how = 'left'
        return self._join_index(other, how, lsuffix, rsuffix)

    def _join_index(self, other, how, lsuffix, rsuffix):
        join_major, join_minor = self._get_join_index(other, how)
        this = self.reindex(major=join_major, minor=join_minor)
        other = other.reindex(major=join_major, minor=join_minor)
        merged_data = this._data.merge(other._data, lsuffix, rsuffix)
        return self._constructor(merged_data)

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

    #----------------------------------------------------------------------
    # Deprecated stuff

    getMinorXS = deprecate('getMinorXS', minor_xs)
    getMajorXS = deprecate('getMajorXS', major_xs)

WidePanel = Panel

#-------------------------------------------------------------------------------
# LongPanel and friends


class LongPanel(DataFrame):
    """
    Represents long or "stacked" format panel data

    Parameters
    ----------
    values : ndarray (N x K)
    items : sequence
    index : MultiIndex

    Note
    ----
    LongPanel will likely disappear in a future release in favor of just using
    DataFrame objects with hierarchical indexes. You should be careful about
    writing production code depending on LongPanel
    """

    @property
    def consistent(self):
        offset = max(len(self.major_axis), len(self.minor_axis))

        major_labels = self.major_labels
        minor_labels = self.minor_labels

        # overflow risk
        if (offset + 1) ** 2 > 2**32:  # pragma: no cover
            major_labels = major_labels.astype(np.int64)
            minor_labels = minor_labels.astype(np.int64)

        keys = major_labels * offset + minor_labels
        unique_keys = np.unique(keys)

        if len(unique_keys) < len(keys):
            return False

        return True

    @property
    def wide_shape(self):
        return (len(self.items), len(self.major_axis), len(self.minor_axis))

    @property
    def items(self):
        return self.columns

    @property
    def _constructor(self):
        return LongPanel

    def __len__(self):
        return len(self.index)

    def __repr__(self):
        return DataFrame.__repr__(self)

    @classmethod
    def fromRecords(cls, data, major_field, minor_field,
                    exclude=None):
        """
        Create LongPanel from DataFrame or record / structured ndarray
        object

        Parameters
        ----------
        data : DataFrame, structured or record array, or dict
        major_field : string
        minor_field : string
            Name of field
        exclude : list-like, default None

        Returns
        -------
        LongPanel
        """
        return cls.from_records(data, [major_field, minor_field],
                                exclude=exclude)

    def toRecords(self):
        major = np.asarray(self.major_axis).take(self.major_labels)
        minor = np.asarray(self.minor_axis).take(self.minor_labels)

        arrays = [major, minor] + list(self.values[:, i]
                                       for i in range(len(self.items)))

        names = ['major', 'minor'] + list(self.items)

        return np.rec.fromarrays(arrays, names=names)

    @property
    def major_axis(self):
        return self.index.levels[0]

    @property
    def minor_axis(self):
        return self.index.levels[1]

    @property
    def major_labels(self):
        return self.index.labels[0]

    @property
    def minor_labels(self):
        return self.index.labels[1]

    def _combine(self, other, func, axis='items'):
        if isinstance(other, LongPanel):
            return self._combine_frame(other, func)
        elif isinstance(other, DataFrame):
            return self._combine_panel_frame(other, func, axis=axis)
        elif isinstance(other, Series):
            return self._combine_series(other, func, axis=axis)
        elif np.isscalar(other):
            return LongPanel(func(self.values, other), columns=self.items,
                             index=self.index)
        else:  # pragma: no cover
            raise Exception('type %s not supported' % type(other))

    def _combine_panel_frame(self, other, func, axis='items'):
        """
        Arithmetic op

        Parameters
        ----------
        other : DataFrame
        func : function
        axis : int / string

        Returns
        -------
        y : LongPanel
        """
        wide = self.to_wide()
        result = wide._combine_frame(other, func, axis=axis)
        return result.to_long()

    add = _panel_arith_method(operator.add, 'add')
    subtract = sub = _panel_arith_method(operator.sub, 'subtract')
    multiply = mul = _panel_arith_method(operator.mul, 'multiply')

    try:
        divide = div = _panel_arith_method(operator.div, 'divide')
    except AttributeError:   # Python 3
        divide = div = _panel_arith_method(operator.truediv, 'divide')

    def to_wide(self):
        """
        Transform long (stacked) format into wide format

        Returns
        -------
        Panel
        """
        assert(self.consistent)
        mask = make_mask(self.index)
        if self._data.is_mixed_dtype():
            return self._to_wide_mixed(mask)
        else:
            return self._to_wide_homogeneous(mask)

    def _to_wide_homogeneous(self, mask):
        values = np.empty(self.wide_shape, dtype=self.values.dtype)

        if not issubclass(self.values.dtype.type, np.integer):
            values.fill(np.nan)

        for i in xrange(len(self.items)):
            values[i].flat[mask] = self.values[:, i]

        return Panel(values, self.items, self.major_axis, self.minor_axis)

    def _to_wide_mixed(self, mask):
        _, N, K = self.wide_shape

        # TODO: make much more efficient

        data = {}
        for i, item in enumerate(self.items):
            item_vals = self[item].values

            values = np.empty((N, K), dtype=item_vals.dtype)
            values.ravel()[mask] = item_vals
            data[item] = DataFrame(values, index=self.major_axis,
                                   columns=self.minor_axis)
        return Panel.from_dict(data)

    toWide = deprecate('toWide', to_wide)

    def toCSV(self, path):
        def format_cols(items):
            cols = ['Major', 'Minor'] + list(items)
            return '"%s"' % '","'.join(cols)

        def format_row(major, minor, values):
            vals = ','.join('%.12f' % val for val in values)
            return '%s,%s,%s' % (major, minor, vals)

        f = open(path, 'w')
        self._textConvert(f, format_cols, format_row)
        f.close()

    def _textConvert(self, buf, format_cols, format_row):
        print >> buf, format_cols(self.items)

        label_pairs = zip(self.major_axis.take(self.major_labels),
                          self.minor_axis.take(self.minor_labels))
        for i, (major, minor) in enumerate(label_pairs):
            row = format_row(major, minor, self.values[i])
            print >> buf, row

    def swapaxes(self):
        """
        Swap major and minor axes and reorder values to be grouped by
        minor axis values

        Returns
        -------
        LongPanel (new object)
        """
        # Order everything by minor labels. Have to use mergesort
        # because NumPy quicksort is not stable. Here of course I'm
        # using the property that the major labels are ordered.
        indexer = self.minor_labels.argsort(kind='mergesort')

        new_major = self.minor_labels.take(indexer)
        new_minor = self.major_labels.take(indexer)
        new_values = self.values.take(indexer, axis=0)

        new_index = MultiIndex(levels=[self.minor_axis, self.major_axis],
                               labels=[new_major, new_minor])

        return LongPanel(new_values, columns=self.items,
                         index=new_index)

    def truncate(self, before=None, after=None):
        """
        Slice panel between two major axis values, return complete LongPanel

        Parameters
        ----------
        before : type of major_axis values or None, default None
            None defaults to start of panel

        after : type of major_axis values or None, default None
            None defaults to end of panel

        Returns
        -------
        LongPanel
        """
        left, right = self.index.slice_locs(before, after)
        new_index = self.index.truncate(before, after)

        return LongPanel(self.values[left : right],
                         columns=self.items, index=new_index)

    def get_axis_dummies(self, axis='minor', transform=None,
                         prefix=None):
        """
        Construct 1-0 dummy variables corresponding to designated axis
        labels

        Parameters
        ----------
        axis : {'major', 'minor'}, default 'minor'
        transform : function, default None

            Function to apply to axis labels first. For example, to
            get "day of week" dummies in a time series regression you might
            call:

                panel.get_axis_dummies(axis='major',
                                       transform=lambda d: d.weekday())
        Returns
        -------
        LongPanel, item names taken from chosen axis
        """
        if axis == 'minor':
            dim = len(self.minor_axis)
            items = self.minor_axis
            labels = self.minor_labels
        elif axis == 'major':
            dim = len(self.major_axis)
            items = self.major_axis
            labels = self.major_labels
        else: # pragma: no cover
            raise ValueError('Do not recognize axis %s' % axis)

        if transform:
            mapped = np.array([transform(val) for val in items])

            items = np.array(sorted(set(mapped)))
            labels = Index(items).get_indexer(mapped[labels])
            dim = len(items)

        values = np.eye(dim, dtype=float)
        values = values.take(labels, axis=0)

        result = LongPanel(values, columns=items, index=self.index)

        if prefix is None:
            prefix = ''

        result = result.add_prefix(prefix)

        return result

    def get_dummies(self, item):
        """
        Use unique values in column of panel to construct LongPanel
        containing dummy

        Parameters
        ----------
        item : object
            Value in panel items Index

        Returns
        -------
        LongPanel
        """
        idx = self.items.indexMap[item]
        values = self.values[:, idx]

        distinct_values = np.array(sorted(set(values)))
        mapping = distinct_values.searchsorted(values)

        values = np.eye(len(distinct_values))

        dummy_mat = values.take(mapping, axis=0)

        return LongPanel(dummy_mat, columns=distinct_values,
                         index=self.index)

    def mean(self, axis='major', broadcast=False):
        return self.apply(lambda x: np.mean(x, axis=0), axis, broadcast)

    def sum(self, axis='major', broadcast=False):
        return self.apply(lambda x: np.sum(x, axis=0), axis, broadcast)

    def apply(self, f, axis='major', broadcast=False):
        """
        Aggregate over a particular axis

        Parameters
        ----------
        f : function
            NumPy function to apply to each group
        axis : {'major', 'minor'}

        broadcast : boolean

        Returns
        -------
        broadcast=True  -> LongPanel
        broadcast=False -> DataFrame
        """
        try:
            return self._apply_level(f, axis=axis, broadcast=broadcast)
        except Exception:
            # ufunc
            new_values = f(self.values)
            return LongPanel(new_values, columns=self.items,
                             index=self.index)

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

    index = _get_combined_index(adj_frames, intersect=intersect)
    columns = _get_combined_columns(adj_frames, intersect=intersect)

    for key, frame in adj_frames.iteritems():
        result[key] = frame.reindex(index=index, columns=columns,
                                    copy=False)

    return result, index, columns

def _get_combined_columns(frames, intersect=False):
    columns = None

    if intersect:
        combine = set.intersection
    else:
        combine = set.union

    for _, frame in frames.iteritems():
        this_cols = set(frame.columns)

        if columns is None:
            columns = this_cols
        else:
            columns = combine(columns, this_cols)

    return Index(sorted(columns))

def _get_combined_index(frames, intersect=False):
    from pandas.core.frame import _union_indexes

    indexes = _get_distinct_indexes([df.index for df in frames.values()])
    if len(indexes) == 1:
        return indexes[0]
    if intersect:
        index = indexes[0]
        for other in indexes[1:]:
            index = index.intersection(other)
        return index
    union =  _union_indexes(indexes)
    return Index(union)

def _get_distinct_indexes(indexes):
    from itertools import groupby
    indexes = sorted(indexes, key=id)
    return [gp.next() for _, gp in groupby(indexes, id)]

def make_mask(index):
    """
    Create observation selection vector using major and minor
    labels, for converting to wide format.
    """
    N, K = index.levshape
    selector = index.labels[1] + K * index.labels[0]
    mask = np.zeros(N * K, dtype=bool)
    mask.put(selector, True)
    return mask

def _monotonic(arr):
    return not (arr[1:] < arr[:-1]).any()
