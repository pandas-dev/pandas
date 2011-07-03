"""
Contains data structures designed for manipulating panel (3-dimensional) data
"""
# pylint: disable=E1103,W0231,W0212,W0621

import operator
import sys
import warnings

import numpy as np

from pandas.core.common import (PandasError, _mut_exclusive,
                                _ensure_index, _pfixed)
from pandas.core.index import Index
from pandas.core.internals import BlockManager, make_block
from pandas.core.frame import DataFrame
from pandas.core.generic import PandasGeneric, Picklable
import pandas.core.common as common
import pandas._tseries as _tseries

class PanelError(Exception):
    pass

def _arith_method(func, name):
    # work only for scalars

    def f(self, other):
        if not np.isscalar(other):
            raise ValueError('Simple arithmetic with WidePanel can only be '
                            'done with scalar values')

        return self._combine(other, func)

    return f

def _long_arith_method(op, name):
    def f(self, other, axis='items'):
        """
        Wrapper method for %s

        Parameters
        ----------
        other : DataFrame or Panel class
        axis : {'items', 'major', 'minor'}

        Returns
        -------
        LongPanel
        """
        return self._combine(other, op, axis=axis)

    f.__name__ = name
    f.__doc__ = f.__doc__ % str(op)

    return f

def _wide_arith_method(op, name):
    def f(self, other, axis='items'):
        """
        Wrapper method for %s

        Parameters
        ----------
        other : DataFrame or Panel class
        axis : {'items', 'major', 'minor'}

        Returns
        -------
        WidePanel
        """
        return self._combine(other, op, axis=axis)

    f.__name__ = name
    f.__doc__ = f.__doc__ % str(op)

    return f

class PanelAxis(object):

    def __init__(self, cache_field):
        self.cache_field = cache_field

    def __get__(self, obj, type=None):
        return getattr(obj, self.cache_field, None)

    def __set__(self, obj, value):
        value = _ensure_index(value)
        setattr(obj, self.cache_field, value)

class Panel(object):
    """
    Abstract superclass for LongPanel and WidePanel data structures
    """
    _values = None
    factors = None

    __add__ = _arith_method(operator.add, '__add__')
    __sub__ = _arith_method(operator.sub, '__sub__')
    __mul__ = _arith_method(operator.mul, '__mul__')
    __div__ = _arith_method(operator.div, '__div__')
    __pow__ = _arith_method(operator.pow, '__pow__')

    __radd__ = _arith_method(operator.add, '__radd__')
    __rmul__ = _arith_method(operator.mul, '__rmul__')
    __rsub__ = _arith_method(lambda x, y: y - x, '__rsub__')
    __rdiv__ = _arith_method(lambda x, y: y / x, '__rdiv__')
    __rpow__ = _arith_method(lambda x, y: y ** x, '__rpow__')

    items = PanelAxis('_items')
    major_axis = PanelAxis('_major_axis')
    minor_axis = PanelAxis('_minor_axis')

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

        if self.factors:
            output += '\nFactors: %s' % ', '.join(self.factors)

        return output

    def __iter__(self):
        return iter(self.items)

    def iteritems(self):
        for item in self.items:
            yield item, self[item]

    @property
    def shape(self):
        return len(self.items), len(self.major_axis), len(self.minor_axis)

    @property
    def dims(self): # pragma: no cover
        warnings.warn("Please change panel.dims to panel.shape, will be removed"
                      " in future release",
                      FutureWarning)
        return self.shape

class AxisProperty(object):

    def __init__(self, axis=0):
        self.axis = axis

    def __get__(self, obj, type=None):
        data = getattr(obj, '_data')
        return data.axes[self.axis]

    def __set__(self, obj, value):
        data = getattr(obj, '_data')
        data.set_axis(self.axis, value)

class WidePanel(Panel, PandasGeneric):
    """
    Represents wide format panel data, stored as 3-dimensional array

    Parameters
    ----------
    values : ndarray (items x major x minor)
    items : sequence
    major_axis : sequence
    minor_axis : sequence
    """
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

    items = AxisProperty(0)
    major_axis = AxisProperty(1)
    minor_axis = AxisProperty(2)

    def __init__(self, data, items=None, major_axis=None, minor_axis=None,
                 copy=False, dtype=None):
        if isinstance(data, BlockManager):
            mgr = data
            if copy and dtype is None:
                mgr = mgr.copy()
            elif dtype is not None:
                # no choice but to copy
                mgr = mgr.cast(dtype)
        elif isinstance(data, np.ndarray):
            mgr = self._init_matrix(data, [items, major_axis, minor_axis],
                                    dtype=dtype, copy=copy)
        else:
            raise PandasError('Panel constructor not properly called!')

        self.factors = {}
        self._data = mgr

    def _consolidate_inplace(self):
        self._data = self._data.consolidate()

    def consolidate(self):
        """
        Compute DataFrame with "consolidated" internals (data of each dtype
        grouped together in a single ndarray). Mainly an internal API function,
        but available here to the savvy user

        Returns
        -------
        consolidated : DataFrame
        """
        cons_data = self._data.consolidate()
        if cons_data is self._data:
            cons_data = cons_data.copy()
        return type(self)(cons_data)

    def _init_matrix(self, data, axes, dtype=None, copy=False):
        values = _prep_ndarray(data, copy=copy)

        if dtype is not None:
            try:
                values = values.astype(dtype)
            except Exception:
                pass

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
        return BlockManager([block], axes)

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

    def copy(self):
        """
        Return a copy of WidePanel (only values ndarray copied)

        Returns
        -------
        y : WidePanel
        """
        return WidePanel(self._data.copy())

    @classmethod
    def from_dict(cls, data, intersect=False, dtype=float):
        """
        Construct WidePanel from dict of DataFrame objects

        Parameters
        ----------
        data : dict
            {field : DataFrame}
        intersect : boolean

        Returns
        -------
        WidePanel
        """
        data, index, columns = _homogenize(data, intersect=intersect)
        items = Index(sorted(data.keys()))
        values = np.array([data[k].values for k in items], dtype=dtype)
        return cls(values, items, index, columns)

    fromDict = from_dict

    def to_sparse(self, fill_value=None, kind='block'):
        """
        Convert to SparseWidePanel

        Parameters
        ----------
        fill_value : float, default NaN
        kind : {'block', 'integer'}

        Returns
        -------
        y : SparseDataFrame
        """
        from pandas.core.sparse import SparseWidePanel
        frames = dict(self.iteritems())
        return SparseWidePanel(frames, items=self.items,
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
        # old WidePanel pickle
        if isinstance(state, BlockManager):
            self._data = state
        elif len(state) == 4: # pragma: no cover
            self._unpickle_panel_compat(state)
        else:
            raise ValueError('unrecognized pickle')

    def _unpickle_panel_compat(self, state): # pragma: no cover
        "Unpickle the panel"
        _unpickle = common._unpickle_array
        vals, items, major, minor = state

        self.items = _unpickle(items)
        self.major_axis = _unpickle(major)
        self.minor_axis = _unpickle(minor)
        self.values = _unpickle(vals)

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
                major_axis=None, minor_axis=None):
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

        Returns
        -------
        WidePanel (new object)
        """
        result = self

        major = _mut_exclusive(major, major_axis)
        minor = _mut_exclusive(minor, minor_axis)

        if major is not None:
            result = result._reindex_axis(major, method, 1)

        if minor is not None:
            result = result._reindex_axis(minor, method, 2)

        if items is not None:
            result = result._reindex_axis(items, method, 0)

        if result is self:
            raise ValueError('Must specify at least one axis')

        return result

    def reindex_like(self, other, method=None):
        """
        Reindex WidePanel to match indices of another WidePanel

        Parameters
        ----------
        other : WidePanel
        method : string or None

        Returns
        -------
        reindexed : WidePanel
        """
        # todo: object columns
        return self.reindex(major=other.major_axis, items=other.items,
                            minor=other.minor_axis, method=method)

    def _reindex_axis(self, new_index, fill_method, axis):
        if axis == 0:
            new_data = self._data.reindex_items(new_index)
        else:
            new_data = self._data.reindex_axis(new_index, axis=axis,
                                               method=fill_method)
        return WidePanel(new_data)

    def _combine(self, other, func, axis=0):
        if isinstance(other, DataFrame):
            return self._combineFrame(other, func, axis=axis)
        elif isinstance(other, Panel):
            return self._combinePanel(other, func)
        elif np.isscalar(other):
            newValues = func(self.values, other)

            return WidePanel(newValues, self.items, self.major_axis,
                             self.minor_axis)

    def __neg__(self):
        return -1 * self

    def _combineFrame(self, other, func, axis=0):
        index, columns = self._get_plane_axes(axis)
        axis = self._get_axis_number(axis)

        other = other.reindex(index=index, columns=columns)

        if axis == 0:
            newValues = func(self.values, other.values)
        elif axis == 1:
            newValues = func(self.values.swapaxes(0, 1), other.values.T)
            newValues = newValues.swapaxes(0, 1)
        elif axis == 2:
            newValues = func(self.values.swapaxes(0, 2), other.values)
            newValues = newValues.swapaxes(0, 2)

        return WidePanel(newValues, self.items, self.major_axis,
                         self.minor_axis)

    def _combinePanel(self, other, func):
        if isinstance(other, LongPanel):
            other = other.to_wide()

        items = self.items + other.items
        major = self.major_axis + other.major_axis
        minor = self.minor_axis + other.minor_axis

        # could check that everything's the same size, but forget it

        this = self.reindex(items=items, major=major, minor=minor)
        other = other.reindex(items=items, major=major, minor=minor)

        result_values = func(this.values, other.values)

        return WidePanel(result_values, items, major, minor)

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
            for col, s in self.iteritems():
                result[col] = s.fillna(method=method, value=value)

            return WidePanel.from_dict(result)
        else:
            # Float type values
            if len(self.items) == 0:
                return self

            new_data = self._data.fillna(value)
            return WidePanel(new_data)

    add = _wide_arith_method(operator.add, 'add')
    subtract = _wide_arith_method(operator.sub, 'subtract')
    divide = _wide_arith_method(operator.div, 'divide')
    multiply = _wide_arith_method(operator.mul, 'multiply')

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
        loc = self.major_axis.get_loc(key)
        mat = np.array(self.values[:, loc, :].T)
        return DataFrame(mat, index=self.minor_axis, columns=self.items)

    def minor_xs(self, key):
        """
        Return slice of panel along minor axis

        Parameters
        ----------
        key : object
            Minor axis label

        Returns
        -------
        y : DataFrame
            index -> major axis, columns -> items
        """
        loc = self.minor_axis.get_loc(key)
        mat = np.array(self.values[:, :, loc].T)
        return DataFrame(mat, index=self.major_axis, columns=self.items)

    def getMinorXS(self, key): # pragma: no cover
        warnings.warn("getMinorXS has been replaced by the minor_xs function "
                      "please modify your code accordingly",
                      FutureWarning)

        return self.minor_xs(key)

    def getMajorXS(self, key): # pragma: no cover
        warnings.warn("getMajorXS has been replaced by the major_xs function "
                      "please modify your code accordingly",
                      FutureWarning)

        return self.major_xs(key)

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
        grouped : WidePanelGroupBy
        """
        from pandas.core.groupby import WidePanelGroupBy
        axis = self._get_axis_number(axis)
        return WidePanelGroupBy(self, function, axis=axis)

    def swapaxes(self, axis1='major', axis2='minor'):
        """
        Interchange axes and swap values axes appropriately

        Returns
        -------
        y : WidePanel (new object)
        """
        i = self._get_axis_number(axis1)
        j = self._get_axis_number(axis2)

        if i == j:
            raise ValueError('Cannot specify the same axis')

        mapping = {i : j, j : i}

        new_axes = (self._get_axis(mapping.get(k, k))
                    for k in range(3))
        new_values = self.values.swapaxes(i, j).copy()

        return WidePanel(new_values, *new_axes)

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
            mask = np.isfinite(self.values).all(axis=0)
            size = mask.sum()
            selector = mask.ravel()
        else:
            size = N * K
            selector = slice(None, None)

        values = np.empty((size, I), dtype=float)

        for i in xrange(len(self.items)):
            values[:, i] = self.values[i].ravel()[selector]

        major_labels = np.arange(N).repeat(K)[selector]

        # Anyone think of a better way to do this? np.repeat does not
        # do what I want
        minor_labels = np.arange(K).reshape(1, K)[np.zeros(N, dtype=int)]
        minor_labels = minor_labels.ravel()[selector]

        if filter_observations:
            mask = selector
        else:
            mask = None

        index = LongPanelIndex(self.major_axis,
                               self.minor_axis,
                               major_labels,
                               minor_labels,
                               mask=mask)

        return LongPanel(values, self.items, index)

    toLong = to_long

    def filter(self, items):
        """
        Restrict items in panel to input list

        Parameters
        ----------
        items : sequence

        Returns
        -------
        y : WidePanel
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
        result : DataFrame or WidePanel
        """
        i = self._get_axis_number(axis)

        result = np.apply_along_axis(func, i, self.values)

        return self._wrap_result(result, axis=axis)

    def _values_aggregate(self, func, axis, fill_value):
        axis = self._get_axis_number(axis)

        values = self.values
        mask = np.isfinite(values)

        if fill_value is not None:
            values = values.copy()
            values[-mask] = fill_value

        result = func(values, axis=axis)
        count = mask.sum(axis=axis)

        result[count == 0] = np.NaN

        return result

    def _values_accum(self, func, axis, fill_value):
        axis = self._get_axis_number(axis)

        values = self.values
        mask = np.isfinite(values)

        if fill_value is not None:
            values = values.copy()
            values[-mask] = fill_value

        result = func(values, axis=axis)

        if fill_value is not None:
            result[-mask] = np.NaN

        return result

    def _array_method(self, func, axis='major', fill_value=None):
        """
        Parameters
        ----------
        func : numpy function
            Signature should match numpy.{sum, mean, var, std} etc.
        axis : {'major', 'minor', 'items'}
        fill_value : boolean, default True
            Replace NaN values with specified first

        Returns
        -------
        y : DataFrame
        """
        result = self._values_aggregate(func, axis, fill_value)
        return self._wrap_result(result, axis=axis)

    def _wrap_result(self, result, axis):
        axis = self._get_axis_name(axis)

        if result.ndim == 2:
            index, columns = self._get_plane_axes(axis)

            if axis != 'items':
                result = result.T

            return DataFrame(result, index=index, columns=columns)
        else:
            return WidePanel(result, self.items, self.major_axis,
                             self.minor_axis)

    def count(self, axis='major'):
        """
        Return DataFrame of observation counts along desired axis

        Returns
        -------
        y : DataFrame
        """
        i = self._get_axis_number(axis)

        values = self.values
        mask = np.isfinite(values)
        result = mask.sum(axis=i)

        return self._wrap_result(result, axis)

    def sum(self, axis='major'):
        """

        Returns
        -------
        y : DataFrame
        """
        return self._array_method(np.sum, axis=axis, fill_value=0)

    def cumsum(self, axis='major'):
        """

        Returns
        -------
        y : WidePanel
        """
        result = self._values_accum(np.cumsum, axis=axis, fill_value=0)
        return self._wrap_result(result, axis)

    def mean(self, axis='major'):
        """

        Returns
        -------
        y : DataFrame
        """
        return self.sum(axis=axis) / self.count(axis=axis)

    def var(self, axis='major'):
        """

        Returns
        -------
        y : DataFrame
        """
        i = self._get_axis_number(axis)
        index, columns = self._get_plane_axes(axis)

        y = np.array(self.values)
        mask = np.isnan(y)

        count = (-mask).sum(axis=i).astype(float)
        y[mask] = 0

        X = y.sum(axis=i)
        XX = (y ** 2).sum(axis=i)

        theVar = (XX - X**2 / count) / (count - 1)

        return self._wrap_result(theVar, axis)

    def std(self, axis='major'):
        """

        Returns
        -------
        y : DataFrame
        """
        return self.var(axis=axis).apply(np.sqrt)

    def skew(self, axis='major'):
        raise NotImplementedError

    def prod(self, axis='major'):
        """

        Returns
        -------
        y : DataFrame
        """
        return self._array_method(np.prod, axis=axis, fill_value=1)

    def compound(self, axis='major'):
        """

        Returns
        -------
        y : DataFrame
        """
        return (1 + self).prod(axis=axis) - 1

    def median(self, axis='major'):
        """

        Returns
        -------
        y : DataFrame
        """
        def f(arr):
            return _tseries.median(arr[common.notnull(arr)])

        return self.apply(f, axis=axis)

    def max(self, axis='major'):
        """

        Returns
        -------
        y : DataFrame
        """
        i = self._get_axis_number(axis)

        y = np.array(self.values)
        mask = np.isfinite(y)

        fill_value = y.flat[mask.ravel()].min() - 1

        y[-mask] = fill_value

        result = y.max(axis=i)
        result[result == fill_value] = np.NaN

        return self._wrap_result(result, axis)

    def min(self, axis='major'):
        """

        Returns
        -------
        y : DataFrame
        """
        i = self._get_axis_number(axis)

        y = np.array(self.values)
        mask = np.isfinite(y)

        fill_value = y.flat[mask.ravel()].max() + 1

        y[-mask] = fill_value

        result = y.min(axis=i)
        result[result == fill_value] = np.NaN

        return self._wrap_result(result, axis)

    def shift(self, lags, axis='major'):
        """

        Returns
        -------
        y : WidePanel
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

        return WidePanel(values, items=items, major_axis=major_axis,
                         minor_axis=minor_axis)

    def truncate(self, before=None, after=None, axis='major'):
        """Function truncates a sorted Panel before and/or after
        some particular dates

        Parameters
        ----------
        before : date
            Left boundary
        after : date
            Right boundary

        Returns
        -------
        WidePanel
        """
        axis = self._get_axis_name(axis)
        index = self._get_axis(axis)

        beg_slice, end_slice = index.slice_locs(before, after)
        new_index = index[beg_slice:end_slice]

        return self.reindex(**{axis : new_index})

    def select(self, crit, axis=0):
        """
        Return data corresponding to axis labels matching criteria

        Parameters
        ----------
        crit : function
            To be called on each index (label). Should return True or False
        axis : {0, 1, 2} or {'items', 'major', 'minor'}, default 'items'
            Axis to select on

        Returns
        -------
        selection : DataFrame
        """
        return self._select_generic(crit, axis=axis)

#-------------------------------------------------------------------------------
# LongPanel and friends


class LongPanel(Panel, Picklable):
    """
    Represents long or "stacked" format panel data

    Parameters
    ----------
    values : ndarray (N x K)
    items : sequence
    index : LongPanelIndex

    Note
    ----
    Constructor should probably not be called directly since it
    requires creating the major and minor axis label vectors for for
    the LongPanelIndex
    """

    def __init__(self, values, items, index, factors=None):
        self.items = items
        self.index = index

        self.values = values

        self.factors = factors or {}

    def __len__(self):
        return len(self.index)

    @classmethod
    def fromRecords(cls, data, major_field, minor_field,
                    factors=None, exclude=None):
        """
        Create LongPanel from DataFrame or record / structured ndarray
        object

        Parameters
        ----------
        data : DataFrame, structured or record array, or dict
        major_field : string
        minor_field : string
            Name of field
        factors : list-like, default None
        exclude : list-like, default None

        Returns
        -------
        LongPanel
        """
        if isinstance(data, np.ndarray):
            # Dtype when you have data
            if not issubclass(data.dtype.type, np.void):
                raise ValueError('Input was not a structured array!')

            columns = data.dtype.names
            data = dict((k, data[k]) for k in columns)
        elif isinstance(data, DataFrame):
            data = data._series.copy()
        elif isinstance(data, dict):
            # otherwise will pop columns out of original
            data = data.copy()

        if exclude is None:
            exclude = set()
        else:
            exclude = set(exclude)

        major_vec = data.pop(major_field)
        minor_vec = data.pop(minor_field)

        major_axis = Index(sorted(set(major_vec)))
        minor_axis = Index(sorted(set(minor_vec)))

        major_labels, _ = _tseries.getMergeVec(major_vec, major_axis.indexMap)
        minor_labels, _ = _tseries.getMergeVec(minor_vec, minor_axis.indexMap)

        for col in exclude:
            del data[col]

        factor_dict = {}
        for col in data.keys():
            series = data[col]

            # Is it a factor?
            if not np.issctype(series.dtype):
                factor_dict[col] = factor = Factor.fromarray(series)
                data[col] = factor.labels

        items = sorted(data)
        values = np.array([data[k] for k in items]).T

        index = LongPanelIndex(major_axis, minor_axis,
                               major_labels, minor_labels)

        return LongPanel(values, items, index, factors=factor_dict)

    def toRecords(self):
        major = np.asarray(self.major_axis).take(self.index.major_labels)
        minor = np.asarray(self.minor_axis).take(self.index.minor_labels)

        arrays = [major, minor] + list(self.values[:, i]
                                       for i in range(len(self.items)))

        names = ['major', 'minor'] + list(self.items)

        return np.rec.fromarrays(arrays, names=names)

    @property
    def columns(self):
        """
        So LongPanel can be DataFrame-like at times
        """
        return self.items

    def copy(self):
        """
        Return copy of LongPanel (copies ndarray)

        Returns
        -------
        y : LongPanel
        """
        return LongPanel(self.values.copy(), self.items, self.index,
                         factors=self.factors)

    @property
    def major_axis(self):
        return self.index.major_axis

    @property
    def minor_axis(self):
        return self.index.minor_axis

    def _get_values(self):
        return self._values

    def _set_values(self, values):
        if not values.flags.contiguous:
            values = values.copy()

        shape = len(self.index.major_labels), len(self.items)

        if values.shape != shape:
            raise ValueError('Values shape %s mismatch to %s' % (values.shape,
                                                                shape))

        self._values = values

    values = property(fget=_get_values, fset=_set_values)

    def __getitem__(self, key):
        "Return column of panel as LongPanel"
        loc = self.items.get_loc(key)
        return LongPanel(self.values[:, loc : loc + 1].copy(),
                        [key], self.index, factors=self.factors)

    def __setitem__(self, key, value):
        if np.isscalar(value):
            mat = np.empty((len(self.values), 1), dtype=float)
            mat.fill(value)
        elif isinstance(value, np.ndarray):
            mat = value
        elif isinstance(value, LongPanel):
            if len(value.items) > 1:
                raise ValueError('input LongPanel must only have one column')

            if value.index is not self.index:
                raise ValueError('Only can set identically-indexed LongPanel '
                                'items for now')

            mat = value.values

        # Insert item at end of items for now
        self.items = Index(list(self.items) + [key])
        self.values = np.column_stack((self.values, mat))

    def __getstate__(self):
        "Returned pickled representation of the panel"

        return (common._pickle_array(self.values),
                common._pickle_array(self.items),
                self.index)

    def __setstate__(self, state):
        "Unpickle the panel"
        (vals, items, index) = state

        self.items = common._unpickle_array(items)
        self.index = index
        self.values = common._unpickle_array(vals)

    def _combine(self, other, func, axis='items'):
        if isinstance(other, DataFrame):
            return self._combineFrame(other, func, axis=axis)
        elif isinstance(other, Panel):
            return self._combinePanel(other, func)
        elif np.isscalar(other):
            return LongPanel(func(self.values, other), self.items,
                             self.index, factors=self.factors)

    def _combineFrame(self, other, func, axis='items'):
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
        result = wide._combineFrame(other, func, axis=axis)
        return result.to_long()

    def _combinePanel(self, other, func):
        """
        Arithmetic operation between panels
        """
        if self.index is not other.index:
            raise ValueError("Can only combine identically-indexed "
                            "panels for now")

        if len(other.items) == 1:
            new_values = func(self.values, other.values)
        else:
            new_values = func(self.values, other.values)

        return LongPanel(new_values, self.items, self.index,
                         factors=self.factors)

    add = _long_arith_method(operator.add, 'add')
    subtract = _long_arith_method(operator.sub, 'subtract')
    divide = _long_arith_method(operator.div, 'divide')
    multiply = _long_arith_method(operator.mul, 'multiply')

    def sort(self, axis='major'):
        """
        Sort value by chosen axis (break ties using other axis)

        Note
        ----
        A LongPanel must be sorted to convert to a WidePanel

        Returns
        -------
        LongPanel (in sorted order)
        """
        if axis == 'major':
            first = self.index.major_labels
            second = self.index.minor_labels

        elif axis == 'minor':
            first = self.index.minor_labels
            second = self.index.major_labels

        # Lexsort starts from END
        indexer = np.lexsort((second, first))

        new_major = self.index.major_labels[indexer]
        new_minor = self.index.minor_labels[indexer]
        new_values = self.values[indexer]

        new_index = LongPanelIndex(self.major_axis, self.minor_axis,
                                   new_major, new_minor)

        new_factors = dict((k, v[indexer])
                           for k, v in self.factors.iteritems())

        return LongPanel(new_values, self.items, new_index,
                         factors=new_factors)

    def to_wide(self):
        """
        Transform long (stacked) format into wide format

        Returns
        -------
        WidePanel
        """
        if not self.index.consistent:
            raise PanelError('Panel has duplicate (major, minor) pairs, '
                             'cannot be reliably converted to wide format.')

        I, N, K = self.shape

        values = np.empty((I, N, K), dtype=self.values.dtype)

        mask = self.index.mask
        notmask = -mask

        for i in xrange(len(self.items)):
            values[i].flat[mask] = self.values[:, i]
            values[i].flat[notmask] = np.NaN

        return WidePanel(values, self.items, self.major_axis, self.minor_axis)

    toWide = to_wide

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

    def toString(self, buf=sys.stdout, col_space=15):
        """
        Output a screen-friendly version of this Panel
        """
        _pf = _pfixed
        major_space = max(max([len(str(idx))
                               for idx in self.major_axis]) + 4, 9)
        minor_space = max(max([len(str(idx))
                               for idx in self.minor_axis]) + 4, 9)

        def format_cols(items):
            return '%s%s%s' % (_pf('Major', major_space),
                               _pf('Minor', minor_space),
                               ''.join(_pf(h, col_space) for h in items))

        def format_row(major, minor, values):
            return '%s%s%s' % (_pf(major, major_space),
                               _pf(minor, minor_space),
                               ''.join(_pf(v, col_space) for v in values))

        self._textConvert(buf, format_cols, format_row)

    def _textConvert(self, buf, format_cols, format_row):
        print >> buf, format_cols(self.items)

        label_pairs = zip(self.index.major_labels,
                          self.index.minor_labels)
        major, minor = self.major_axis, self.minor_axis
        for i, (major_i, minor_i) in enumerate(label_pairs):
            row = format_row(major[major_i], minor[minor_i], self.values[i])
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
        indexer = self.index.minor_labels.argsort(kind='mergesort')

        new_major = self.index.minor_labels.take(indexer)
        new_minor = self.index.major_labels.take(indexer)

        new_values = self.values.take(indexer, axis=0)

        new_index = LongPanelIndex(self.minor_axis,
                                   self.major_axis,
                                   new_major,
                                   new_minor,
                                   mask=self.index.mask)

        return LongPanel(new_values, self.items, new_index)

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
        left, right = self.index.get_major_bounds(before, after)
        new_index = self.index.truncate(before, after)

        return LongPanel(self.values[left : right],
                         self.items, new_index)

    def filter(self, items):
        """
        Restrict items in panel to input list

        Parameters
        ----------
        items : sequence

        Returns
        -------
        WidePanel
        """
        intersection = self.items.intersection(items)
        indexer = [self.items.indexMap[col] for col in intersection]

        new_values = self.values.take(indexer, axis=1)
        return LongPanel(new_values, intersection, self.index)

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
            labels = self.index.minor_labels
        elif axis == 'major':
            dim = len(self.major_axis)
            items = self.major_axis
            labels = self.index.major_labels
        else: # pragma: no cover
            raise ValueError('Do not recognize axis %s' % axis)

        if transform:
            mapped = np.array([transform(val) for val in items])

            items = np.array(sorted(set(mapped)))
            labels = items.searchsorted(mapped[labels])
            dim = len(items)

        values = np.eye(dim, dtype=float)
        values = values.take(labels, axis=0)

        result = LongPanel(values, items, self.index)

        if prefix is None:
            prefix = ''

        result = result.addPrefix(prefix)

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

        return LongPanel(dummy_mat, distinct_values, self.index)

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
            return self._apply_axis(f, axis=axis, broadcast=broadcast)
        except Exception:
            # ufunc
            new_values = f(self.values)
            return LongPanel(new_values, self.items, self.index)

    def _apply_axis(self, f, axis='major', broadcast=False):
        if axis == 'major':
            panel = self.swapaxes()
            result = panel._apply_axis(f, axis='minor', broadcast=broadcast)
            if broadcast:
                result = result.swapaxes()

            return result

        bounds = self.index._bounds
        values = self.values
        N, _ = values.shape
        result = group_agg(values, bounds, f)

        if broadcast:
            repeater = np.concatenate((np.diff(bounds), [N - bounds[-1]]))
            panel = LongPanel(result.repeat(repeater, axis=0),
                              self.items, self.index)
        else:
            panel = DataFrame(result, index=self.major_axis,
                               columns=self.items)

        return panel

    def count(self, axis='major'):
        """
        Compute observation counts within each group

        Parameters
        ----------
        axis : {'major', 'minor'}
            major: compute minor_axis obs for each major axis value
            minor: same but for each minor axis value

        Returns
        -------
        counts : ndarray (1d)
            Length will be length of input axis
        """
        if axis == 'major':
            lp = self
        elif axis == 'minor':
            lp = self.swapaxes()
        else: # pragma: no cover
            raise ValueError('invalid axis')

        N = len(lp.values)
        bounds = lp.index._bounds

        return np.concatenate((np.diff(bounds), [N - bounds[-1]]))

    def leftJoin(self, other):
        """

        Parameters
        ----------
        other : LongPanel
        """
        assert(self.index is other.index)

        values = np.concatenate((self.values, other.values), axis=1).copy()
        items = self.items.tolist() + other.items.tolist()

        return LongPanel(values, items, self.index)

    def addPrefix(self, prefix=None):
        """
        Concatenate prefix string with panel items names.

        Parameters
        ----------
        prefix : string

        Returns
        -------
        LongPanel

        Note
        ----
        does *not* copy values matrix
        """
        new_items = [_prefix_item(item, prefix) for item in self.items]

        return LongPanel(self.values, new_items, self.index)


class LongPanelIndex(object):
    """
    Holds axis indexing information for a LongPanel instance

    Parameters
    ----------
    major_axis : Index-like
    minor_axis : Index-like
    major_labels : ndarray
    minor_labels : ndarray
    mask : ndarray (bool), optional
        observation selection vector using major and minor labels, for
        converting to wide format.
    """
    def __init__(self, major_axis, minor_axis, major_labels,
                 minor_labels, mask=None):

        self.major_axis = major_axis
        self.minor_axis = minor_axis

        assert(len(minor_labels) == len(major_labels))

        self.major_labels = major_labels
        self.minor_labels = minor_labels

        self._mask = mask

    def __len__(self):
        return len(self.major_labels)

    def __getstate__(self):
        _pickle = common._pickle_array
        return (_pickle(self.major_axis),
                _pickle(self.minor_axis),
                _pickle(self.major_labels),
                _pickle(self.minor_labels))

    def __setstate__(self, state):
        _unpickle = common._unpickle_array

        major, minor, major_labels, minor_labels = state

        self.major_axis = _unpickle(major)
        self.minor_axis = _unpickle(minor)

        self.major_labels = _unpickle(major_labels)
        self.minor_labels = _unpickle(minor_labels)

    @property
    def consistent(self):
        offset = max(len(self.major_axis), len(self.minor_axis))

        # overflow risk
        if (offset + 1) ** 2 > 2**32:
            keys = (self.major_labels.astype(np.int64) * offset +
                    self.minor_labels.astype(np.int64))
        else:
            keys = self.major_labels * offset + self.minor_labels

        unique_keys = np.unique(keys)

        if len(unique_keys) < len(keys):
            return False

        return True

    def truncate(self, before=None, after=None):
        """
        Slice index between two major axis values, return new
        LongPanelIndex

        Parameters
        ----------
        before : type of major_axis values or None, default None
            None defaults to start of panel

        after : type of major_axis values or None, default None
            None defaults to after of panel

        Returns
        -------
        LongPanelIndex
        """
        i, j = self._get_axis_bounds(before, after)
        left, right = self._get_label_bounds(i, j)

        return LongPanelIndex(self.major_axis[i : j],
                              self.minor_axis,
                              self.major_labels[left : right] - i,
                              self.minor_labels[left : right])

    def get_major_bounds(self, begin=None, end=None):
        """
        Return index bounds for slicing LongPanel labels and / or
        values

        Parameters
        ----------
        begin : axis value or None
        end : axis value or None

        Returns
        -------
        y : tuple
            (left, right) absolute bounds on LongPanel values
        """
        i, j = self._get_axis_bounds(begin, end)
        left, right = self._get_label_bounds(i, j)

        return left, right

    def _get_axis_bounds(self, begin, end):
        """
        Return major axis locations corresponding to interval values
        """
        if begin is not None:
            i = self.major_axis.indexMap.get(begin)
            if i is None:
                i = self.major_axis.searchsorted(begin, side='right')
        else:
            i = 0

        if end is not None:
            j = self.major_axis.indexMap.get(end)
            if j is None:
                j = self.major_axis.searchsorted(end)
            else:
                j = j + 1
        else:
            j = len(self.major_axis)

        if i > j:
            raise ValueError('Must have begin <= end!')

        return i, j

    def _get_label_bounds(self, i, j):
        "Return slice points between two major axis locations"

        left = self._bounds[i]

        if j >= len(self.major_axis):
            right = len(self.major_labels)
        else:
            right = self._bounds[j]

        return left, right

    __bounds = None
    @property
    def _bounds(self):
        "Return or compute and return slice points for major axis"
        if self.__bounds is None:
            inds = np.arange(len(self.major_axis))
            self.__bounds = self.major_labels.searchsorted(inds)

        return self.__bounds

    @property
    def mask(self):
        if self._mask is None:
            self._mask = self._make_mask()

        return self._mask

    def _make_mask(self):
        """
        Create observation selection vector using major and minor
        labels, for converting to wide format.
        """
        N, K = self.shape
        selector = self.minor_labels + K * self.major_labels

        mask = np.zeros(N * K, dtype=bool)
        mask[selector] = True

        return mask

    @property
    def shape(self):
        return len(self.major_axis), len(self.minor_axis)

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

def _infer_dtype(value):
    if isinstance(value, (float, np.floating)):
        return float
    elif isinstance(value, (int, np.integer)):
        return int
    elif isinstance(value, (bool, np.bool_)):
        return bool
    else:
        return object

class Factor(object):
    """
    Represents a categorical variable in classic R / S-plus fashion
    """
    def __init__(self, labels, levels):
        self.labels = labels
        self.levels = levels

    @classmethod
    def fromarray(cls, values):
        levels = np.array(sorted(set(values)), dtype=object)
        labels = levels.searchsorted(values)

        return Factor(labels, levels=levels)

    def asarray(self):
        return self.levels[self.labels]

    def __len__(self):
        return len(self.labels)

    def __repr__(self):
        temp = 'Factor:\n%s\nLevels (%d): %s'
        values = self.asarray()
        return temp % (repr(values), len(self.levels), self.levels)

    def __getitem__(self, key):
        if isinstance(key, (int, np.integer)):
            i = self.labels[key]
            return self.levels[i]
        else:
            new_labels = self.labels[key]
            return Factor(new_labels, self.levels)

def factor_agg(factor, vec, func):
    """
    Aggregate array based on Factor

    Parameters
    ----------
    factor : Factor
        length n
    vec : sequence
        length n
    func : function
        1D array aggregation function

    Returns
    -------
    ndarray corresponding to Factor levels
    """
    indexer = np.argsort(factor.labels)
    unique_labels = np.arange(len(factor.levels))

    ordered_labels = factor.labels.take(indexer)
    ordered_vec = np.asarray(vec).take(indexer)
    bounds = ordered_labels.searchsorted(unique_labels)

    return group_agg(ordered_vec, bounds, func)

def group_agg(values, bounds, f):
    """
    R-style aggregator

    Parameters
    ----------
    values : N-length or N x K ndarray
    bounds : B-length ndarray
    f : ndarray aggregation function

    Returns
    -------
    ndarray with same length as bounds array
    """
    if values.ndim == 1:
        N = len(values)
        result = np.empty(len(bounds), dtype=float)
    elif values.ndim == 2:
        N, K = values.shape
        result = np.empty((len(bounds), K), dtype=float)

    testagg = f(values[:min(1, len(values))])
    if isinstance(testagg, np.ndarray) and testagg.ndim == 2:
        raise Exception('Passed function does not aggregate!')

    for i, left_bound in enumerate(bounds):
        if i == len(bounds) - 1:
            right_bound = N
        else:
            right_bound = bounds[i + 1]

        result[i] = f(values[left_bound : right_bound])

    return result

def _prefix_item(item, prefix=None):
    if prefix is None:
        return item

    template = '%s%s'
    return template % (prefix, item)

def _homogenize(frames, intersect=True):
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
        result[key] = frame.reindex(index=index, columns=columns)

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
    index = None

    if intersect:
        combine = Index.intersection
    else:
        combine = Index.union

    for _, frame in frames.iteritems():
        if index is None:
            index = frame.index
        elif index is not frame.index:
            index = combine(index, frame.index)

    return index

def pivot(index, columns, values):
    """
    Produce 'pivot' table based on 3 columns of this DataFrame.
    Uses unique values from index / columns and fills with values.

    Parameters
    ----------
    index : ndarray
        Labels to use to make new frame's index
    columns : ndarray
        Labels to use to make new frame's columns
    values : ndarray
        Values to use for populating new frame's values

    Note
    ----
    Obviously, all 3 of the input arguments must have the same length

    Returns
    -------
    DataFrame
    """
    assert(len(index) == len(columns) == len(values))

    if len(index) == 0:
        return DataFrame(index=[])

    try:
        longIndex = _make_long_index(index, columns)
        valueMat = values.view(np.ndarray).reshape(len(values), 1)
        longPanel = LongPanel(valueMat, ['foo'], longIndex)
        longPanel = longPanel.sort()
        return longPanel.to_wide()['foo']
    except PanelError:
        return _slow_pivot(index, columns, values)

def _make_long_index(major_values, minor_values):
    major_axis = Index(sorted(set(major_values)))
    minor_axis = Index(sorted(set(minor_values)))

    major_labels, _ = _tseries.getMergeVec(major_values, major_axis.indexMap)
    minor_labels, _ = _tseries.getMergeVec(minor_values, minor_axis.indexMap)

    long_index = LongPanelIndex(major_axis, minor_axis,
                               major_labels, minor_labels)
    return long_index

def _slow_pivot(index, columns, values):
    """
    Produce 'pivot' table based on 3 columns of this DataFrame.
    Uses unique values from index / columns and fills with values.

    Parameters
    ----------
    index : string or object
        Column name to use to make new frame's index
    columns : string or object
        Column name to use to make new frame's columns
    values : string or object
        Column name to use for populating new frame's values

    Could benefit from some Cython here.
    """
    from itertools import izip
    tree = {}
    for i, (idx, col) in enumerate(izip(index, columns)):
        if col not in tree:
            tree[col] = {}
        branch = tree[col]
        branch[idx] = values[i]

    return DataFrame(tree)

def _monotonic(arr):
    return not (arr[1:] < arr[:-1]).any()
