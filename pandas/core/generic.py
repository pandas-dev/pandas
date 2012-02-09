# pylint: disable=W0231

import numpy as np

from pandas.core.common import save, load
from pandas.core.index import MultiIndex
import pandas.core.datetools as datetools

#-------------------------------------------------------------------------------
# Picklable mixin

class Picklable(object):

    def save(self, path):
        save(self, path)

    @classmethod
    def load(cls, path):
        return load(path)

class PandasError(Exception):
    pass

class PandasObject(Picklable):

    _AXIS_NUMBERS = {
        'index' : 0,
        'columns' : 1
    }

    _AXIS_ALIASES = {}
    _AXIS_NAMES = dict((v, k) for k, v in _AXIS_NUMBERS.iteritems())

    #----------------------------------------------------------------------
    # Axis name business

    @classmethod
    def _get_axis_number(cls, axis):
        axis = cls._AXIS_ALIASES.get(axis, axis)

        if isinstance(axis, int):
            if axis in cls._AXIS_NAMES:
                return axis
            else:
                raise Exception('No %d axis' % axis)
        else:
            return cls._AXIS_NUMBERS[axis]

    @classmethod
    def _get_axis_name(cls, axis):
        axis = cls._AXIS_ALIASES.get(axis, axis)
        if isinstance(axis, basestring):
            if axis in cls._AXIS_NUMBERS:
                return axis
            else:
                raise Exception('No axis named %s' % axis)
        else:
            return cls._AXIS_NAMES[axis]

    def _get_axis(self, axis):
        name = self._get_axis_name(axis)
        return getattr(self, name)

    def abs(self):
        """
        Return an object with absolute value taken. Only applicable to objects
        that are all numeric

        Returns
        -------
        abs: type of caller
        """
        return np.abs(self)

    def get(self, key, default=None):
        """
        Get item from object for given key (DataFrame column, Panel slice,
        etc.). Returns default value if not found

        Parameters
        ----------
        key : object

        Returns
        -------
        value : type of items contained in object
        """
        try:
            return self[key]
        except KeyError:
            return default

    def groupby(self, by=None, axis=0, level=None, as_index=True, sort=True):
        """
        Group series using mapper (dict or key function, apply given function
        to group, return result as series) or by a series of columns

        Parameters
        ----------
        by : mapping function / list of functions, dict, Series, or tuple /
            list of column names.
            Called on each element of the object index to determine the groups.
            If a dict or Series is passed, the Series or dict VALUES will be
            used to determine the groups
        axis : int, default 0
        level : int, level name, or sequence of such, default None
            If the axis is a MultiIndex (hierarchical), group by a particular
            level or levels
        as_index : boolean, default True
            For aggregated output, return object with group labels as the
            index. Only relevant for DataFrame input. as_index=False is
            effectively "SQL-style" grouped output
        sort : boolean, default True
            Sort group keys. Get better performance by turning this off

        Examples
        --------
        # DataFrame result
        >>> data.groupby(func, axis=0).mean()

        # DataFrame result
        >>> data.groupby(['col1', 'col2'])['col3'].mean()

        # DataFrame with hierarchical index
        >>> data.groupby(['col1', 'col2']).mean()

        Returns
        -------
        GroupBy object
        """
        from pandas.core.groupby import groupby
        return groupby(self, by, axis=axis, level=level, as_index=as_index,
                       sort=sort)

    def select(self, crit, axis=0):
        """
        Return data corresponding to axis labels matching criteria

        Parameters
        ----------
        crit : function
            To be called on each index (label). Should return True or False
        axis : int

        Returns
        -------
        selection : type of caller
        """
        axis_name = self._get_axis_name(axis)
        axis = self._get_axis(axis)

        if len(axis) > 0:
            new_axis = axis[np.asarray([crit(label) for label in axis])]
        else:
            new_axis = axis

        return self.reindex(**{axis_name : new_axis})

    def drop(self, labels, axis=0):
        """
        Return new object with labels in requested axis removed

        Parameters
        ----------
        labels : array-like
        axis : int

        Returns
        -------
        dropped : type of caller
        """
        axis_name = self._get_axis_name(axis)
        axis = self._get_axis(axis)
        new_axis = axis.drop(labels)
        return self.reindex(**{axis_name : new_axis})

    def sort_index(self, axis=0, ascending=True):
        """
        Sort object by labels (along an axis)

        Parameters
        ----------
        axis : {0, 1}
            Sort index/rows versus columns
        ascending : boolean, default True
            Sort ascending vs. descending

        Returns
        -------
        sorted_obj : type of caller
        """
        axis = self._get_axis_number(axis)
        axis_name = self._get_axis_name(axis)
        labels = self._get_axis(axis)

        sort_index = labels.argsort()
        if not ascending:
            sort_index = sort_index[::-1]

        new_axis = labels.take(sort_index)
        return self.reindex(**{axis_name : new_axis})

    @property
    def ix(self):
        raise NotImplementedError

    def reindex(self, *args, **kwds):
        raise NotImplementedError

class NDFrame(PandasObject):
    """
    N-dimensional analogue of DataFrame. Store multi-dimensional in a
    size-mutable, labeled data structure

    Parameters
    ----------
    data : BlockManager
    axes : list
    copy : boolean, default False
    """
    # kludge
    _default_stat_axis = 0

    def __init__(self, data, axes=None, copy=False, dtype=None):
        if dtype is not None:
            data = data.astype(dtype)
        elif copy:
            data = data.copy()

        if axes is not None:
            for i, ax in enumerate(axes):
                data = data.reindex_axis(ax, axis=i)

        self._data = data
        self._item_cache = {}

    def astype(self, dtype):
        """
        Cast object to input numpy.dtype

        Parameters
        ----------
        dtype : numpy.dtype or Python type

        Returns
        -------
        casted : type of caller
        """
        return self._constructor(self._data, dtype=dtype)

    @property
    def _constructor(self):
        return NDFrame

    @property
    def axes(self):
        return self._data.axes

    def __repr__(self):
        return 'NDFrame'

    @property
    def values(self):
        return self._data.as_matrix()

    @property
    def ndim(self):
        return self._data.ndim

    def _set_axis(self, axis, labels):
        self._data.set_axis(axis, labels)
        self._clear_item_cache()

    def __getitem__(self, item):
        return self._get_item_cache(item)

    def _get_item_cache(self, item):
        cache = self._item_cache
        try:
            return cache[item]
        except Exception:
            values = self._data.get(item)
            res = self._box_item_values(item, values)
            cache[item] = res
            return res

    def _box_item_values(self, key, values):
        raise NotImplementedError

    def _clear_item_cache(self):
        self._item_cache.clear()

    def _set_item(self, key, value):
        if hasattr(self,'columns') and isinstance(self.columns, MultiIndex):
            # Pad the key with empty strings if lower levels of the key
            # aren't specified:
            if not isinstance(key, tuple):
                key = (key,)
            if len(key) != self.columns.nlevels:
                key += ('',)*(self.columns.nlevels - len(key))
        self._data.set(key, value)

        try:
            del self._item_cache[key]
        except KeyError:
            pass

    def __delitem__(self, key):
        """
        Delete item
        """
        deleted = False
        if (hasattr(self,'columns') and 
                isinstance(self.columns, MultiIndex)
                and key not in self.columns):
            # Allow shorthand to delete all columns whose first len(key)
            # elements match key:
            if not isinstance(key,tuple):
                key = (key,)
            for col in self.columns:
                if isinstance(col,tuple) and col[:len(key)] == key:
                    del self[col]
                    deleted = True
        if not deleted:
            # If the above loop ran and didn't delete anything because
            # there was no match, this call should raise the appropriate
            # exception:
            self._data.delete(key)

        try:
            del self._item_cache[key]
        except KeyError:
            pass

    def pop(self, item):
        """
        Return item and drop from frame. Raise KeyError if not found.
        """
        result = self[item]
        del self[item]
        return result

    def _expand_axes(self, key):
        new_axes = []
        for k, ax in zip(key, self.axes):
            if k not in ax:
                new_axes.append(np.concatenate([ax, [k]]))
            else:
                new_axes.append(ax)

        return new_axes

    #----------------------------------------------------------------------
    # Consolidation of internals

    def _consolidate_inplace(self):
        self._clear_item_cache()
        self._data = self._data.consolidate()

    def consolidate(self, inplace=False):
        """
        Compute NDFrame with "consolidated" internals (data of each dtype
        grouped together in a single ndarray). Mainly an internal API function,
        but available here to the savvy user

        Parameters
        ----------
        inplace : boolean, default False
            If False return new object, otherwise modify existing object

        Returns
        -------
        consolidated : type of caller
        """
        if inplace:
            self._consolidate_inplace()
            return self
        else:
            cons_data = self._data.consolidate()
            if cons_data is self._data:
                cons_data = cons_data.copy()
            return self._constructor(cons_data)

    @property
    def _is_mixed_type(self):
        self._consolidate_inplace()
        return len(self._data.blocks) > 1

    def _reindex_axis(self, new_index, fill_method, axis, copy):
        new_data = self._data.reindex_axis(new_index, axis=axis,
                                           method=fill_method, copy=copy)

        if new_data is self._data and not copy:
            return self
        else:
            return self._constructor(new_data)

    def cumsum(self, axis=None, skipna=True):
        """
        Return DataFrame of cumulative sums over requested axis.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise
        skipna : boolean, default True
            Exclude NA/null values. If an entire row/column is NA, the result
            will be NA

        Returns
        -------
        y : DataFrame
        """
        if axis is None:
            axis = self._default_stat_axis
        else:
            axis = self._get_axis_number(axis)

        y = self.values.copy()
        if not issubclass(y.dtype.type, np.integer):
            mask = np.isnan(self.values)

            if skipna:
                np.putmask(y, mask, 0.)

            result = y.cumsum(axis)

            if skipna:
                np.putmask(result, mask, np.nan)
        else:
            result = y.cumsum(axis)
        return self._wrap_array(result, self.axes, copy=False)

    def _wrap_array(self, array, axes, copy=False):
        raise NotImplementedError

    def cumprod(self, axis=None, skipna=True):
        """
        Return cumulative product over requested axis as DataFrame

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise
        skipna : boolean, default True
            Exclude NA/null values. If an entire row/column is NA, the result
            will be NA

        Returns
        -------
        y : DataFrame
        """
        if axis is None:
            axis = self._default_stat_axis
        else:
            axis = self._get_axis_number(axis)

        y = self.values.copy()
        if not issubclass(y.dtype.type, np.integer):
            mask = np.isnan(self.values)

            if skipna:
                np.putmask(y, mask, 1.)
            result = y.cumprod(axis)

            if skipna:
                np.putmask(result, mask, np.nan)
        else:
            result = y.cumprod(axis)
        return self._wrap_array(result, self.axes, copy=False)

    def cummax(self, axis=None, skipna=True):
        """
        Return DataFrame of cumulative max over requested axis.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise
        skipna : boolean, default True
            Exclude NA/null values. If an entire row/column is NA, the result
            will be NA

        Returns
        -------
        y : DataFrame
        """
        if axis is None:
            axis = self._default_stat_axis
        else:
            axis = self._get_axis_number(axis)

        y = self.values.copy()
        if not issubclass(y.dtype.type, np.integer):
            mask = np.isnan(self.values)

            if skipna:
                np.putmask(y, mask, -np.inf)

            result = np.maximum.accumulate(y, axis)

            if skipna:
                np.putmask(result, mask, np.nan)
        else:
            result = np.maximum.accumulate(y,axis)
        return self._wrap_array(result, self.axes, copy=False)

    def cummin(self, axis=None, skipna=True):
        """
        Return DataFrame of cumulative min over requested axis.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise
        skipna : boolean, default True
            Exclude NA/null values. If an entire row/column is NA, the result
            will be NA

        Returns
        -------
        y : DataFrame
        """
        if axis is None:
            axis = self._default_stat_axis
        else:
            axis = self._get_axis_number(axis)

        y = self.values.copy()
        if not issubclass(y.dtype.type, np.integer):
            mask = np.isnan(self.values)

            if skipna:
                np.putmask(y, mask, np.inf)

            result = np.minimum.accumulate(y, axis)

            if skipna:
                np.putmask(result, mask, np.nan)
        else:
            result = np.minimum.accumulate(y,axis)
        return self._wrap_array(result, self.axes, copy=False)

    def copy(self, deep=True):
        """
        Make a copy of this object

        Parameters
        ----------
        deep : boolean, default True
            Make a deep copy, i.e. also copy data

        Returns
        -------
        copy : type of caller
        """
        data = self._data
        if deep:
            data = data.copy()
        return self._constructor(data)

    def swaplevel(self, i, j, axis=0):
        """
        Swap levels i and j in a MultiIndex on a particular axis

        Returns
        -------
        swapped : type of caller (new object)
        """
        axis = self._get_axis_number(axis)
        result = self.copy()
        labels = result._data.axes[axis]
        result._data.set_axis(axis, labels.swaplevel(i, j))
        return result

    def add_prefix(self, prefix):
        """
        Concatenate prefix string with panel items names.

        Parameters
        ----------
        prefix : string

        Returns
        -------
        with_prefix : type of caller
        """
        new_data = self._data.add_prefix(prefix)
        return self._constructor(new_data)

    def add_suffix(self, suffix):
        """
        Concatenate suffix string with panel items names

        Parameters
        ----------
        suffix : string

        Returns
        -------
        with_suffix : type of caller
        """
        new_data = self._data.add_suffix(suffix)
        return self._constructor(new_data)

    def rename_axis(self, mapper, axis=0, copy=True):
        """
        Alter index and / or columns using input function or functions.
        Function / dict values must be unique (1-to-1). Labels not contained in
        a dict / Series will be left as-is.

        Parameters
        ----------
        mapper : dict-like or function, optional
        axis : int, default 0
        copy : boolean, default True
            Also copy underlying data

        See also
        --------
        DataFrame.rename

        Returns
        -------
        renamed : type of caller
        """
        # should move this at some point
        from pandas.core.series import _get_rename_function

        mapper_f = _get_rename_function(mapper)

        if axis == 0:
            new_data = self._data.rename_items(mapper_f, copydata=copy)
        else:
            new_data = self._data.rename_axis(mapper_f, axis=axis)
            if copy:
                new_data = new_data.copy()

        return self._constructor(new_data)

    def take(self, indices, axis=0):
        """
        Analogous to ndarray.take

        Parameters
        ----------
        indices : list / array of ints
        axis : int, default 0

        Returns
        -------
        taken : type of caller
        """
        if axis == 0:
            labels = self._get_axis(axis)
            new_items = labels.take(indices)
            new_data = self._data.reindex_axis(new_items, axis=0)
        else:
            new_data = self._data.take(indices, axis=axis)
        return self._constructor(new_data)

# Good for either Series or DataFrame

def truncate(self, before=None, after=None, copy=True):
    """Function truncate a sorted DataFrame / Series before and/or after
    some particular dates.

    Parameters
    ----------
    before : date
        Truncate before date
    after : date
        Truncate after date

    Returns
    -------
    truncated : type of caller
    """
    before = datetools.to_datetime(before)
    after = datetools.to_datetime(after)

    if before is not None and after is not None:
        assert(before <= after)

    result = self.ix[before:after]

    if isinstance(self.index, MultiIndex):
        result.index = self.index.truncate(before, after)

    if copy:
        result = result.copy()

    return result

