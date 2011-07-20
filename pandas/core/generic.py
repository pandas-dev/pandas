import numpy as np
import cPickle

import pandas.core.datetools as datetools

#-------------------------------------------------------------------------------
# Picklable mixin

class Picklable(object):

    def save(self, fileName):
        f = open(fileName, 'wb')
        try:
            cPickle.dump(self, f, protocol=cPickle.HIGHEST_PROTOCOL)
        finally:
            f.close()

    @classmethod
    def load(cls, fileName):
        f = open(fileName, 'rb')
        try:
            return cPickle.load(f)
        finally:
            f.close()

class PandasError(Exception):
    pass

class AxisProperty(object):

    def __init__(self, axis=0):
        self.axis = axis

    def __get__(self, obj, type=None):
        data = getattr(obj, '_data')
        return data.axes[self.axis]

    def __set__(self, obj, value):
        data = getattr(obj, '_data')
        data.set_axis(self.axis, value)

class NDFrame(Picklable):

    _AXIS_NUMBERS = {
        'index' : 0,
        'columns' : 1
    }

    _AXIS_ALIASES = {}
    _AXIS_NAMES = dict((v, k) for k, v in _AXIS_NUMBERS.iteritems())

    _default_stat_axis = 0

    def __init__(self, data, axes=None, copy=False):
        self._data = data

    @property
    def _constructor(self):
        return NDFrame

    @property
    def axes(self):
        return self._data.axes

    def __repr__(self):
        # TODO
        return 'NDFrame'

    @property
    def values(self):
        return self._data.as_matrix()

    @property
    def ndim(self):
        return self._data.ndim

    #----------------------------------------------------------------------
    # Consolidation of internals

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
        return self._constructor(cons_data)

    @property
    def _is_mixed_type(self):
        self._consolidate_inplace()
        return len(self._data.blocks) > 1

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

    @property
    def axes(self):
        return self._data.axes

    def groupby(self, mapper):
        """
        Goup series using mapper (dict or key function, apply given
        function to group, return result as series).

        Parameters
        ----------
        mapper: function, dict or Series
            Called on each element of the object index to determine
            the groups.  If a dict or Series is passed, the Series or
            dict VALUES will be used to determine the groups

        Returns
        -------
        GroupBy object
        """
        from pandas.core.groupby import groupby
        return groupby(self, mapper)

    def _select_generic(self, crit, axis=0):
        """
        Return data corresponding to axis labels matching criteria

        Parameters
        ----------
        crit : function
            To be called on each index (label). Should return True or False
        axis : {0, 1}

        Returns
        -------
        selection : type of caller
        """
        axis_name = self._get_axis_name(axis)
        axis = self._get_axis(axis)
        new_axis = axis[np.asarray([crit(label) for label in axis])]
        return self.reindex(**{axis_name : new_axis})

    def _reindex_axis(self, new_index, fill_method, axis):
        if axis == 0:
            new_data = self._data.reindex_items(new_index)
        else:
            new_data = self._data.reindex_axis(new_index, axis=axis,
                                               method=fill_method)
        return self._constructor(new_data)

    def truncate(self, before=None, after=None):
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
        # returns view, want to copy
        return self.ix[before:after].copy()

    def cumsum(self, axis=None):
        """
        Return DataFrame of cumulative sums over requested axis.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise

        Returns
        -------
        y : DataFrame
        """
        if axis is None:
            axis = self._default_stat_axis
        else:
            axis = self._get_axis_number(axis)

        y = self.values.copy()
        if not issubclass(y.dtype.type, np.int_):
            mask = np.isnan(self.values)
            np.putmask(y, mask, 0.)
            result = y.cumsum(axis)
            np.putmask(result, mask, np.nan)
        else:
            result = y.cumsum(axis)
        return self._wrap_array(result, self.axes, copy=False)

    def _wrap_array(self, array, axes, copy=False):
        raise NotImplementedError

    def cumprod(self, axis=None):
        """
        Return cumulative product over requested axis as DataFrame

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise

        Returns
        -------
        y : DataFrame
        """
        if axis is None:
            axis = self._default_stat_axis
        else:
            axis = self._get_axis_number(axis)

        y = self.values.copy()
        if not issubclass(y.dtype.type, np.int_):
            mask = np.isnan(self.values)
            np.putmask(y, mask, 1.)
            result = y.cumprod(axis)
            np.putmask(result, mask, np.nan)
        else:
            result = y.cumprod(axis)
        return self._wrap_array(result, self.axes, copy=False)

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
