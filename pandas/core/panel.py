"""
Contains data structures designed for manipulating panel (3-dimensional) data
"""
# pylint: disable-msg=E1103
# pylint: disable-msg=W0231
# pylint: disable-msg=W0212
# pylint: disable-msg=W0621

from cStringIO import StringIO
from functools import partial
import operator

import numpy as np
from numpy.lib.format import write_array, read_array

from pandas.core.groupby import GroupBy
from pandas.core.index import Index
from pandas.core.frame import DataFrame
from pandas.core.matrix import DataMatrix
from pandas.core.mixins import Picklable

class PanelError(Exception):
    pass


def _pickle(arr):
    "Render text representation of array"

    io = StringIO()
    write_array(io, arr)

    return io.getvalue()


def _interpret(s):
    "Read text representation of ndarray"
    arr = read_array(StringIO(s))
    return arr

class Panel(Picklable):
    """
    Abstract superclass for LongPanel and WidePanel data structures
    """
    _items = None
    _major_axis = None
    _minor_axis = None
    _values = None

    def __repr__(self):
        class_name = str(self.__class__)

        I, N, K = len(self.items), len(self.major_axis), len(self.minor_axis)

        dims = 'Dimensions: %d (items) x %d (major) x %d (minor)' % (I, N, K)

        major = 'Major axis: %s to %s' % (self.major_axis[0],
                                          self.major_axis[-1])

        minor = 'Minor axis: %s to %s' % (self.minor_axis[0],
                                          self.minor_axis[-1])

        items = 'Items: %s to %s' % (self.items[0], self.items[-1])

        output = '%s\n%s\n%s\n%s\n%s' % (class_name, dims, items, major, minor)

        if self.factors:
            output += '\nFactors: %s' % ', '.join(self.factors)

        return output

    def _get_items(self):
        return self._items

    def _set_items(self, items):
        if not isinstance(items, Index):
            items = Index(items)

        self._items = items

    items = property(fget=_get_items, fset=_set_items)

    def _get_major_axis(self):
        return self._major_axis

    def _set_major_axis(self, major_axis):
        if not isinstance(major_axis, Index):
            major_axis = Index(major_axis)

        self._major_axis = major_axis

    major_axis = property(fget=_get_major_axis, fset=_set_major_axis)

    def _get_minor_axis(self):
        return self._minor_axis

    def _set_minor_axis(self, minor_axis):
        if not isinstance(minor_axis, Index):
            minor_axis = Index(minor_axis)

        self._minor_axis = minor_axis

    minor_axis = property(fget=_get_minor_axis, fset=_set_minor_axis)

    def _get_values(self):
        return self._values

    def _set_values(self, values):
        if not values.flags.contiguous:
            values = values.copy()

        self._values = values

    values = property(fget=_get_values, fset=_set_values)

    @property
    def dims(self):
        return len(self.items), len(self.major_axis), len(self.minor_axis)

_WIDE_AXIS_NUMBERS = {
    'items' : 0,
    'major' : 1,
    'minor' : 2
}
_WIDE_AXIS_NAMES = dict((v, k) for k, v in _WIDE_AXIS_NUMBERS.iteritems())


class WidePanel(Panel):
    """
    Represents wide format panel data, stored as 3-dimensional array

    Parameters
    ----------
    values: ndarray (items x major x minor)
    items: sequence
    major_axis: sequence
    minor_axis: sequence
    """
    def __init__(self, values, items, major_axis, minor_axis):
        self.items = items
        self.major_axis = major_axis
        self.minor_axis = minor_axis

#        self.factors = factors or {}
        self.factors = {}
        self.values = values

    @classmethod
    def _wide_axis_number(cls, axis):
        if axis in (0, 1, 2):
            return axis
        else:
            return _WIDE_AXIS_NUMBERS[axis]

    @classmethod
    def _wide_axis_name(cls, axis):
        if axis in _WIDE_AXIS_NUMBERS:
            return axis
        else:
            return _WIDE_AXIS_NAMES[axis]

    def _get_axis(self, axis):
        results = {
            0 : self.items,
            1 : self.major_axis,
            2 : self.minor_axis
        }

        return results[self._wide_axis_number(axis)]

    def _get_plane_axes(self, axis):
        """

        """
        axis = self._wide_axis_name(axis)

        if axis == 'major':
            index = self.minor_axis
            columns = self.items
        if axis == 'minor':
            index = self.major_axis
            columns = self.items
        elif axis == 'items':
            index = self.major_axis
            columns = self.minor_axis

        return index, columns

    @classmethod
    def fromDict(cls, data, intersect=True):
        """
        Construct WidePanel from dict of DataFrame objects

        Parameters
        ----------
        data: dict
            {field : DataFrame}
        intersect: boolean

        Returns
        -------
        WidePanel
        """
        data, index, columns = _homogenize(data, intersect=intersect)
        items = Index(sorted(data.keys()))

        values = np.array([data[k].values for k in items], dtype=float)

        return cls(values, items, index, columns)

    def keys(self):
        return list(self.items)

    def iteritems(self):
        for item in self.items:
            yield item, self[item]

    def _get_values(self):
        return self._values

    def _set_values(self, values):
        if not values.flags.contiguous:
            values = values.copy()

        if self.dims != values.shape:
            raise PanelError('Values shape %s did not match axes / items %s' %
                             (values.shape, self.dims))

        self._values = values

    values = property(fget=_get_values, fset=_set_values)

    def __getitem__(self, key):
        try:
            loc = self.items.indexMap[key]
        except KeyError:
            raise KeyError('%s not contained in panel data items!' % key)

        mat = self.values[loc]

        return DataMatrix(mat, index=self.major_axis, columns=self.minor_axis)

    def __setitem__(self, key, value):
        """
        Insert item at end of items for now
        """
        _, N, K = self.dims

        # XXX
        if isinstance(value, LongPanel):
            if len(value.items) != 1:
                raise Exception('Input panel must have only one item!')

            value = value.toWide()[value.items[0]]

        if isinstance(value, DataFrame):
            value = value.reindex(self.major_axis)
            value = value._withColumns(self.minor_axis)

            mat = value.values.reshape((1, N, K))

        elif np.isscalar(value):
            mat = np.empty((1, N, K), dtype=float)
            mat.fill(value)

        self.items = Index(list(self.items) + [key])
        self.values = np.row_stack((self.values, mat))

    def __getstate__(self):
        "Returned pickled representation of the panel"

        return (_pickle(self.values),
                _pickle(self.items),
                _pickle(self.major_axis),
                _pickle(self.minor_axis))

    def __setstate__(self, state):
        "Unpickle the panel"
        vals, items, major, minor = state

        self.items = _interpret(items)
        self.major_axis = _interpret(major)
        self.minor_axis = _interpret(minor)
        self.values = _interpret(vals)

    def conform(self, frame, axis='items'):
        """
        Conform input DataFrame to align with chosen axis pair.

        Parameters
        ----------
        frame: DataFrame
        axis: {'items', 'major', 'minor'}
            Axis the input corresponds to. E.g., if axis='major', then
            the frame's columns would be items, and the index would be
            values of the minor axis

        Returns
        -------
        DataFrame (or DataMatrix)
        """
        index, columns = self._get_plane_axes(axis)

        return frame.reindex(index)._withColumns(columns)

    def reindex(self, new_index, axis='major', fill_method=None):
        """
        Conform

        Parameters
        ----------
        new_index: Index or sequence
        axis: {'items', 'major', 'minor'}
            Axis to reindex
        fill_method: {'backfill', 'pad', 'interpolate', None}
            Method to use for filling holes in reindexed panel

        Returns
        -------
        WidePanel (new object)
        """
        import pandas.lib.tseries as tseries

        axis_i = self._wide_axis_number(axis)
        current_axis = self._get_axis(axis)

        if new_index is current_axis:
            return self.copy()

        if not isinstance(new_index, Index):
            new_index = Index(new_index)

        if not fill_method:
            fill_method = ''

        fill_method = fill_method.upper()

        if fill_method not in ['BACKFILL', 'PAD', '']:
            raise Exception("Don't recognize fill_method: %s" % fill_method)

        indexer, mask = tseries.getFillVec(current_axis, new_index,
                                           current_axis.indexMap,
                                           new_index.indexMap, fill_method)

        new_values = self.values.take(indexer, axis=axis_i)

        new_items = self.items
        new_major = self.major_axis
        new_minor = self.minor_axis
#        new_factors = dict((k, v.take(indexer))
#                           for k, v in self.factors.iteritems())

        if axis_i == 0:
            new_values[-mask] = np.NaN
            new_items = new_index
        elif axis_i == 1:
            new_values[:, -mask, :] = np.NaN
            new_major = new_index
        else:
            new_values[:, :, -mask] = np.NaN
            new_minor = new_index

        return WidePanel(new_values, new_items, new_major, new_minor)

    def _combine(self, other, func, axis=0):
        if isinstance(other, DataFrame):
            return self._combineFrame(other, func, axis=axis)
        elif isinstance(other, Panel):
            return self._combinePanel(other, func)
        elif np.isscalar(other):
            pass

    def _combineFrame(self, other, func, axis=0):
        index, columns = self._get_plane_axes(axis)
        axis = self._wide_axis_number(axis)

        other = other.reindex(index)._withColumns(columns)

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
        pass

    def add(self, other, axis='major'):
        """

        """
        return self._combine(other, operator.add, axis=axis)

    def subtract(self, other, axis='major'):
        """

        """
        return self._combine(other, operator.sub, axis=axis)

    def multiply(self, other, axis='major'):
        """

        """
        return self._combine(other, operator.mul, axis=axis)

    def divide(self, other, axis='major'):
        """

        """
        return self._combine(other, operator.div, axis=axis)

    def getMajorXS(self, key):
        """
        Parameters
        ----------

        Returns
        -------
        DataMatrix: index -> minor axis, columns -> items
        """
        try:
            loc = self.major_axis.indexMap[key]
        except KeyError:
            raise KeyError('%s not contained in major axis!' % key)

        mat = np.array(self.values[:, loc, :].T)
        return DataMatrix(mat, index=self.minor_axis, columns=self.items)

    def getMinorXS(self, key):
        """
        Parameters
        ----------

        Returns
        -------
        DataMatrix: index -> major axis, columns -> items
        """
        try:
            loc = self.minor_axis.indexMap[key]
        except KeyError:
            raise KeyError('%s not contained in minor axis!' % key)

        mat = np.array(self.values[:, :, loc].T)
        return DataMatrix(mat, index=self.major_axis, columns=self.items)

    def groupby(self, function, axis='major'):
        """
        Parameters
        ----------
        function: callable
            Mapping function for chosen access
        axis: {'major', 'minor', 'items'}, default 'major'

        Returns
        -------
        WidePanelGroupBy
        """
        return WidePanelGroupBy(self, function, axis=axis)

    def swapaxes(self):
        """
        Switch minor and major axes (and transpose values to reflect
        the change)

        Returns
        -------
        WidePanel (new object)
        """
        new_values = self.values.swapaxes(1, 2)

        return WidePanel(new_values, self.items,
                         self.minor_axis, self.major_axis)

    def toLong(self, filter_observations=True):
        """
        Transform wide format into long (stacked) format

        Parameters
        ----------
        filter_observations: boolean, default True
            Drop (major, minor) pairs without a complete set of observations
            across all the items

        Returns
        -------
        LongPanel
        """
        I, N, K = self.dims

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

    def filterItems(self, items):
        """
        Restrict items in panel to input list

        Parameters
        ----------
        items: sequence

        Returns
        -------
        WidePanel
        """
        intersection = self.items.intersection(items)
        indexer = [self.items.indexMap[col] for col in intersection]

        new_values = self.values.take(indexer, axis=0)
        return WidePanel(new_values, intersection, self.major_axis,
                         self.minor_axis)

    def _apply(self, func, axis='major', fill_na=True):
        """
        Parameters
        ----------
        func: numpy function
            Signature should match numpy.{sum, mean, var, std} etc.
        axis: {'major', 'minor', 'items'}
        fill_na: boolean, default True
            Replace NaN values with 0 first

        Returns
        -------
        DataMatrix
        """

        i = self._wide_axis_number(axis)
        index, columns = self._get_plane_axes(axis)

        values = self.values
        if fill_na:
            values = values.copy()
            values[-np.isfinite(values)] = 0

        result = func(values, axis=i)

        if axis != 'items':
            result = result.T

        if not result.ndim == 2:
            raise Exception('function %s incompatible' % func)

        return DataMatrix(result, index=index, columns=columns)

    def sum(self, axis='major'):
        return self._apply(np.sum, axis=axis)

    def mean(self, axis='major'):
        return self._apply(np.mean, axis=axis)

    def var(self, axis='major'):
        def _var(arr, axis=0):
            return np.std(arr, axis=axis, ddof=1)

        return self._apply(_var, axis=axis)

    def std(self, axis='major'):
        def _std(arr, axis=0):
            return np.std(arr, axis=axis, ddof=1)

        return self._apply(_std, axis=axis)

class LongPanelIndex(object):
    """
    Parameters
    ----------

    """
    def __init__(self, major_axis, minor_axis, major_labels,
                 minor_labels, mask=None):

        self.major_axis = major_axis
        self.minor_axis = minor_axis

        self.major_labels = major_labels
        self.minor_labels = minor_labels

        self._mask = mask

    def __getstate__(self):
        return (_pickle(self.major_axis),
                _pickle(self.minor_axis),
                _pickle(self.major_labels),
                _pickle(self.minor_labels))

    def __setstate__(self, state):
        major, minor, major_labels, minor_labels = state

        self.major_axis = _interpret(major)
        self.minor_axis = _interpret(minor)

        self.major_axis = _interpret(major_labels)
        self.minor_axis = _interpret(minor_labels)

    def isConsistent(self):
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
        Slice index between two major axis values, return complete LongPanel

        Parameters
        ----------
        before: type of major_axis values or None, default None
            None defaults to start of panel

        after: type of major_axis values or None, default None
            None defaults to after of panel

        Returns
        -------
        LongPanel
        """
        i, j = self._getAxisBounds(before, after)
        left, right = self._getLabelBounds(i, j)

        return LongPanelIndex(self.major_axis[i : j],
                              self.minor_axis,
                              self.major_labels[left : right] - i,
                              self.minor_labels[left : right])

    def getMajorBounds(self, begin=None, end=None):
        """

        """
        i, j = self._getAxisBounds(begin, end)
        left, right = self._getLabelBounds(i, j)

        return left, right

    def _getAxisBounds(self, begin, end):
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
            raise Exception('Must have begin <= end!')

        return i, j

    def _getLabelBounds(self, i, j):
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
        """

        """
        if self._mask is None:
            self._mask = self._makeMask()

        return self._mask

    def _makeMask(self):
        """
        Create observation selection vector using major and minor
        labels, for converting to wide format.
        """
        N, K = self.dims
        selector = self.minor_labels + K * self.major_labels

        mask = np.zeros(N * K, dtype=bool)
        mask[selector] = True

        return mask

    @property
    def dims(self):
        return len(self.major_axis), len(self.minor_axis)


class LongPanel(Panel):
    """
    Represents long or "stacked" format panel data
    """

    def __init__(self, values, items, index, factors=None):
        self.items = items
        self.index = index

        self.values = values

        self.factors = factors or {}

    @classmethod
    def fromRecords(cls, data, major_field, minor_field,
                    factors=None, exclude=None):
        """
        Convert LongPanel

        Parameters
        ----------
        data: DataFrame, structured or record array, or dict
        major_field: string
        minor_field: string
            Name of field
        factors: list-like, default None
        exclude: list-like, default None

        Returns
        -------
        LongPanel
        """
        from pandas.lib.tseries import getMergeVec

        if isinstance(data, np.ndarray):
            # Dtype when you have data
            if data.dtype.type != np.void:
                raise Exception('Input was not a structured array!')

            columns = data.dtype.names
            data = dict((k, data[k]) for k in columns)
        elif isinstance(data, DataFrame):
            data = data._series.copy()

        exclude = set(exclude) if exclude is not None else set()

        if major_field in data:
            major_vec = data.pop(major_field)
        else:
            raise Exception('No field named %s' % major_field)

        if minor_field in data:
            minor_vec = data.pop(minor_field)
        else:
            raise Exception('No field named %s' % minor_field)

        major_axis = Index(sorted(set(major_vec)))
        minor_axis = Index(sorted(set(minor_vec)))

        major_labels, _ = getMergeVec(major_vec, major_axis.indexMap)
        minor_labels, _ = getMergeVec(minor_vec, minor_axis.indexMap)

        factor_dict = {}
        for col in data.keys():
            series = data[col]

            # Is it a factor?
            if not np.issctype(series.dtype):
                factor_dict[col] = fac = Factor.fromarray(series)
                del data[col]

        items = sorted(data)
        values = np.array([data[k] for k in items]).T

        index = LongPanelIndex(major_axis, minor_axis,
                               major_labels, minor_labels)

        return LongPanel(values, items, index, factors=factor_dict)

    @property
    def columns(self):
        """
        So LongPanel can be DataMatrix-like at times
        """
        return self.items

    def cols(self):
        "DataMatrix compatibility"
        return self.columns

    def copy(self):
        values = self.values.copy()
        items = self.items
        index = self.index
        return LongPanel(values, items, index, factors=self.factors)

    def _get_major_axis(self):
        return self.index.major_axis

    major_axis = property(fget=_get_major_axis)

    def _get_minor_axis(self):
        return self.index.minor_axis

    minor_axis = property(fget=_get_minor_axis)

    def _get_values(self):
        return self._values

    def _set_values(self, values):
        if not values.flags.contiguous:
            values = values.copy()

        shape = len(self.index.major_labels), len(self.items)

        if values.shape != shape:
            raise Exception('Values shape %s mismatch to %s' % (values.shape,
                                                                shape))

        self._values = values

    values = property(fget=_get_values, fset=_set_values)

    def __getitem__(self, key):
        "Return column of panel as LongPanel"

        loc = self.items.indexMap[key]

        return LongPanel(self.values[:, loc : loc + 1].copy(),
                        [key], self.index, factors=self.factors)

    def __setitem__(self, key, value):
        """
        Insert item at end of items for now
        """
        if np.isscalar(value):
            mat = np.empty((len(self.values), 1), dtype=float)
            mat.fill(value)

        self.items = Index(list(self.items) + [key])
        self.values = np.column_stack((self.values, mat))

    def __getstate__(self):
        "Returned pickled representation of the panel"

        return (_pickle(self.values),
                _pickle(self.items),
                self.index)

    def __setstate__(self, state):
        "Unpickle the panel"
        (vals, items, index) = state

        self.items = _interpret(items)
        self.index = index
        self.values = _interpret(vals)

    def _combine(self, other, func, axis=0):
        if isinstance(other, DataFrame):
            return self._combineFrame(other, func, axis=axis)
        elif isinstance(other, Panel):
            return self._combinePanel(other, func)
        elif np.isscalar(other):
            pass

    def _combineFrame(self, other, axis=0):
        pass

    def _combinePanel(self, other, func):
        """
        Arithmetic operation between panels
        """
        if self.index is not other.index:
            raise Exception("Can only combine identically-indexed "
                            "panels for now")

        if len(other.items) == 1:
            new_values = func(self.values, other.values)
        else:
            new_values = func(self.values, other.values)

        return LongPanel(new_values, self.items, self.index,
                         factors=self.factors)

    def add(self, other, axis='major'):
        """

        """
        return self._combine(other, operator.add, axis=axis)

    def subtract(self, other, axis='major'):
        """

        """
        return self._combine(other, operator.sub, axis=axis)

    def multiply(self, other, axis='major'):
        """

        """
        return self._combine(other, operator.mul, axis=axis)

    def divide(self, other, axis='major'):
        """

        """
        return self._combine(other, operator.div, axis=axis)

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

    def toWide(self):
        """
        Transform long (stacked) format into wide format

        Returns
        -------
        WidePanel
        """
        if not self.index.isConsistent():
            raise PanelError('Panel has duplicate (major, minor) pairs, '
                             'cannot be reliably converted to wide format.')

        I, N, K = self.dims

        values = np.empty((I, N, K), dtype=float)

        mask = self.index.mask
        notmask = -mask

        for i in xrange(len(self.items)):
            values[i].flat[mask] = self.values[:, i]
            values[i].flat[notmask] = np.NaN

        return WidePanel(values, self.items, self.major_axis, self.minor_axis)

    def toCSV(self, path):
        def format_cols(items):
            cols = ['Major', 'Minor'] + list(items)
            return '"%s"' % '","'.join(cols)

        def format_row(major, minor, values):
            vals = ','.join('%.12f' % val for val in values)
            return '%s,%s,%s' % (major, minor, vals)

        output = self._textConvert(format_cols, format_row)

        f = open(path, 'w')
        f.write(output)
        f.close()

    def toString(self, col_space=15, return_=False):
        """
        Output a screen-friendly version of this Panel
        """
        from pandas.core.frame import _pfixed

        major_space = max(max([len(str(idx))
                               for idx in self.major_axis]) + 4, 9)
        minor_space = max(max([len(str(idx))
                               for idx in self.minor_axis]) + 4, 9)

        def format_cols(items):
            return '%s%s%s' % (_pfixed('Major', major_space),
                               _pfixed('Minor', minor_space),
                               ''.join(_pfixed(h, col_space) for h in items))

        def format_row(major, minor, values):
            return '%s%s%s' % (_pfixed(major, major_space),
                               _pfixed(minor, minor_space),
                               ''.join(_pfixed(v, col_space) for v in values))

        output = self._textConvert(format_cols, format_row)

        if return_:
            return output
        else:
            print output

    def _textConvert(self, format_cols, format_row):
        output = StringIO()
        print >> output, format_cols(self.items)

        label_pairs = zip(self.index.major_labels,
                          self.index.minor_labels)
        major, minor = self.major_axis, self.minor_axis
        for i, (major_i, minor_i) in enumerate(label_pairs):
            row = format_row(major[major_i], minor[minor_i], self.values[i])
            print >> output, row

        return output.getvalue()

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
        # using the invariant that the major labels are ordered.
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
        before: type of major_axis values or None, default None
            None defaults to start of panel

        after: type of major_axis values or None, default None
            None defaults to end of panel

        Returns
        -------
        LongPanel
        """
        left, right = self.index.getMajorBounds(before, after)
        new_index = self.index.truncate(before, after)

        return LongPanel(self.values[left : right],
                         self.items, new_index)

    def filterItems(self, items):
        """
        Restrict items in panel to input list

        Parameters
        ----------
        items: sequence

        Returns
        -------
        WidePanel
        """
        intersection = self.items.intersection(items)
        indexer = [self.items.indexMap[col] for col in intersection]

        new_values = self.values.take(indexer, axis=1)
        return LongPanel(new_values, intersection, self.index)

    def getAxisDummies(self, axis='minor'):
        """
        Construct 1-0 dummy variables corresponding to designated axis
        labels

        Parameters
        ----------
        axis: {'major', 'minor'}, default 'minor'

        Returns
        -------
        LongPanel, item names taken from chosen axis
        """
        if axis == 'minor':
            dim = len(self.minor_axis)
            items = self.minor_axis
        elif axis == 'major':
            dim = len(self.major_axis)
            items = self.major_axis
        else:
            raise Exception('Do not recognize axis %s' % axis)

        vals = np.eye(dim, dtype=float)

        return self._makeDummyPanel(vals, items, axis=axis)

    def getFrameDummies(self, dataFrame, axis='minor', prefix=None):
        if axis == 'minor':
            dataFrame = dataFrame.reindex(self.minor_axis)
        elif axis == 'major':
            dataFrame = dataFrame.reindex(self.major_axis)

        items = dataFrame.columns

        return self._makeDummyPanel(dataFrame.values, items, axis=axis,
                                    prefix=prefix)

    def getItemDummies(self, item):
        """
        Use unique values in column of panel to construct LongPanel
        containing dummy

        Parameters
        ----------
        item: object
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

    def _makeDummyPanel(self, values, items, axis='minor'):
        """
        Construct 1-0 dummy variables corresponding to designated axis
        labels

        Parameters
        ----------
        axis: {'major', 'minor'}, default 'minor'

        Returns
        -------
        LongPanel, item names taken from chosen axis
        """

        N, K = values.shape

        if len(items) != K:
            raise Exception('items length does not match values matrix')

        if axis == 'minor':
            if len(self.minor_axis) != N:
                raise Exception('Axis length does not match values matrix')
            dummy_mat = values.take(self.index.minor_labels, axis=0)

        elif axis == 'major':
            if len(self.major_axis) != N:
                raise Exception('Axis length does not match values matrix')
            dummy_mat = values.take(self.index.major_labels, axis=0)
        else:
            raise Exception('Do not recognize axis %s' % axis)

        return LongPanel(dummy_mat, items, self.index)

    def applyToAxis(self, f, axis='major', broadcast=False):
        """
        Aggregate over a particular axis

        Parameters
        ----------
        f: function
            NumPy function to apply to each group
        axis: {'major', 'minor'}

        broadcast: boolean

        Returns
        -------
        broadcast=True  -> LongPanel
        broadcast=False -> DataMatrix
        """
        if axis == 'minor':
            panel = self.swapaxes()
            result = panel.applyToAxis(f, axis='major', broadcast=broadcast)
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
            panel = DataMatrix(result, index=self.major_axis,
                               columns=self.items)

        return panel

    def mean(self, axis='major', broadcast=False):
        return self.applyToAxis(partial(np.mean, axis=0), axis, broadcast)

    def sum(self, axis='major', broadcast=False):
        return self.applyToAxis(partial(np.sum, axis=0), axis, broadcast)

    def apply(self, f):
        return LongPanel(f(self.values), self.items, self.index)

    def square(self):
        return self.apply(np.square)

    def count(self, axis=0):
        if axis == 0:
            lp = self
        else:
            lp = self.swapaxes()

        N = len(lp.values)
        bounds = lp.index._bounds

        return np.concatenate((np.diff(bounds), [N - bounds[-1]]))

    def leftJoin(self, other):
        """

        Parameters
        ----------
        other: LongPanel
        """
        if other.index is self.index:
            pass

    def merge(self, other):
        """

        Parameters
        ----------
        other: LongPanel

        Returns
        -------
        LongPanel
        """
        assert(self.index is other.index)

        values = np.concatenate((self.values, other.values), axis=1).copy()
        items = self.items.tolist() + other.items.tolist()

        return LongPanel(values, items, self.index)

    def addPrefix(self, prefix):
        """
        Concatenate prefix string with panel items names.

        Parameters
        ----------
        prefix: string

        Returns
        -------
        LongPanel

        Note: does *not* copy values matrix
        """
        new_items = [_makeItemName(item, prefix) for item in self.items]

        return LongPanel(self.values, new_items, self.index)


class Factor(object):
    """
    Represents a categorical variable in classic R / S+ fashion
    """
    def __init__(self, labels, levels):
        self.labels = labels
        self.levels = levels

    @classmethod
    def fromarray(cls, values):
        levels = np.array(sorted(set(values)), dtype=object)
        labels = levels.searchsorted(values)

        return Factor(labels, levels=levels)

    def __repr__(self):
        temp = 'Factor:\n%s\nLevels (%d): %s'

        values = self.levels[self.labels]
        return temp % (repr(values), len(self.levels), self.levels)

    def __getitem__(self, key):
        if key is None and key not in self.index:
            raise Exception('None/Null object requested of Series!')

        if isinstance(key, int):
            i = self.labels[key]
            return self.levels[i]
        else:
            new_labels = self.labels[key]
            return Factor(new_labels, self.levels)


def _makeItemName(item, prefix=None):
    if prefix is None:
        return item

    template = '%g%s' if isinstance(item, float) else '%s%s'
    return template % (prefix, item)

def _makePrefixedLongPanel(values, items, index, prefix):
    items = [_makeItemName(item, prefix) for item in items]

    return LongPanel(values, items, index)

def _convert(data, order, factors=None):
    """

    Parameters
    ----------

    Returns
    -------

    """

def _homogenize(frames, intersect=True):
    """
    Conform set of DataFrame-like objects to either an intersection
    of indices / columns or a union.

    Parameters
    ----------
    frames: dict
    intersect: boolean, default True

    Returns
    -------
    dict of aligned frames, index, columns
    """
    result = {}

    index = None
    columns = None

    if intersect:
        for key, frame in frames.iteritems():
            if index is None:
                index = frame.index
            elif index is not frame.index:
                index = index.intersection(frame.index)

            if columns is None:
                columns = set(frame.cols())
            else:
                columns &= set(frame.cols())
    else:
        for key, frame in frames.iteritems():
            if index is None:
                index = frame.index
            elif index is not frame.index:
                index = index.union(frame.index)

            if columns is None:
                columns = set(frame.cols())
            else:
                columns |= set(frame.cols())

    columns = sorted(columns)

    if intersect:
        for key, frame in frames.iteritems():
            result[key] = frame.filterItems(columns).reindex(index)
    else:
        for key, frame in frames.iteritems():
            if not isinstance(frame, DataMatrix):
                frame = frame.toDataMatrix()

            result[key] = frame._withColumns(columns).reindex(index)

    return result, index, columns

def pivot(index, columns, values):
    """
    Produce 'pivot' table based on 3 columns of this DataFrame.
    Uses unique values from index / columns and fills with values.

    Parameters
    ----------
    index: ndarray
        Labels to use to make new frame's index
    columns: ndarray
        Labels to use to make new frame's columns
    values: ndarray
        Values to use for populating new frame's values

    Note
    ----
    Obviously, all 3 of the input arguments must have the same length

    Returns
    -------
    DataMatrix
    """
    import pandas.lib.tseries as tseries

    if not (len(index) == len(columns) == len(values)):
        raise Exception('Pivot inputs must all be same length!')

    major_axis = Index(sorted(set(index)))
    minor_axis = Index(sorted(set(columns)))

    major_labels, _ = tseries.getMergeVec(index, major_axis.indexMap)
    minor_labels, _ = tseries.getMergeVec(columns, minor_axis.indexMap)

    valueMat = values.view(np.ndarray).reshape(len(values), 1)

    longIndex = LongPanelIndex(major_axis, minor_axis,
                               major_labels, minor_labels)

    longPanel = LongPanel(valueMat, ['foo'], longIndex)
    longPanel = longPanel.sort()

    try:
        return longPanel.toWide()['foo']
    except PanelError:
        return _slow_pivot(index, columns, values)

def _slow_pivot(index, columns, values):
    """
    Produce 'pivot' table based on 3 columns of this DataFrame.
    Uses unique values from index / columns and fills with values.

    Parameters
    ----------
    index: string or object
        Column name to use to make new frame's index
    columns: string or object
        Column name to use to make new frame's columns
    values: string or object
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

    return DataFrame.fromDict(tree)

def test():
    return pivot(np.array([1, 2, 3, 4, 4]),
                 np.array(['a', 'a', 'a', 'a', 'a']),
                 np.array([1, 2, 3, 5, 4]))

def _monotonic(arr):
    return not (arr[1:] < arr[:-1]).any()

def group_agg(values, bounds, f):
    """
    R-style aggregator
    """
    N, K = values.shape
    result = np.empty((len(bounds), K), dtype=float)

    for i, left_bound in enumerate(bounds):
        if i == len(bounds) - 1:
            right_bound = N
        else:
            right_bound = bounds[i + 1]

        result[i] = f(values[left_bound : right_bound])

    return result

class WidePanelGroupBy(GroupBy):
    pass

class LongPanelGroupBy(GroupBy):
    pass

if __name__ == '__main__':
    from datetime import datetime
    import string

    import numpy as np

    from pandas.core.api import DataMatrix, DateRange

    N = 50
    K = 4

    start = datetime(2009, 9, 2)
    dateRange = DateRange(start, periods=N)

    cols = ['Col' + c for c in string.ascii_uppercase[:K]]

    def makeDataMatrix():
        data = DataMatrix(np.random.randn(N, K),
                          columns=cols,
                          index=dateRange)

        return data

    def makeDataMatrixForWeekday():
        values = [d.weekday() for d in dateRange]
        data = DataMatrix(dict((k, values) for k in cols),
                          index=dateRange)

        return data

    data = {
        'ItemA' : makeDataMatrix(),
        'ItemB' : makeDataMatrix(),
        'ItemC' : makeDataMatrix(),
#        'ItemD' : makeDataMatrixForWeekday(),
    }

    Y = makeDataMatrix()

    data['ItemA']['ColA'][:10] = np.NaN

    panel = WidePanel.fromDict(data)

    longPanel = panel.toLong(filter_observations=True)
    widePanel = longPanel.toWide()
