"""
Contains data structures designed for manipulating panel (3-dimensional) data
"""

from pandas.core.index import Index
from pandas.core.frame import _pfixed
from pandas.core.matrix import DataMatrix

class Panel(object):
    """
    Abstract superclass for LongPanel and WidePanel data structures
    """
    _items = None
    _major_axis = None
    _minor_axis = None
    _values = None

    def __init__(self, *args, **kwargs):
        raise NotImplementedError()
    
    def __repr__(self):
        class_name = str(self.__class__)

        I, N, K = len(self.items), len(self.major_axis), len(self.minor_axis)
        
        dims = 'Dimensions: %d (items) x %d (major) x %d (minor)' % (I, N, K)

        major = 'Major axis: %s to %s' % (self.major_axis[0],
                                          self.major_axis[-1])

        minor = 'Minor axis: %s to %s' % (self.minor_axis[0],
                                          self.minor_axis[-1])
        
        items = 'Items: %s to %s' % (self.items[0], self.items[-1]) 
        
        return ('%(class_name)s\n%(dims)s\n%(items)s\n'
                '%(major)s\n%(minor)s' % locals())

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
            raise Exception('Values must be C-contiguous!')

        self._values = values
        
    values = property(fget=_get_values, fset=_set_values)

    @property
    def dims(self):
        return len(self.items), len(self.major_axis), len(self.minor_axis)
    

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
        self.values = values
        self.major_axis = major_axis
        self.minor_axis = minor_axis
    
    def __getitem__(self, key):
        loc = self.items.indexMap[key]

        mat = self.values[loc]

        return DataMatrix(mat, index=self.major_axis, columns=self.minor_axis)
        
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
        
    def swapaxes(self):
        """
        Switch minor and major axes (and transpose values to reflect
        the change)
        
        Returns
        -------
        WidePanel (new object)
        """
        new_values = self.values.swapaxes(1, 2)
        
        return cls(new_values, self.items, self.minor_axis, self.major_axis)
        
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
            mask = np.isfinite(self.values).sum(axis=0) == I
            size = mask.sum()
            selector = mask.ravel()
        else:
            size = N * K
            selector = slice(None, None)
            
        values = np.empty((size, I), dtype=float)        
            
        for i, field in enumerate(self.items):
            values[:, i] = self.values[i].ravel()[selector]
            
        major_labels = np.arange(N).repeat(K)[selector]

        # Anyone think of a better way to do this? np.repeat does not
        # do what I want
        minor_labels = np.arange(K).reshape(1, K)[np.zeros(N, dtype=int)]
        minor_labels = minor_labels.ravel()[selector]

        if filter_observations:
            return LongPanel(values, self.items, self.major_axis,
                             self.minor_axis, major_labels, minor_labels,
                             mask=selector)
        else:
            return LongPanel(values, self.items, self.major_axis,
                             self.minor_axis, major_labels, minor_labels)

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

class LongPanel(Panel):
    """
    Represents long or "stacked" format panel data
    """

    def __init__(self, values, items, major_axis, minor_axis, major_labels,
                 minor_labels, mask=None):

        self.items = items
        self.values = values
        self.major_axis = major_axis
        self.minor_axis = minor_axis

        self.major_labels = major_labels
        self.minor_labels = minor_labels

        self.__mask = mask

    def toWide(self):
        """
        Transform long (stacked) format into wide format

        Returns
        -------
        WidePanel
        """
        I, N, K = self.dims

        values = np.empty((I, N, K), dtype=float)

        mask = self._mask
        notmask = -mask
        
        for i in xrange(len(self.items)):
            values[i].flat[mask] = self.values[:, i]
            values[i].flat[notmask] = np.NaN
            
        return WidePanel(values, self.items, self.major_axis, self.minor_axis)

    def toString(self, col_space=15, return_=False):
        """
        Output a screen-friendly version of this Panel
        """
        from cStringIO import StringIO

        output = StringIO()

        major_space = max([len(str(idx)) for idx in self.major_axis]) + 4
        minor_space = max([len(str(idx)) for idx in self.minor_axis]) + 4
        
        for h in ['Major', 'Minor'] + list(self.items):
             output.write(_pfixed(h, col_space))

        output.write('\n')

        label_pairs = zip(self.major_labels, self.minor_labels)
        for i, (major_i, minor_i) in enumerate(label_pairs):
            vals = ''.join(_pfixed(v, col_space) for v in self.values[i])

            row = '%s%s%s\n' % (_pfixed(self.major_axis[major_i], col_space),
                                _pfixed(self.minor_axis[minor_i], col_space),
                                vals)

            output.write(row)
        
        if return_:
            return output.getvalue()
        else:
            print output.getvalue()

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
        indexer = self.minor_labels.argsort(kind='mergesort')

        new_major = self.minor_labels[indexer]
        new_minor = self.major_labels[indexer]
        new_values = self.values[indexer]
        
        return LongPanel(new_values, self.items, self.minor_axis,
                         self.major_axis, new_major, new_minor,
                         mask=self._mask)

    def getSlice(self, begin=None, end=None):
        """
        Slice panel between two major axis values, return complete LongPanel
        
        Parameters
        ----------
        begin: type of major_axis values or None, default None
            None defaults to start of panel

        end: type of major_axis values or None, default None
            None defaults to end of panel

        Returns
        -------
        LongPanel
        """
        i, j = self._getAxisBounds(begin, end)
        left, right = self._getLabelBounds(i, j)
        
        return LongPanel(self.values[left : right],
                         self.items,
                         self.major_axis[i : j],
                         self.minor_axis,
                         self.major_labels[left : right] - i,
                         self.minor_labels[left : right])

    def getValueSlice(self, begin=None, end=None):
        """
        Slice panel between two major axis values and return only
        values array
        
        Parameters
        ----------
        begin: type of major_axis values or None, default None
            None defaults to start of panel

        end: type of major_axis values or None, default None
            None defaults to end of panel

        Returns
        -------
        ndarray
        """
        i, j = self._getAxisBounds(begin, end)
        left, right = self._getLabelBounds(i, j)

        return self.values[left : right]
    
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
        return LongPanel(new_values, intersection, self.major_axis,
                         self.minor_axis, self.major_labels,
                         self.minor_labels, mask=self._mask)

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
            right = len(self.values)
        else:
            right = self._bounds[j]
        
        return left, right

    __bounds = None
    @property
    def _bounds(self):
        if self.__bounds is None:
            inds = np.arange(len(self.major_axis))
            self.__bounds = self.major_labels.searchsorted(inds)
            
        return self.__bounds
    
    @property
    def _mask(self):
        """
        
        """
        if self.__mask is None:
            self.__mask = self._makeMask()

        return self.__mask
        
    def _makeMask(self):
        """
        Create observation selection vector using major and minor
        labels, for converting to wide format.
        """
        _, N, K = self.dims
        selector = self.minor_labels + K * self.major_labels
        
        mask = np.zeros(N * K, dtype=bool)
        mask[selector] = True
        
        return mask


def _homogenize(frames, intersect=True):
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


if __name__ == '__main__':
    from datetime import datetime
    import string

    import numpy as np

    from pandas.core.api import DataMatrix, DateRange
    from pandas.stats.linmodel import LinearModel, XSLinearModel

    N = 5000
    K = 4

    start = datetime(2009, 9, 2)
    dateRange = DateRange(start, periods=N)

    def makeDataMatrix():
        data = DataMatrix(np.random.randn(N, K),
                          columns=list(string.ascii_uppercase[:K]),
                          index=dateRange)

        return data

    data = {
        'A' : makeDataMatrix(),
        'B' : makeDataMatrix(),
        'C' : makeDataMatrix()
    }

    data['A']['A'][:10] = np.NaN
    
    panel = WidePanel.fromDict(data)

    long = panel.toLong(filter_observations=True)
    wide = long.toWide()
