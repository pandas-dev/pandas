"Contains data structures designed for manipulating panel (3-dimensional) data"

from pandas.core.index import Index
from pandas.core.frame import _pfixed
from pandas.core.matrix import DataMatrix

class Panel(object):
    """
    Abstract superclass for LongPanel and WidePanel data structures
    """
    _fields = None
    _major_axis = None
    _minor_axis = None

    def __init__(self, *args, **kwargs):
        raise NotImplementedError()
    
    def __repr__(self):
        class_name = str(self.__class__)

        F, N, K = len(self.fields), len(self.major_axis), len(self.minor_axis)
        
        dims = 'Dimensions: %d (fields) x %d (major) x %d (minor)' % (F, N, K)

        major = 'Major axis: %s to %s' % (self.major_axis[0],
                                          self.major_axis[-1])

        minor = 'Minor axis: %s to %s' % (self.minor_axis[0],
                                          self.minor_axis[-1])
        
        fields = 'Fields: %s to %s' % (self.fields[0], self.fields[-1]) 
        
        return ('%(class_name)s\n%(dims)s\n%(fields)s\n'
                '%(major)s\n%(minor)s' % locals())

    def _get_fields(self):
        return self._fields

    def _set_fields(self, fields):
        if not isinstance(fields, Index):
            fields = Index(fields)

        self._fields = fields
        
    fields = property(fget=_get_fields, fset=_set_fields)

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

    @property
    def dims(self):
        return len(self.fields), len(self.major_axis), len(self.minor_axis)
    

class WidePanel(Panel):
    """

    Parameters
    ----------
    
    """
    def __init__(self, values, fields, major_axis, minor_axis):
        self.fields = fields
        self.values = values
        self.major_axis = major_axis
        self.minor_axis = minor_axis
    
    def __getitem__(self, key):
        loc = self.fields.indexMap[key]

        mat = self.values[loc]

        return DataMatrix(mat, index=self.major_axis, columns=self.minor_axis)
        
    @classmethod
    def fromDict(cls, data):
        data, index, columns = _homogenize(data)
        fields = Index(sorted(data.keys()))

        values = np.array([data[k].values for k in fields], dtype=float)
        
        return cls(values, fields, index, columns)
        
    def swapAxes(self):
        new_values = self.values.swapaxes(1, 2)
        
        return cls(new_values, self.fields, self.minor_axis, self.major_axis)
        
    def toLong(self):
        """
        Transform wide format into long (stacked) format

        Returns
        -------
        LongPanel
        """
        F, N, K = self.dims
        
        mask = np.isfinite(self.values).sum(axis=0) == F

        size = mask.sum()
        values = np.empty((size, F), dtype=float)
        
        selector = mask.ravel()
        for i, field in enumerate(self.fields):
            values[:, i] = self.values[i].ravel()[selector]
            
        major_labels = np.arange(N).repeat(K)[selector]

        # Anyone think of a better way to do this? np.repeat does not
        # do what I want
        minor_labels = np.arange(K).reshape(1, K)[np.zeros(N, dtype=int)]
        minor_labels = minor_labels.ravel()[selector]
        
        return LongPanel(values, self.fields, self.major_axis, self.minor_axis,
                         major_labels, minor_labels, mask=selector)


class LongPanel(Panel):
    """

    """
    def __init__(self, values, fields, major_axis, minor_axis, major_labels,
                 minor_labels, mask=None):

        self.fields = fields
        self.values = values
        self.major_axis = major_axis
        self.minor_axis = minor_axis

        self.major_labels = major_labels
        self.minor_labels = minor_labels

        self._mask = mask
        
    def _makeMask(self):
        """
        Seperate method for testing purposes
        """
        _, N, K = self.dims
        selector = self.minor_labels + K * self.major_labels
        
        mask = np.zeros(N * K, dtype=bool)
        mask[selector] = True
        
        return mask

    @property
    def mask(self):
        if self._mask is None:
            self._mask = self._maskMask()

        return self._makeMask()

    def toString(self, col_space=15, return_=False):
        """
        Output a screen-friendly version of this Panel
        """
        from cStringIO import StringIO

        output = StringIO()

        major_space = max([len(str(idx)) for idx in self.major_axis]) + 4
        minor_space = max([len(str(idx)) for idx in self.minor_axis]) + 4
        
        for h in ['Major', 'Minor'] + list(self.fields):
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

    def swapAxes(self):
        """
        Swap major and minor axes and reorder values to be grouped by
        minor axis values
        """
        # Order everything by minor labels. Have to use mergesort
        # because NumPy quicksort is not stable. Here of course I'm
        # using the invariant that the major labels are ordered.
        indexer = self.minor_labels.argsort(kind='mergesort')

        new_major = self.minor_labels[indexer]
        new_minor = self.major_labels[indexer]
        new_values = self.values[indexer]
        
        return LongPanel(new_values, self.fields, self.minor_axis,
                         self.major_axis, new_major, new_minor,
                         mask=self._mask)
    
    def toWide(self):
        """
        Transform long (stacked) format into wide format

        Returns
        -------
        WidePanel
        """
        F, N, K = self.dims

        values = np.empty((F, N, K), dtype=float)

        mask = self.mask
        notmask = -self.mask
        
        for i in xrange(len(self.fields)):
            values[i].flat[mask] = self.values[:, i]
            values[i].flat[notmask] = np.NaN
            
        return WidePanel(values, self.fields, self.major_axis, self.minor_axis)


def _homogenize(frames):
    result = {}

    index = None
    columns = None

    for key, frame in frames.iteritems():
        if index is None:
            index = frame.index
        elif index is not frame.index:
            index = index.intersection(frame.index)

        if columns is None:
            columns = set(frame.cols())
        else:
            columns &= set(frame.cols())

    columns = sorted(columns)
            
    for key, frame in frames.iteritems():
        result[key] = frame.filterItems(columns).reindex(index)

    return result, index, columns


if __name__ == '__main__':
    from datetime import datetime
    import string

    import numpy as np

    from pandas.core.api import DataMatrix, DateRange
    from pandas.stats.linmodel import LinearModel, XSLinearModel

    N = 20
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

    long = panel.toLong()
    wide = long.toWide()
