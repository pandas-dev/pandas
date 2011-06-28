import numpy as np
import cPickle

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

class NDFrame(object):
    """
    N-dimensional labeled array data structure with potentially heterogenous
    dtypes along one axis
    """
    def __init__(self, data):
        pass

class PandasGeneric(Picklable):

    _AXIS_NUMBERS = {
        'index' : 0,
        'columns' : 1
    }

    _AXIS_ALIASES = {}

    _AXIS_NAMES = dict((v, k) for k, v in _AXIS_NUMBERS.iteritems())

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

