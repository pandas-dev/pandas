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

#-------------------------------------------------------------------------------
# Groupable mixin


class Groupable(object):

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
