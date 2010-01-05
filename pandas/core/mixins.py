from __future__ import with_statement

import cPickle

#-------------------------------------------------------------------------------
# Picklable mixin

class Picklable(object):
    def save(self, fileName):
        with open(fileName, 'wb') as f:
            cPickle.dump(self, f, protocol=cPickle.HIGHEST_PROTOCOL)

    @classmethod
    def load(cls, fileName):
        with open(fileName, 'rb') as f:
            obj = cPickle.load(f)
            return obj
        raise Exception("Error trying to unpickle %s" % fileName)

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

        from pandas.core.groupby import GroupBy
        return GroupBy(self, mapper)


