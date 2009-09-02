from __future__ import with_statement

import cPickle

#-------------------------------------------------------------------------------
# Picklable mixin

class Picklable(object):
    def save(self, fileName):
        with open(fileName, 'wb') as f:
            cPickle.dump(self, f)

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
        from pandas.core.groupby import GroupBy
        return GroupBy(self, mapper)

