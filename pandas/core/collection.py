from collections import defaultdict
from pandas.core.mixins import Picklable
from pandas.core.index import Index
from pandas.core.pytools import rands, adjoin, groupby
import cPickle
import os

__all__ = ['PickleContainer']

class PickleContainer(Picklable):
    """
    Store collection of objects on disk with this dict-like object.

    Parameters
    ----------
    dirPath: string
       Directory where to store the objects
    lruSize: int
       Number of objects to keep in memory (not implemented yet)
    """
    def __init__(self, dirPath, lruSize=5):
        self.dirPath = dirPath
        if not os.path.exists(dirPath):
            os.mkdir(dirPath)

        self._lruSize = lruSize

        self._paths = {}
        self._classes = {}
        self._lru = {}

    def __repr__(self):
        output = str(self.__class__) + '\n'
        keys, values = zip(*self._classes.iteritems())
        output += adjoin(5, map(str, keys), map(str, values))
        return output

    def __setitem__(self, key, value):
        theKey = rands(10)
        filePath = self.dirPath + '/' + theKey

        self._paths[key] = filePath

        if isinstance(value, Picklable):
            value.save(filePath)
        else:
            f = open(filePath, 'w')
            try:
                cPickle.dump(value, f)
            finally:
                f.close()

        self._paths[key] = filePath
        self._classes[key] = value.__class__

    def __getitem__(self, key):
        if key not in self._paths:
            raise Exception('Requested key not in this container!')

        thePath = self._paths[key]
        theClass = self._classes[key]

        if issubclass(theClass, Picklable):
            obj = theClass.load(thePath)
        else:
            f = open(thePath, 'rb')
            try:
                obj = cPickle.load(f)
            finally:
                f.close()

        return obj

    def __delitem__(self, key):
        del self._paths[key]
        del self._classes[key]

    def __iter__(self):
        return iter(self._paths)

    def iteritems(self):
        for key, path in self._paths.iteritems():
            yield key, self[key]

    def keys(self):
        return self._paths.keys()

    def values(self):
        result = []
        for key in self._paths:
            result.append(self[key])
        return result
