import numpy as np
import pandas.core.common as com
import pandas._tseries as lib


class Factor(np.ndarray):
    """
    Represents a categorical variable in classic R / S-plus fashion

    Parameters
    ----------
    data : array-like

    Returns
    -------
    **Attributes**
      * labels : ndarray
      * levels : ndarray
    """
    def __new__(cls, data):
        data = np.asarray(data, dtype=object)
        levels, factor = unique_with_labels(data)
        factor = factor.view(Factor)
        factor.levels = levels
        return factor

    levels = None

    def __array_finalize__(self, obj):
        self.levels = getattr(obj, 'levels', None)

    @property
    def labels(self):
        return self.view(np.ndarray)

    def asarray(self):
        return np.asarray(self.levels).take(self.labels)

    def __len__(self):
        return len(self.labels)

    def __repr__(self):
        temp = 'Factor:\n%s\nLevels (%d): %s'
        values = self.asarray()
        return temp % (repr(values), len(self.levels), self.levels)

    def __getitem__(self, key):
        if isinstance(key, (int, np.integer)):
            i = self.labels[key]
            return self.levels[i]
        else:
            return np.ndarray.__getitem__(self, key)


def unique_with_labels(values):
    from pandas.core.index import Index
    rizer = lib.Factorizer(len(values))
    labels, _ = rizer.factorize(values, sort=False)
    uniques = Index(rizer.uniques)
    labels = com._ensure_platform_int(labels)
    try:
        sorter = uniques.argsort()
        reverse_indexer = np.empty(len(sorter), dtype=np.int_)
        reverse_indexer.put(sorter, np.arange(len(sorter)))
        labels = reverse_indexer.take(labels)
        uniques = uniques.take(sorter)
    except TypeError:
        pass

    return uniques, labels

