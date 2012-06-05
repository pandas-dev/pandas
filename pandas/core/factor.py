# pylint: disable=E1101,W0232

import numpy as np
import pandas.core.common as com


class Factor(object):
    """
    Represents a categorical variable in classic R / S-plus fashion

    Parameters
    ----------
    labels : ndarray of integers
    levels : Index-like (unique)

    data : array-like

    Returns
    -------
    **Attributes**
      * labels : ndarray
      * levels : ndarray
    """
    def __init__(self, labels, levels, name=None):
        from pandas.core.index import _ensure_index

        levels = _ensure_index(levels)
        if not levels.is_unique:
            raise ValueError('Factor levels must be unique')

        self.labels = labels
        self.levels = levels
        self.name = name

    @classmethod
    def from_array(cls, data):
        from pandas.core.algorithms import factorize

        try:
            labels, levels, _ = factorize(data, sort=True)
        except TypeError:
            labels, levels, _ = factorize(data, sort=False)

        return Factor(labels, levels)

    levels = None

    def __array__(self, dtype=None):
        return com.take_1d(self.levels, self.labels)

    def __len__(self):
        return len(self.labels)

    def __repr__(self):
        temp = 'Factor:%s\n%s\nLevels (%d): %s'
        values = np.asarray(self)
        return temp % ('' if self.name is None else self.name,
                       repr(values), len(self.levels), self.levels)

    def __getitem__(self, key):
        if isinstance(key, (int, np.integer)):
            i = self.labels[key]
            if i == -1:
                return np.nan
            else:
                return self.levels[i]
        else:
            return Factor(self.labels[key], self.levels)


