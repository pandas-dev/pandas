# pylint: disable=E1101,W0232

import numpy as np

from pandas.core.algorithms import factorize
from pandas.core.index import Index
import pandas.core.common as com


def _cat_compare_op(op):
    def f(self, other):
        if isinstance(other, (Categorical, np.ndarray)):
            values = np.asarray(self)
            f = getattr(values, op)
            return f(np.asarray(other))
        else:
            if other in self.levels:
                i = self.levels.get_loc(other)
                return getattr(self.labels, op)(i)
            else:
                return np.repeat(False, len(self))

    f.__name__ = op

    return f


class Categorical(object):
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
        self.labels = labels
        self.levels = levels
        self.name = name

    @classmethod
    def from_array(cls, data):
        if isinstance(data, Index) and hasattr(data, 'factorize'):
            labels, levels = data.factorize()
        else:
            try:
                labels, levels = factorize(data, sort=True)
            except TypeError:
                labels, levels = factorize(data, sort=False)

        return Categorical(labels, levels,
                           name=getattr(data, 'name', None))

    _levels = None

    def _set_levels(self, levels):
        from pandas.core.index import _ensure_index

        levels = _ensure_index(levels)
        if not levels.is_unique:
            raise ValueError('Categorical levels must be unique')
        self._levels = levels

    def _get_levels(self):
        return self._levels

    levels = property(fget=_get_levels, fset=_set_levels)

    __eq__ = _cat_compare_op('__eq__')
    __ne__ = _cat_compare_op('__ne__')
    __lt__ = _cat_compare_op('__lt__')
    __gt__ = _cat_compare_op('__gt__')
    __le__ = _cat_compare_op('__le__')
    __ge__ = _cat_compare_op('__ge__')

    def __array__(self, dtype=None):
        return com.take_1d(self.levels.values, self.labels)

    def __len__(self):
        return len(self.labels)

    def __repr__(self):
        temp = 'Categorical: %s\n%s\n%s'
        values = np.asarray(self)
        levheader = 'Levels (%d): ' % len(self.levels)
        levstring = np.array_repr(self.levels,
                                  max_line_width=60)

        indent = ' ' * (levstring.find('[') + len(levheader) + 1)
        lines = levstring.split('\n')
        levstring = '\n'.join([lines[0]] +
                              [indent + x.lstrip() for x in lines[1:]])

        return temp % ('' if self.name is None else self.name,
                       repr(values), levheader + levstring)

    def __getitem__(self, key):
        if isinstance(key, (int, np.integer)):
            i = self.labels[key]
            if i == -1:
                return np.nan
            else:
                return self.levels[i]
        else:
            return Categorical(self.labels[key], self.levels)

    def equals(self, other):
        """
        Returns True if categorical arrays are equal

        Parameters
        ----------
        other : Categorical

        Returns
        -------
        are_equal : boolean
        """
        if not isinstance(other, Categorical):
            return False

        return (self.levels.equals(other.levels) and
                np.array_equal(self.labels, other.labels))

Factor = Categorical
