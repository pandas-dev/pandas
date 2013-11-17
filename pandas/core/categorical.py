# pylint: disable=E1101,W0232

import numpy as np

from pandas import compat
from pandas.compat import u

from pandas.core.algorithms import factorize
from pandas.core.base import PandasObject
from pandas.core.index import Index
import pandas.core.common as com
from pandas.util.terminal import get_terminal_size
from pandas.core.config import get_option
from pandas.core import format as fmt


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


class Categorical(PandasObject):

    """
    Represents a categorical variable in classic R / S-plus fashion

    Parameters
    ----------
    labels : ndarray of integers
        If levels is given, the integer at label `i` is the index of the level
        for that label. I.e., the level at labels[i] is levels[labels[i]].
        Otherwise, if levels is None, these are just the labels and the levels
        are assumed to be the unique labels. See from_array.
    levels : Index-like (unique), optional
        The unique levels for each label. If not given, the levels are assumed
        to be the unique values of labels.
    name : str, optional
        Name for the Categorical variable. If levels is None, will attempt
        to infer from labels.

    Returns
    -------
    **Attributes**
      * labels : ndarray
      * levels : ndarray

    Examples
    --------
    >>> from pandas import Categorical
    >>> Categorical([0, 1, 2, 0, 1, 2], [1, 2, 3])
    Categorical:
    array([1, 2, 3, 1, 2, 3])
    Levels (3): Int64Index([1, 2, 3])

    >>> Categorical([0,1,2,0,1,2], ['a', 'b', 'c'])
    Categorical:
    array(['a', 'b', 'c', 'a', 'b', 'c'], dtype=object)
    Levels (3): Index(['a', 'b', 'c'], dtype=object)

    >>> Categorical(['a', 'b', 'c', 'a', 'b', 'c'])
    Categorical:
    array(['a', 'b', 'c', 'a', 'b', 'c'], dtype=object)
    Levels (3): Index(['a', 'b', 'c'], dtype=object)
    """

    def __init__(self, labels, levels=None, name=None):
        if levels is None:
            if name is None:
                name = getattr(labels, 'name', None)
            if isinstance(labels, Index) and hasattr(labels, 'factorize'):
                labels, levels = labels.factorize()
            else:
                try:
                    labels, levels = factorize(labels, sort=True)
                except TypeError:
                    labels, levels = factorize(labels, sort=False)

        self.labels = labels
        self.levels = levels
        self.name = name

    @classmethod
    def from_array(cls, data):
        """
        Make a Categorical type from a single array-like object.

        Parameters
        ----------
        data : array-like
            Can be an Index or array-like. The levels are assumed to be
            the unique values of `data`.
        """
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

    def _tidy_repr(self, max_vals=20):
        num = max_vals // 2
        head = self[:num]._get_repr(length=False, name=False, footer=False)
        tail = self[-(max_vals - num):]._get_repr(length=False,
                                                  name=False,
                                                  footer=False)

        result = '%s\n...\n%s' % (head, tail)
        # TODO: tidy_repr for footer since there may be a ton of levels?
        result = '%s\n%s' % (result, self._repr_footer())

        return compat.text_type(result)

    def _repr_footer(self):
        levheader = 'Levels (%d): ' % len(self.levels)
        # TODO: should max_line_width respect a setting?
        levstring = np.array_repr(self.levels, max_line_width=60)
        indent = ' ' * (levstring.find('[') + len(levheader) + 1)
        lines = levstring.split('\n')
        levstring = '\n'.join([lines[0]] +
                              [indent + x.lstrip() for x in lines[1:]])

        namestr = "Name: %s, " % self.name if self.name is not None else ""
        return u('%s\n%sLength: %d' % (levheader + levstring, namestr,
                                       len(self)))

    def _get_repr(self, name=False, length=True, na_rep='NaN', footer=True):
        formatter = fmt.CategoricalFormatter(self, name=name,
                                             length=length, na_rep=na_rep,
                                             footer=footer)
        result = formatter.to_string()
        return compat.text_type(result)

    def __unicode__(self):
        width, height = get_terminal_size()
        max_rows = (height if get_option("display.max_rows") == 0
                    else get_option("display.max_rows"))
        if len(self.labels) > (max_rows or 1000):
            result = self._tidy_repr(min(30, max_rows) - 4)
        elif len(self.labels) > 0:
            result = self._get_repr(length=len(self) > 50,
                                    name=True)
        else:
            result = 'Categorical([], %s' % self._get_repr(name=True,
                                                           length=False,
                                                           footer=True,
                                                           )

        return result

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

    def describe(self):
        """
        Returns a dataframe with frequency and counts by level.
        """
        # Hack?
        from pandas.core.frame import DataFrame
        grouped = DataFrame(self.labels).groupby(0)
        counts = grouped.count().values.squeeze()
        freqs = counts / float(counts.sum())
        return DataFrame.from_dict({
            'counts': counts,
            'freqs': freqs,
            'levels': self.levels
        }).set_index('levels')
