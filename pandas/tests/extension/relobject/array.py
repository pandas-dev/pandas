import sys
import random
import numbers


import numpy as np

import pandas as pd
from pandas.core.arrays import ExtensionArray
from pandas.core.dtypes.base import ExtensionDtype
from pandas.core.dtypes.common import _ensure_platform_int

# Idea here is to test if we have an array of objects where
# relational operators are defined that return strings rather
# than booleans, but not for all of the operators. Also test
# if things work when _can_hold_na == False.
# NOTE: It is not expected that sorting and groupby() will
# work on these objects !!


def relation(left, op, right):
    return str(left) + op + str(right)


class RelObj(object):
    def __init__(self, value):
        if type(value) == int:
            self._value = value
        else:
            raise TypeError("Value should be int")

    def __eq__(self, other):
        if isinstance(other, RelObj):
            return relation(self, " == ", other)
        else:
            raise TypeError("Cannot compare RelObj to " +
                            "object of different type")

    def __le__(self, other):
        if isinstance(other, RelObj):
            return relation(self, " <= ", other)
        else:
            raise TypeError("Cannot compare RelObj to " +
                            "object of different type")

    def __ge__(self, other):
        if isinstance(other, RelObj):
            return relation(self, " >= ", other)
        else:
            raise TypeError("Cannot compare RelObj to " +
                            "object of different type")

    def __lt__(self, other):
        raise Exception('lt not supported')

    def __gt__(self, other):
        raise Exception('gt not supported')

    def __ne__(self, other):
        raise Exception('ne not supported')

    @property
    def value(self):
        return self._value

    def __hash__(self):
        return id(self)

    def __repr__(self):
        return "RelObjValue: " + str(self._value)


class RelObjectDtype(ExtensionDtype):
    type = RelObj
    name = 'relobj'
    kind = 'O'

    @classmethod
    def construct_from_string(cls, string):
        if string == cls.name:
            return cls()
        else:
            raise TypeError("Cannot construct a '{}' from "
                            "'{}'".format(cls, string))


class RelObjectArray(ExtensionArray):
    dtype = RelObjectDtype()

    def __init__(self, values):
        values = np.asarray(values, dtype=object)

        self.values = values

    @classmethod
    def _from_sequence(cls, scalars):
        return cls(scalars)

    def __getitem__(self, item):
        if isinstance(item, numbers.Integral):
            return self.values[item]
        else:
            return type(self)(self.values[item])

    def copy(self, deep=False):
        if deep:
            return type(self)(self.values.copy())
        return type(self)(self)

    def __setitem__(self, key, value):
        if pd.api.types.is_list_like(value):
            value = [int(v) for v in value]
        else:
            value = int(value)
        self.values[key] = value

    def __len__(self):
        return len(self.values)

    def __repr__(self):
        return 'RelObjArray({!r})'.format([i for i in self.values])

    @property
    def nbytes(self):
        n = len(self)
        if n:
            return n * sys.getsizeof(self[0])
        return 0

    def isna(self):
        return np.array([False]*len(self.values))

    def _values_for_argsort(self):
        # type: () -> ndarray
        """Return values for sorting.

        Returns
        -------
        ndarray
            The transformed values should maintain the ordering between values
            within the array.

        See Also
        --------
        ExtensionArray.argsort
        """
        # Note: this is used in `ExtensionArray.argsort`.
        # We will sort based on the value of the object
        return np.array([ro.value for ro in self.values])

    def unique(self):
        """Compute unique values using the ID of the objects

        Cannot use pandas.unique() because it requires __eq__()
        to return a boolean

        Returns
        -------
        uniques : RelObjArray
        """
        seen = set()
        uniques = [x for x in self.values if x not in seen and not seen.add(x)]
        return self._from_sequence(uniques)

    def factorize(self, na_sentinel=-1):
        # type: (int) -> Tuple[ndarray, ExtensionArray]
        """Encode the extension array as an enumerated type.

        Parameters
        ----------
        na_sentinel : int, default -1
            Value to use in the `labels` array to indicate missing values.

        Returns
        -------
        labels : ndarray
            An integer NumPy array that's an indexer into the original
            ExtensionArray.
        uniques : ExtensionArray
            An ExtensionArray containing the unique values of `self`.

            .. note::

               uniques will *not* contain an entry for the NA value of
               the ExtensionArray if there are any missing values present
               in `self`.
        """
        uniques = self.unique()
        dun = {id(v): i for i, v in enumerate(uniques)}
        labels = np.array([dun[id(v)] for v in self.values], dtype=np.intp)
        return labels, uniques

    def take(self, indexer, allow_fill=True, fill_value=None):
        indexer = np.asarray(indexer)
        mask = indexer == -1

        indexer = _ensure_platform_int(indexer)
        out = self.values.take(indexer)
        out[mask] = np.nan

        return type(self)(out)

    @classmethod
    def _concat_same_type(cls, to_concat):
        return cls(np.concatenate([x.values for x in to_concat]))

    _can_hold_na = False


def make_data():
    return [RelObj(random.randint(0, 2000)) for _ in range(100)]
