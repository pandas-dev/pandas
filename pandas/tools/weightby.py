"""
functions to compute weighting Indexes
"""

import numpy as np
from pandas.types.generic import ABCSeries, ABCDataFrame
from pandas.compat import string_types
from pandas.util.decorators import Substitution

_shared_docs = {}
_shared_docs['weights'] = """weights : str or ndarray-like, optional
    Default 'None' results in equal probability weighting.

    If passed a Series, will align with target object on index.
    Index values in weights not found in the target object
    will be ignored and index values in the target object
    not in weights will be assigned weights of zero.

    If called on a DataFrame, will accept the name of a column
    when axis = 0.

    Unless weights are a Series, weights must be same length
    as axis of the target object.

    If weights do not sum to 1, they will be normalized to sum to 1.

    Missing values in the weights column will be treated as zero.
    inf and -inf values not allowed."""


@Substitution(weights=_shared_docs['weights'])
def weightby(obj, weights=None, axis=0):
    """returns a weights Series for the specified weights

Paramaters
----------
obj : Series/DataFrame
%(weights)s
axis : {0 (index), 1 (columns)}
    axis to compute weights on obj

Returns
-------
tuple of (obj, ndarray of weights, like indexed to obj)"""

    # If a series, align with frame
    if isinstance(weights, ABCSeries):
        weights = weights.reindex(obj.axes[axis])

    # Strings acceptable if a dataframe and axis = 0
    if isinstance(weights, string_types):

        # we use self.obj as we may have a selection here
        if isinstance(obj, ABCDataFrame):
            if axis == 0:
                try:
                    w, weights = weights, obj[weights]

                    # remove the weights column from obj
                    obj = obj.drop([w], axis=1)
                except KeyError:
                    raise KeyError("String passed to weights is not a "
                                   "valid column")
            else:
                raise ValueError("Strings can only be passed to "
                                 "weights when weighting by the rows on "
                                 "a DataFrame")
        else:
            raise ValueError("Strings cannot be passed as weights "
                             "when weighting from a Series or Panel.")

    from pandas import Series
    weights = Series(weights, dtype='float64')

    if len(weights) != len(obj.axes[axis]):
        raise ValueError("Weights and axis to be must be of "
                         "same length")

    if (weights == np.inf).any() or (weights == -np.inf).any():
        raise ValueError("weight vector may not include `inf` values")

    if (weights < 0).any():
        raise ValueError("weight vector many not include negative "
                         "values")

    # If has nan, set to zero.
    weights = weights.fillna(0)

    # Renormalize if don't sum to 1
    if weights.sum() != 1:
        if weights.sum() != 0:
            weights = weights / weights.sum()
        else:
            raise ValueError("Invalid weights: weights sum to zero")

    return obj, weights.values


def weight(values, weights):
    """
    Return the values * weights, broadcasting if needed

    Parameters
    ----------
    values : ndarray
    weights : 1d-ndarray

    Returns
    -------
    values shaped ndarray
    """

    if weights is None:
        return values

    if values.ndim == 1:
        return values * weights

    elif values.ndim == 2:

        return values * weights

    raise NotImplementedError
