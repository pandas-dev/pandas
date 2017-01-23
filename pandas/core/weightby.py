"""

provide a lazy structure to support
weights for calculation
similar to how we have a Groupby object


"""

import numpy as np
import pandas as pd

from pandas import Series, DataFrame
from pandas.compat import string_types, set_function_name
from pandas.types.generic import ABCSeries, ABCDataFrame
from pandas.types.common import is_scalar, is_list_like
from pandas.core.base import PandasObject, SelectionMixin
from pandas.util.decorators import Appender, Substitution


class Weightby(PandasObject, SelectionMixin):
    _attributes = ['weights', 'axis']

    def __init__(self, obj, weights=None, axis=0):

        self.exclusions = set()
        self._weights = None
        self.weights = weights
        self.axis = axis
        self.obj = obj

    def _gotitem(self, key, ndim, subset=None):
        """
        sub-classes to define
        return a sliced object

        Parameters
        ----------
        key : string / list of selections
        ndim : 1,2
            requested ndim of result
        subset : object, default None
            subset to act on
        """

        # create a new object to prevent aliasing
        if subset is None:
            subset = self.obj

        newself = self._shallow_copy(subset, obj_type=type(self))
        newself._reset_cache()
        if subset.ndim == 2:
            if is_scalar(key) and key in subset or is_list_like(key):
                newself._selection = key
        return newself

    def __getattr__(self, attr):
        if attr in self._internal_names_set:
            return object.__getattribute__(self, attr)
        if attr in self.obj:
            return self[attr]

        raise AttributeError("%r object has no attribute %r" %
                             (type(self).__name__, attr))

    def _compute_weights(self):
        """
        compute our _weights
        """
        if self._weights is not None:
            return self._weights

        obj = self._selected_obj

        weights = self.weights
        axis = self.axis

        # If a series, align with frame
        if isinstance(weights, Series):
            weights = weights.reindex(obj.axes[axis])

        # Strings acceptable if a dataframe and axis = 0
        if isinstance(weights, string_types):

            # we use self.obj as we may have a selection here
            if isinstance(self.obj, pd.DataFrame):
                if axis == 0:
                    try:

                        # exclude this as an aggregator
                        self.exclusions.add(weights)

                        weights = self.obj[weights]

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

        self._weights = weights.values
        return self._weights

    def _apply(self, func, *args, **kwargs):
        """
        Apply the function with weights

        Parameters
        ----------
        func : string/callable to apply

        Returns
        -------
        y : type of input
        """

        weights = self._compute_weights()

        # we may need to drop the dim
        # before operations
        obj = self._obj_with_exclusions
        if self._selection is not None:
            obj = obj[self._selection]

        f = getattr(obj, func)

        kwargs['axis'] = self.axis
        kwargs['_weights'] = weights

        result = f(*args, **kwargs)
        result = self._wrap_results(result)
        return result

    def _wrap_results(self, result):
        return result


class SeriesWeightBy(Weightby):

    @property
    def _constructor(self):
        return Series

    @Substitution(name='weightby')
    @Appender(SelectionMixin._see_also_template)
    @Appender(SelectionMixin._agg_doc)
    def aggregate(self, arg, *args, **kwargs):
        return super(SeriesWeightBy, self).aggregate(arg, *args, **kwargs)

    agg = aggregate

    @Appender(Series.sample.__doc__)
    def sample(self, n=None, frac=None, replace=False,
               random_state=None):
        return self._apply('sample', n=n, frac=frac, replace=replace,
                           random_state=random_state)


class DataFrameWeightBy(Weightby):

    @property
    def _constructor(self):
        return DataFrame

    @Substitution(name='weightby')
    @Appender(SelectionMixin._see_also_template)
    @Appender(SelectionMixin._agg_doc)
    def aggregate(self, arg, *args, **kwargs):
        return super(DataFrameWeightBy, self).aggregate(arg, *args, **kwargs)

    agg = aggregate

    @Appender(DataFrame.sample.__doc__)
    def sample(self, n=None, frac=None, replace=False,
               random_state=None):
        return self._apply('sample', n=n, frac=frac, replace=replace,
                           random_state=random_state)


def _add_stat_function(cls, ref_obj, name):

    @Appender(getattr(ref_obj, name).__doc__)
    def stat_func(self, axis=None, skipna=None, level=None, numeric_only=None,
                  **kwargs):
        return self._apply(name, axis=axis, skipna=skipna, level=level,
                           numeric_only=numeric_only, **kwargs)

    setattr(cls, name, set_function_name(stat_func, name, cls))


# add in stat methods
for method in ['sum', 'mean', 'std', 'var',
               'sem', 'kurt', 'skew', 'sem']:

    _add_stat_function(SeriesWeightBy, Series, method)
    _add_stat_function(DataFrameWeightBy, DataFrame, method)


# Top-level exports
def weightby(obj, *args, **kwds):
    if isinstance(obj, ABCSeries):
        klass = SeriesWeightBy
    elif isinstance(obj, ABCDataFrame):
        klass = DataFrameWeightBy
    else:
        raise TypeError('invalid type: %s' % type(obj))

    return klass(obj, *args, **kwds)


weightby.__doc__ = Weightby.__doc__
