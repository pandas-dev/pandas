"""Common utility functions for rolling operations"""
from collections import defaultdict
from typing import Callable, Optional
import warnings

import numpy as np

from pandas.core.dtypes.generic import ABCDataFrame, ABCSeries

from pandas.core.generic import _shared_docs
from pandas.core.groupby.base import GroupByMixin
from pandas.core.indexes.api import MultiIndex

_shared_docs = dict(**_shared_docs)
_doc_template = """
        Returns
        -------
        Series or DataFrame
            Return type is determined by the caller.

        See Also
        --------
        pandas.Series.%(name)s : Calling object with Series data.
        pandas.DataFrame.%(name)s : Calling object with DataFrame data.
        pandas.Series.%(func_name)s : Similar method for Series.
        pandas.DataFrame.%(func_name)s : Similar method for DataFrame.
"""


def _dispatch(name: str, *args, **kwargs):
    """
    Dispatch to apply.
    """

    def outer(self, *args, **kwargs):
        def f(x):
            x = self._shallow_copy(x, groupby=self._groupby)
            return getattr(x, name)(*args, **kwargs)

        return self._groupby.apply(f)

    outer.__name__ = name
    return outer


class WindowGroupByMixin(GroupByMixin):
    """
    Provide the groupby facilities.
    """

    def __init__(self, obj, *args, **kwargs):
        kwargs.pop("parent", None)
        groupby = kwargs.pop("groupby", None)
        if groupby is None:
            groupby, obj = obj, obj._selected_obj
        self._groupby = groupby
        self._groupby.mutated = True
        self._groupby.grouper.mutated = True
        super().__init__(obj, *args, **kwargs)

    count = _dispatch("count")
    corr = _dispatch("corr", other=None, pairwise=None)
    cov = _dispatch("cov", other=None, pairwise=None)

    def _apply(
        self,
        func: Callable,
        center: bool,
        require_min_periods: int = 0,
        floor: int = 1,
        is_weighted: bool = False,
        name: Optional[str] = None,
        use_numba_cache: bool = False,
        **kwargs,
    ):
        """
        Dispatch to apply; we are stripping all of the _apply kwargs and
        performing the original function call on the grouped object.
        """
        kwargs.pop("floor", None)
        kwargs.pop("original_func", None)

        # TODO: can we de-duplicate with _dispatch?
        def f(x, name=name, *args):
            x = self._shallow_copy(x)

            if isinstance(name, str):
                return getattr(x, name)(*args, **kwargs)

            return x.apply(name, *args, **kwargs)

        return self._groupby.apply(f)


def _flex_binary_moment(arg1, arg2, f, pairwise=False):

    if not (
        isinstance(arg1, (np.ndarray, ABCSeries, ABCDataFrame))
        and isinstance(arg2, (np.ndarray, ABCSeries, ABCDataFrame))
    ):
        raise TypeError(
            "arguments to moment function must be of type np.ndarray/Series/DataFrame"
        )

    if isinstance(arg1, (np.ndarray, ABCSeries)) and isinstance(
        arg2, (np.ndarray, ABCSeries)
    ):
        X, Y = prep_binary(arg1, arg2)
        return f(X, Y)

    elif isinstance(arg1, ABCDataFrame):
        from pandas import DataFrame

        def dataframe_from_int_dict(data, frame_template):
            result = DataFrame(data, index=frame_template.index)
            if len(result.columns) > 0:
                result.columns = frame_template.columns[result.columns]
            return result

        results = {}
        if isinstance(arg2, ABCDataFrame):
            if pairwise is False:
                if arg1 is arg2:
                    # special case in order to handle duplicate column names
                    for i, col in enumerate(arg1.columns):
                        results[i] = f(arg1.iloc[:, i], arg2.iloc[:, i])
                    return dataframe_from_int_dict(results, arg1)
                else:
                    if not arg1.columns.is_unique:
                        raise ValueError("'arg1' columns are not unique")
                    if not arg2.columns.is_unique:
                        raise ValueError("'arg2' columns are not unique")
                    with warnings.catch_warnings(record=True):
                        warnings.simplefilter("ignore", RuntimeWarning)
                        X, Y = arg1.align(arg2, join="outer")
                    X = X + 0 * Y
                    Y = Y + 0 * X

                    with warnings.catch_warnings(record=True):
                        warnings.simplefilter("ignore", RuntimeWarning)
                        res_columns = arg1.columns.union(arg2.columns)
                    for col in res_columns:
                        if col in X and col in Y:
                            results[col] = f(X[col], Y[col])
                    return DataFrame(results, index=X.index, columns=res_columns)
            elif pairwise is True:
                results = defaultdict(dict)
                for i, k1 in enumerate(arg1.columns):
                    for j, k2 in enumerate(arg2.columns):
                        if j < i and arg2 is arg1:
                            # Symmetric case
                            results[i][j] = results[j][i]
                        else:
                            results[i][j] = f(
                                *prep_binary(arg1.iloc[:, i], arg2.iloc[:, j])
                            )

                from pandas import concat

                result_index = arg1.index.union(arg2.index)
                if len(result_index):

                    # construct result frame
                    result = concat(
                        [
                            concat(
                                [results[i][j] for j, c in enumerate(arg2.columns)],
                                ignore_index=True,
                            )
                            for i, c in enumerate(arg1.columns)
                        ],
                        ignore_index=True,
                        axis=1,
                    )
                    result.columns = arg1.columns

                    # set the index and reorder
                    if arg2.columns.nlevels > 1:
                        result.index = MultiIndex.from_product(
                            arg2.columns.levels + [result_index]
                        )
                        # GH 34440
                        num_levels = len(result.index.levels)
                        new_order = [num_levels - 1] + list(range(num_levels - 1))
                        result = result.reorder_levels(new_order).sort_index()
                    else:
                        result.index = MultiIndex.from_product(
                            [range(len(arg2.columns)), range(len(result_index))]
                        )
                        result = result.swaplevel(1, 0).sort_index()
                        result.index = MultiIndex.from_product(
                            [result_index] + [arg2.columns]
                        )
                else:

                    # empty result
                    result = DataFrame(
                        index=MultiIndex(
                            levels=[arg1.index, arg2.columns], codes=[[], []]
                        ),
                        columns=arg2.columns,
                        dtype="float64",
                    )

                # reset our index names to arg1 names
                # reset our column names to arg2 names
                # careful not to mutate the original names
                result.columns = result.columns.set_names(arg1.columns.names)
                result.index = result.index.set_names(
                    result_index.names + arg2.columns.names
                )

                return result

            else:
                raise ValueError("'pairwise' is not True/False")
        else:
            results = {
                i: f(*prep_binary(arg1.iloc[:, i], arg2))
                for i, col in enumerate(arg1.columns)
            }
            return dataframe_from_int_dict(results, arg1)

    else:
        return _flex_binary_moment(arg2, arg1, f)


def zsqrt(x):
    with np.errstate(all="ignore"):
        result = np.sqrt(x)
        mask = x < 0

    if isinstance(x, ABCDataFrame):
        if mask._values.any():
            result[mask] = 0
    else:
        if mask.any():
            result[mask] = 0

    return result


def prep_binary(arg1, arg2):
    if not isinstance(arg2, type(arg1)):
        raise Exception("Input arrays must be of the same type!")

    # mask out values, this also makes a common index...
    X = arg1 + 0 * arg2
    Y = arg2 + 0 * arg1

    return X, Y
