from __future__ import annotations

from typing import Any
import warnings

from pandas.util._exceptions import find_stack_level

from pandas.core.dtypes.common import is_list_like

import pandas as pd
import pandas.core.common as com


def warn_and_override_index(series, series_type, index):
    warnings.warn(
        f"Series {series_type} will be reindexed to match obj index.",
        UserWarning,
        stacklevel=find_stack_level(),
    )
    return pd.Series(series.values, index=index)


def case_when(obj: pd.DataFrame | pd.Series, *args, default: Any) -> pd.Series:
    """
    Returns a Series based on multiple conditions assignment.

    This is useful when you want to assign a column based on multiple conditions.
    Uses `Series.mask` to perform the assignment.

    The returned Series have the same index as `obj`.

    Parameters
    ----------
    obj : Dataframe or Series on which the conditions will be applied.
    args : Variable argument of conditions and expected values.
        Takes the form:
            `condition0`, `value0`, `condition1`, `value1`, ...
        `condition` can be a 1-D boolean array/series or a callable
        that evaluate to a 1-D boolean array/series. See examples below.
    default : Any
        The default value to be used if all conditions evaluate False. This value
        will be used to create the `Series` on which `Series.mask` will be called.
        If this value is not already an array like (i.e. it is not of type `Series`,
        `np.array` or `list`) it will be repeated `obj.shape[0]` times in order to
        create an array like object from it and then apply the `Series.mask`.

    Returns
    -------
    Series
        Series with the corresponding values based on the conditions, their values
        and the default value.


    Examples
    --------
    >>> df = pd.DataFrame(dict(a=[1, 2, 3], b=[4, 5, 6]))
    >>> df
       a  b
    0  1  4
    1  2  5
    2  3  6

    >>> pd.case_when(
    ...     df,
    ...     lambda x: x.a == 1,
    ...     'first',
    ...     lambda x: (x.a == 2) & (x.b == 5),
    ...     'second',
    ...     default='default',
    ... )
    0      first
    1     second
    2    default
    dtype: object

    >>> pd.case_when(
    ...     df,
    ...     lambda x: (x.a == 1) & (x.b == 4),
    ...     df.b,
    ...     default=0,
    ... )
    0    4
    1    0
    2    0
    dtype: int64

    >>> pd.case_when(
    ...     df,
    ...     lambda x: (x.a > 1) & (x.b > 1),
    ...     -1,
    ...     default=df.a,
    ... )
    0    1
    1   -1
    2   -1
    Name: a, dtype: int64

    >>> pd.case_when(
    ...     df.a,
    ...     lambda x: x == 1,
    ...     -1,
    ...     default=df.a,
    ... )
    0   -1
    1    2
    2    3
    Name: a, dtype: int64

    >>> pd.case_when(
    ...     df.a,
    ...     df.a  > 1,
    ...     -1,
    ...     default=df.a,
    ... )
    0    1
    1   -1
    2   -1
    Name: a, dtype: int64

    The index will always follow that of `obj`. For example:
    >>> df = pd.DataFrame(
    ...     dict(a=[1, 2, 3], b=[4, 5, 6]),
    ...     index=['index 1', 'index 2', 'index 3']
    ... )
    >>> df
             a  b
    index 1  1  4
    index 2  2  5
    index 3  3  6

    >>> pd.case_when(
    ...     df,
    ...     lambda x: (x.a == 1) & (x.b == 4),
    ...     df.b,
    ...     default=0,
    ... )
    index 1    4
    index 2    0
    index 3    0
    dtype: int64
    """
    len_args = len(args)

    if len_args < 2:
        raise ValueError("At least two arguments are required for `case_when`")
    if len_args % 2:
        raise ValueError(
            "The number of conditions and values do not match. "
            f"There are {len_args - len_args//2} conditions "
            f"and {len_args//2} values."
        )

    # construct series on which we will apply `Series.mask`
    if is_list_like(default):
        series = pd.Series(default.values, index=obj.index)
    else:
        series = pd.Series([default] * obj.shape[0], index=obj.index)

    for i in range(0, len_args, 2):
        # get conditions
        if callable(args[i]):
            conditions = com.apply_if_callable(args[i], obj)
        else:
            conditions = args[i]

        # get replacements
        replacements = args[i + 1]

        if isinstance(replacements, pd.Series) and not replacements.index.equals(
            obj.index
        ):
            replacements = warn_and_override_index(
                replacements, f"(in args[{i+1}])", obj.index
            )

        if isinstance(conditions, pd.Series) and not conditions.index.equals(obj.index):
            conditions = warn_and_override_index(
                conditions, f"(in args[{i}])", obj.index
            )

        # `Series.mask` call
        series = series.mask(conditions, replacements)

    return series
