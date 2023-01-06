from __future__ import annotations

from typing import (
    Any,
    Callable,
)

import numpy as np

from pandas._libs import lib

import pandas as pd
import pandas.core.common as com


def case_when(*args, default: Any = lib.no_default) -> Callable:
    """
    Create a callable for assignment based on a condition or multiple conditions.

    This is useful when you want to assign a column based on multiple conditions.

    Parameters
    ----------
    args : Variable argument of conditions and expected values.
        Takes the form:
            `condition0`, `value0`, `condition1`, `value1`, ...
        `condition` can be a 1-D boolean array/series or a callable
        that evaluate to a 1-D boolean array/series.
    default : Any, default is `None`.
        The default value to be used if all conditions evaluate False.

    Returns
    -------
    Callable
        The Callable returned in `case_when` can be used with `df.assign(...)`
        for multi-condition assignment. See examples below for more info.

    See Also
    --------
    DataFrame.assign: Assign new columns to a DataFrame.

    Examples
    --------
    >>> df = pd.DataFrame(dict(a=[1, 2, 3], b=[4, 5, 6]))
    >>> df
       a  b
    0  1  4
    1  2  5
    2  3  6

    >>> df.assign(
    ...     new_column = pd.case_when(
    ...         lambda x: x.a == 1, 'first',
    ...         lambda x: (x.a > 1) & (x.b == 5), 'second',
    ...         default='default'
    ...     )
    ... )
       a  b new_column
    0  1  4      first
    1  2  5     second
    2  3  6    default
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

    if default is lib.no_default:
        default = None

    def _eval(df: pd.DataFrame) -> np.ndarray:
        booleans = []
        replacements = []

        for index, value in enumerate(args):
            if not index % 2:
                if callable(value):
                    value = com.apply_if_callable(value, df)
                booleans.append(value)
            else:
                replacements.append(value)

        return np.select(booleans, replacements, default=default)

    return _eval
