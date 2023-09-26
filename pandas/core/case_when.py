from __future__ import annotations

import numpy as np

from pandas.core.dtypes.cast import (
    construct_1d_arraylike_from_scalar,
    find_common_type,
)
from pandas.core.dtypes.common import (
    is_bool_dtype,
    is_scalar,
)
from pandas.core.dtypes.generic import ABCSeries
from pandas.core.dtypes.missing import na_value_for_dtype

from pandas.core.common import convert_to_list_like
from pandas.core.construction import array as pd_array


def case_when(
    *args,
    default=None,
    level=None,
):
    """
    Replace values where the conditions are True.

    Parameters
    ----------
    args : Variable argument of conditions and expected replacements.
        Takes the form:
            `condition0`, `replacement0`,
            `condition1`, `replacement1`, ... .
        `condition` should be a 1-D boolean array.
         The provided `conditions` should have the same size.
        `replacement` should be a 1-D array or a scalar.
        If `replacement` is a 1-D array, it should have the same
        shape as the paired `condition`.
        When multiple conditions are satisfied, the first one is used.

    default: scalar, one-dimensional array, default None
        If provided, it is the replacement value to use
        if all conditions evaluate to False.
        If default is a 1-D array, it should have the same shape as
        `condition` and `replacement`.

    level : int, default None
        Alignment level if needed.

        .. versionchanged:: 2.2.0

    Returns
    -------
    Series

    See Also
    --------
    Series.mask : Replace values where the condition is True.

    Examples
    --------
    >>> s = pd.Series([2, 0, 4, 8, np.nan])

    Boundary values are included by default:

    >>> s.between(1, 4)
    0     True
    1    False
    2     True
    3    False
    4    False
    dtype: bool

    With `inclusive` set to ``"neither"`` boundary values are excluded:

    >>> s.between(1, 4, inclusive="neither")
    0     True
    1    False
    2    False
    3    False
    4    False
    dtype: bool

    `left` and `right` can be any scalar value:

    >>> s = pd.Series(['Alice', 'Bob', 'Carol', 'Eve'])
    >>> s.between('Anna', 'Daniel')
    0    False
    1     True
    2     True
    3    False
    dtype: bool
    """
    from pandas import Series

    len_args = len(args)
    if not len_args:
        raise ValueError(
            "Kindly provide at least one boolean condition, "
            "with a corresponding replacement."
        )
    if len_args % 2:
        raise ValueError(
            "The number of boolean conditions should be equal "
            "to the number of replacements. "
            "However, the total number of conditions and replacements "
            f"is {len(args)}, which is an odd number."
        )

    counter = len_args // 2 - 1
    counter = range(counter, -1, -1)
    conditions = []
    for num, condition in zip(counter, args[-2::-2]):
        if not hasattr(condition, "shape"):
            condition = np.asanyarray(condition)
        if condition.ndim > 1:
            raise ValueError(f"condition{num} is not a one dimensional array.")
        if not is_bool_dtype(condition):
            raise ValueError(f"condition{num} is not a boolean array.")
        conditions.append(condition)
    bool_length = {condition.size for condition in conditions}
    if len(bool_length) > 1:
        raise ValueError("All boolean conditions should have the same length.")
    bool_length = conditions[0].size
    if default is not None:
        if is_scalar(default):
            default = construct_1d_arraylike_from_scalar(
                default, length=bool_length, dtype=None
            )
        default = convert_to_list_like(default)
        if not hasattr(default, "shape"):
            default = pd_array(default, copy=False)
        if default.ndim > 1:
            raise ValueError(
                "The provided default argument should "
                "either be a scalar or a 1-D array."
            )
        if default.size != conditions[0].size:
            raise ValueError(
                "The length of the default argument should "
                "be the same as the length of any "
                "of the boolean conditions."
            )
    replacements = []
    for num, replacement in zip(counter, args[-1::-2]):
        if is_scalar(replacement):
            replacement = construct_1d_arraylike_from_scalar(
                replacement, length=bool_length, dtype=None
            )
        replacement = convert_to_list_like(replacement)
        if not hasattr(replacement, "shape"):
            replacement = pd_array(replacement, copy=False)
        if replacement.ndim > 1:
            raise ValueError(f"replacement{num} should be a 1-D array.")
        if replacement.size != bool_length:
            raise ValueError(
                f"The size of replacement{num} array"
                f"does not match the size of condition{num} array."
            )
        replacements.append(replacement)
    common_dtype = [arr.dtype for arr in replacements]
    if default is not None:
        common_dtype.append(default.dtype)
    if len(set(common_dtype)) > 1:
        common_dtype = find_common_type(common_dtype)
        replacements = [
            arr.astype(common_dtype, copy=False)
            if isinstance(arr, ABCSeries)
            else pd_array(arr, dtype=common_dtype, copy=False)
            for arr in replacements
        ]
        if default is not None:
            if isinstance(default, ABCSeries):
                default = default.astype(common_dtype, copy=False)
            else:
                default = pd_array(default, dtype=common_dtype, copy=False)
    else:
        common_dtype = common_dtype[0]

    if default is not None:
        if not isinstance(default, ABCSeries):
            ser = Series(default)
        else:
            ser = default[:]
    else:
        ser = construct_1d_arraylike_from_scalar(
            na_value_for_dtype(common_dtype, compat=False),
            length=bool_length,
            dtype=common_dtype,
        )
        ser = Series(ser)

    for position, condition, replacement in zip(counter, conditions, replacements):
        try:
            ser = ser.mask(
                condition, other=replacement, axis=0, inplace=False, level=level
            )
        except Exception as error:
            raise ValueError(
                f"condition{position} and replacement{position} failed to evaluate. "
                f"Original error message: {error}"
            ) from error
    return ser
