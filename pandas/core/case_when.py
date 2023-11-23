from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from pandas._libs import lib

from pandas.core.dtypes.cast import (
    construct_1d_arraylike_from_scalar,
    find_common_type,
    infer_dtype_from,
)
from pandas.core.dtypes.common import is_scalar
from pandas.core.dtypes.generic import ABCSeries

from pandas.core.construction import array as pd_array

if TYPE_CHECKING:
    from pandas._typing import Series


def case_when(
    *args: tuple[tuple],
    default=lib.no_default,
) -> Series:
    """
    Replace values where the conditions are True.

    Parameters
    ----------
    *args : tuple(s) of array-like, scalar
        Variable argument of tuples of conditions and expected replacements.
        Takes the form:  ``(condition0, replacement0)``,
        ``(condition1, replacement1)``, ... .
        ``condition`` should be a 1-D boolean array.
        When multiple boolean conditions are satisfied,
        the first replacement is used.
        If ``condition`` is a Series, and the equivalent ``replacement``
        is a Series, they must have the same index.
        If there are multiple replacement options,
        and they are Series, they must have the same index.

    default : scalar, array-like, default None
        If provided, it is the replacement value to use
        if all conditions evaluate to False.
        If not specified, entries will be filled with the
        corresponding NULL value.

        .. versionadded:: 2.2.0

    Returns
    -------
    Series

    See Also
    --------
    Series.mask : Replace values where the condition is True.

    Examples
    --------
    >>> df = pd.DataFrame({
    ...     "a": [0,0,1,2],
    ...     "b": [0,3,4,5],
    ...     "c": [6,7,8,9]
    ... })
    >>> df
       a  b  c
    0  0  0  6
    1  0  3  7
    2  1  4  8
    3  2  5  9

    >>> pd.case_when((df.a.gt(0), df.a),   # condition, replacement
    ...              (df.b.gt(0), df.b),   # condition, replacement
    ...              default=df.c)        # optional
    0    6
    1    3
    2    1
    3    2
    Name: c, dtype: int64
    """
    from pandas import Series

    validate_case_when(args=args)

    conditions, replacements = zip(*args)
    common_dtypes = [infer_dtype_from(replacement)[0] for replacement in replacements]

    if default is not lib.no_default:
        arg_dtype, _ = infer_dtype_from(default)
        common_dtypes.append(arg_dtype)
    else:
        default = None
    if len(set(common_dtypes)) > 1:
        common_dtypes = find_common_type(common_dtypes)
        updated_replacements = []
        for condition, replacement in zip(conditions, replacements):
            if is_scalar(replacement):
                replacement = construct_1d_arraylike_from_scalar(
                    value=replacement, length=len(condition), dtype=common_dtypes
                )
            elif isinstance(replacement, ABCSeries):
                replacement = replacement.astype(common_dtypes)
            else:
                replacement = pd_array(replacement, dtype=common_dtypes)
            updated_replacements.append(replacement)
        replacements = updated_replacements
        if (default is not None) and isinstance(default, ABCSeries):
            default = default.astype(common_dtypes)
    else:
        common_dtypes = common_dtypes[0]
    if not isinstance(default, ABCSeries):
        cond_indices = [cond for cond in conditions if isinstance(cond, ABCSeries)]
        replacement_indices = [
            replacement
            for replacement in replacements
            if isinstance(replacement, ABCSeries)
        ]
        cond_length = None
        if replacement_indices:
            for left, right in zip(replacement_indices, replacement_indices[1:]):
                if not left.index.equals(right.index):
                    raise AssertionError(
                        "All replacement objects must have the same index."
                    )
        if cond_indices:
            for left, right in zip(cond_indices, cond_indices[1:]):
                if not left.index.equals(right.index):
                    raise AssertionError(
                        "All condition objects must have the same index."
                    )
            if replacement_indices:
                if not replacement_indices[0].index.equals(cond_indices[0].index):
                    raise AssertionError(
                        "All replacement objects and condition objects "
                        "should have the same index."
                    )
        else:
            conditions = [
                np.asanyarray(cond) if not hasattr(cond, "shape") else cond
                for cond in conditions
            ]
            cond_length = {len(cond) for cond in conditions}
            if len(cond_length) > 1:
                raise ValueError("The boolean conditions should have the same length.")
            cond_length = len(conditions[0])
            if not is_scalar(default):
                if len(default) != cond_length:
                    raise ValueError(
                        "length of `default` does not match the length "
                        "of any of the conditions."
                    )
            if not replacement_indices:
                for num, replacement in enumerate(replacements):
                    if is_scalar(replacement):
                        continue
                    if not hasattr(replacement, "shape"):
                        replacement = np.asanyarray(replacement)
                    if len(replacement) != cond_length:
                        raise ValueError(
                            f"Length of condition{num} does not match "
                            f"the length of replacement{num}; "
                            f"{cond_length} != {len(replacement)}"
                        )
        if cond_indices:
            default_index = cond_indices[0].index
        elif replacement_indices:
            default_index = replacement_indices[0].index
        else:
            default_index = range(cond_length)
        default = Series(default, index=default_index, dtype=common_dtypes)
    counter = reversed(range(len(conditions)))
    for position, condition, replacement in zip(
        counter, conditions[::-1], replacements[::-1]
    ):
        try:
            default = default.mask(
                condition, other=replacement, axis=0, inplace=False, level=None
            )
        except Exception as error:
            raise ValueError(
                f"Failed to apply condition{position} and replacement{position}."
            ) from error
    return default


def validate_case_when(args: tuple) -> tuple:
    """
    Validates the variable arguments for the case_when function.
    """
    if not len(args):
        raise ValueError(
            "provide at least one boolean condition, "
            "with a corresponding replacement."
        )

    for num, entry in enumerate(args):
        if not isinstance(entry, tuple):
            raise TypeError(f"Argument {num} must be a tuple; got {type(entry)}.")
        if len(entry) != 2:
            raise ValueError(
                f"Argument {num} must have length 2; "
                "a condition and replacement; "
                f"got length {len(entry)}."
            )
    return None
