"""
Methods used by Block.replace and related methods.
"""
import operator
import re
from typing import Pattern, Union

import numpy as np

from pandas._typing import ArrayLike, Scalar

from pandas.core.dtypes.common import (
    is_datetimelike_v_numeric,
    is_numeric_v_string_like,
    is_scalar,
)


def compare_or_regex_search(
    a: ArrayLike, b: Union[Scalar, Pattern], regex: bool, mask: ArrayLike
) -> Union[ArrayLike, bool]:
    """
    Compare two array_like inputs of the same shape or two scalar values

    Calls operator.eq or re.search, depending on regex argument. If regex is
    True, perform an element-wise regex matching.

    Parameters
    ----------
    a : array_like
    b : scalar or regex pattern
    regex : bool
    mask : array_like

    Returns
    -------
    mask : array_like of bool
    """

    def _check_comparison_types(
        result: Union[ArrayLike, bool], a: ArrayLike, b: Union[Scalar, Pattern]
    ):
        """
        Raises an error if the two arrays (a,b) cannot be compared.
        Otherwise, returns the comparison result as expected.
        """
        if is_scalar(result) and isinstance(a, np.ndarray):
            type_names = [type(a).__name__, type(b).__name__]

            if isinstance(a, np.ndarray):
                type_names[0] = f"ndarray(dtype={a.dtype})"

            raise TypeError(
                f"Cannot compare types {repr(type_names[0])} and {repr(type_names[1])}"
            )

    if not regex:
        op = lambda x: operator.eq(x, b)
    else:
        op = np.vectorize(
            lambda x: bool(re.search(b, x))
            if isinstance(x, str) and isinstance(b, (str, Pattern))
            else False
        )

    # GH#32621 use mask to avoid comparing to NAs
    if isinstance(a, np.ndarray):
        a = a[mask]

    if is_numeric_v_string_like(a, b):
        # GH#29553 avoid deprecation warnings from numpy
        return np.zeros(a.shape, dtype=bool)

    elif is_datetimelike_v_numeric(a, b):
        # GH#29553 avoid deprecation warnings from numpy
        _check_comparison_types(False, a, b)
        return False

    result = op(a)

    if isinstance(result, np.ndarray) and mask is not None:
        # The shape of the mask can differ to that of the result
        # since we may compare only a subset of a's or b's elements
        tmp = np.zeros(mask.shape, dtype=np.bool_)
        tmp[mask] = result
        result = tmp

    _check_comparison_types(result, a, b)
    return result
