"""
Data structure for 1-dimensional cross-sectional and time series data
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from pandas.core.dtypes.missing import isna

from pandas.core import ops

if TYPE_CHECKING:
    from pandas.core.base import Axis


def _flex_method(
    self,
    other,
    op,
    *,
    level=None,
    fill_value=None,
    axis: Axis = 0,
):
    if axis is not None:
        self._get_axis_number(axis)
    res_name = ops.get_op_result_name(self, other)
    if isinstance(other, type(self)):
        return self._binop(other, op, level=level, fill_value=fill_value)
    elif isinstance(other, (np.ndarray, list, tuple)):
        # Fix for issue #59053: Handle 0-dimensional numpy arrays as scalars
        # Skip length check for scalar numpy arrays (ndim == 0) to avoid TypeError
        if hasattr(other, "ndim") and other.ndim == 0:
            # Treat 0-dimensional arrays as scalars, skip len() check
            pass
        else:
            if len(other) != len(self):
                raise ValueError("Lengths must be equal")
        other = self._constructor(other, self.index, copy=False)
        result = self._binop(other, op, level=level, fill_value=fill_value)
        result._name = res_name
        return result
    else:
        if fill_value is not None:
            if isna(other):
                return op(self, fill_value)
            self = self.fillna(fill_value)
        return op(self, other)
