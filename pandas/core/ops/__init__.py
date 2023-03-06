"""
Arithmetic operations for PandasObjects

This is not a public API.
"""
from __future__ import annotations

from typing import cast

import numpy as np

from pandas._libs.ops_dispatch import maybe_dispatch_ufunc_to_dunder_op
from pandas._typing import Axis
from pandas.util._decorators import Appender

from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCSeries,
)
from pandas.core.dtypes.missing import isna

from pandas.core.ops.array_ops import (
    arithmetic_op,
    comp_method_OBJECT_ARRAY,
    comparison_op,
    get_array_op,
    logical_op,
    maybe_prepare_scalar_for_op,
)
from pandas.core.ops.common import (
    get_op_result_name,
    unpack_zerodim_and_defer,
)
from pandas.core.ops.docstrings import (
    _flex_comp_doc_FRAME,
    _op_descriptions,
    make_flex_doc,
)
from pandas.core.ops.invalid import invalid_comparison
from pandas.core.ops.mask_ops import (
    kleene_and,
    kleene_or,
    kleene_xor,
)
from pandas.core.ops.methods import add_flex_arithmetic_methods
from pandas.core.roperator import (
    radd,
    rand_,
    rdiv,
    rdivmod,
    rfloordiv,
    rmod,
    rmul,
    ror_,
    rpow,
    rsub,
    rtruediv,
    rxor,
)

# -----------------------------------------------------------------------------
# constants
ARITHMETIC_BINOPS: set[str] = {
    "add",
    "sub",
    "mul",
    "pow",
    "mod",
    "floordiv",
    "truediv",
    "divmod",
    "radd",
    "rsub",
    "rmul",
    "rpow",
    "rmod",
    "rfloordiv",
    "rtruediv",
    "rdivmod",
}


COMPARISON_BINOPS: set[str] = {"eq", "ne", "lt", "gt", "le", "ge"}


# -----------------------------------------------------------------------------
# Masking NA values and fallbacks for operations numpy does not support


def fill_binop(left, right, fill_value):
    """
    If a non-None fill_value is given, replace null entries in left and right
    with this value, but only in positions where _one_ of left/right is null,
    not both.

    Parameters
    ----------
    left : array-like
    right : array-like
    fill_value : object

    Returns
    -------
    left : array-like
    right : array-like

    Notes
    -----
    Makes copies if fill_value is not None and NAs are present.
    """
    if fill_value is not None:
        left_mask = isna(left)
        right_mask = isna(right)

        # one but not both
        mask = left_mask ^ right_mask

        if left_mask.any():
            # Avoid making a copy if we can
            left = left.copy()
            left[left_mask & mask] = fill_value

        if right_mask.any():
            # Avoid making a copy if we can
            right = right.copy()
            right[right_mask & mask] = fill_value

    return left, right


# -----------------------------------------------------------------------------
# Series


def flex_method_SERIES(op):
    name = op.__name__.strip("_")
    doc = make_flex_doc(name, "series")

    @Appender(doc)
    def flex_wrapper(self, other, level=None, fill_value=None, axis: Axis = 0):
        # validate axis
        if axis is not None:
            self._get_axis_number(axis)

        res_name = get_op_result_name(self, other)

        if isinstance(other, ABCSeries):
            return self._binop(other, op, level=level, fill_value=fill_value)
        elif isinstance(other, (np.ndarray, list, tuple)):
            if len(other) != len(self):
                raise ValueError("Lengths must be equal")
            other = self._constructor(other, self.index)
            result = self._binop(other, op, level=level, fill_value=fill_value)
            result.name = res_name
            return result
        else:
            if fill_value is not None:
                self = self.fillna(fill_value)

            return op(self, other)

    flex_wrapper.__name__ = name
    return flex_wrapper


# -----------------------------------------------------------------------------
# DataFrame


def flex_arith_method_FRAME(op):
    op_name = op.__name__.strip("_")

    na_op = get_array_op(op)
    doc = make_flex_doc(op_name, "dataframe")

    @Appender(doc)
    def f(self, other, axis: Axis = "columns", level=None, fill_value=None):
        axis = self._get_axis_number(axis) if axis is not None else 1
        axis = cast(int, axis)

        if self._should_reindex_frame_op(other, op, axis, fill_value, level):
            return self._arith_method_with_reindex(other, op)

        if isinstance(other, ABCSeries) and fill_value is not None:
            # TODO: We could allow this in cases where we end up going
            #  through the DataFrame path
            raise NotImplementedError(f"fill_value {fill_value} not supported.")

        other = maybe_prepare_scalar_for_op(other, self.shape)
        self, other = self._align_for_op(other, axis, flex=True, level=level)

        if isinstance(other, ABCDataFrame):
            # Another DataFrame
            new_data = self._combine_frame(other, na_op, fill_value)

        elif isinstance(other, ABCSeries):
            new_data = self._dispatch_frame_op(other, op, axis=axis)
        else:
            # in this case we always have `np.ndim(other) == 0`
            if fill_value is not None:
                self = self.fillna(fill_value)

            new_data = self._dispatch_frame_op(other, op)

        return self._construct_result(new_data)

    f.__name__ = op_name

    return f


def flex_comp_method_FRAME(op):
    op_name = op.__name__.strip("_")

    doc = _flex_comp_doc_FRAME.format(
        op_name=op_name, desc=_op_descriptions[op_name]["desc"]
    )

    @Appender(doc)
    def f(self, other, axis: Axis = "columns", level=None):
        axis = self._get_axis_number(axis) if axis is not None else 1

        self, other = self._align_for_op(other, axis, flex=True, level=level)

        new_data = self._dispatch_frame_op(other, op, axis=axis)
        return self._construct_result(new_data)

    f.__name__ = op_name

    return f


__all__ = [
    "add_flex_arithmetic_methods",
    "ARITHMETIC_BINOPS",
    "arithmetic_op",
    "COMPARISON_BINOPS",
    "comparison_op",
    "comp_method_OBJECT_ARRAY",
    "fill_binop",
    "flex_arith_method_FRAME",
    "flex_comp_method_FRAME",
    "flex_method_SERIES",
    "invalid_comparison",
    "kleene_and",
    "kleene_or",
    "kleene_xor",
    "logical_op",
    "maybe_dispatch_ufunc_to_dunder_op",
    "radd",
    "rand_",
    "rdiv",
    "rdivmod",
    "rfloordiv",
    "rmod",
    "rmul",
    "ror_",
    "rpow",
    "rsub",
    "rtruediv",
    "rxor",
    "unpack_zerodim_and_defer",
]
