"""
Arithmetic operations for PandasObjects

This is not a public API.
"""
from __future__ import annotations

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
from pandas.core.ops.invalid import invalid_comparison
from pandas.core.ops.mask_ops import (
    kleene_and,
    kleene_or,
    kleene_xor,
)
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


__all__ = [
    "ARITHMETIC_BINOPS",
    "arithmetic_op",
    "COMPARISON_BINOPS",
    "comparison_op",
    "comp_method_OBJECT_ARRAY",
    "fill_binop",
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
    "get_op_result_name",
    "maybe_prepare_scalar_for_op",
    "get_array_op",
]
