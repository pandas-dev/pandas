"""
Module for correlation related implementation
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from pandas.core.dtypes.dtypes import CategoricalDtype

if TYPE_CHECKING:
    from pandas import DataFrame


def transform_ord_cat_cols_to_coded_cols(df: DataFrame) -> DataFrame:
    """
    Replace ordered categoricals with their codes, making a shallow copy if necessary.
    """

    result = df
    made_copy = False
    for idx, dtype in enumerate(df.dtypes):
        if not isinstance(dtype, CategoricalDtype) or not dtype.ordered:
            continue
        col = result._ixs(idx, axis=1)
        if not made_copy:
            made_copy = True
            result = result.copy(deep=False)
        result._iset_item(idx, col.cat.codes.replace(-1, np.nan))
    return result
