import numpy as np

from pandas import DataFrame


def transform_ord_cat_cols_to_coded_cols(df: DataFrame) -> DataFrame:
    """
    any ordered categorical columns are transformed to the respective
    categorical codes while other columns remain untouched
    """

    result = df
    made_copy = False
    for idx, dtype in enumerate(df.dtypes):
        if not dtype == "category" or not dtype.ordered:
            continue
        col = result._ixs(idx, axis=1)
        if not made_copy:
            made_copy = True
            result = result.copy(deep=False)
        result._iset_item(idx, col.cat.codes.replace(-1, np.nan))
    return result
