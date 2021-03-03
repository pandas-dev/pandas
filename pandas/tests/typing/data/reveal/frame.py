# flake8: noqa

import numpy as np

import pandas as pd

empty_df = pd.DataFrame()
empty_ser = pd.Series()
# TODO: np.array resolves to Any
empty_arr = np.array([])
# TODO: Index.__getitem__ resolves to Any
empty_idx: pd.Index = pd.Index([1, 2, 3])[:0]


reveal_type(empty_df.dot(empty_df))  # E: DataFrame
reveal_type(empty_df.dot(empty_ser))  # E: Series
reveal_type(empty_df.dot(empty_arr))  # E: Any
reveal_type(empty_df.dot(empty_idx))  # E: DataFrame

reveal_type(empty_df @ empty_df)  # E: Union[{DataFrame}, {Series}]
reveal_type(empty_df @ empty_ser)  # E: Series
reveal_type(empty_df @ empty_arr)  # E: Any
reveal_type(empty_df @ empty_idx)  # E: Union[{DataFrame}, {Series}]
