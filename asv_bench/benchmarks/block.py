from datetime import datetime

import numpy as np

import pandas as pd


class Block:
    params = [
        (True, "True"),
        (np.array(True), "np.array(True)"),
    ]

    def setup(self, true_value):
        self.df = pd.DataFrame(
            False,
            columns=np.arange(500).astype(str),
            index=pd.date_range("2010-01-01", "2011-01-01"),
        )

        self.true_value = true_value

    def time_test(self, true_value):
        """Test time for assigning a slice `True` and `np.array(True)`"""
        tmp_df = self.df.copy()

        start = datetime(2010, 5, 1)
        end = datetime(2010, 9, 1)
        tmp_df.loc[start:end, :] = true_value
