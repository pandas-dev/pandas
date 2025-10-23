import datetime

import numpy as np
import pytest

import pandas as pd


def test_multiindex_date_npdatetime_mismatch_raises():
    idx = pd.MultiIndex.from_product(
        [[datetime.date(2023, 11, 1)], ["A"], ["C"]], names=["date", "t1", "t2"]
    )
    df = pd.DataFrame({"vals": [1]}, index=idx)
    with pytest.raises(TypeError, match="Type mismatch"):
        df.loc[(np.datetime64("2023-11-01"), "A", "C")]
