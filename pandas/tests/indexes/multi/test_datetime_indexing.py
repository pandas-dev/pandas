import datetime as dt

import numpy as np
import pytest

import pandas as pd


def test_multiindex_date_npdatetime_mismatch_raises():
    dates = [dt.date(2023, 11, 1), dt.date(2023, 11, 1), dt.date(2023, 11, 2)]
    t1 = ["A", "B", "C"]
    t2 = ["C", "D", "E"]
    vals = [10, 20, 30]

    df = pd.DataFrame(
        data=np.array([dates, t1, t2, vals]).T, columns=["dates", "t1", "t2", "vals"]
    )
    df.set_index(["dates", "t1", "t2"], inplace=True)

    # Exact type match
    result = df.loc[(dt.date(2023, 11, 1), "A", "C")]
    assert result["vals"] == 10

    # TypeError
    with pytest.raises(KeyError):
        df.loc[(np.datetime64("2023-11-01"), "A", "C")]
