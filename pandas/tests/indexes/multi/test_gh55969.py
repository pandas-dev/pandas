import numpy as np
import pytest

from pandas import (
    DataFrame,
    MultiIndex,
    Timestamp,
)
import pandas._testing as tm


def test_mixed_datetime_types_lookup():

    import datetime as dt

    dates = [dt.date(2023, 11, 1), dt.date(2023, 11, 1), dt.date(2023, 11, 2)]
    t1 = ["A", "B", "C"]
    t2 = ["C", "D", "E"]
    vals = [10, 20, 30]

    df = DataFrame({"dates": dates, "t1": t1, "t2": t2, "vals": vals}).set_index(
        ["dates", "t1", "t2"]
    )

    date_np = np.datetime64("2023-11-01")

    result = df.loc[(date_np, "A")]
    expected_val = 10
    assert len(result) == 1
    assert result["vals"].iloc[0] == expected_val

    msg = "'C'"
    with pytest.raises(KeyError, match=msg):
        df.loc[(date_np, "C")]
