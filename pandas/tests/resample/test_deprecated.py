import pandas as pd
from pandas import Series
import pandas._testing as tm


def test_interpolate_posargs_deprecation():
    # GH 41485
    idx = pd.to_datetime(["1992-08-27 07:46:48", "1992-08-27 07:46:59"])
    s = Series([1, 4], index=idx)

    msg = (
        r"In a future version of pandas all arguments of Resampler\.interpolate "
        r"except for the argument 'method' will be keyword-only"
    )

    with tm.assert_produces_warning(FutureWarning, match=msg):
        result = s.resample("3s").interpolate("linear", 0)

    idx = pd.to_datetime(
        [
            "1992-08-27 07:46:48",
            "1992-08-27 07:46:51",
            "1992-08-27 07:46:54",
            "1992-08-27 07:46:57",
        ]
    )
    expected = Series([1.0, 1.0, 1.0, 1.0], index=idx)

    expected.index._data.freq = "3s"
    tm.assert_series_equal(result, expected)
