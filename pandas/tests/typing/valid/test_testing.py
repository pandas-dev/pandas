import pandas as pd
from pandas._testing import assert_series_equal


def test_types_assert_series_equal() -> None:
    s1 = pd.Series([0, 1, 1, 0])
    s2 = pd.Series([0, 1, 1, 0])
    assert_series_equal(left=s1, right=s2)
    assert_series_equal(
        s1,
        s2,
        check_freq=False,
        check_categorical=True,
        check_flags=True,
        check_datetimelike_compat=True,
    )
    assert_series_equal(
        s1, s2, check_dtype=True, check_less_precise=True, check_names=True
    )
