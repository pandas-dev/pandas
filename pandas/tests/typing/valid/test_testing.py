"""
Copyright (c) Virtus Lab sp. z o.o. (Ltd.)

Distributed under the terms of the MIT license.

The full license is in the STUBS_LICENSE file, distributed with this software.
"""
# pyright: reportGeneralTypeIssues = true

import pandas as pd
import pandas._testing as tm


def test_types_assert_series_equal() -> None:
    s1 = pd.Series([0, 1, 1, 0])
    s2 = pd.Series([0, 1, 1, 0])
    tm.assert_series_equal(left=s1, right=s2)
    tm.assert_series_equal(
        s1,
        s2,
        check_freq=False,
        check_categorical=True,
        check_flags=True,
        check_datetimelike_compat=True,
    )
    tm.assert_series_equal(s1, s2, check_dtype=True, check_names=True)
