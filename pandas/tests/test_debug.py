from pandas._libs.tslibs.timedeltas import (
    debug_2,
    debug_divmod_bug,
)

import pandas as pd


def non_buggy_divmod(delta):
    conv = 1000 * 1
    n = (
        delta.days * 24 * 3600 * 1_000_000
        + delta.seconds * 1_000_000
        + delta.microseconds
    )
    div, mod = divmod(n, conv)
    return div


def test_debug_divmod_bug():
    td = pd.Timedelta(minutes=-7).to_pytimedelta()
    result = debug_divmod_bug(td)
    expected = non_buggy_divmod(td)
    assert result == expected


def test_debug_2():
    assert False, debug_2()
