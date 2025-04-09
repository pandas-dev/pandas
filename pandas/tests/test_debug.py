from datetime import timedelta

from pandas._libs.tslibs.timedeltas import debug_divmod_bug


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
    td = timedelta(minutes=-7)
    result = debug_divmod_bug(td)
    expected = non_buggy_divmod(td)
    assert result == expected
