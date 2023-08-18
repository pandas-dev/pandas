from pytest import mark

from pandas.util._exceptions import find_stack_level


@mark.parametrize("expected", [1])
def test_find_stack_level(expected):
    """
    Note that this test would not be expected to pass on CPython implementations which
    don't support getting the frame with currentframe (which would always return None).
    """

    top_lvl = find_stack_level()
    assert top_lvl == expected, f"Expected stack level {expected} but got {top_lvl}"
