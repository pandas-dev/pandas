from pytest import (
    fixture,
    mark,
)

from pandas.util._exceptions import find_stack_level


@mark.parametrize("above,below", [(0, 1)])
def test_find_stack_level(above, below):
    """
    Note that this test would not be expected to pass on CPython implementations which
    don't support getting the frame with currentframe (which would always return None).
    """

    def test_stack_level() -> int:
        return find_stack_level()

    top_lvl = find_stack_level()
    assert top_lvl == above, f"Expected stack level {above} but got {top_lvl}"
    sub_lvl = nested_call()
    assert sub_lvl == nest_lvl, f"Expected stack level {below} but got {sub_lvl}"
