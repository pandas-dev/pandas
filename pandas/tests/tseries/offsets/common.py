"""
Assertion helpers for offsets tests
"""


def assert_offset_equal(offset, base, expected):
    actual = offset + base
    actual_swapped = base + offset
    actual_apply = offset.apply(base)
    try:
        assert actual == expected
        assert actual_swapped == expected
        assert actual_apply == expected
    except AssertionError as err:
        raise AssertionError(
            f"\nExpected: {expected}\nActual: {actual}\nFor Offset: {offset})"
            f"\nAt Date: {base}"
        ) from err


def assert_is_on_offset(offset, date, expected):
    actual = offset.is_on_offset(date)
    assert actual == expected, (
        f"\nExpected: {expected}\nActual: {actual}\nFor Offset: {offset})"
        f"\nAt Date: {date}"
    )
