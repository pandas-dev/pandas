from scripts.no_bool_in_generic import check_for_bool_in_generic

BAD_FILE = "def foo(a: bool) -> bool:\n    return bool(0)"
GOOD_FILE = "def foo(a: bool_t) -> bool_t:\n    return bool(0)"


def test_bad_file_with_replace():
    content = BAD_FILE
    _, result = check_for_bool_in_generic(content)
    expected = GOOD_FILE
    assert result == expected


def test_good_file_with_replace():
    content = GOOD_FILE
    _, result = check_for_bool_in_generic(content)
    expected = content
    assert result == expected
