import pytest

from pandas import Series


@pytest.mark.parametrize("dtype", [str])
def test_string_array(dtype):
    test_series = Series(["asdf", "as"], dtype=dtype)
    regex = r"((as)|(as))"
    regex2 = r"(as)|(as)"
    assert list(test_series.str.fullmatch(regex)) == [False, True]
    assert list(test_series.str.fullmatch(regex2)) == [False, True]


@pytest.mark.parametrize(
    "data, pattern, expected",
    [
        (["cat", "duck", "dove"], r"d.+", [False, True, True]),
    ],
)
def test_string_match(data, pattern, expected):
    ser = Series(data)
    assert list(ser.str.fullmatch(pattern)) == expected


@pytest.mark.parametrize("dtype", [str])
@pytest.mark.parametrize(
    "pattern, expected",
    [
        (r"(foo)|((as)(df)?)", [True, True, True]),
        ("foo|as", [False, True, True]),
    ],
)
def test_string_alternation_patterns(dtype, pattern, expected):
    ser = Series(["asdf", "foo", "as"], dtype=dtype)
    assert list(ser.str.fullmatch(pattern)) == expected
