import pandas as pd
import pandas._testing as tm


def test_astype_str_from_bytes():
    # https://github.com/pandas-dev/pandas/issues/38607
    # GH#49658 pre-2.0 Index called .values.astype(str) here, which effectively
    #  did a .decode() on the bytes object.  In 2.0 we go through
    #  ensure_string_array which does f"{val}"
    idx = pd.Index(["あ", b"a"], dtype="object")
    result = idx.astype(str)
    expected = pd.Index(["あ", "a"], dtype="str")
    tm.assert_index_equal(result, expected)

    # while we're here, check that Series.astype behaves the same
    result = pd.Series(idx).astype(str)
    expected = pd.Series(expected, dtype="str")
    tm.assert_series_equal(result, expected)
