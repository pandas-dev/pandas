import numpy as np
import pytest

from pandas import Categorical, Series
import pandas._testing as tm


@pytest.mark.parametrize(
    "keep, expected",
    [
        ("first", Series([False, False, False, False, True, True, False])),
        ("last", Series([False, True, True, False, False, False, False])),
        (False, Series([False, True, True, False, True, True, False])),
    ],
)
def test_drop_duplicates(any_numpy_dtype, keep, expected):
    tc = Series([1, 0, 3, 5, 3, 0, 4], dtype=np.dtype(any_numpy_dtype))

    if tc.dtype == "bool":
        pytest.skip("tested separately in test_drop_duplicates_bool")

    tm.assert_series_equal(tc.duplicated(keep=keep), expected)
    tm.assert_series_equal(tc.drop_duplicates(keep=keep), tc[~expected])
    sc = tc.copy()
    sc.drop_duplicates(keep=keep, inplace=True)
    tm.assert_series_equal(sc, tc[~expected])


@pytest.mark.parametrize(
    "keep, expected",
    [
        ("first", Series([False, False, True, True])),
        ("last", Series([True, True, False, False])),
        (False, Series([True, True, True, True])),
    ],
)
def test_drop_duplicates_bool(keep, expected):
    tc = Series([True, False, True, False])

    tm.assert_series_equal(tc.duplicated(keep=keep), expected)
    tm.assert_series_equal(tc.drop_duplicates(keep=keep), tc[~expected])
    sc = tc.copy()
    sc.drop_duplicates(keep=keep, inplace=True)
    tm.assert_series_equal(sc, tc[~expected])


@pytest.mark.parametrize("values", [[], list(range(5))])
def test_drop_duplicates_no_duplicates(any_numpy_dtype, keep, values):
    tc = Series(values, dtype=np.dtype(any_numpy_dtype))
    expected = Series([False] * len(tc), dtype="bool")

    if tc.dtype == "bool":
        # 0 -> False and 1-> True
        # any other value would be duplicated
        tc = tc[:2]
        expected = expected[:2]

    tm.assert_series_equal(tc.duplicated(keep=keep), expected)

    result_dropped = tc.drop_duplicates(keep=keep)
    tm.assert_series_equal(result_dropped, tc)

    # validate shallow copy
    assert result_dropped is not tc


class TestSeriesDropDuplicates:
    @pytest.mark.parametrize(
        "dtype",
        ["int_", "uint", "float_", "unicode_", "timedelta64[h]", "datetime64[D]"],
    )
    def test_drop_duplicates_categorical_non_bool(self, dtype, ordered_fixture):
        cat_array = np.array([1, 2, 3, 4, 5], dtype=np.dtype(dtype))

        # Test case 1
        input1 = np.array([1, 2, 3, 3], dtype=np.dtype(dtype))
        tc1 = Series(Categorical(input1, categories=cat_array, ordered=ordered_fixture))
        if dtype == "datetime64[D]":
            # pre-empty flaky xfail, tc1 values are seemingly-random
            if not (np.array(tc1) == input1).all():
                pytest.xfail(reason="GH#7996")

        expected = Series([False, False, False, True])
        tm.assert_series_equal(tc1.duplicated(), expected)
        tm.assert_series_equal(tc1.drop_duplicates(), tc1[~expected])
        sc = tc1.copy()
        sc.drop_duplicates(inplace=True)
        tm.assert_series_equal(sc, tc1[~expected])

        expected = Series([False, False, True, False])
        tm.assert_series_equal(tc1.duplicated(keep="last"), expected)
        tm.assert_series_equal(tc1.drop_duplicates(keep="last"), tc1[~expected])
        sc = tc1.copy()
        sc.drop_duplicates(keep="last", inplace=True)
        tm.assert_series_equal(sc, tc1[~expected])

        expected = Series([False, False, True, True])
        tm.assert_series_equal(tc1.duplicated(keep=False), expected)
        tm.assert_series_equal(tc1.drop_duplicates(keep=False), tc1[~expected])
        sc = tc1.copy()
        sc.drop_duplicates(keep=False, inplace=True)
        tm.assert_series_equal(sc, tc1[~expected])

        # Test case 2
        input2 = np.array([1, 2, 3, 5, 3, 2, 4], dtype=np.dtype(dtype))
        tc2 = Series(Categorical(input2, categories=cat_array, ordered=ordered_fixture))
        if dtype == "datetime64[D]":
            # pre-empty flaky xfail, tc2 values are seemingly-random
            if not (np.array(tc2) == input2).all():
                pytest.xfail(reason="GH#7996")

        expected = Series([False, False, False, False, True, True, False])
        tm.assert_series_equal(tc2.duplicated(), expected)
        tm.assert_series_equal(tc2.drop_duplicates(), tc2[~expected])
        sc = tc2.copy()
        sc.drop_duplicates(inplace=True)
        tm.assert_series_equal(sc, tc2[~expected])

        expected = Series([False, True, True, False, False, False, False])
        tm.assert_series_equal(tc2.duplicated(keep="last"), expected)
        tm.assert_series_equal(tc2.drop_duplicates(keep="last"), tc2[~expected])
        sc = tc2.copy()
        sc.drop_duplicates(keep="last", inplace=True)
        tm.assert_series_equal(sc, tc2[~expected])

        expected = Series([False, True, True, False, True, True, False])
        tm.assert_series_equal(tc2.duplicated(keep=False), expected)
        tm.assert_series_equal(tc2.drop_duplicates(keep=False), tc2[~expected])
        sc = tc2.copy()
        sc.drop_duplicates(keep=False, inplace=True)
        tm.assert_series_equal(sc, tc2[~expected])

    def test_drop_duplicates_categorical_bool(self, ordered_fixture):
        tc = Series(
            Categorical(
                [True, False, True, False],
                categories=[True, False],
                ordered=ordered_fixture,
            )
        )

        expected = Series([False, False, True, True])
        tm.assert_series_equal(tc.duplicated(), expected)
        tm.assert_series_equal(tc.drop_duplicates(), tc[~expected])
        sc = tc.copy()
        sc.drop_duplicates(inplace=True)
        tm.assert_series_equal(sc, tc[~expected])

        expected = Series([True, True, False, False])
        tm.assert_series_equal(tc.duplicated(keep="last"), expected)
        tm.assert_series_equal(tc.drop_duplicates(keep="last"), tc[~expected])
        sc = tc.copy()
        sc.drop_duplicates(keep="last", inplace=True)
        tm.assert_series_equal(sc, tc[~expected])

        expected = Series([True, True, True, True])
        tm.assert_series_equal(tc.duplicated(keep=False), expected)
        tm.assert_series_equal(tc.drop_duplicates(keep=False), tc[~expected])
        sc = tc.copy()
        sc.drop_duplicates(keep=False, inplace=True)
        tm.assert_series_equal(sc, tc[~expected])
