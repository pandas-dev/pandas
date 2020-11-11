import numpy as np
import pytest

from pandas.compat import is_platform_windows

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
    return_value = sc.drop_duplicates(keep=keep, inplace=True)
    assert return_value is None
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
    return_value = sc.drop_duplicates(keep=keep, inplace=True)
    tm.assert_series_equal(sc, tc[~expected])
    assert return_value is None


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
    @pytest.fixture(
        params=["int_", "uint", "float_", "unicode_", "timedelta64[h]", "datetime64[D]"]
    )
    def dtype(self, request):
        return request.param

    @pytest.fixture
    def tc1(self, dtype, ordered):
        # Test case 1
        cat_array = np.array([1, 2, 3, 4, 5], dtype=np.dtype(dtype))

        input1 = np.array([1, 2, 3, 3], dtype=np.dtype(dtype))

        tc1 = Series(Categorical(input1, categories=cat_array, ordered=ordered))
        return tc1

    def _maybe_xfail_tc(self, tc, request):
        if tc.cat.categories.dtype.kind == "M":
            if len(tc) == 4:
                # This is tc1
                input_arr = np.array([1, 2, 3, 3], dtype=np.dtype("datetime64[D]"))
            else:
                # This is tc2
                input_arr = np.array(
                    [1, 2, 3, 5, 3, 2, 4], dtype=np.dtype("datetime64[D]")
                )

            if not (np.array(tc) == input_arr).all() and is_platform_windows():
                mark = pytest.mark.xfail(
                    reason="GH#7996 tc1/tc2 values are seemingly-random",
                    raises=AssertionError,
                )
                request.node.add_marker(mark)

    def _check_drop_duplicates_vs_duplicated(self, tc, keep, expected):
        result = tc.duplicated(keep=keep)
        tm.assert_series_equal(result, expected)

        result = tc.drop_duplicates(keep=keep)
        tm.assert_series_equal(result, tc[~expected])

        sc = tc.copy()
        return_value = sc.drop_duplicates(keep=keep, inplace=True)
        assert return_value is None
        tm.assert_series_equal(sc, tc[~expected])

    def test_drop_duplicates_categorical_non_bool(self, tc1, request):
        tc = tc1
        self._maybe_xfail_tc(tc, request)

        keep = "first"
        expected = Series([False, False, False, True])
        self._check_drop_duplicates_vs_duplicated(tc, keep, expected)

    def test_drop_duplicates_categorical_non_bool_keeplast(self, tc1, request):
        tc = tc1
        self._maybe_xfail_tc(tc, request)

        keep = "last"
        expected = Series([False, False, True, False])
        self._check_drop_duplicates_vs_duplicated(tc, keep, expected)

    def test_drop_duplicates_categorical_non_bool_keepfalse(self, tc1, request):
        tc = tc1
        self._maybe_xfail_tc(tc, request)

        keep = False
        expected = Series([False, False, True, True])
        self._check_drop_duplicates_vs_duplicated(tc, keep, expected)

    @pytest.fixture
    def tc2(self, dtype, ordered):
        # Test case 2; TODO: better name
        cat_array = np.array([1, 2, 3, 4, 5], dtype=np.dtype(dtype))

        input2 = np.array([1, 2, 3, 5, 3, 2, 4], dtype=np.dtype(dtype))
        tc2 = Series(Categorical(input2, categories=cat_array, ordered=ordered))
        return tc2

    def test_drop_duplicates_categorical_non_bool2(self, tc2, request):
        # Test case 2; TODO: better name
        tc = tc2
        self._maybe_xfail_tc(tc, request)

        keep = "first"
        expected = Series([False, False, False, False, True, True, False])
        self._check_drop_duplicates_vs_duplicated(tc, keep, expected)

    def test_drop_duplicates_categorical_non_bool2_keeplast(self, tc2, request):
        tc = tc2
        self._maybe_xfail_tc(tc, request)

        keep = "last"
        expected = Series([False, True, True, False, False, False, False])
        self._check_drop_duplicates_vs_duplicated(tc, keep, expected)

    def test_drop_duplicates_categorical_non_bool2_keepfalse(self, tc2, request):
        tc = tc2
        self._maybe_xfail_tc(tc, request)

        keep = False
        expected = Series([False, True, True, False, True, True, False])
        self._check_drop_duplicates_vs_duplicated(tc, keep, expected)

    def test_drop_duplicates_categorical_bool(self, ordered):
        tc = Series(
            Categorical(
                [True, False, True, False], categories=[True, False], ordered=ordered
            )
        )

        expected = Series([False, False, True, True])
        self._check_drop_duplicates_vs_duplicated(tc, "first", expected)

        expected = Series([True, True, False, False])
        self._check_drop_duplicates_vs_duplicated(tc, "last", expected)

        expected = Series([True, True, True, True])
        self._check_drop_duplicates_vs_duplicated(tc, False, expected)
