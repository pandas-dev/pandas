import pytest

from pandas import DataFrame, Series
import pandas._testing as tm


@pytest.mark.parametrize("n, frac", [(2, None), (None, 0.2)])
def test_groupby_sample_balanced_groups_shape(n, frac):
    df = DataFrame({"a": [1] * 10 + [2] * 10, "b": [1] * 20})

    result = df.groupby("a").sample(n=n, frac=frac)
    expected = DataFrame({"a": [1] * 2 + [2] * 2, "b": [1] * 4}, index=result.index)
    tm.assert_frame_equal(result, expected)

    result = df.groupby("a")["b"].sample(n=n, frac=frac)
    expected = Series([1] * 4, name="b", index=result.index)
    tm.assert_series_equal(result, expected)


def test_groupby_sample_unbalanced_groups_shape():
    df = DataFrame({"a": [1] * 10 + [2] * 20, "b": [1] * 30})

    result = df.groupby("a").sample(n=5)
    expected = DataFrame({"a": [1] * 5 + [2] * 5, "b": [1] * 10}, index=result.index)
    tm.assert_frame_equal(result, expected)

    result = df.groupby("a")["b"].sample(n=5)
    expected = Series([1] * 10, name="b", index=result.index)
    tm.assert_series_equal(result, expected)


def test_groupby_sample_n_and_frac_raises():
    df = DataFrame({"a": [1] * 10 + [2] * 10, "b": [1] * 20})
    msg = "Please enter a value for frac or n but not both"

    with pytest.raises(ValueError, match=msg):
        df.groupby("a").sample(n=1, frac=1.0)

    with pytest.raises(ValueError, match=msg):
        df.groupby("a")["b"].sample(n=1, frac=1.0)


def test_groupby_sample_frac_gt_one_without_replacement_raises():
    df = DataFrame({"a": [1] * 10 + [2] * 10, "b": [1] * 20})
    msg = "replace must be set to True when frac > 1"

    with pytest.raises(ValueError, match=msg):
        df.groupby("a").sample(frac=1.5, replace=False)

    with pytest.raises(ValueError, match=msg):
        df.groupby("a")["b"].sample(frac=1.5, replace=False)


@pytest.mark.parametrize("n", [-1, 1.5])
def test_groupby_sample_invalid_n(n):
    df = DataFrame({"a": [1] * 10 + [2] * 10, "b": [1] * 20})
    msg = "Only non-negative integers accepted as n values"

    with pytest.raises(ValueError, match=msg):
        df.groupby("a").sample(n=n)

    with pytest.raises(ValueError, match=msg):
        df.groupby("a")["b"].sample(n=n)


def test_groupby_sample_oversample():
    df = DataFrame({"a": [1] * 10 + [2] * 10, "b": [1] * 20})

    result = df.groupby("a").sample(frac=2.0, replace=True)
    expected = DataFrame({"a": [1] * 20 + [2] * 20, "b": [1] * 40}, index=result.index)
    tm.assert_frame_equal(result, expected)

    result = df.groupby("a")["b"].sample(frac=2.0, replace=True)
    expected = Series([1] * 40, name="b", index=result.index)
    tm.assert_series_equal(result, expected)


def test_groupby_sample_without_n_or_frac():
    df = DataFrame({"a": [1] * 10 + [2] * 10, "b": [1] * 20})

    result = df.groupby("a").sample(n=None, frac=None)
    expected = DataFrame({"a": [1, 2], "b": [1, 1]}, index=result.index)
    tm.assert_frame_equal(result, expected)

    result = df.groupby("a")["b"].sample(n=None, frac=None)
    expected = Series([1, 1], name="b", index=result.index)
    tm.assert_series_equal(result, expected)
