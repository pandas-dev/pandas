import pytest

from pandas import DataFrame, Index, Series
import pandas._testing as tm


@pytest.mark.parametrize(
    "before, after, indices",
    [(1, 2, [2, 1]), (None, 2, [2, 1, 0]), (1, None, [3, 2, 1])],
)
@pytest.mark.parametrize("klass", [Series, DataFrame])
def test_truncate_decreasing_index(before, after, indices, klass):
    # https://github.com/pandas-dev/pandas/issues/33756
    idx = Index([3, 2, 1, 0])
    values = klass(range(len(idx)), index=idx)
    result = values.truncate(before=before, after=after)
    expected = values.loc[indices]
    if klass is Series:
        tm.assert_series_equal(result, expected)
    else:
        tm.assert_frame_equal(result, expected)
