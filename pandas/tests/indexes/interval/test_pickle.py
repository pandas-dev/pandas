from pandas import IntervalIndex
import pandas._testing as tm


class TestPickle:
    def test_pickle_round_trip_closed(self, closed, temp_file):
        # https://github.com/pandas-dev/pandas/issues/35658
        idx = IntervalIndex.from_tuples([(1, 2), (2, 3)], closed=closed)
        result = tm.round_trip_pickle(idx, temp_file)
        tm.assert_index_equal(result, idx)
