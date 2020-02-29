from pandas import MultiIndex, Timestamp, date_range


class TestGetLevelValues:
    def test_get_level_values_box_datetime64(self):

        dates = date_range("1/1/2000", periods=4)
        levels = [dates, [0, 1]]
        codes = [[0, 0, 1, 1, 2, 2, 3, 3], [0, 1, 0, 1, 0, 1, 0, 1]]

        index = MultiIndex(levels=levels, codes=codes)

        assert isinstance(index.get_level_values(0)[0], Timestamp)
