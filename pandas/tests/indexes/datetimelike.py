""" generic datetimelike tests """
import pandas as pd
from .common import Base
import pandas.util.testing as tm


class DatetimeLike(Base):

    def test_shift_identity(self):

        idx = self.create_index()
        tm.assert_index_equal(idx, idx.shift(0))

    def test_str(self):

        # test the string repr
        idx = self.create_index()
        idx.name = 'foo'
        assert not "length=%s" % len(idx) in str(idx)
        assert "'foo'" in str(idx)
        assert idx.__class__.__name__ in str(idx)

        if hasattr(idx, 'tz'):
            if idx.tz is not None:
                assert idx.tz in str(idx)
        if hasattr(idx, 'freq'):
            assert "freq='%s'" % idx.freqstr in str(idx)

    def test_view(self, indices):
        super(DatetimeLike, self).test_view(indices)

        i = self.create_index()

        i_view = i.view('i8')
        result = self._holder(i)
        tm.assert_index_equal(result, i)

        i_view = i.view(self._holder)
        result = self._holder(i)
        tm.assert_index_equal(result, i_view)

    def test_map(self):
        expected = self.index + 1
        tm.assert_index_equal(self.index.map(lambda x: x + 1), expected)

        series_map = pd.Series(expected, self.index)
        tm.assert_index_equal(self.index.map(series_map), expected)

        dict_map = {i: e for e, i in zip(expected, self.index)}
        tm.assert_index_equal(self.index.map(dict_map), expected)

        # map to NaT
        result = self.index.map(lambda x: pd.NaT if x == self.index[0] else x)
        expected = pd.Index([pd.NaT] + self.index[1:].tolist())
        tm.assert_index_equal(result, expected)

        series_map = pd.Series(expected, self.index)
        tm.assert_index_equal(self.index.map(series_map), expected)

        dict_map = {i: e for e, i in zip(expected, self.index)}
        tm.assert_index_equal(self.index.map(dict_map), expected)
