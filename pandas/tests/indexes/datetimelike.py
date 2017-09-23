""" generic datetimelike tests """

import numpy as np
import pandas as pd
import pandas.util.testing as tm

from .common import Base


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

    def test_view(self):
        super(DatetimeLike, self).test_view()

        i = self.create_index()

        i_view = i.view('i8')
        result = self._holder(i)
        tm.assert_index_equal(result, i)

        i_view = i.view(self._holder)
        result = self._holder(i)
        tm.assert_index_equal(result, i_view)

    def test_add_timedelta(self):
        # GH 17558
        # Check that tz-aware DatetimeIndex + np.array(dtype="timedelta64")
        # and DatetimeIndex + TimedeltaIndex work as expected
        idx = self.create_index()
        idx.name = "x"
        if isinstance(idx, pd.DatetimeIndex):
            idx = idx.tz_localize("US/Eastern")

        expected = idx + np.timedelta64(1, 'D')
        tm.assert_index_equal(idx, expected - np.timedelta64(1, 'D'))

        deltas = np.array([np.timedelta64(1, 'D')] * len(idx),
                          dtype="timedelta64[ns]")
        results = [idx + deltas,     # add numpy array
                   idx + deltas.astype(dtype="timedelta64[m]"),
                   idx + pd.TimedeltaIndex(deltas, name=idx.name),
                   idx + pd.to_timedelta(deltas[0]),
                   ]
        for actual in results:
            tm.assert_index_equal(actual, expected)

        errmsg = (r"cannot add {cls} and np.ndarray\[float64\]"
                  .format(cls=idx.__class__.__name__))
        with tm.assert_raises_regex(TypeError, errmsg):
            idx + np.array([0.1], dtype=np.float64)
