""" generic datetimelike tests """
import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm

from .common import Base


class DatetimeLike(Base):
    def test_argmax_axis_invalid(self):
        # GH#23081
        msg = r"`axis` must be fewer than the number of dimensions \(1\)"
        rng = self.create_index()
        with pytest.raises(ValueError, match=msg):
            rng.argmax(axis=1)
        with pytest.raises(ValueError, match=msg):
            rng.argmin(axis=2)
        with pytest.raises(ValueError, match=msg):
            rng.min(axis=-2)
        with pytest.raises(ValueError, match=msg):
            rng.max(axis=-3)

    def test_can_hold_identifiers(self):
        idx = self.create_index()
        key = idx[0]
        assert idx._can_hold_identifiers_and_holds_name(key) is False

    def test_shift_identity(self):

        idx = self.create_index()
        tm.assert_index_equal(idx, idx.shift(0))

    def test_str(self):

        # test the string repr
        idx = self.create_index()
        idx.name = "foo"
        assert not (f"length={len(idx)}" in str(idx))
        assert "'foo'" in str(idx)
        assert type(idx).__name__ in str(idx)

        if hasattr(idx, "tz"):
            if idx.tz is not None:
                assert idx.tz in str(idx)
        if hasattr(idx, "freq"):
            assert f"freq='{idx.freqstr}'" in str(idx)

    def test_view(self):
        i = self.create_index()

        i_view = i.view("i8")
        result = self._holder(i)
        tm.assert_index_equal(result, i)

        i_view = i.view(self._holder)
        result = self._holder(i)
        tm.assert_index_equal(result, i_view)

    def test_map_callable(self):
        index = self.create_index()
        expected = index + index.freq
        result = index.map(lambda x: x + x.freq)
        tm.assert_index_equal(result, expected)

        # map to NaT
        result = index.map(lambda x: pd.NaT if x == index[0] else x)
        expected = pd.Index([pd.NaT] + index[1:].tolist())
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize(
        "mapper",
        [
            lambda values, index: {i: e for e, i in zip(values, index)},
            lambda values, index: pd.Series(values, index, dtype=object),
        ],
    )
    def test_map_dictlike(self, mapper):
        index = self.create_index()
        expected = index + index.freq

        # don't compare the freqs
        if isinstance(expected, (pd.DatetimeIndex, pd.TimedeltaIndex)):
            expected = expected._with_freq(None)

        result = index.map(mapper(expected, index))
        tm.assert_index_equal(result, expected)

        expected = pd.Index([pd.NaT] + index[1:].tolist())
        result = index.map(mapper(expected, index))
        tm.assert_index_equal(result, expected)

        # empty map; these map to np.nan because we cannot know
        # to re-infer things
        expected = pd.Index([np.nan] * len(index))
        result = index.map(mapper([], []))
        tm.assert_index_equal(result, expected)

    def test_getitem_preserves_freq(self):
        index = self.create_index()
        assert index.freq is not None

        result = index[:]
        assert result.freq == index.freq
