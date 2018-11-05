"""
Tests for TimedeltaArray
"""

import numpy as np
import pytest

from pandas.core.arrays import TimedeltaArrayMixin as TimedeltaArray

import pandas as pd
import pandas.util.testing as tm


# TODO: Many of these tests are mirrored in test_datetimes; see if these
#  can be shared
class TestTimedeltaArrayConstructors(object):
    def test_scalar_raises_type_error(self):
        # GH#23493
        with pytest.raises(TypeError):
            TimedeltaArray(2)

        with pytest.raises(TypeError):
            pd.TimedeltaIndex(pd.Timedelta(days=4))

    def test_from_sequence_requires_1dim(self):
        arr2d = np.arange(10).view('m8[s]').astype(object).reshape(2, 5)
        with pytest.raises(TypeError):
            TimedeltaArray(arr2d)

        with pytest.raises(TypeError):
            pd.TimedeltaIndex(arr2d)

        arr0d = np.array(pd.Timedelta('49 days'))
        with pytest.raises(TypeError):
            TimedeltaArray(arr0d)

        with pytest.raises(TypeError):
            pd.TimedeltaIndex(arr0d)

    def test_init_from_object_dtype(self):
        # GH#23493

        # arbitrary TimedeltaIndex; this should work for any TimedeltaIndex
        #  with non-None freq
        tdi = pd.timedelta_range('3 Days', freq='ms', periods=9)

        # Fails because np.array(tdi, dtype=object) incorrectly returns Longs
        result = TimedeltaArray(np.array(tdi, dtype=object), freq='infer')
        tm.assert_equal(pd.TimedeltaIndex(result), tdi)

        # Fails because `pd.Index(tdi, dtype=object) raises incorrectly
        result = TimedeltaArray(pd.Index(tdi, dtype=object), freq='infer')
        tm.assert_equal(pd.TimedeltaIndex(result), tdi)

    # NB: for now we re-wrap in TimedeltaIndex to use assert_index_equal
    #  once assert_timedelta_array_equal is in place, this will be changed
    def test_init_only_freq_infer(self):
        # GH#23493
        # just pass data and freq='infer' if relevant; no other kwargs

        # arbitrary TimedeltaIndex; this should work for any TimedeltaIndex
        #  with non-None freq
        tdi = pd.timedelta_range('3 Days', freq='H', periods=9)
        expected = TimedeltaArray(tdi)
        assert expected.freq == tdi.freq

        assert (tdi == expected).all()
        assert (expected == tdi).all()

        result = TimedeltaArray(expected)
        tm.assert_equal(pd.TimedeltaIndex(result), tdi)

        result = TimedeltaArray(list(tdi), freq='infer')
        tm.assert_equal(pd.TimedeltaIndex(result), tdi)

        result = TimedeltaArray(tuple(tdi), freq='infer')
        tm.assert_equal(pd.TimedeltaIndex(result), tdi)

        result = TimedeltaArray(np.array(tdi), freq='infer')
        tm.assert_equal(pd.TimedeltaIndex(result), tdi)

        result = TimedeltaArray(np.array(tdi).astype('m8[s]'), freq='infer')
        tm.assert_equal(pd.TimedeltaIndex(result), tdi)

        result = TimedeltaArray(pd.Series(tdi), freq='infer')
        tm.assert_equal(pd.TimedeltaIndex(result), tdi)

        result = TimedeltaArray(pd.Series(tdi, dtype=object), freq='infer')
        tm.assert_equal(pd.TimedeltaIndex(result), tdi)
