"""
Tests for DatetimeArray
"""
import pytest
import numpy as np

from pandas.core.arrays import DatetimeArrayMixin as DatetimeArray

import pandas as pd
import pandas.util.testing as tm


class TestDatetimeArrayConstructors(object):

    def test_init_from_object_dtype(self, tz_naive_fixture):
        tz = tz_naive_fixture
        if tz is not None:
            pytest.xfail(reason="Casting DatetimeIndex to object-dtype raises "
                                "for pd.Index and is incorrect for np.array; "
                                "GH#24391")

        # arbitrary DatetimeIndex; this should work for any DatetimeIndex
        #  with non-None freq
        dti = pd.date_range('2016-01-1', freq='MS', periods=9, tz=tz)
        expected = DatetimeArray(dti)

        # Fails because np.array(dti, dtype=object) incorrectly returns Longs
        result = DatetimeArray(np.array(dti, dtype=object), freq='infer')
        tm.assert_equal(pd.DatetimeIndex(result), dti)

        # Fails because `pd.Index(dti, dtype=object) raises incorrectly
        result = DatetimeArray(pd.Index(dti, dtype=object), freq='infer')
        tm.assert_equal(pd.DatetimeIndex(result), dti)

    # NB: for now we re-wrap in DatetimeIndex to use assert_index_equal
    #  once assert_datetime_array_equal is in place, this will be changed
    def test_init_only_freq_infer(self, tz_naive_fixture):
        # just pass data and freq='infer' if relevant; no other kwargs
        tz = tz_naive_fixture

        # arbitrary DatetimeIndex; this should work for any DatetimeIndex
        #  with non-None freq
        dti = pd.date_range('2016-01-1', freq='MS', periods=9, tz=tz)
        expected = DatetimeArray(dti)
        assert expected.freq == dti.freq
        assert expected.tz == dti.tz

        # broken until ABCDatetimeArray and isna is fixed
        # assert (dti == expected).all()
        # assert (expected == dti).all()

        result = DatetimeArray(expected)
        tm.assert_equal(pd.DatetimeIndex(result), dti)

        result = DatetimeArray(list(dti), freq='infer')
        tm.assert_equal(pd.DatetimeIndex(result), dti)

        result = DatetimeArray(tuple(dti), freq='infer')
        tm.assert_equal(pd.DatetimeIndex(result), dti)

        if tz is None:
            result = DatetimeArray(np.array(dti), freq='infer')
            tm.assert_equal(pd.DatetimeIndex(result), dti)

            result = DatetimeArray(np.array(dti).astype('M8[s]'), freq='infer')
            tm.assert_equal(pd.DatetimeIndex(result), dti)

        result = DatetimeArray(pd.Series(dti), freq='infer')
        tm.assert_equal(pd.DatetimeIndex(result), dti)

        result = DatetimeArray(pd.Series(dti, dtype=object), freq='infer')
        tm.assert_equal(pd.DatetimeIndex(result), dti)
