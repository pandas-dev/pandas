""" generic tests from the Datetimelike class """
import pytest

from pandas import DatetimeIndex, date_range
import pandas._testing as tm

from ..datetimelike import DatetimeLike


class TestDatetimeIndex(DatetimeLike):
    _holder = DatetimeIndex

    @pytest.fixture(
        params=[tm.makeDateIndex(10), date_range("20130110", periods=10, freq="-1D")],
        ids=["index_inc", "index_dec"],
    )
    def indices(self, request):
        return request.param

    def create_index(self):
        return date_range("20130101", periods=5)

    def test_shift(self):
        pass  # handled in test_ops

    def test_pickle_compat_construction(self):
        pass

    def test_intersection(self):
        pass  # handled in test_setops

    def test_union(self):
        pass  # handled in test_setops
