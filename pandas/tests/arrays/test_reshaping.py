import pytest

import pandas as pd
from pandas.core.arrays import ReshapeableArray


class TestReshapeableArray:
    def test_repr(self):
        dti = pd.date_range('2016-01-01', periods=3, tz='US/Pacific')
        ea = dti._data
        ra = ReshapeableArray(ea, shape=ea.shape)

        result = repr(ra)
        expected = (
            "<ReshapeableArray> shape=(3,) Wrapping:\n"
            "<DatetimeArray>\n"
            "['2016-01-01 00:00:00-08:00', '2016-01-02 00:00:00-08:00',\n"
            " '2016-01-03 00:00:00-08:00']\n"
            "Length: 3, dtype: datetime64[ns, US/Pacific]"
        )
        assert result == expected

    def test_reshape(self):
        dti = pd.date_range('2016-01-01', periods=3, tz='US/Pacific')
        ea = dti._data
        ra = ReshapeableArray(ea, shape=ea.shape)
        assert ra.shape == (3,)

        result = ra.reshape(1, -1)
        assert result.shape == (1, 3)

        result = ra.reshape(-1, 1)
        assert result.shape == (3, 1)

        with pytest.raises(ValueError, match="Product of shape"):
            # must match original size
            ra.reshape(2, 2)
        with pytest.raises(ValueError, match="Invalid shape"):
            # No more than 1 "-1"
            ra.reshape(-1, -1)
        with pytest.raises(ValueError, match="Invalid shape"):
            # Nothing less than -1
            ra.reshape(-2, 3)

    def test_ravel(self):
        dti = pd.date_range('2016-01-01', periods=4, tz='US/Pacific')
        ea = dti._data
        ra = ReshapeableArray(ea, shape=(1, 4))
        # TODO: case with e.g. (2, 2) with potential ravel ambiguity

        result = ra.ravel()
        assert result.shape == (4,)
        assert list(result) == list(dti)

    def test_transpose(self):
        dti = pd.date_range('2016-01-01', periods=4, tz='US/Pacific')
        ea = dti._data
        ra = ReshapeableArray(ea, shape=(1, 4))

        result = ra.T
        assert result.shape == (4, 1)

    def test_getitem(self):
        dti = pd.date_range('2016-01-01', periods=4, tz='US/Pacific')
        ea = dti._data

        flat = ReshapeableArray(ea, shape=ea.shape)
        collike = ReshapeableArray(ea, shape=(4, 1))
        rowlike = ReshapeableArray(ea, shape=(1, 4))
        square = ReshapeableArray(ea, shape=(2, 2))

        assert flat[0] == ea[0]
        result = flat[:2]
        assert isinstance(result, ReshapeableArray)
        assert list(flat[:2]) == list(ea[:2])

        result = rowlike[0]
        assert isinstance(result, ReshapeableArray)
        assert result.shape == (4,)
        assert list(result) == list(ea)

        result = rowlike[:]
        assert result.shape == rowlike.shape
        assert result._1dvalues is ea

        # TODO: many more untested cases
