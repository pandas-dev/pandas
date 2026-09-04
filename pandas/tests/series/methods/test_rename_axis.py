import pytest

from pandas.errors import Pandas4Warning

import pandas as pd
import pandas._testing as tm


class TestSeriesRenameAxis:
    def test_rename_axis_mapper(self):
        # GH 19978
        mi = pd.MultiIndex.from_product([["a", "b", "c"], [1, 2]], names=["ll", "nn"])
        ser = pd.Series(list(range(len(mi))), index=mi)

        result = ser.rename_axis(index={"ll": "foo"})
        assert result.index.names == ["foo", "nn"]

        result = ser.rename_axis(index=str.upper, axis=0)
        assert result.index.names == ["LL", "NN"]

        result = ser.rename_axis(index=["foo", "goo"])
        assert result.index.names == ["foo", "goo"]

        with pytest.raises(TypeError, match="unexpected"):
            ser.rename_axis(columns="wrong")

    @pytest.mark.filterwarnings("ignore:The inplace keyword in Series.rename_axis")
    def test_rename_axis_inplace(self, datetime_series):
        # GH 15704
        expected = datetime_series.rename_axis("foo")
        result = datetime_series
        no_return = result.rename_axis("foo", inplace=True)

        assert no_return is None
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("kwargs", [{"mapper": None}, {"index": None}, {}])
    def test_rename_axis_none(self, kwargs):
        # GH 25034
        index = pd.Index(list("abc"), name="foo")
        ser = pd.Series([1, 2, 3], index=index)

        result = ser.rename_axis(**kwargs)
        expected_index = index.rename(None) if kwargs else index
        expected = pd.Series([1, 2, 3], index=expected_index)
        tm.assert_series_equal(result, expected)


def test_rename_axis_inplace_depr():
    msg = "The inplace keyword in Series.rename_axis is deprecated"

    ser = pd.Series([1, 2, 3], index=pd.Index([0, 1, 2], name="idx"))
    ser_orig = ser.copy()
    expected = ser.rename_axis("foo")

    # does not use keyword, no warning
    with tm.assert_produces_warning(False):
        result = ser.rename_axis("foo")
    tm.assert_series_equal(result, expected)
    tm.assert_series_equal(ser, ser_orig)

    # uses keyword, set to false, warning
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = ser.rename_axis("foo", inplace=False)
    tm.assert_series_equal(result, expected)
    tm.assert_series_equal(ser, ser_orig)

    # uses keyword, set to true, warning
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        ser.rename_axis("foo", inplace=True)
    tm.assert_series_equal(ser, expected)
