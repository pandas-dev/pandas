from pandas import (
    DataFrame,
    Series,
)
import pandas._testing as tm


class TestToFrame:
    def test_to_frame(self, datetime_series):
        datetime_series.name = None
        rs = datetime_series.to_frame()
        xp = DataFrame(datetime_series.values, index=datetime_series.index)
        tm.assert_frame_equal(rs, xp)

        datetime_series.name = "testname"
        rs = datetime_series.to_frame()
        xp = DataFrame(
            {"testname": datetime_series.values}, index=datetime_series.index
        )
        tm.assert_frame_equal(rs, xp)

        rs = datetime_series.to_frame(name="testdifferent")
        xp = DataFrame(
            {"testdifferent": datetime_series.values}, index=datetime_series.index
        )
        tm.assert_frame_equal(rs, xp)

    def test_to_frame_expanddim(self):
        # GH#9762

        class SubclassedSeries(Series):
            @property
            def _constructor_expanddim(self):
                return SubclassedFrame

        class SubclassedFrame(DataFrame):
            pass

        ser = SubclassedSeries([1, 2, 3], name="X")
        result = ser.to_frame()
        assert isinstance(result, SubclassedFrame)
        expected = SubclassedFrame({"X": [1, 2, 3]})
        tm.assert_frame_equal(result, expected)

    def test_to_frame_finalize(self, datetime_series):
        # GH#28283 Call __finalize__
        expected_attrs = {"a": 1}

        datetime_series.attrs.update(expected_attrs)

        df = datetime_series.to_frame()
        assert df.attrs == expected_attrs

        df = datetime_series.to_frame(name="X")
        assert df.attrs == expected_attrs
