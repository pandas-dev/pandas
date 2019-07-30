""" test scalar indexing, including at and iat """

import numpy as np
import pytest

from pandas import DataFrame, Series, Timedelta, Timestamp, date_range
from pandas.tests.indexing.common import Base
from pandas.util import testing as tm


class TestScalar(Base):
    def test_at_and_iat_get(self):
        def _check(f, func, values=False):

            if f is not None:
                indicies = self.generate_indices(f, values)
                for i in indicies:
                    result = getattr(f, func)[i]
                    expected = self.get_value(f, i, values)
                    tm.assert_almost_equal(result, expected)

        for o in self._objs:

            d = getattr(self, o)

            # iat
            for f in [d["ints"], d["uints"]]:
                _check(f, "iat", values=True)

            for f in [d["labels"], d["ts"], d["floats"]]:
                if f is not None:
                    msg = "iAt based indexing can only have integer indexers"
                    with pytest.raises(ValueError, match=msg):
                        self.check_values(f, "iat")

            # at
            for f in [d["ints"], d["uints"], d["labels"], d["ts"], d["floats"]]:
                _check(f, "at")

    def test_at_and_iat_set(self):
        def _check(f, func, values=False):

            if f is not None:
                indicies = self.generate_indices(f, values)
                for i in indicies:
                    getattr(f, func)[i] = 1
                    expected = self.get_value(f, i, values)
                    tm.assert_almost_equal(expected, 1)

        for t in self._objs:

            d = getattr(self, t)

            # iat
            for f in [d["ints"], d["uints"]]:
                _check(f, "iat", values=True)

            for f in [d["labels"], d["ts"], d["floats"]]:
                if f is not None:
                    msg = "iAt based indexing can only have integer indexers"
                    with pytest.raises(ValueError, match=msg):
                        _check(f, "iat")

            # at
            for f in [d["ints"], d["uints"], d["labels"], d["ts"], d["floats"]]:
                _check(f, "at")

    def test_at_iat_coercion(self):

        # as timestamp is not a tuple!
        dates = date_range("1/1/2000", periods=8)
        df = DataFrame(np.random.randn(8, 4), index=dates, columns=["A", "B", "C", "D"])
        s = df["A"]

        result = s.at[dates[5]]
        xp = s.values[5]
        assert result == xp

        # GH 7729
        # make sure we are boxing the returns
        s = Series(["2014-01-01", "2014-02-02"], dtype="datetime64[ns]")
        expected = Timestamp("2014-02-02")

        for r in [lambda: s.iat[1], lambda: s.iloc[1]]:
            result = r()
            assert result == expected

        s = Series(["1 days", "2 days"], dtype="timedelta64[ns]")
        expected = Timedelta("2 days")

        for r in [lambda: s.iat[1], lambda: s.iloc[1]]:
            result = r()
            assert result == expected

    def test_iat_invalid_args(self):
        pass

    def test_imethods_with_dups(self):

        # GH6493
        # iat/iloc with dups

        s = Series(range(5), index=[1, 1, 2, 2, 3], dtype="int64")
        result = s.iloc[2]
        assert result == 2
        result = s.iat[2]
        assert result == 2

        msg = "index 10 is out of bounds for axis 0 with size 5"
        with pytest.raises(IndexError, match=msg):
            s.iat[10]
        msg = "index -10 is out of bounds for axis 0 with size 5"
        with pytest.raises(IndexError, match=msg):
            s.iat[-10]

        result = s.iloc[[2, 3]]
        expected = Series([2, 3], [2, 2], dtype="int64")
        tm.assert_series_equal(result, expected)

        df = s.to_frame()
        result = df.iloc[2]
        expected = Series(2, index=[0], name=2)
        tm.assert_series_equal(result, expected)

        result = df.iat[2, 0]
        assert result == 2

    def test_at_to_fail(self):
        # at should not fallback
        # GH 7814
        s = Series([1, 2, 3], index=list("abc"))
        result = s.at["a"]
        assert result == 1
        msg = (
            "At based indexing on an non-integer index can only have"
            " non-integer indexers"
        )
        with pytest.raises(ValueError, match=msg):
            s.at[0]

        df = DataFrame({"A": [1, 2, 3]}, index=list("abc"))
        result = df.at["a", "A"]
        assert result == 1
        with pytest.raises(ValueError, match=msg):
            df.at["a", 0]

        s = Series([1, 2, 3], index=[3, 2, 1])
        result = s.at[1]
        assert result == 3
        msg = "At based indexing on an integer index can only have integer indexers"
        with pytest.raises(ValueError, match=msg):
            s.at["a"]

        df = DataFrame({0: [1, 2, 3]}, index=[3, 2, 1])
        result = df.at[1, 0]
        assert result == 3
        with pytest.raises(ValueError, match=msg):
            df.at["a", 0]

        # GH 13822, incorrect error string with non-unique columns when missing
        # column is accessed
        df = DataFrame({"x": [1.0], "y": [2.0], "z": [3.0]})
        df.columns = ["x", "x", "z"]

        # Check that we get the correct value in the KeyError
        with pytest.raises(KeyError, match=r"\['y'\] not in index"):
            df[["x", "y", "z"]]

    def test_at_with_tz(self):
        # gh-15822
        df = DataFrame(
            {
                "name": ["John", "Anderson"],
                "date": [
                    Timestamp(2017, 3, 13, 13, 32, 56),
                    Timestamp(2017, 2, 16, 12, 10, 3),
                ],
            }
        )
        df["date"] = df["date"].dt.tz_localize("Asia/Shanghai")

        expected = Timestamp("2017-03-13 13:32:56+0800", tz="Asia/Shanghai")

        result = df.loc[0, "date"]
        assert result == expected

        result = df.at[0, "date"]
        assert result == expected

    def test_series_set_tz_timestamp(self, tz_naive_fixture):
        # GH 25506
        ts = Timestamp("2017-08-05 00:00:00+0100", tz=tz_naive_fixture)
        result = Series(ts)
        result.at[1] = ts
        expected = Series([ts, ts])
        tm.assert_series_equal(result, expected)

    def test_mixed_index_at_iat_loc_iloc_series(self):
        # GH 19860
        s = Series([1, 2, 3, 4, 5], index=["a", "b", "c", 1, 2])
        for el, item in s.items():
            assert s.at[el] == s.loc[el] == item
        for i in range(len(s)):
            assert s.iat[i] == s.iloc[i] == i + 1

        with pytest.raises(KeyError, match="^4$"):
            s.at[4]
        with pytest.raises(KeyError, match="^4$"):
            s.loc[4]

    def test_mixed_index_at_iat_loc_iloc_dataframe(self):
        # GH 19860
        df = DataFrame(
            [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9]], columns=["a", "b", "c", 1, 2]
        )
        for rowIdx, row in df.iterrows():
            for el, item in row.items():
                assert df.at[rowIdx, el] == df.loc[rowIdx, el] == item

        for row in range(2):
            for i in range(5):
                assert df.iat[row, i] == df.iloc[row, i] == row * 5 + i

        with pytest.raises(KeyError, match="^3$"):
            df.at[0, 3]
        with pytest.raises(KeyError, match="^3$"):
            df.loc[0, 3]

    def test_iat_setter_incompatible_assignment(self):
        # GH 23236
        result = DataFrame({"a": [0, 1], "b": [4, 5]})
        result.iat[0, 0] = None
        expected = DataFrame({"a": [None, 1], "b": [4, 5]})
        tm.assert_frame_equal(result, expected)

    def test_getitem_zerodim_np_array(self):
        # GH24924
        # dataframe __getitem__
        df = DataFrame([[1, 2], [3, 4]])
        result = df[np.array(0)]
        expected = Series([1, 3], name=0)
        tm.assert_series_equal(result, expected)

        # series __getitem__
        s = Series([1, 2])
        result = s[np.array(0)]
        assert result == 1
