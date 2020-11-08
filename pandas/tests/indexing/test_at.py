from datetime import datetime, timezone

import numpy as np
import pytest

from pandas import DataFrame, Series
import pandas._testing as tm


def test_at_timezone():
    # https://github.com/pandas-dev/pandas/issues/33544
    result = DataFrame({"foo": [datetime(2000, 1, 1)]})
    result.at[0, "foo"] = datetime(2000, 1, 2, tzinfo=timezone.utc)
    expected = DataFrame(
        {"foo": [datetime(2000, 1, 2, tzinfo=timezone.utc)]}, dtype=object
    )
    tm.assert_frame_equal(result, expected)


class TestAtWithDuplicates:
    def test_at_with_duplicate_axes_requires_scalar_lookup(self):
        # GH#33041 check that falling back to loc doesn't allow non-scalar
        #  args to slip in

        arr = np.random.randn(6).reshape(3, 2)
        df = DataFrame(arr, columns=["A", "A"])

        msg = "Invalid call for scalar access"
        with pytest.raises(ValueError, match=msg):
            df.at[[1, 2]]
        with pytest.raises(ValueError, match=msg):
            df.at[1, ["A"]]
        with pytest.raises(ValueError, match=msg):
            df.at[:, "A"]

        with pytest.raises(ValueError, match=msg):
            df.at[[1, 2]] = 1
        with pytest.raises(ValueError, match=msg):
            df.at[1, ["A"]] = 1
        with pytest.raises(ValueError, match=msg):
            df.at[:, "A"] = 1


class TestAtErrors:
    # TODO: De-duplicate/parametrize
    #  test_at_series_raises_key_error, test_at_frame_raises_key_error,
    #  test_at_series_raises_key_error2, test_at_frame_raises_key_error2

    def test_at_series_raises_key_error(self):
        # GH#31724 .at should match .loc

        ser = Series([1, 2, 3], index=[3, 2, 1])
        result = ser.at[1]
        assert result == 3
        result = ser.loc[1]
        assert result == 3

        with pytest.raises(KeyError, match="a"):
            ser.at["a"]
        with pytest.raises(KeyError, match="a"):
            # .at should match .loc
            ser.loc["a"]

    def test_at_frame_raises_key_error(self):
        # GH#31724 .at should match .loc

        df = DataFrame({0: [1, 2, 3]}, index=[3, 2, 1])

        result = df.at[1, 0]
        assert result == 3
        result = df.loc[1, 0]
        assert result == 3

        with pytest.raises(KeyError, match="a"):
            df.at["a", 0]
        with pytest.raises(KeyError, match="a"):
            df.loc["a", 0]

        with pytest.raises(KeyError, match="a"):
            df.at[1, "a"]
        with pytest.raises(KeyError, match="a"):
            df.loc[1, "a"]

    def test_at_series_raises_key_error2(self):
        # at should not fallback
        # GH#7814
        # GH#31724 .at should match .loc
        ser = Series([1, 2, 3], index=list("abc"))
        result = ser.at["a"]
        assert result == 1
        result = ser.loc["a"]
        assert result == 1

        with pytest.raises(KeyError, match="^0$"):
            ser.at[0]
        with pytest.raises(KeyError, match="^0$"):
            ser.loc[0]

    def test_at_frame_raises_key_error2(self):
        # GH#31724 .at should match .loc
        df = DataFrame({"A": [1, 2, 3]}, index=list("abc"))
        result = df.at["a", "A"]
        assert result == 1
        result = df.loc["a", "A"]
        assert result == 1

        with pytest.raises(KeyError, match="^0$"):
            df.at["a", 0]
        with pytest.raises(KeyError, match="^0$"):
            df.loc["a", 0]


@pytest.mark.parametrize("func", ["at", "loc"])
def test_at_assign_float_to_int_frame(func):
    # GH: 26395
    obj = DataFrame([0, 0, 0], index=["A", "B", "C"], columns=["D"])
    getattr(obj, func)["C", "D"] = 44.5
    expected = DataFrame([0, 0, 44.5], index=["A", "B", "C"], columns=["D"])
    tm.assert_frame_equal(obj, expected)


@pytest.mark.parametrize("func", ["at", "loc"])
def test_at_assign_float_to_int_series(func):
    # GH: 26395
    obj = Series([0, 0, 0], index=["A", "B", "C"])
    getattr(obj, func)["C"] = 44.5
    expected = Series([0, 0, 44.5], index=["A", "B", "C"])
    tm.assert_series_equal(obj, expected)
