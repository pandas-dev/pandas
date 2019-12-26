import pytest

import pandas as pd
import pandas.util.testing as tm


@pytest.mark.parametrize("dtype", [int, float])
def test_read(dtype):
    ser = pd.Series(range(5))
    view = memoryview(ser)
    assert list(view) == ser.tolist()


@pytest.mark.parametrize("dtype", [int, float])
def test_write(dtype):
    ser = pd.Series(range(5))
    expected = ser.copy()
    view = memoryview(ser)

    expected.iloc[2] = 200
    view[2] = 200
    tm.assert_series_equal(ser, expected)


def test_object_raises():
    ser = pd.Series(range(5), dtype=object)

    with pytest.raises(NotImplementedError, match="format O not supported"):
        # Note that this only raises on indexing...
        memoryview(ser)[0]
