import pandas as pd


def test_basic_read():
    ser = pd.Series(range(5))
    view = memoryview(ser)
    assert list(view) == ser.tolist()
