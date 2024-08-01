from io import StringIO
import pytest
import pandas as pd

def test_large_number():
    assert pd.read_json(StringIO('["9999999999999999"]'),orient="values",typ="series",convert_dates=False)[0] == 9999999999999999
