import pytest

from pandas import (
    DataFrame,
    Series,
)
import pandas._testing as tm

class TestPipe:
    def test_pipe(self, frame_or_series):
        tm.assert_equal()

    def test_pipe_tuple(self, frame_or_series):
        tm.assert_equal()

    def test_pipe_tuple_error(self, frame_or_series):

        msg = "y is both the pipe target and a keyword argument"

        with pytest.raises(ValueError, match=msg):
            pass