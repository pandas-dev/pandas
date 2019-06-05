import pandas as pd
import pytest
import operator


class Pipe:
    def __init__(self, function):
        self.function = function

    def __pandas_ufunc__(self, func, method, *args, **kwargs):
        if func == operator.__or__:
            return self.function(args[0])
        return NotImplemented


def test_ufunc_implemented():
    df = pd.DataFrame({'x': [1, 2]})
    result = df | Pipe(lambda x: x + 1)
    assert result.equals(df + 1)


def test_ufunc_not_implemented():
    df = pd.DataFrame({'x': [1, 2]})
    with pytest.raises(AssertionError):
        df * Pipe(lambda x: x + 1)
