import pandas as pd
import pytest
import operator


class CustomOperatorOverload(object):
    def __pandas_ufunc__(self, func, method, *args, **kwargs):
        if func == operator.add:
            return 5
        return NotImplemented


def test_ufunc_implemented():
    assert pd.DataFrame() + CustomOperatorOverload() == 5


def test_ufunc_not_implemented():
    with pytest.raises(AssertionError):
        pd.DataFrame() + CustomOperatorOverload()
