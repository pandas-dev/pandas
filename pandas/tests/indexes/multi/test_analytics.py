import pytest


def test_shift(idx):

    # GH8083 test the base class for shift
    pytest.raises(NotImplementedError, idx.shift, 1)
    pytest.raises(NotImplementedError, idx.shift, 1, 2)
