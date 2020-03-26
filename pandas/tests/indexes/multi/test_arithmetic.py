"""
Tests for arithmetic being disabled on MultiIndex.
"""
import pytest


def test_sub(idx):

    first = idx

    # - now raises (previously was set op difference)
    msg = "cannot perform __sub__ with this index type: MultiIndex"
    with pytest.raises(TypeError, match=msg):
        first - idx[-3:]
    with pytest.raises(TypeError, match=msg):
        idx[-3:] - first
    with pytest.raises(TypeError, match=msg):
        idx[-3:] - first.tolist()
    msg = "cannot perform __rsub__ with this index type: MultiIndex"
    with pytest.raises(TypeError, match=msg):
        first.tolist() - idx[-3:]


def test_numeric_compat(idx):
    with pytest.raises(TypeError, match="cannot perform __mul__"):
        idx * 1

    with pytest.raises(TypeError, match="cannot perform __rmul__"):
        1 * idx

    div_err = "cannot perform __truediv__"
    with pytest.raises(TypeError, match=div_err):
        idx / 1

    div_err = div_err.replace(" __", " __r")
    with pytest.raises(TypeError, match=div_err):
        1 / idx

    with pytest.raises(TypeError, match="cannot perform __floordiv__"):
        idx // 1

    with pytest.raises(TypeError, match="cannot perform __rfloordiv__"):
        1 // idx
