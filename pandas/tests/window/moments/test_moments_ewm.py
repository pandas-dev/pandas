import pytest

from pandas import (
    DataFrame,
    Series,
)
import pandas._testing as tm


@pytest.mark.parametrize("name", ["var", "std", "mean"])
def test_ewma_series(series, name):
    series_result = getattr(series.ewm(com=10), name)()
    assert isinstance(series_result, Series)


@pytest.mark.parametrize("name", ["var", "std", "mean"])
def test_ewma_frame(frame, name):
    frame_result = getattr(frame.ewm(com=10), name)()
    assert isinstance(frame_result, DataFrame)


def test_ewma_span_com_args(series):
    A = series.ewm(com=9.5).mean()
    B = series.ewm(span=20).mean()
    tm.assert_almost_equal(A, B)
    msg = "comass, span, halflife, and alpha are mutually exclusive"
    with pytest.raises(ValueError, match=msg):
        series.ewm(com=9.5, span=20)

    msg = "Must pass one of comass, span, halflife, or alpha"
    with pytest.raises(ValueError, match=msg):
        series.ewm().mean()


def test_ewma_halflife_arg(series):
    A = series.ewm(com=13.932726172912965).mean()
    B = series.ewm(halflife=10.0).mean()
    tm.assert_almost_equal(A, B)
    msg = "comass, span, halflife, and alpha are mutually exclusive"
    with pytest.raises(ValueError, match=msg):
        series.ewm(span=20, halflife=50)
    with pytest.raises(ValueError, match=msg):
        series.ewm(com=9.5, halflife=50)
    with pytest.raises(ValueError, match=msg):
        series.ewm(com=9.5, span=20, halflife=50)
    msg = "Must pass one of comass, span, halflife, or alpha"
    with pytest.raises(ValueError, match=msg):
        series.ewm()


def test_ewm_alpha_arg(series):
    # GH 10789
    s = series
    msg = "Must pass one of comass, span, halflife, or alpha"
    with pytest.raises(ValueError, match=msg):
        s.ewm()

    msg = "comass, span, halflife, and alpha are mutually exclusive"
    with pytest.raises(ValueError, match=msg):
        s.ewm(com=10.0, alpha=0.5)
    with pytest.raises(ValueError, match=msg):
        s.ewm(span=10.0, alpha=0.5)
    with pytest.raises(ValueError, match=msg):
        s.ewm(halflife=10.0, alpha=0.5)
