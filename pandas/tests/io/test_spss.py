import pytest
import numpy as np
import pandas as pd
from pandas.io.spss import read_spss
from pandas.util import testing as tm


try:
    import pyreadstat  # noqa
except ImportError:
    _HAVE_PYREADSTAT = False
else:
    _HAVE_PYREADSTAT = True


@pytest.mark.skipif(not _HAVE_PYREADSTAT)
def test_spss_labelled_num():
    fname = "data/labelled-num.sav"

    df = read_spss(fname, categorical=True)
    expected = pd.DataFrame({"VAR00002": "This is one"}, index=[0])
    expected["VAR00002"] = pd.Categorical(expected["VAR00002"])
    tm.assert_frame_equal(df, expected)

    df = read_spss(fname, categorical=False)
    expected = pd.DataFrame({"VAR00002": 1.0}, index=[0])
    tm.assert_frame_equal(df, expected)


@pytest.mark.skipif(not _HAVE_PYREADSTAT)
def test_spss_labelled_num_na():
    fname = "data/labelled-num-na.sav"

    df = read_spss(fname, categorical=True)
    expected = pd.DataFrame({"VAR00002": ["This is one", None]})
    expected["VAR00002"] = pd.Categorical(expected["VAR00002"])
    tm.assert_frame_equal(df, expected)

    df = read_spss(fname, categorical=False)
    expected = pd.DataFrame({"VAR00002": [1.0, np.nan]})
    tm.assert_frame_equal(df, expected)


@pytest.mark.skipif(not _HAVE_PYREADSTAT)
def test_spss_labelled_str():
    fname = "data/labelled-str.sav"

    df = read_spss(fname, categorical=True)
    expected = pd.DataFrame({"gender": ["Male", "Female"]})
    expected["gender"] = pd.Categorical(expected["gender"])
    tm.assert_frame_equal(df, expected)

    df = read_spss(fname, categorical=False)
    expected = pd.DataFrame({"gender": ["M", "F"]})
    tm.assert_frame_equal(df, expected)


@pytest.mark.skipif(not _HAVE_PYREADSTAT)
def test_spss_umlauts():
    fname = "data/umlauts.sav"

    df = read_spss(fname, categorical=True)
    expected = pd.DataFrame({"var1": ["the ä umlaut",
                                      "the ü umlaut",
                                      "the ä umlaut",
                                      "the ö umlaut"]})
    expected["var1"] = pd.Categorical(expected["var1"])
    tm.assert_frame_equal(df, expected)

    df = read_spss(fname, categorical=False)
    expected = pd.DataFrame({"var1": [1.0, 2.0, 1.0, 3.0]})
    tm.assert_frame_equal(df, expected)
