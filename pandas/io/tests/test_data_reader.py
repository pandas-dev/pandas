import unittest

from pandas.core.generic import PandasObject
from pandas.io.data import DataReader
from pandas.util.testing import network


class TestDataReader(unittest.TestCase):
    @network
    def test_read_yahoo(self):
        gs = DataReader("GS", "yahoo")
        assert isinstance(gs, PandasObject)

    @network
    def test_read_google(self):
        gs = DataReader("GS", "google")
        assert isinstance(gs, PandasObject)

    @network
    def test_read_fred(self):
        vix = DataReader("VIXCLS", "fred")
        assert isinstance(vix, PandasObject)

    @network
    def test_read_famafrench(self):
        for name in ("F-F_Research_Data_Factors",
                     "F-F_Research_Data_Factors_weekly", "6_Portfolios_2x3",
                     "F-F_ST_Reversal_Factor"):
            ff = DataReader(name, "famafrench")
            assert isinstance(ff, dict)
