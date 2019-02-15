import os

import pytest

import pandas.util._test_decorators as td

from pandas import DataFrame
import pandas.util.testing as tm

from pandas.io.excel import ExcelFile, read_excel
from pandas.io.parsers import read_csv

_seriesd = tm.getSeriesData()
_tsd = tm.getTimeSeriesData()
_frame = DataFrame(_seriesd)[:10]
_frame2 = DataFrame(_seriesd, columns=['D', 'C', 'B', 'A'])[:10]
_tsframe = tm.makeTimeDataFrame()[:5]
_mixed_frame = _frame.copy()
_mixed_frame['foo'] = 'bar'


@td.skip_if_no('xlrd', '1.0.0')
class SharedItems(object):

    @pytest.fixture(autouse=True)
    def setup_method(self, datapath):
        self.dirpath = datapath("io", "data")
        self.frame = _frame.copy()
        self.frame2 = _frame2.copy()
        self.tsframe = _tsframe.copy()
        self.mixed_frame = _mixed_frame.copy()

    def get_csv_refdf(self, basename):
        """
        Obtain the reference data from read_csv with the Python engine.

        Parameters
        ----------

        basename : str
            File base name, excluding file extension.

        Returns
        -------

        dfref : DataFrame
        """
        pref = os.path.join(self.dirpath, basename + '.csv')
        dfref = read_csv(pref, index_col=0, parse_dates=True, engine='python')
        return dfref

    def get_excelfile(self, basename, ext):
        """
        Return test data ExcelFile instance.

        Parameters
        ----------

        basename : str
            File base name, excluding file extension.

        Returns
        -------

        excel : io.excel.ExcelFile
        """
        return ExcelFile(os.path.join(self.dirpath, basename + ext))

    def get_exceldf(self, basename, ext, *args, **kwds):
        """
        Return test data DataFrame.

        Parameters
        ----------

        basename : str
            File base name, excluding file extension.

        Returns
        -------

        df : DataFrame
        """
        pth = os.path.join(self.dirpath, basename + ext)
        return read_excel(pth, *args, **kwds)
