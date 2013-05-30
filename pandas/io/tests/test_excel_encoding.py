# pylint: disable=E1101
# -*- coding: utf-8 -*-


from pandas.util.py3compat import StringIO, BytesIO, PY3
from datetime import datetime
from os.path import split as psplit
import csv
import os
import sys
import re
import unittest

import nose

from numpy import nan
import numpy as np

from pandas import DataFrame, Series, Index, MultiIndex, DatetimeIndex
import pandas.io.parsers as parsers
from pandas.io.parsers import (read_csv, read_table, read_fwf,
                               ExcelFile, TextFileReader, TextParser)
from pandas.util.testing import (assert_almost_equal,
                                 assert_series_equal, 
                                 network,
                                 ensure_clean)
import pandas.util.testing as tm
import pandas as pd

import pandas.lib as lib
from pandas.util import py3compat
from pandas.lib import Timestamp
from pandas.tseries.index import date_range
import pandas.tseries.tools as tools

from numpy.testing.decorators import slow

from pandas._parser import OverflowError

from pandas.io.parsers import (ExcelFile, ExcelWriter, read_csv)


def _skip_if_no_xlrd():
    try:
        import xlrd
        ver = tuple(map(int, xlrd.__VERSION__.split(".")[:2]))
        if ver < (0, 9):
            raise nose.SkipTest('xlrd not installed, skipping')
    except ImportError:
        raise nose.SkipTest('xlrd not installed, skipping')


def _skip_if_no_xlwt():
    try:
        import xlwt
    except ImportError:
        raise nose.SkipTest('xlwt not installed, skipping')


def _skip_if_no_openpyxl():
    try:
        import openpyxl
    except ImportError:
        raise nose.SkipTest('openpyxl not installed, skipping')


def _skip_if_no_excelsuite():
    _skip_if_no_xlrd()
    _skip_if_no_xlwt()
    _skip_if_no_openpyxl()


_seriesd = tm.getSeriesData()
_tsd = tm.getTimeSeriesData()
_frame = DataFrame(_seriesd)[:10]
_frame2 = DataFrame(_seriesd, columns=['D', 'C', 'B', 'A'])[:10]
_tsframe = tm.makeTimeDataFrame()[:5]
_mixed_frame = _frame.copy()
_mixed_frame['foo'] = 'bar'


class ExcelTests(unittest.TestCase):

    def setUp(self):
        self.dirpath = tm.get_data_path()
        self.xls_ta = os.path.join(self.dirpath, 'excel_test_ascii.xls')
        self.xls_tna = os.path.join(self.dirpath, 'excel_test_noascii.xls')
        self.xls_wa = os.path.join(self.dirpath, 'excel_writer_ascii.xls')
        self.xls_wna = os.path.join(self.dirpath, 'excel_writer_noascii.xls')
        
    def test_excel_output_encoding(self):
        _skip_if_no_xlrd()
        _skip_if_no_xlwt()

        # TESTS IF DataFrame.to_excel() WORKS WITH ENCODING PARAMETER MAKING POSSIBLE TO
        # WORK WITH ENCODINGS OTHER TAN ASCII
        
        #FIRST WITH ONLY ASCII 
        
        data_ascii = {
        'index' : ['A', 'B', 'C', 'C', 'B', 'A'],
        'columns' : ['One', 'One', 'One', 'Two', 'Two', 'Two'],
        'values' : [1., 2., 3., 3., 2., 1.]
        }

        original_ascii = DataFrame(data_ascii)
        
        original_ascii.to_excel(self.xls_ta, sheet_name='DataFrame_TEST')
        
        get_xls_ascii = ExcelFile(self.xls_ta)
        
        saved_ascii = get_xls_ascii.parse('DataFrame_TEST', index_col=None, na_values=['NA'])
        
        # NOW WITH NON-ASCII CHARS AND SUPPLYING THE PARAMETER encoding TO DataFrame.to_excel()
        
        data_noascii = {
            'index' : ['Año', 'Baldío', 'Trócola', 'Mínimo', 'Barça', 'Cigüeña'],
            'columns' : ['Año', 'Narices', 'Búlgaro', 'Libélula', 'Cínico', '1º'],
            'values' : ['Céfiro', 'Tarugo', 'Déspota', 'Camión', 'Añejo', 'º']
        }

        original_noascii = DataFrame(data_noascii)
                    
        original_noascii.to_excel(self.xls_tna, sheet_name='DataFrame_TEST', encoding='utf8')
        
        get_xls_noascii = ExcelFile(self.xls_tna, encoding = 'uft8')
        
        #saved_noascii = get_xls_noascii.parse('DataFrame_TEST', index_col=None, na_values=['NA'])
        
        saved_noascii = get_xls_noascii.parse('DataFrame_TEST', index_col=None, na_values=['NA'])
        
        print original_noascii,saved_noascii
        
        tm.assert_frame_equal(original_ascii, saved_ascii)
        tm.assert_frame_equal(original_noascii, saved_noascii)
        
        
        # TESTS IF CLASS ExcelWriter WORKS WITH ENCODING PARAMETER MAKING POSSIBLE TO
        # WORK WITH ENCODINGS OTHER TAN ASCII
        
        #FIRST WITH ONLY ASCII 
        
        data_ascii_1 = {
        'index' : ['A', 'B', 'C', 'C', 'B', 'A'],
        'columns' : ['One', 'One', 'One', 'Two', 'Two', 'Two'],
        'values' : [1., 2., 3., 3., 2., 1.]
        }
        
        data_ascii_2 = {
        'index' : ['A', 'B', 'C', 'C', 'B', 'A'],
        'columns' : ['One', 'One', 'One', 'Two', 'Two', 'Two'],
        'values' : [1., 2., 3., 3., 2., 1.]
        }
        
        excel_writer_ascii=ExcelWriter(self.xls_wa)
        
        original_ascii_1 = DataFrame(data_ascii_1)
        
        original_ascii_2 = DataFrame(data_ascii_2)
        
        original_ascii_1.to_excel(excel_writer_ascii, sheet_name = 'DataFrame_TEST')
        
        original_ascii_2.to_excel(excel_writer_ascii, sheet_name = 'DataFrame_TEST_2')
        
        excel_writer_ascii.save()
        
        get_xls_writer_ascii = ExcelFile(self.xls_wa)
        
        saved_ascii_1 = get_xls_writer_ascii.parse('DataFrame_TEST', index_col = None, na_values = ['NA'])
        
        saved_ascii_2 = get_xls_writer_ascii.parse('DataFrame_TEST_2', index_col = None, na_values = ['NA'])
        
        # NOW WITH NON-ASCII CHARS AND SUPPLYING THE PARAMETER encoding TO class ExcelWriter
        
        data_noascii_1 = {
            'index' : ['Puño', 'Mísero', 'Brújula', 'Pájaro', 'Barça', 'Cigüeña'],
            'columns' : ['Años', 'Nariz', 'Bígaro', 'Céfiro', '2º', '2€'],
            'values' : ['Tímido', 'Variado', 'Efímero', 'Trágico', 'Compañero', '5º']
        }
        
        data_noascii_2 = {
            'index' : ['Año', 'Baldío', 'Trócola', 'Mínimo', 'Barça', 'Cigüeña'],
            'columns' : ['Año', 'Narices', 'Búlgaro', 'Libélula', 'Cínico', '1º'],
            'values' : ['Céfiro', 'Tarugo', 'Déspota', 'Camión', 'Añejo', 'º']
        }

        excel_writer_noascii=ExcelWriter(self.xls_wna,encoding = 'utf8')
        
        original_noascii_1 = DataFrame(data_noascii_1)
        
        original_noascii_2 = DataFrame(data_noascii_2)
        
        original_noascii_1.to_excel(excel_writer_noascii, sheet_name = 'DataFrame_TEST')
        
        original_noascii_2.to_excel(excel_writer_noascii, sheet_name = 'DataFrame_TEST_2')
        
        excel_writer_noascii.save()
        
        get_xls_writer_noascii = ExcelFile(self.xls_wna,encoding = 'uft8')
        
        saved_noascii_1 = get_xls_writer_noascii.parse('DataFrame_TEST', index_col = None, na_values = ['NA'])
        
        saved_noascii_2 = get_xls_writer_noascii.parse('DataFrame_TEST_2', index_col = None, na_values = ['NA'])
        
        tm.assert_frame_equal(original_ascii_1, saved_ascii_1)
        tm.assert_frame_equal(original_ascii_2, saved_ascii_2)
        
        tm.assert_frame_equal(original_noascii_1, saved_noascii_1)
        tm.assert_frame_equal(original_noascii_2, saved_noascii_2)

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)