from .pandas_vb_common import *
from numpy.random import randint
import pandas as pd
from collections import OrderedDict
from pandas.compat import BytesIO
import os
import numpy as np
from random import randrange


class _Packers(object):
    goal_time = 0.2

    def _setup(self):
        self.f = '__test__.msg'
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame({'float{0}'.format(i): randn(self.N) for i in range(self.C)}, index=self.index)
        self.df2 = self.df.copy()
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        self.remove(self.f)

    def remove(self, f):
        try:
            os.remove(f)
        except:
            pass

    def teardown(self):
        self.remove(self.f)


class Packers(_Packers):

    def setup(self):
        self._setup()
        self.df.to_csv(self.f)

    def time_packers_read_csv(self):
        pd.read_csv(self.f)


class packers_read_excel(_Packers):

    def setup(self):
        self._setup()
        self.bio = BytesIO()
        self.writer = pd.io.excel.ExcelWriter(self.bio, engine='xlsxwriter')
        self.df[:2000].to_excel(self.writer)
        self.writer.save()

    def time_packers_read_excel(self):
        self.bio.seek(0)
        pd.read_excel(self.bio)


class Excel(_Packers):

    def setup(self):
        self._setup()
        self.bio = BytesIO()

    def time_write_excel_openpyxl(self):
        self.bio.seek(0)
        self.writer = pd.io.excel.ExcelWriter(self.bio, engine='openpyxl')
        self.df[:2000].to_excel(self.writer)
        self.writer.save()

    def time_write_excel_xlsxwriter(self):
        self.bio.seek(0)
        self.writer = pd.io.excel.ExcelWriter(self.bio, engine='xlsxwriter')
        self.df[:2000].to_excel(self.writer)
        self.writer.save()

    def time_write_excel_xlwt(self):
        self.bio.seek(0)
        self.writer = pd.io.excel.ExcelWriter(self.bio, engine='xlwt')
        self.df[:2000].to_excel(self.writer)
        self.writer.save()
