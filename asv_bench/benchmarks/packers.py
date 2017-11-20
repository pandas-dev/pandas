from .pandas_vb_common import *
from numpy.random import randint
import pandas as pd
from collections import OrderedDict
from pandas.compat import BytesIO
import sqlite3
import os
from sqlalchemy import create_engine
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


class packers_read_hdf_store(_Packers):

    def setup(self):
        self._setup()
        self.df2.to_hdf(self.f, 'df')

    def time_packers_read_hdf_store(self):
        pd.read_hdf(self.f, 'df')


class packers_read_hdf_table(_Packers):

    def setup(self):
        self._setup()
        self.df2.to_hdf(self.f, 'df', format='table')

    def time_packers_read_hdf_table(self):
        pd.read_hdf(self.f, 'df')


class packers_read_json(_Packers):

    def setup(self):
        self._setup()
        self.df.to_json(self.f, orient='split')
        self.df.index = np.arange(self.N)

    def time_packers_read_json(self):
        pd.read_json(self.f, orient='split')


class packers_read_json_date_index(_Packers):

    def setup(self):
        self._setup()
        self.remove(self.f)
        self.df.to_json(self.f, orient='split')

    def time_packers_read_json_date_index(self):
        pd.read_json(self.f, orient='split')


class packers_read_pack(_Packers):

    def setup(self):
        self._setup()
        self.df2.to_msgpack(self.f)

    def time_packers_read_pack(self):
        pd.read_msgpack(self.f)


class packers_read_pickle(_Packers):

    def setup(self):
        self._setup()
        self.df2.to_pickle(self.f)

    def time_packers_read_pickle(self):
        pd.read_pickle(self.f)


class packers_read_sql(_Packers):

    def setup(self):
        self._setup()
        self.engine = create_engine('sqlite:///:memory:')
        self.df2.to_sql('table', self.engine, if_exists='replace')

    def time_packers_read_sql(self):
        pd.read_sql_table('table', self.engine)


class packers_read_stata(_Packers):

    def setup(self):
        self._setup()
        self.df.to_stata(self.f, {'index': 'tc', })

    def time_packers_read_stata(self):
        pd.read_stata(self.f)


class packers_read_stata_with_validation(_Packers):

    def setup(self):
        self._setup()
        self.df['int8_'] = [randint(np.iinfo(np.int8).min, (np.iinfo(np.int8).max - 27)) for _ in range(self.N)]
        self.df['int16_'] = [randint(np.iinfo(np.int16).min, (np.iinfo(np.int16).max - 27)) for _ in range(self.N)]
        self.df['int32_'] = [randint(np.iinfo(np.int32).min, (np.iinfo(np.int32).max - 27)) for _ in range(self.N)]
        self.df['float32_'] = np.array(randn(self.N), dtype=np.float32)
        self.df.to_stata(self.f, {'index': 'tc', })

    def time_packers_read_stata_with_validation(self):
        pd.read_stata(self.f)


class packers_read_sas(_Packers):

    def setup(self):

        testdir = os.path.join(os.path.dirname(__file__), '..', '..',
                               'pandas', 'tests', 'io', 'sas')
        if not os.path.exists(testdir):
            testdir = os.path.join(os.path.dirname(__file__), '..', '..',
                                   'pandas', 'io', 'tests', 'sas')
        self.f = os.path.join(testdir, 'data', 'test1.sas7bdat')
        self.f2 = os.path.join(testdir, 'data', 'paxraw_d_short.xpt')

    def time_read_sas7bdat(self):
        pd.read_sas(self.f, format='sas7bdat')

    def time_read_xport(self):
        pd.read_sas(self.f2, format='xport')


class CSV(_Packers):

    def setup(self):
        self._setup()

    def time_write_csv(self):
        self.df.to_csv(self.f)


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


class HDF(_Packers):

    def setup(self):
        self._setup()

    def time_write_hdf_store(self):
        self.df2.to_hdf(self.f, 'df')

    def time_write_hdf_table(self):
        self.df2.to_hdf(self.f, 'df', table=True)


class JSON(_Packers):

    def setup(self):
        self._setup()
        self.df_date = self.df.copy()
        self.df.index = np.arange(self.N)
        self.cols = [(lambda i: ('{0}_timedelta'.format(i), [pd.Timedelta(('%d seconds' % randrange(1000000.0))) for _ in range(self.N)])), (lambda i: ('{0}_int'.format(i), randint(100000000.0, size=self.N))), (lambda i: ('{0}_timestamp'.format(i), [pd.Timestamp((1418842918083256000 + randrange(1000000000.0, 1e+18, 200))) for _ in range(self.N)]))]
        self.df_mixed = DataFrame(OrderedDict([self.cols[(i % len(self.cols))](i) for i in range(self.C)]), index=self.index)

        self.cols = [(lambda i: ('{0}_float'.format(i), randn(self.N))), (lambda i: ('{0}_int'.format(i), randint(100000000.0, size=self.N)))]
        self.df_mixed2 = DataFrame(OrderedDict([self.cols[(i % len(self.cols))](i) for i in range(self.C)]), index=self.index)

        self.cols = [(lambda i: ('{0}_float'.format(i), randn(self.N))), (lambda i: ('{0}_int'.format(i), randint(100000000.0, size=self.N))), (lambda i: ('{0}_str'.format(i), [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]))]
        self.df_mixed3 = DataFrame(OrderedDict([self.cols[(i % len(self.cols))](i) for i in range(self.C)]), index=self.index)

    def time_write_json(self):
        self.df.to_json(self.f, orient='split')

    def time_write_json_T(self):
        self.df.to_json(self.f, orient='columns')

    def time_write_json_date_index(self):
        self.df_date.to_json(self.f, orient='split')

    def time_write_json_mixed_delta_int_tstamp(self):
        self.df_mixed.to_json(self.f, orient='split')

    def time_write_json_mixed_float_int(self):
        self.df_mixed2.to_json(self.f, orient='index')

    def time_write_json_mixed_float_int_T(self):
        self.df_mixed2.to_json(self.f, orient='columns')

    def time_write_json_mixed_float_int_str(self):
        self.df_mixed3.to_json(self.f, orient='split')

    def time_write_json_lines(self):
        self.df.to_json(self.f, orient="records", lines=True)


class MsgPack(_Packers):

    def setup(self):
        self._setup()

    def time_write_msgpack(self):
        self.df2.to_msgpack(self.f)


class Pickle(_Packers):

    def setup(self):
        self._setup()

    def time_write_pickle(self):
        self.df2.to_pickle(self.f)


class SQL(_Packers):

    def setup(self):
        self._setup()
        self.engine = create_engine('sqlite:///:memory:')

    def time_write_sql(self):
        self.df2.to_sql('table', self.engine, if_exists='replace')


class STATA(_Packers):

    def setup(self):
        self._setup()

        self.df3=self.df.copy()
        self.df3['int8_'] = [randint(np.iinfo(np.int8).min, (np.iinfo(np.int8).max - 27)) for _ in range(self.N)]
        self.df3['int16_'] = [randint(np.iinfo(np.int16).min, (np.iinfo(np.int16).max - 27)) for _ in range(self.N)]
        self.df3['int32_'] = [randint(np.iinfo(np.int32).min, (np.iinfo(np.int32).max - 27)) for _ in range(self.N)]
        self.df3['float32_'] = np.array(randn(self.N), dtype=np.float32)

    def time_write_stata(self):
        self.df.to_stata(self.f, {'index': 'tc', })

    def time_write_stata_with_validation(self):
        self.df3.to_stata(self.f, {'index': 'tc', })
