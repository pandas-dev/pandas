from numpy.random import randint
import pandas as pd
from collections import OrderedDict
from pandas.compat import BytesIO
import sqlite3
from pandas_vb_common import *
import os
from sqlalchemy import create_engine
import numpy as np
from random import randrange
from pandas.core import common as com


class packers_read_csv(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.df.to_csv(self.f)

    def time_packers_read_csv(self):
        pd.read_csv(self.f)


class packers_read_excel(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.bio = BytesIO()
        self.writer = pd.io.excel.ExcelWriter(self.bio, engine='xlsxwriter')
        self.df[:2000].to_excel(self.writer)
        self.writer.save()

    def time_packers_read_excel(self):
        self.bio.seek(0)
        pd.read_excel(self.bio)


class packers_read_hdf_store(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.df2.to_hdf(self.f, 'df')

    def time_packers_read_hdf_store(self):
        pd.read_hdf(self.f, 'df')


class packers_read_hdf_table(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.df2.to_hdf(self.f, 'df', format='table')

    def time_packers_read_hdf_table(self):
        pd.read_hdf(self.f, 'df')


class packers_read_json(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.df.to_json(self.f, orient='split')
        self.df.index = np.arange(self.N)

    def time_packers_read_json(self):
        pd.read_json(self.f, orient='split')


class packers_read_json_date_index(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.df.to_json(self.f, orient='split')

    def time_packers_read_json_date_index(self):
        pd.read_json(self.f, orient='split')


class packers_read_pack(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.df2.to_msgpack(self.f)

    def time_packers_read_pack(self):
        pd.read_msgpack(self.f)


class packers_read_pickle(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.df2.to_pickle(self.f)

    def time_packers_read_pickle(self):
        pd.read_pickle(self.f)


class packers_read_sql(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.engine = create_engine('sqlite:///:memory:')
        self.df2.to_sql('table', self.engine, if_exists='replace')

    def time_packers_read_sql(self):
        pd.read_sql_table('table', self.engine)


class packers_read_stata(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.df.to_stata(self.f, {'index': 'tc', })

    def time_packers_read_stata(self):
        pd.read_stata(self.f)


class packers_read_stata_with_validation(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.df['int8_'] = [randint(np.iinfo(np.int8).min, (np.iinfo(np.int8).max - 27)) for _ in range(self.N)]
        self.df['int16_'] = [randint(np.iinfo(np.int16).min, (np.iinfo(np.int16).max - 27)) for _ in range(self.N)]
        self.df['int32_'] = [randint(np.iinfo(np.int32).min, (np.iinfo(np.int32).max - 27)) for _ in range(self.N)]
        self.df['float32_'] = np.array(randn(self.N), dtype=np.float32)
        self.df.to_stata(self.f, {'index': 'tc', })

    def time_packers_read_stata_with_validation(self):
        pd.read_stata(self.f)


class packers_write_csv(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)

    def time_packers_write_csv(self):
        self.df.to_csv(self.f)

    def teardown(self):
        remove(self.f)


class packers_write_excel_openpyxl(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.bio = BytesIO()

    def time_packers_write_excel_openpyxl(self):
        self.bio.seek(0)
        self.writer = pd.io.excel.ExcelWriter(self.bio, engine='openpyxl')
        self.df[:2000].to_excel(self.writer)
        self.writer.save()


class packers_write_excel_xlsxwriter(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.bio = BytesIO()

    def time_packers_write_excel_xlsxwriter(self):
        self.bio.seek(0)
        self.writer = pd.io.excel.ExcelWriter(self.bio, engine='xlsxwriter')
        self.df[:2000].to_excel(self.writer)
        self.writer.save()


class packers_write_excel_xlwt(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.bio = BytesIO()

    def time_packers_write_excel_xlwt(self):
        self.bio.seek(0)
        self.writer = pd.io.excel.ExcelWriter(self.bio, engine='xlwt')
        self.df[:2000].to_excel(self.writer)
        self.writer.save()


class packers_write_hdf_store(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)

    def time_packers_write_hdf_store(self):
        self.df2.to_hdf(self.f, 'df')

    def teardown(self):
        remove(self.f)


class packers_write_hdf_table(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)

    def time_packers_write_hdf_table(self):
        self.df2.to_hdf(self.f, 'df', table=True)

    def teardown(self):
        remove(self.f)


class packers_write_json(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.df.index = np.arange(self.N)

    def time_packers_write_json(self):
        self.df.to_json(self.f, orient='split')

    def teardown(self):
        remove(self.f)


class packers_write_json_T(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.df.index = np.arange(self.N)

    def time_packers_write_json_T(self):
        self.df.to_json(self.f, orient='columns')

    def teardown(self):
        remove(self.f)


class packers_write_json_date_index(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)

    def time_packers_write_json_date_index(self):
        self.df.to_json(self.f, orient='split')

    def teardown(self):
        remove(self.f)


class packers_write_json_mixed_delta_int_tstamp(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.cols = [(lambda i: ('{0}_timedelta'.format(i), [pd.Timedelta(('%d seconds' % randrange(1000000.0))) for _ in range(self.N)])), (lambda i: ('{0}_int'.format(i), randint(100000000.0, size=self.N))), (lambda i: ('{0}_timestamp'.format(i), [pd.Timestamp((1418842918083256000 + randrange(1000000000.0, 1e+18, 200))) for _ in range(self.N)]))]
        self.df_mixed = DataFrame(OrderedDict([self.cols[(i % len(self.cols))](i) for i in range(self.C)]), index=self.index)

    def time_packers_write_json_mixed_delta_int_tstamp(self):
        self.df_mixed.to_json(self.f, orient='split')

    def teardown(self):
        remove(self.f)


class packers_write_json_mixed_float_int(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.cols = [(lambda i: ('{0}_float'.format(i), randn(self.N))), (lambda i: ('{0}_int'.format(i), randint(100000000.0, size=self.N)))]
        self.df_mixed = DataFrame(OrderedDict([self.cols[(i % len(self.cols))](i) for i in range(self.C)]), index=self.index)

    def time_packers_write_json_mixed_float_int(self):
        self.df_mixed.to_json(self.f, orient='index')

    def teardown(self):
        remove(self.f)


class packers_write_json_mixed_float_int_T(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.cols = [(lambda i: ('{0}_float'.format(i), randn(self.N))), (lambda i: ('{0}_int'.format(i), randint(100000000.0, size=self.N)))]
        self.df_mixed = DataFrame(OrderedDict([self.cols[(i % len(self.cols))](i) for i in range(self.C)]), index=self.index)

    def time_packers_write_json_mixed_float_int_T(self):
        self.df_mixed.to_json(self.f, orient='columns')

    def teardown(self):
        remove(self.f)


class packers_write_json_mixed_float_int_str(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.cols = [(lambda i: ('{0}_float'.format(i), randn(self.N))), (lambda i: ('{0}_int'.format(i), randint(100000000.0, size=self.N))), (lambda i: ('{0}_str'.format(i), [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]))]
        self.df_mixed = DataFrame(OrderedDict([self.cols[(i % len(self.cols))](i) for i in range(self.C)]), index=self.index)

    def time_packers_write_json_mixed_float_int_str(self):
        self.df_mixed.to_json(self.f, orient='split')

    def teardown(self):
        remove(self.f)


class packers_write_pack(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)

    def time_packers_write_pack(self):
        self.df2.to_msgpack(self.f)

    def teardown(self):
        remove(self.f)


class packers_write_pickle(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)

    def time_packers_write_pickle(self):
        self.df2.to_pickle(self.f)

    def teardown(self):
        remove(self.f)


class packers_write_sql(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.engine = create_engine('sqlite:///:memory:')

    def time_packers_write_sql(self):
        self.df2.to_sql('table', self.engine, if_exists='replace')


class packers_write_stata(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.df.to_stata(self.f, {'index': 'tc', })

    def time_packers_write_stata(self):
        self.df.to_stata(self.f, {'index': 'tc', })

    def teardown(self):
        remove(self.f)


class packers_write_stata_with_validation(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.msg'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.N = 100000
        self.C = 5
        self.index = date_range('20000101', periods=self.N, freq='H')
        self.df2 = DataFrame(dict([('float{0}'.format(i), randn(self.N)) for i in range(self.C)]), index=self.index)
        self.df2['object'] = [('%08x' % randrange((16 ** 8))) for _ in range(self.N)]
        remove(self.f)
        self.df['int8_'] = [randint(np.iinfo(np.int8).min, (np.iinfo(np.int8).max - 27)) for _ in range(self.N)]
        self.df['int16_'] = [randint(np.iinfo(np.int16).min, (np.iinfo(np.int16).max - 27)) for _ in range(self.N)]
        self.df['int32_'] = [randint(np.iinfo(np.int32).min, (np.iinfo(np.int32).max - 27)) for _ in range(self.N)]
        self.df['float32_'] = np.array(randn(self.N), dtype=np.float32)
        self.df.to_stata(self.f, {'index': 'tc', })

    def time_packers_write_stata_with_validation(self):
        self.df.to_stata(self.f, {'index': 'tc', })

    def teardown(self):
        remove(self.f)