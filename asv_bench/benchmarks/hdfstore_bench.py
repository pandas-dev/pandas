from .pandas_vb_common import *
import os


class HDF5(object):
    goal_time = 0.2

    def setup(self):
        self.index = tm.makeStringIndex(25000)
        self.df = DataFrame({'float1': randn(25000), 'float2': randn(25000),},
                            index=self.index)

        self.df_mixed = DataFrame(
            {'float1': randn(25000), 'float2': randn(25000),
             'string1': (['foo'] * 25000),
             'bool1': ([True] * 25000),
             'int1': np.random.randint(0, 250000, size=25000),},
            index=self.index)

        self.df_wide = DataFrame(np.random.randn(25000, 100))

        self.df2 = DataFrame({'float1': randn(25000), 'float2': randn(25000)},
                             index=date_range('1/1/2000', periods=25000))
        self.df_wide2 = DataFrame(np.random.randn(25000, 100),
                                  index=date_range('1/1/2000', periods=25000))

        self.df_dc = DataFrame(np.random.randn(10000, 10),
                               columns=[('C%03d' % i) for i in range(10)])

        self.f = '__test__.h5'
        self.remove(self.f)

        self.store = HDFStore(self.f)
        self.store.put('fixed', self.df)
        self.store.put('fixed_mixed', self.df_mixed)
        self.store.append('table', self.df2)
        self.store.append('table_mixed', self.df_mixed)
        self.store.append('table_wide', self.df_wide)
        self.store.append('table_wide2', self.df_wide2)

    def teardown(self):
        self.store.close()

    def remove(self, f):
        try:
            os.remove(self.f)
        except:
            pass

    def time_read_store(self):
        self.store.get('fixed')

    def time_read_store_mixed(self):
        self.store.get('fixed_mixed')

    def time_write_store(self):
        self.store.put('fixed_write', self.df)

    def time_write_store_mixed(self):
        self.store.put('fixed_mixed_write', self.df_mixed)

    def time_read_store_table_mixed(self):
        self.store.select('table_mixed')

    def time_write_store_table_mixed(self):
        self.store.append('table_mixed_write', self.df_mixed)

    def time_read_store_table(self):
        self.store.select('table')

    def time_write_store_table(self):
        self.store.append('table_write', self.df)

    def time_read_store_table_wide(self):
        self.store.select('table_wide')

    def time_write_store_table_wide(self):
        self.store.append('table_wide_write', self.df_wide)

    def time_write_store_table_dc(self):
        self.store.append('table_dc_write', self.df_dc, data_columns=True)

    def time_query_store_table_wide(self):
        start = self.df_wide2.index[10000]
        stop = self.df_wide2.index[15000]
        self.store.select('table_wide', where="index > start and index < stop")

    def time_query_store_table(self):
        start = self.df2.index[10000]
        stop = self.df2.index[15000]
        self.store.select('table', where="index > start and index < stop")

    def time_store_repr(self):
        repr(self.store)

    def time_store_str(self):
        str(self.store)

    def time_store_info(self):
        self.store.info()


class HDF5Panel(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.h5'
        self.p = Panel(randn(20, 1000, 25),
                       items=[('Item%03d' % i) for i in range(20)],
                       major_axis=date_range('1/1/2000', periods=1000),
                       minor_axis=[('E%03d' % i) for i in range(25)])
        self.remove(self.f)
        self.store = HDFStore(self.f)
        self.store.append('p1', self.p)

    def teardown(self):
        self.store.close()

    def remove(self, f):
        try:
            os.remove(self.f)
        except:
            pass

    def time_read_store_table_panel(self):
        self.store.select('p1')

    def time_write_store_table_panel(self):
        self.store.append('p2', self.p)
