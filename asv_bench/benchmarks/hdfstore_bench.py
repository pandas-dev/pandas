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
        self.store.put('df1', self.df)
        self.store.put('df_mixed', self.df_mixed)

        self.store.append('df5', self.df_mixed)
        self.store.append('df7', self.df)

        self.store.append('df9', self.df_wide)

        self.store.append('df11', self.df_wide2)
        self.store.append('df12', self.df2)

    def teardown(self):
        self.store.close()

    def remove(self, f):
        try:
            os.remove(self.f)
        except:
            pass

    def time_read_store(self):
        self.store.get('df1')

    def time_read_store_mixed(self):
        self.store.get('df_mixed')

    def time_write_store(self):
        self.store.put('df2', self.df)

    def time_write_store_mixed(self):
        self.store.put('df_mixed2', self.df_mixed)

    def time_read_store_table_mixed(self):
        self.store.select('df5')

    def time_write_store_table_mixed(self):
        self.store.append('df6', self.df_mixed)

    def time_read_store_table(self):
        self.store.select('df7')

    def time_write_store_table(self):
        self.store.append('df8', self.df)

    def time_read_store_table_wide(self):
        self.store.select('df9')

    def time_write_store_table_wide(self):
        self.store.append('df10', self.df_wide)

    def time_write_store_table_dc(self):
        self.store.append('df15', self.df, data_columns=True)

    def time_query_store_table_wide(self):
        self.store.select('df11', [('index', '>', self.df_wide2.index[10000]),
                                   ('index', '<', self.df_wide2.index[15000])])

    def time_query_store_table(self):
        self.store.select('df12', [('index', '>', self.df2.index[10000]),
                                   ('index', '<', self.df2.index[15000])])


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
