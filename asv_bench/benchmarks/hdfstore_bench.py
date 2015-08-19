from pandas_vb_common import *
import os


class query_store_table(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.h5'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.index = date_range('1/1/2000', periods=25000)
        self.df = DataFrame({'float1': randn(25000), 'float2': randn(25000), }, index=self.index)
        remove(self.f)
        self.store = HDFStore(self.f)
        self.store.append('df12', self.df)

    def time_query_store_table(self):
        self.store.select('df12', [('index', '>', self.df.index[10000]), ('index', '<', self.df.index[15000])])

    def teardown(self):
        self.store.close()


class query_store_table_wide(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.h5'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.index = date_range('1/1/2000', periods=25000)
        self.df = DataFrame(np.random.randn(25000, 100), index=self.index)
        remove(self.f)
        self.store = HDFStore(self.f)
        self.store.append('df11', self.df)

    def time_query_store_table_wide(self):
        self.store.select('df11', [('index', '>', self.df.index[10000]), ('index', '<', self.df.index[15000])])

    def teardown(self):
        self.store.close()


class read_store(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.h5'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.index = tm.makeStringIndex(25000)
        self.df = DataFrame({'float1': randn(25000), 'float2': randn(25000), }, index=self.index)
        remove(self.f)
        self.store = HDFStore(self.f)
        self.store.put('df1', self.df)

    def time_read_store(self):
        self.store.get('df1')

    def teardown(self):
        self.store.close()


class read_store_mixed(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.h5'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.index = tm.makeStringIndex(25000)
        self.df = DataFrame({'float1': randn(25000), 'float2': randn(25000), 'string1': (['foo'] * 25000), 'bool1': ([True] * 25000), 'int1': np.random.randint(0, 250000, size=25000), }, index=self.index)
        remove(self.f)
        self.store = HDFStore(self.f)
        self.store.put('df3', self.df)

    def time_read_store_mixed(self):
        self.store.get('df3')

    def teardown(self):
        self.store.close()


class read_store_table(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.h5'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.index = tm.makeStringIndex(25000)
        self.df = DataFrame({'float1': randn(25000), 'float2': randn(25000), }, index=self.index)
        remove(self.f)
        self.store = HDFStore(self.f)
        self.store.append('df7', self.df)

    def time_read_store_table(self):
        self.store.select('df7')

    def teardown(self):
        self.store.close()


class read_store_table_mixed(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.h5'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.N = 10000
        self.index = tm.makeStringIndex(self.N)
        self.df = DataFrame({'float1': randn(self.N), 'float2': randn(self.N), 'string1': (['foo'] * self.N), 'bool1': ([True] * self.N), 'int1': np.random.randint(0, self.N, size=self.N), }, index=self.index)
        remove(self.f)
        self.store = HDFStore(self.f)
        self.store.append('df5', self.df)

    def time_read_store_table_mixed(self):
        self.store.select('df5')

    def teardown(self):
        self.store.close()


class read_store_table_panel(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.h5'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.p = Panel(randn(20, 1000, 25), items=[('Item%03d' % i) for i in xrange(20)], major_axis=date_range('1/1/2000', periods=1000), minor_axis=[('E%03d' % i) for i in xrange(25)])
        remove(self.f)
        self.store = HDFStore(self.f)
        self.store.append('p1', self.p)

    def time_read_store_table_panel(self):
        self.store.select('p1')

    def teardown(self):
        self.store.close()


class read_store_table_wide(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.h5'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.df = DataFrame(np.random.randn(25000, 100))
        remove(self.f)
        self.store = HDFStore(self.f)
        self.store.append('df9', self.df)

    def time_read_store_table_wide(self):
        self.store.select('df9')

    def teardown(self):
        self.store.close()


class write_store(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.h5'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.index = tm.makeStringIndex(25000)
        self.df = DataFrame({'float1': randn(25000), 'float2': randn(25000), }, index=self.index)
        remove(self.f)
        self.store = HDFStore(self.f)

    def time_write_store(self):
        self.store.put('df2', self.df)

    def teardown(self):
        self.store.close()


class write_store_mixed(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.h5'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.index = tm.makeStringIndex(25000)
        self.df = DataFrame({'float1': randn(25000), 'float2': randn(25000), 'string1': (['foo'] * 25000), 'bool1': ([True] * 25000), 'int1': np.random.randint(0, 250000, size=25000), }, index=self.index)
        remove(self.f)
        self.store = HDFStore(self.f)

    def time_write_store_mixed(self):
        self.store.put('df4', self.df)

    def teardown(self):
        self.store.close()


class write_store_table(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.h5'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.index = tm.makeStringIndex(25000)
        self.df = DataFrame({'float1': randn(25000), 'float2': randn(25000), }, index=self.index)
        remove(self.f)
        self.store = HDFStore(self.f)

    def time_write_store_table(self):
        self.store.append('df8', self.df)

    def teardown(self):
        self.store.close()


class write_store_table_dc(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.h5'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.df = DataFrame(np.random.randn(10000, 10), columns=[('C%03d' % i) for i in xrange(10)])
        remove(self.f)
        self.store = HDFStore(self.f)

    def time_write_store_table_dc(self):
        self.store.append('df15', self.df, data_columns=True)

    def teardown(self):
        self.store.close()


class write_store_table_mixed(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.h5'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.index = tm.makeStringIndex(25000)
        self.df = DataFrame({'float1': randn(25000), 'float2': randn(25000), 'string1': (['foo'] * 25000), 'bool1': ([True] * 25000), 'int1': np.random.randint(0, 25000, size=25000), }, index=self.index)
        remove(self.f)
        self.store = HDFStore(self.f)

    def time_write_store_table_mixed(self):
        self.store.append('df6', self.df)

    def teardown(self):
        self.store.close()


class write_store_table_panel(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.h5'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.p = Panel(randn(20, 1000, 25), items=[('Item%03d' % i) for i in xrange(20)], major_axis=date_range('1/1/2000', periods=1000), minor_axis=[('E%03d' % i) for i in xrange(25)])
        remove(self.f)
        self.store = HDFStore(self.f)

    def time_write_store_table_panel(self):
        self.store.append('p2', self.p)

    def teardown(self):
        self.store.close()


class write_store_table_wide(object):
    goal_time = 0.2

    def setup(self):
        self.f = '__test__.h5'

        def remove(f):
            try:
                os.remove(self.f)
            except:
                pass
        self.df = DataFrame(np.random.randn(25000, 100))
        remove(self.f)
        self.store = HDFStore(self.f)

    def time_write_store_table_wide(self):
        self.store.append('df10', self.df)

    def teardown(self):
        self.store.close()