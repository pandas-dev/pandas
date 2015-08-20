from pandas_vb_common import *
from pandas import concat, Timestamp
from StringIO import StringIO


class frame_to_csv(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(3000, 30))

    def time_frame_to_csv(self):
        self.df.to_csv('__test__.csv')


class frame_to_csv2(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame({'A': range(50000), })
        self.df['B'] = (self.df.A + 1.0)
        self.df['C'] = (self.df.A + 2.0)
        self.df['D'] = (self.df.A + 3.0)

    def time_frame_to_csv2(self):
        self.df.to_csv('__test__.csv')


class frame_to_csv_date_formatting(object):
    goal_time = 0.2

    def setup(self):
        self.rng = date_range('1/1/2000', periods=1000)
        self.data = DataFrame(self.rng, index=self.rng)

    def time_frame_to_csv_date_formatting(self):
        self.data.to_csv('__test__.csv', date_format='%Y%m%d')


class frame_to_csv_mixed(object):
    goal_time = 0.2

    def setup(self):

        def create_cols(name):
            return [('%s%03d' % (name, i)) for i in xrange(5)]
        self.df_float = DataFrame(np.random.randn(5000, 5), dtype='float64', columns=create_cols('float'))
        self.df_int = DataFrame(np.random.randn(5000, 5), dtype='int64', columns=create_cols('int'))
        self.df_bool = DataFrame(True, index=self.df_float.index, columns=create_cols('bool'))
        self.df_object = DataFrame('foo', index=self.df_float.index, columns=create_cols('object'))
        self.df_dt = DataFrame(Timestamp('20010101'), index=self.df_float.index, columns=create_cols('date'))
        self.df_float.ix[30:500, 1:3] = np.nan
        self.df = concat([self.df_float, self.df_int, self.df_bool, self.df_object, self.df_dt], axis=1)

    def time_frame_to_csv_mixed(self):
        self.df.to_csv('__test__.csv')


class read_csv_infer_datetime_format_custom(object):
    goal_time = 0.2

    def setup(self):
        self.rng = date_range('1/1/2000', periods=1000)
        self.data = '\n'.join(self.rng.map((lambda x: x.strftime('%m/%d/%Y %H:%M:%S.%f'))))

    def time_read_csv_infer_datetime_format_custom(self):
        read_csv(StringIO(self.data), header=None, names=['foo'], parse_dates=['foo'], infer_datetime_format=True)


class read_csv_infer_datetime_format_iso8601(object):
    goal_time = 0.2

    def setup(self):
        self.rng = date_range('1/1/2000', periods=1000)
        self.data = '\n'.join(self.rng.map((lambda x: x.strftime('%Y-%m-%d %H:%M:%S'))))

    def time_read_csv_infer_datetime_format_iso8601(self):
        read_csv(StringIO(self.data), header=None, names=['foo'], parse_dates=['foo'], infer_datetime_format=True)


class read_csv_infer_datetime_format_ymd(object):
    goal_time = 0.2

    def setup(self):
        self.rng = date_range('1/1/2000', periods=1000)
        self.data = '\n'.join(self.rng.map((lambda x: x.strftime('%Y%m%d'))))

    def time_read_csv_infer_datetime_format_ymd(self):
        read_csv(StringIO(self.data), header=None, names=['foo'], parse_dates=['foo'], infer_datetime_format=True)


class read_csv_skiprows(object):
    goal_time = 0.2

    def setup(self):
        self.index = tm.makeStringIndex(20000)
        self.df = DataFrame({'float1': randn(20000), 'float2': randn(20000), 'string1': (['foo'] * 20000), 'bool1': ([True] * 20000), 'int1': np.random.randint(0, 200000, size=20000), }, index=self.index)
        self.df.to_csv('__test__.csv')

    def time_read_csv_skiprows(self):
        read_csv('__test__.csv', skiprows=10000)


class read_csv_standard(object):
    goal_time = 0.2

    def setup(self):
        self.index = tm.makeStringIndex(10000)
        self.df = DataFrame({'float1': randn(10000), 'float2': randn(10000), 'string1': (['foo'] * 10000), 'bool1': ([True] * 10000), 'int1': np.random.randint(0, 100000, size=10000), }, index=self.index)
        self.df.to_csv('__test__.csv')

    def time_read_csv_standard(self):
        read_csv('__test__.csv')


class read_parse_dates_iso8601(object):
    goal_time = 0.2

    def setup(self):
        self.rng = date_range('1/1/2000', periods=1000)
        self.data = '\n'.join(self.rng.map((lambda x: x.strftime('%Y-%m-%d %H:%M:%S'))))

    def time_read_parse_dates_iso8601(self):
        read_csv(StringIO(self.data), header=None, names=['foo'], parse_dates=['foo'])


class write_csv_standard(object):
    goal_time = 0.2

    def setup(self):
        self.index = tm.makeStringIndex(10000)
        self.df = DataFrame({'float1': randn(10000), 'float2': randn(10000), 'string1': (['foo'] * 10000), 'bool1': ([True] * 10000), 'int1': np.random.randint(0, 100000, size=10000), }, index=self.index)

    def time_write_csv_standard(self):
        self.df.to_csv('__test__.csv')