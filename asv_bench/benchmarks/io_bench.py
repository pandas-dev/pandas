import os
from .pandas_vb_common import *
from pandas import concat, Timestamp, compat
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import timeit


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
        self.df_float = DataFrame(np.random.randn(5000, 5), dtype='float64', columns=self.create_cols('float'))
        self.df_int = DataFrame(np.random.randn(5000, 5), dtype='int64', columns=self.create_cols('int'))
        self.df_bool = DataFrame(True, index=self.df_float.index, columns=self.create_cols('bool'))
        self.df_object = DataFrame('foo', index=self.df_float.index, columns=self.create_cols('object'))
        self.df_dt = DataFrame(Timestamp('20010101'), index=self.df_float.index, columns=self.create_cols('date'))
        self.df_float.ix[30:500, 1:3] = np.nan
        self.df = concat([self.df_float, self.df_int, self.df_bool, self.df_object, self.df_dt], axis=1)

    def time_frame_to_csv_mixed(self):
        self.df.to_csv('__test__.csv')

    def create_cols(self, name):
        return [('%s%03d' % (name, i)) for i in range(5)]


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


class read_uint64_integers(object):
    goal_time = 0.2

    def setup(self):
        self.na_values = [2**63 + 500]

        self.arr1 = np.arange(10000).astype('uint64') + 2**63
        self.data1 = '\n'.join(map(lambda x: str(x), self.arr1))

        self.arr2 = self.arr1.copy().astype(object)
        self.arr2[500] = -1
        self.data2 = '\n'.join(map(lambda x: str(x), self.arr2))

    def time_read_uint64(self):
        read_csv(StringIO(self.data1), header=None)

    def time_read_uint64_neg_values(self):
        read_csv(StringIO(self.data2), header=None)

    def time_read_uint64_na_values(self):
        read_csv(StringIO(self.data1), header=None, na_values=self.na_values)


class write_csv_standard(object):
    goal_time = 0.2

    def setup(self):
        self.index = tm.makeStringIndex(10000)
        self.df = DataFrame({'float1': randn(10000), 'float2': randn(10000), 'string1': (['foo'] * 10000), 'bool1': ([True] * 10000), 'int1': np.random.randint(0, 100000, size=10000), }, index=self.index)

    def time_write_csv_standard(self):
        self.df.to_csv('__test__.csv')


class read_csv_from_s3(object):
    # Make sure that we can read part of a file from S3 without
    # needing to download the entire thing. Use the timeit.default_timer
    # to measure wall time instead of CPU time -- we want to see
    # how long it takes to download the data.
    timer = timeit.default_timer
    params = ([None, "gzip", "bz2"], ["python", "c"])
    param_names = ["compression", "engine"]

    def setup(self, compression, engine):
        if compression == "bz2" and engine == "c" and compat.PY2:
            # The Python 2 C parser can't read bz2 from open files.
            raise NotImplementedError
        try:
            import s3fs
        except ImportError:
            # Skip these benchmarks if `boto` is not installed.
            raise NotImplementedError

        self.big_fname = "s3://pandas-test/large_random.csv"

    def time_read_nrows(self, compression, engine):
        # Read a small number of rows from a huge (100,000 x 50) table.
        ext = ""
        if compression == "gzip":
            ext = ".gz"
        elif compression == "bz2":
            ext = ".bz2"
        pd.read_csv(self.big_fname + ext, nrows=10,
                    compression=compression, engine=engine)


class read_json_lines(object):
    goal_time = 0.2
    fname = "__test__.json"

    def setup(self):
        self.N = 100000
        self.C = 5
        self.df = DataFrame({('float{0}'.format(i), randn(self.N)) for i in range(self.C)})
        self.df.to_json(self.fname,orient="records",lines=True)

    def teardown(self):
        try:
            os.remove(self.fname)
        except:
            pass

    def time_read_json_lines(self):
        pd.read_json(self.fname, lines=True)

    def time_read_json_lines_chunk(self):
        pd.concat(pd.read_json(self.fname, lines=True, chunksize=self.N//4))

    def peakmem_read_json_lines(self):
        pd.read_json(self.fname, lines=True)

    def peakmem_read_json_lines_chunk(self):
        pd.concat(pd.read_json(self.fname, lines=True, chunksize=self.N//4))
