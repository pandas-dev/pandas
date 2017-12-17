import timeit

import numpy as np
import pandas.util.testing as tm
from pandas import DataFrame, date_range, read_csv
from pandas.compat import PY2, StringIO

from .pandas_vb_common import setup, BaseIO  # noqa


class ToCSV(BaseIO):

    goal_time = 0.2
    fname = '__test__.csv'
    params = ['wide', 'long', 'mixed']
    param_names = ['kind']

    def setup(self, kind):
        wide_frame = DataFrame(np.random.randn(3000, 30))
        long_frame = DataFrame({'A': np.arange(50000),
                                'B': np.arange(50000) + 1.,
                                'C': np.arange(50000) + 2.,
                                'D': np.arange(50000) + 3.})
        mixed_frame = DataFrame({'float': np.random.randn(5000),
                                 'int': np.random.randn(5000).astype(int),
                                 'bool': (np.arange(5000) % 2) == 0,
                                 'datetime': date_range('2001',
                                                        freq='s',
                                                        periods=5000),
                                 'object': ['foo'] * 5000})
        mixed_frame.loc[30:500, 'float'] = np.nan
        data = {'wide': wide_frame,
                'long': long_frame,
                'mixed': mixed_frame}
        self.df = data[kind]

    def time_frame(self, kind):
        self.df.to_csv(self.fname)


class ToCSVDatetime(BaseIO):

    goal_time = 0.2
    fname = '__test__.csv'

    def setup(self):
        rng = date_range('1/1/2000', periods=1000)
        self.data = DataFrame(rng, index=rng)

    def time_frame_date_formatting(self):
        self.data.to_csv(self.fname, date_format='%Y%m%d')


class ReadCSVDInferDatetimeFormat(object):

    goal_time = 0.2
    params = ([True, False], ['custom', 'iso8601', 'ymd'])
    param_names = ['infer_datetime_format', 'format']

    def setup(self, infer_datetime_format, format):
        rng = date_range('1/1/2000', periods=1000)
        formats = {'custom': '%m/%d/%Y %H:%M:%S.%f',
                   'iso8601': '%Y-%m-%d %H:%M:%S',
                   'ymd': '%Y%m%d'}
        dt_format = formats[format]
        self.data = StringIO('\n'.join(rng.strftime(dt_format).tolist()))

    def time_read_csv(self, infer_datetime_format, format):
        read_csv(self.data, header=None, names=['foo'], parse_dates=['foo'],
                 infer_datetime_format=infer_datetime_format)


class ReadCSVSkipRows(BaseIO):

    goal_time = 0.2
    fname = '__test__.csv'
    params = [None, 10000]
    param_names = ['skiprows']

    def setup(self, skiprows):
        N = 20000
        index = tm.makeStringIndex(N)
        df = DataFrame({'float1': np.random.randn(N),
                        'float2': np.random.randn(N),
                        'string1': ['foo'] * N,
                        'bool1': [True] * N,
                        'int1': np.random.randint(0, N, size=N)},
                       index=index)
        df.to_csv(self.fname)

    def time_skipprows(self, skiprows):
        read_csv(self.fname, skiprows=skiprows)


class ReadUint64Integers(object):

    goal_time = 0.2

    def setup(self):
        self.na_values = [2**63 + 500]
        arr = np.arange(10000).astype('uint64') + 2**63
        self.data1 = StringIO('\n'.join(arr.astype(str).tolist()))
        arr = arr.astype(object)
        arr[500] = -1
        self.data2 = StringIO('\n'.join(arr.astype(str).tolist()))

    def time_read_uint64(self):
        read_csv(self.data1, header=None, names=['foo'])

    def time_read_uint64_neg_values(self):
        read_csv(self.data2, header=None, names=['foo'])

    def time_read_uint64_na_values(self):
        read_csv(self.data1, header=None, names=['foo'],
                 na_values=self.na_values)


class S3(object):
    # Make sure that we can read part of a file from S3 without
    # needing to download the entire thing. Use the timeit.default_timer
    # to measure wall time instead of CPU time -- we want to see
    # how long it takes to download the data.
    timer = timeit.default_timer
    params = ([None, "gzip", "bz2"], ["python", "c"])
    param_names = ["compression", "engine"]

    def setup(self, compression, engine):
        if compression == "bz2" and engine == "c" and PY2:
            # The Python 2 C parser can't read bz2 from open files.
            raise NotImplementedError
        try:
            import s3fs
        except ImportError:
            # Skip these benchmarks if `boto` is not installed.
            raise NotImplementedError

        ext = ""
        if compression == "gzip":
            ext = ".gz"
        elif compression == "bz2":
            ext = ".bz2"
        self.big_fname = "s3://pandas-test/large_random.csv" + ext

    def time_read_csv_10_rows(self, compression, engine):
        # Read a small number of rows from a huge (100,000 x 50) table.
        read_csv(self.big_fname, nrows=10, compression=compression,
                 engine=engine)
