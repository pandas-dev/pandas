import random
import string

import numpy as np
import pandas.util.testing as tm
from pandas import DataFrame, Categorical, date_range, read_csv
from pandas.compat import cStringIO as StringIO

from ..pandas_vb_common import BaseIO


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


class StringIORewind(object):

    def data(self, stringio_object):
        stringio_object.seek(0)
        return stringio_object


class ReadCSVDInferDatetimeFormat(StringIORewind):

    goal_time = 0.2
    params = ([True, False], ['custom', 'iso8601', 'ymd'])
    param_names = ['infer_datetime_format', 'format']

    def setup(self, infer_datetime_format, format):
        rng = date_range('1/1/2000', periods=1000)
        formats = {'custom': '%m/%d/%Y %H:%M:%S.%f',
                   'iso8601': '%Y-%m-%d %H:%M:%S',
                   'ymd': '%Y%m%d'}
        dt_format = formats[format]
        self.StringIO_input = StringIO('\n'.join(
                                       rng.strftime(dt_format).tolist()))

    def time_read_csv(self, infer_datetime_format, format):
        read_csv(self.data(self.StringIO_input),
                 header=None, names=['foo'], parse_dates=['foo'],
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


class ReadUint64Integers(StringIORewind):

    goal_time = 0.2

    def setup(self):
        self.na_values = [2**63 + 500]
        arr = np.arange(10000).astype('uint64') + 2**63
        self.data1 = StringIO('\n'.join(arr.astype(str).tolist()))
        arr = arr.astype(object)
        arr[500] = -1
        self.data2 = StringIO('\n'.join(arr.astype(str).tolist()))

    def time_read_uint64(self):
        read_csv(self.data(self.data1), header=None, names=['foo'])

    def time_read_uint64_neg_values(self):
        read_csv(self.data(self.data2), header=None, names=['foo'])

    def time_read_uint64_na_values(self):
        read_csv(self.data(self.data1), header=None, names=['foo'],
                 na_values=self.na_values)


class ReadCSVThousands(BaseIO):

    goal_time = 0.2
    fname = '__test__.csv'
    params = ([',', '|'], [None, ','])
    param_names = ['sep', 'thousands']

    def setup(self, sep, thousands):
        N = 10000
        K = 8
        data = np.random.randn(N, K) * np.random.randint(100, 10000, (N, K))
        df = DataFrame(data)
        if thousands is not None:
            fmt = ':{}'.format(thousands)
            fmt = '{' + fmt + '}'
            df = df.applymap(lambda x: fmt.format(x))
        df.to_csv(self.fname, sep=sep)

    def time_thousands(self, sep, thousands):
        read_csv(self.fname, sep=sep, thousands=thousands)


class ReadCSVComment(StringIORewind):

    goal_time = 0.2

    def setup(self):
        data = ['A,B,C'] + (['1,2,3 # comment'] * 100000)
        self.StringIO_input = StringIO('\n'.join(data))

    def time_comment(self):
        read_csv(self.data(self.StringIO_input), comment='#',
                 header=None, names=list('abc'))


class ReadCSVFloatPrecision(StringIORewind):

    goal_time = 0.2
    params = ([',', ';'], ['.', '_'], [None, 'high', 'round_trip'])
    param_names = ['sep', 'decimal', 'float_precision']

    def setup(self, sep, decimal, float_precision):
        floats = [''.join(random.choice(string.digits) for _ in range(28))
                  for _ in range(15)]
        rows = sep.join(['0{}'.format(decimal) + '{}'] * 3) + '\n'
        data = rows * 5
        data = data.format(*floats) * 200  # 1000 x 3 strings csv
        self.StringIO_input = StringIO(data)

    def time_read_csv(self, sep, decimal, float_precision):
        read_csv(self.data(self.StringIO_input), sep=sep, header=None,
                 names=list('abc'), float_precision=float_precision)

    def time_read_csv_python_engine(self, sep, decimal, float_precision):
        read_csv(self.data(self.StringIO_input), sep=sep, header=None,
                 engine='python', float_precision=None, names=list('abc'))


class ReadCSVCategorical(BaseIO):

    goal_time = 0.2
    fname = '__test__.csv'

    def setup(self):
        N = 100000
        group1 = ['aaaaaaaa', 'bbbbbbb', 'cccccccc', 'dddddddd', 'eeeeeeee']
        df = DataFrame(np.random.choice(group1, (N, 3)), columns=list('abc'))
        df.to_csv(self.fname, index=False)

    def time_convert_post(self):
        read_csv(self.fname).apply(Categorical)

    def time_convert_direct(self):
        read_csv(self.fname, dtype='category')


class ReadCSVParseDates(StringIORewind):

    goal_time = 0.2

    def setup(self):
        data = """{},19:00:00,18:56:00,0.8100,2.8100,7.2000,0.0000,280.0000\n
                  {},20:00:00,19:56:00,0.0100,2.2100,7.2000,0.0000,260.0000\n
                  {},21:00:00,20:56:00,-0.5900,2.2100,5.7000,0.0000,280.0000\n
                  {},21:00:00,21:18:00,-0.9900,2.0100,3.6000,0.0000,270.0000\n
                  {},22:00:00,21:56:00,-0.5900,1.7100,5.1000,0.0000,290.0000\n
               """
        two_cols = ['KORD,19990127'] * 5
        data = data.format(*two_cols)
        self.StringIO_input = StringIO(data)

    def time_multiple_date(self):
        read_csv(self.data(self.StringIO_input), sep=',', header=None,
                 names=list(string.digits[:9]),
                 parse_dates=[[1, 2], [1, 3]])

    def time_baseline(self):
        read_csv(self.data(self.StringIO_input), sep=',', header=None,
                 parse_dates=[1],
                 names=list(string.digits[:9]))


from ..pandas_vb_common import setup  # noqa: F401
