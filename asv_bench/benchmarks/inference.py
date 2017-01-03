from .pandas_vb_common import *
import pandas as pd


class DtypeInfer(object):
    goal_time = 0.2

    # from GH 7332

    def setup(self):
        self.N = 500000
        self.df_int64 = DataFrame(dict(A=np.arange(self.N, dtype='int64'),
                                       B=np.arange(self.N, dtype='int64')))
        self.df_int32 = DataFrame(dict(A=np.arange(self.N, dtype='int32'),
                                       B=np.arange(self.N, dtype='int32')))
        self.df_uint32 = DataFrame(dict(A=np.arange(self.N, dtype='uint32'),
                                        B=np.arange(self.N, dtype='uint32')))
        self.df_float64 = DataFrame(dict(A=np.arange(self.N, dtype='float64'),
                                         B=np.arange(self.N, dtype='float64')))
        self.df_float32 = DataFrame(dict(A=np.arange(self.N, dtype='float32'),
                                         B=np.arange(self.N, dtype='float32')))
        self.df_datetime64 = DataFrame(dict(A=pd.to_datetime(np.arange(self.N, dtype='int64'), unit='ms'),
                                            B=pd.to_datetime(np.arange(self.N, dtype='int64'), unit='ms')))
        self.df_timedelta64 = DataFrame(dict(A=(self.df_datetime64['A'] - self.df_datetime64['B']),
                                             B=self.df_datetime64['B']))

    def time_int64(self):
        (self.df_int64['A'] + self.df_int64['B'])

    def time_int32(self):
        (self.df_int32['A'] + self.df_int32['B'])

    def time_uint32(self):
        (self.df_uint32['A'] + self.df_uint32['B'])

    def time_float64(self):
        (self.df_float64['A'] + self.df_float64['B'])

    def time_float32(self):
        (self.df_float32['A'] + self.df_float32['B'])

    def time_datetime64(self):
        (self.df_datetime64['A'] - self.df_datetime64['B'])

    def time_timedelta64_1(self):
        (self.df_timedelta64['A'] + self.df_timedelta64['B'])

    def time_timedelta64_2(self):
        (self.df_timedelta64['A'] + self.df_timedelta64['A'])


class to_numeric(object):
    goal_time = 0.2

    def setup(self):
        self.n = 10000
        self.float = Series(np.random.randn(self.n * 100))
        self.numstr = self.float.astype('str')
        self.str = Series(tm.makeStringIndex(self.n))

    def time_from_float(self):
        pd.to_numeric(self.float)

    def time_from_numeric_str(self):
        pd.to_numeric(self.numstr)

    def time_from_str_ignore(self):
        pd.to_numeric(self.str, errors='ignore')

    def time_from_str_coerce(self):
        pd.to_numeric(self.str, errors='coerce')


class to_numeric_downcast(object):

    param_names = ['dtype', 'downcast']
    params = [['string-float', 'string-int', 'string-nint', 'datetime64',
               'int-list', 'int32'],
              [None, 'integer', 'signed', 'unsigned', 'float']]

    N = 500000
    N2 = int(N / 2)

    data_dict = {
        'string-int': (['1'] * N2) + ([2] * N2),
        'string-nint': (['-1'] * N2) + ([2] * N2),
        'datetime64': np.repeat(np.array(['1970-01-01', '1970-01-02'],
                                         dtype='datetime64[D]'), N),
        'string-float': (['1.1'] * N2) + ([2] * N2),
        'int-list': ([1] * N2) + ([2] * N2),
        'int32': np.repeat(np.int32(1), N)
        }

    def setup(self, dtype, downcast):
        self.data = self.data_dict[dtype]

    def time_downcast(self, dtype, downcast):
        pd.to_numeric(self.data, downcast=downcast)


class MaybeConvertNumeric(object):

    def setup(self):
        n = 1000000
        arr = np.repeat([2**63], n)
        arr = arr + np.arange(n).astype('uint64')
        arr = np.array([arr[i] if i%2 == 0 else
                        str(arr[i]) for i in range(n)],
                       dtype=object)

        arr[-1] = -1
        self.data = arr
        self.na_values = set()

    def time_convert(self):
        pd.lib.maybe_convert_numeric(self.data, self.na_values,
                                     coerce_numeric=False)
