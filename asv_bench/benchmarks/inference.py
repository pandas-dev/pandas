from .pandas_vb_common import *
import pandas as pd


class dtype_infer_datetime64(object):
    goal_time = 0.2

    def setup(self):
        self.N = 500000
        self.df_int64 = DataFrame(dict(A=np.arange(self.N, dtype='int64'), B=np.arange(self.N, dtype='int64')))
        self.df_int32 = DataFrame(dict(A=np.arange(self.N, dtype='int32'), B=np.arange(self.N, dtype='int32')))
        self.df_uint32 = DataFrame(dict(A=np.arange(self.N, dtype='uint32'), B=np.arange(self.N, dtype='uint32')))
        self.df_float64 = DataFrame(dict(A=np.arange(self.N, dtype='float64'), B=np.arange(self.N, dtype='float64')))
        self.df_float32 = DataFrame(dict(A=np.arange(self.N, dtype='float32'), B=np.arange(self.N, dtype='float32')))
        self.df_datetime64 = DataFrame(dict(A=pd.to_datetime(np.arange(self.N, dtype='int64'), unit='ms'), B=pd.to_datetime(np.arange(self.N, dtype='int64'), unit='ms')))
        self.df_timedelta64 = DataFrame(dict(A=(self.df_datetime64['A'] - self.df_datetime64['B']), B=self.df_datetime64['B']))

    def time_dtype_infer_datetime64(self):
        (self.df_datetime64['A'] - self.df_datetime64['B'])


class dtype_infer_float32(object):
    goal_time = 0.2

    def setup(self):
        self.N = 500000
        self.df_int64 = DataFrame(dict(A=np.arange(self.N, dtype='int64'), B=np.arange(self.N, dtype='int64')))
        self.df_int32 = DataFrame(dict(A=np.arange(self.N, dtype='int32'), B=np.arange(self.N, dtype='int32')))
        self.df_uint32 = DataFrame(dict(A=np.arange(self.N, dtype='uint32'), B=np.arange(self.N, dtype='uint32')))
        self.df_float64 = DataFrame(dict(A=np.arange(self.N, dtype='float64'), B=np.arange(self.N, dtype='float64')))
        self.df_float32 = DataFrame(dict(A=np.arange(self.N, dtype='float32'), B=np.arange(self.N, dtype='float32')))
        self.df_datetime64 = DataFrame(dict(A=pd.to_datetime(np.arange(self.N, dtype='int64'), unit='ms'), B=pd.to_datetime(np.arange(self.N, dtype='int64'), unit='ms')))
        self.df_timedelta64 = DataFrame(dict(A=(self.df_datetime64['A'] - self.df_datetime64['B']), B=self.df_datetime64['B']))

    def time_dtype_infer_float32(self):
        (self.df_float32['A'] + self.df_float32['B'])


class dtype_infer_float64(object):
    goal_time = 0.2

    def setup(self):
        self.N = 500000
        self.df_int64 = DataFrame(dict(A=np.arange(self.N, dtype='int64'), B=np.arange(self.N, dtype='int64')))
        self.df_int32 = DataFrame(dict(A=np.arange(self.N, dtype='int32'), B=np.arange(self.N, dtype='int32')))
        self.df_uint32 = DataFrame(dict(A=np.arange(self.N, dtype='uint32'), B=np.arange(self.N, dtype='uint32')))
        self.df_float64 = DataFrame(dict(A=np.arange(self.N, dtype='float64'), B=np.arange(self.N, dtype='float64')))
        self.df_float32 = DataFrame(dict(A=np.arange(self.N, dtype='float32'), B=np.arange(self.N, dtype='float32')))
        self.df_datetime64 = DataFrame(dict(A=pd.to_datetime(np.arange(self.N, dtype='int64'), unit='ms'), B=pd.to_datetime(np.arange(self.N, dtype='int64'), unit='ms')))
        self.df_timedelta64 = DataFrame(dict(A=(self.df_datetime64['A'] - self.df_datetime64['B']), B=self.df_datetime64['B']))

    def time_dtype_infer_float64(self):
        (self.df_float64['A'] + self.df_float64['B'])


class dtype_infer_int32(object):
    goal_time = 0.2

    def setup(self):
        self.N = 500000
        self.df_int64 = DataFrame(dict(A=np.arange(self.N, dtype='int64'), B=np.arange(self.N, dtype='int64')))
        self.df_int32 = DataFrame(dict(A=np.arange(self.N, dtype='int32'), B=np.arange(self.N, dtype='int32')))
        self.df_uint32 = DataFrame(dict(A=np.arange(self.N, dtype='uint32'), B=np.arange(self.N, dtype='uint32')))
        self.df_float64 = DataFrame(dict(A=np.arange(self.N, dtype='float64'), B=np.arange(self.N, dtype='float64')))
        self.df_float32 = DataFrame(dict(A=np.arange(self.N, dtype='float32'), B=np.arange(self.N, dtype='float32')))
        self.df_datetime64 = DataFrame(dict(A=pd.to_datetime(np.arange(self.N, dtype='int64'), unit='ms'), B=pd.to_datetime(np.arange(self.N, dtype='int64'), unit='ms')))
        self.df_timedelta64 = DataFrame(dict(A=(self.df_datetime64['A'] - self.df_datetime64['B']), B=self.df_datetime64['B']))

    def time_dtype_infer_int32(self):
        (self.df_int32['A'] + self.df_int32['B'])


class dtype_infer_int64(object):
    goal_time = 0.2

    def setup(self):
        self.N = 500000
        self.df_int64 = DataFrame(dict(A=np.arange(self.N, dtype='int64'), B=np.arange(self.N, dtype='int64')))
        self.df_int32 = DataFrame(dict(A=np.arange(self.N, dtype='int32'), B=np.arange(self.N, dtype='int32')))
        self.df_uint32 = DataFrame(dict(A=np.arange(self.N, dtype='uint32'), B=np.arange(self.N, dtype='uint32')))
        self.df_float64 = DataFrame(dict(A=np.arange(self.N, dtype='float64'), B=np.arange(self.N, dtype='float64')))
        self.df_float32 = DataFrame(dict(A=np.arange(self.N, dtype='float32'), B=np.arange(self.N, dtype='float32')))
        self.df_datetime64 = DataFrame(dict(A=pd.to_datetime(np.arange(self.N, dtype='int64'), unit='ms'), B=pd.to_datetime(np.arange(self.N, dtype='int64'), unit='ms')))
        self.df_timedelta64 = DataFrame(dict(A=(self.df_datetime64['A'] - self.df_datetime64['B']), B=self.df_datetime64['B']))

    def time_dtype_infer_int64(self):
        (self.df_int64['A'] + self.df_int64['B'])


class dtype_infer_timedelta64_1(object):
    goal_time = 0.2

    def setup(self):
        self.N = 500000
        self.df_int64 = DataFrame(dict(A=np.arange(self.N, dtype='int64'), B=np.arange(self.N, dtype='int64')))
        self.df_int32 = DataFrame(dict(A=np.arange(self.N, dtype='int32'), B=np.arange(self.N, dtype='int32')))
        self.df_uint32 = DataFrame(dict(A=np.arange(self.N, dtype='uint32'), B=np.arange(self.N, dtype='uint32')))
        self.df_float64 = DataFrame(dict(A=np.arange(self.N, dtype='float64'), B=np.arange(self.N, dtype='float64')))
        self.df_float32 = DataFrame(dict(A=np.arange(self.N, dtype='float32'), B=np.arange(self.N, dtype='float32')))
        self.df_datetime64 = DataFrame(dict(A=pd.to_datetime(np.arange(self.N, dtype='int64'), unit='ms'), B=pd.to_datetime(np.arange(self.N, dtype='int64'), unit='ms')))
        self.df_timedelta64 = DataFrame(dict(A=(self.df_datetime64['A'] - self.df_datetime64['B']), B=self.df_datetime64['B']))

    def time_dtype_infer_timedelta64_1(self):
        (self.df_timedelta64['A'] + self.df_timedelta64['B'])


class dtype_infer_timedelta64_2(object):
    goal_time = 0.2

    def setup(self):
        self.N = 500000
        self.df_int64 = DataFrame(dict(A=np.arange(self.N, dtype='int64'), B=np.arange(self.N, dtype='int64')))
        self.df_int32 = DataFrame(dict(A=np.arange(self.N, dtype='int32'), B=np.arange(self.N, dtype='int32')))
        self.df_uint32 = DataFrame(dict(A=np.arange(self.N, dtype='uint32'), B=np.arange(self.N, dtype='uint32')))
        self.df_float64 = DataFrame(dict(A=np.arange(self.N, dtype='float64'), B=np.arange(self.N, dtype='float64')))
        self.df_float32 = DataFrame(dict(A=np.arange(self.N, dtype='float32'), B=np.arange(self.N, dtype='float32')))
        self.df_datetime64 = DataFrame(dict(A=pd.to_datetime(np.arange(self.N, dtype='int64'), unit='ms'), B=pd.to_datetime(np.arange(self.N, dtype='int64'), unit='ms')))
        self.df_timedelta64 = DataFrame(dict(A=(self.df_datetime64['A'] - self.df_datetime64['B']), B=self.df_datetime64['B']))

    def time_dtype_infer_timedelta64_2(self):
        (self.df_timedelta64['A'] + self.df_timedelta64['A'])


class dtype_infer_uint32(object):
    goal_time = 0.2

    def setup(self):
        self.N = 500000
        self.df_int64 = DataFrame(dict(A=np.arange(self.N, dtype='int64'), B=np.arange(self.N, dtype='int64')))
        self.df_int32 = DataFrame(dict(A=np.arange(self.N, dtype='int32'), B=np.arange(self.N, dtype='int32')))
        self.df_uint32 = DataFrame(dict(A=np.arange(self.N, dtype='uint32'), B=np.arange(self.N, dtype='uint32')))
        self.df_float64 = DataFrame(dict(A=np.arange(self.N, dtype='float64'), B=np.arange(self.N, dtype='float64')))
        self.df_float32 = DataFrame(dict(A=np.arange(self.N, dtype='float32'), B=np.arange(self.N, dtype='float32')))
        self.df_datetime64 = DataFrame(dict(A=pd.to_datetime(np.arange(self.N, dtype='int64'), unit='ms'), B=pd.to_datetime(np.arange(self.N, dtype='int64'), unit='ms')))
        self.df_timedelta64 = DataFrame(dict(A=(self.df_datetime64['A'] - self.df_datetime64['B']), B=self.df_datetime64['B']))

    def time_dtype_infer_uint32(self):
        (self.df_uint32['A'] + self.df_uint32['B'])


class to_numeric(object):

    param_names = ['dtype', 'downcast']
    params = [['string-float', 'string-int', 'string-nint', 'datetime64',
               'int-list', 'int32'],
              [None, 'integer', 'signed', 'unsigned', 'float']]

    N = 500000

    data_dict = {
        'string-int': (['1'] * (N // 2)) + ([2] * (N // 2)),
        'string-nint': (['-1'] * (N // 2)) + ([2] * (N // 2)),
        'datetime64': np.repeat(np.array(['1970-01-01', '1970-01-02'],
                                         dtype='datetime64[D]'), N),
        'string-float': (['1.1'] * (N // 2)) + ([2] * (N // 2)),
        'int-list': ([1] * (N // 2)) + ([2] * (N // 2)),
        'int32': np.repeat(np.int32(1), N)
        }

    def setup(self, dtype, downcast):
        self.data = self.data_dict[dtype]

    def time_downcast(self, dtype, downcast):
        pd.to_numeric(self.data, downcast=downcast)
