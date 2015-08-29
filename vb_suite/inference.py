from vbench.api import Benchmark
from datetime import datetime
import sys

# from GH 7332

setup = """from .pandas_vb_common import *
import pandas as pd
N = 500000
df_int64 = DataFrame(dict(A = np.arange(N,dtype='int64'), B = np.arange(N,dtype='int64')))
df_int32 = DataFrame(dict(A = np.arange(N,dtype='int32'), B = np.arange(N,dtype='int32')))
df_uint32 = DataFrame(dict(A = np.arange(N,dtype='uint32'), B = np.arange(N,dtype='uint32')))
df_float64 = DataFrame(dict(A = np.arange(N,dtype='float64'), B = np.arange(N,dtype='float64')))
df_float32 = DataFrame(dict(A = np.arange(N,dtype='float32'), B = np.arange(N,dtype='float32')))
df_datetime64 = DataFrame(dict(A = pd.to_datetime(np.arange(N,dtype='int64'),unit='ms'),
                               B = pd.to_datetime(np.arange(N,dtype='int64'),unit='ms')))
df_timedelta64 = DataFrame(dict(A = df_datetime64['A']-df_datetime64['B'],
                                B = df_datetime64['B']))
"""

dtype_infer_int64 = Benchmark('df_int64["A"] + df_int64["B"]', setup,
                               start_date=datetime(2014, 1, 1))
dtype_infer_int32 = Benchmark('df_int32["A"] + df_int32["B"]', setup,
                               start_date=datetime(2014, 1, 1))
dtype_infer_uint32 = Benchmark('df_uint32["A"] + df_uint32["B"]', setup,
                               start_date=datetime(2014, 1, 1))
dtype_infer_float64 = Benchmark('df_float64["A"] + df_float64["B"]', setup,
                               start_date=datetime(2014, 1, 1))
dtype_infer_float32 = Benchmark('df_float32["A"] + df_float32["B"]', setup,
                               start_date=datetime(2014, 1, 1))
dtype_infer_datetime64 = Benchmark('df_datetime64["A"] - df_datetime64["B"]', setup,
                               start_date=datetime(2014, 1, 1))
dtype_infer_timedelta64_1 = Benchmark('df_timedelta64["A"] + df_timedelta64["B"]', setup,
                               start_date=datetime(2014, 1, 1))
dtype_infer_timedelta64_2 = Benchmark('df_timedelta64["A"] + df_timedelta64["A"]', setup,
                               start_date=datetime(2014, 1, 1))
