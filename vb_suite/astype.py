from vbench.api import Benchmark

common_setup = """from pandas_vb_common import *
from datetime import timedelta
import pandas as pd
import numpy as np

N = 1000000
df = pd.DataFrame({'a': 1.,
                   'b': 2,
                   'c': 'foo',
                   'float32' : np.array([1.]*N,dtype='float32'),
                   'int32' : np.array([1]*N,dtype='int32'),
                   },
                   index=np.arange(N))

mn = df._get_numeric_data()
mn['little_float'] = np.array(12345.,dtype='float16')
mn['big_float']    = np.array(123456789101112.,dtype='float64')
"""

astype_test = Benchmark('s.astype(np.int64)',
                        common_setup,
                        name='astype_test')
