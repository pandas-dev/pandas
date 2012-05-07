from vbench.api import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
"""

setup = common_setup + """
from pandas import read_csv
import os
N = 10000
K = 8
df = DataFrame(np.random.randn(N, K) * np.random.randint(100, 10000, (N, K)))
df.to_csv('test.csv')
"""

read_csv_vb = Benchmark("read_csv('test.csv')", setup,
                        cleanup="os.remove('test.csv')",
                        start_date=datetime(2012, 5, 7))


setup = common_setup + """
from pandas import read_csv
import os
N = 10000
K = 8
format = lambda x: '{:,}'.format(x)
df = DataFrame(np.random.randn(N, K) * np.random.randint(100, 10000, (N, K)))
df = df.applymap(format)
df.to_csv('test.csv', sep='|')
"""

read_csv_thou_vb = Benchmark("read_csv('test.csv')", setup,
                             cleanup="os.remove('test.csv')",
                             start_date=datetime(2012, 5, 7))
