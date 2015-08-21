from vbench.api import Benchmark
from datetime import datetime

common_setup = """from .pandas_vb_common import *
from pandas import read_csv, read_table
"""

setup = common_setup + """
import os
N = 10000
K = 8
df = DataFrame(np.random.randn(N, K) * np.random.randint(100, 10000, (N, K)))
df.to_csv('test.csv', sep='|')
"""

read_csv_vb = Benchmark("read_csv('test.csv', sep='|')", setup,
                        cleanup="os.remove('test.csv')",
                        start_date=datetime(2012, 5, 7))


setup = common_setup + """
import os
N = 10000
K = 8
format = lambda x: '{:,}'.format(x)
df = DataFrame(np.random.randn(N, K) * np.random.randint(100, 10000, (N, K)))
df = df.applymap(format)
df.to_csv('test.csv', sep='|')
"""

read_csv_thou_vb = Benchmark("read_csv('test.csv', sep='|', thousands=',')",
                             setup,
                             cleanup="os.remove('test.csv')",
                             start_date=datetime(2012, 5, 7))

setup = common_setup + """
data = ['A,B,C']
data = data + ['1,2,3 # comment'] * 100000
data = '\\n'.join(data)
"""

stmt = "read_csv(StringIO(data), comment='#')"
read_csv_comment2 = Benchmark(stmt, setup,
                              start_date=datetime(2011, 11, 1))

setup = common_setup + """
try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

import os
N = 10000
K = 8
data = '''\
KORD,19990127, 19:00:00, 18:56:00, 0.8100, 2.8100, 7.2000, 0.0000, 280.0000
KORD,19990127, 20:00:00, 19:56:00, 0.0100, 2.2100, 7.2000, 0.0000, 260.0000
KORD,19990127, 21:00:00, 20:56:00, -0.5900, 2.2100, 5.7000, 0.0000, 280.0000
KORD,19990127, 21:00:00, 21:18:00, -0.9900, 2.0100, 3.6000, 0.0000, 270.0000
KORD,19990127, 22:00:00, 21:56:00, -0.5900, 1.7100, 5.1000, 0.0000, 290.0000
'''
data = data * 200
"""
cmd = ("read_table(StringIO(data), sep=',', header=None, "
       "parse_dates=[[1,2], [1,3]])")
sdate = datetime(2012, 5, 7)
read_table_multiple_date = Benchmark(cmd, setup, start_date=sdate)

setup = common_setup + """
try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

import os
N = 10000
K = 8
data = '''\
KORD,19990127 19:00:00, 18:56:00, 0.8100, 2.8100, 7.2000, 0.0000, 280.0000
KORD,19990127 20:00:00, 19:56:00, 0.0100, 2.2100, 7.2000, 0.0000, 260.0000
KORD,19990127 21:00:00, 20:56:00, -0.5900, 2.2100, 5.7000, 0.0000, 280.0000
KORD,19990127 21:00:00, 21:18:00, -0.9900, 2.0100, 3.6000, 0.0000, 270.0000
KORD,19990127 22:00:00, 21:56:00, -0.5900, 1.7100, 5.1000, 0.0000, 290.0000
'''
data = data * 200
"""
cmd = "read_table(StringIO(data), sep=',', header=None, parse_dates=[1])"
sdate = datetime(2012, 5, 7)
read_table_multiple_date_baseline = Benchmark(cmd, setup, start_date=sdate)

setup = common_setup + """
try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

data = '''\
0.1213700904466425978256438611,0.0525708283766902484401839501,0.4174092731488769913994474336
0.4096341697147408700274695547,0.1587830198973579909349496119,0.1292545832485494372576795285
0.8323255650024565799327547210,0.9694902427379478160318626578,0.6295047811546814475747169126
0.4679375305798131323697930383,0.2963942381834381301075609371,0.5268936082160610157032465394
0.6685382761849776311890991564,0.6721207066140679753374342908,0.6519975277021627935170045020
'''
data = data * 200
"""
cmd = "read_csv(StringIO(data), sep=',', header=None, float_precision=None)"
sdate = datetime(2014, 8, 20)
read_csv_default_converter = Benchmark(cmd, setup, start_date=sdate)
cmd = "read_csv(StringIO(data), sep=',', header=None, float_precision='high')"
read_csv_precise_converter = Benchmark(cmd, setup, start_date=sdate)
cmd = "read_csv(StringIO(data), sep=',', header=None, float_precision='round_trip')"
read_csv_roundtrip_converter = Benchmark(cmd, setup, start_date=sdate)
