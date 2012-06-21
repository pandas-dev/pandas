from vbench.api import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
index = MultiIndex.from_arrays([np.arange(100).repeat(100),
                               np.roll(np.tile(np.arange(100), 100), 25)])
df = DataFrame(np.random.randn(10000, 4), index=index)
"""

reshape_unstack_simple = Benchmark('df.unstack(1)', common_setup,
                                   start_date=datetime(2011, 10, 1))

setup = common_setup + """
udf = df.unstack(1)
"""

reshape_stack_simple = Benchmark('udf.stack()', setup,
                                 start_date=datetime(2011, 10, 1))

setup = common_setup + """
def unpivot(frame):
    N, K = frame.shape
    data = {'value' : frame.values.ravel('F'),
            'variable' : np.asarray(frame.columns).repeat(N),
            'date' : np.tile(np.asarray(frame.index), K)}
    return DataFrame(data, columns=['date', 'variable', 'value'])
index = date_range('1/1/2000', periods=10000, freq='h')
df = DataFrame(randn(10000, 50), index=index, columns=range(50))
pdf = unpivot(df)
f = lambda: pdf.pivot('date', 'variable', 'value')
"""

reshape_pivot_time_series = Benchmark('f()', setup,
                                      start_date=datetime(2012, 5, 1))
