from vbench.api import Benchmark
from datetime import datetime

common_setup = """from .pandas_vb_common import *
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

# Sparse key space, re: #2278

setup = common_setup + """
NUM_ROWS = 1000
for iter in range(10):
    df = DataFrame({'A' : np.random.randint(50, size=NUM_ROWS),
                    'B' : np.random.randint(50, size=NUM_ROWS),
                    'C' : np.random.randint(-10,10, size=NUM_ROWS),
                    'D' : np.random.randint(-10,10, size=NUM_ROWS),
                    'E' : np.random.randint(10, size=NUM_ROWS),
                    'F' : np.random.randn(NUM_ROWS)})
    idf = df.set_index(['A', 'B', 'C', 'D', 'E'])
    if len(idf.index.unique()) == NUM_ROWS:
        break
"""

unstack_sparse_keyspace = Benchmark('idf.unstack()', setup,
                                    start_date=datetime(2011, 10, 1))

# Melt

setup = common_setup + """
from pandas.core.reshape import melt
df = DataFrame(np.random.randn(10000, 3), columns=['A', 'B', 'C'])
df['id1'] = np.random.randint(0, 10, 10000)
df['id2'] = np.random.randint(100, 1000, 10000)
"""

melt_dataframe = Benchmark("melt(df, id_vars=['id1', 'id2'])", setup,
                           start_date=datetime(2012, 8, 1))
