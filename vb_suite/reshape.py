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
