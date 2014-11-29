from vbench.api import Benchmark

common_setup = """from pandas_vb_common import *
from datetime import timedelta
import numpy as np

N = 1000000
arr = np.random.randint(1,10,size=1000000)
s = pd.Series(arr)
"""

astype_test = Benchmark('s.astype(np.str)',
                        common_setup,
                        name='astype_test')
