from pandas import *

import numpy as np
import os

from vbench.api import Benchmark
from pandas.util.testing import rands
import pandas.lib as lib
import pandas._sandbox as sbx
import time

import psutil

pid = os.getpid()
proc = psutil.Process(pid)

lst = SparseList()
lst.append([5] * 10000)
lst.append(np.repeat(np.nan, 1000000))

for _ in xrange(10000):
    print proc.get_memory_info()
    sdf = SparseDataFrame({'A' : lst.to_array()})
    chunk = sdf[sdf['A'] == 5]
