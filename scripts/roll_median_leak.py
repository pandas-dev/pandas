from pandas import *

import numpy as np
import os

from vbench.api import Benchmark
from pandas.util.testing import rands
import pandas._tseries as lib
import pandas._sandbox as sbx
import time

import psutil

pid = os.getpid()
proc = psutil.Process(pid)


s = Series(np.random.randn(10000))

for _ in xrange(1000):
    # print proc.get_memory_info()
    result = rolling_median(s, 1000)
