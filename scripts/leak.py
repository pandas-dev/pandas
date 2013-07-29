from pandas import *
from pandas.compat import range
import numpy as np
import pandas.util.testing as tm
import os
import psutil

pid = os.getpid()
proc = psutil.Process(pid)

df = DataFrame(index=np.arange(100))
for i in range(5000):
    df[i] = 5
