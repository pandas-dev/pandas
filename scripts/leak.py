from pandas import *
import numpy as np

import pandas.util.testing as tm

tm.N = 2000
tm.K = 25

for i in xrange(100):
    print i
    df = tm.makeTimeDataFrame()
    y = df.pop('A')
    model = ols(y=y, x=df, window=1999).beta

