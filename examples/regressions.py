from datetime import datetime
import string

import numpy as np

from pandas.core.api import Series, DataFrame, DateRange
from pandas.stats.api import ols

N = 100

start = datetime(2009, 9, 2)
dateRange = DateRange(start, periods=N)

def makeDataFrame():
    data = DataFrame(np.random.randn(N, 7),
                      columns=list(string.ascii_uppercase[:7]),
                      index=dateRange)

    return data

def makeSeries():
    return Series(np.random.randn(N), index=dateRange)

#-------------------------------------------------------------------------------
# Standard rolling linear regression

X = makeDataFrame()
Y =  makeSeries()

model = ols(y=Y, x=X)

print model

#-------------------------------------------------------------------------------
# Panel regression

data = {
    'A' : makeDataFrame(),
    'B' : makeDataFrame(),
    'C' : makeDataFrame()
}

Y = makeDataFrame()

panelModel = ols(y=Y, x=data, window=50)

model = ols(y=Y, x=data)

print panelModel
