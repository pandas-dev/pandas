"""
Some examples playing around with yahoo finance data
"""

from datetime import datetime

import matplotlib.finance as fin
import numpy as np
from pylab import show


from pandas import Index, DataMatrix
from pandas.core.datetools import BMonthEnd

startDate = datetime(2009, 1, 1)
endDate = datetime(2009, 9, 1)

def getQuotes(symbol, start, end):
    quotes = fin.quotes_historical_yahoo(symbol, start, end)
    dates, open, close, high, low, volume = zip(*quotes)

    data = {
        'open' : open,
        'close' : close,
        'high' : high,
        'low' : low,
        'volume' : volume
    }

    dates = Index([datetime.fromordinal(int(d)) for d in dates])
    return DataMatrix(data, index=dates)

msft = getQuotes('MSFT', startDate, endDate)
ibm = getQuotes('IBM', startDate, endDate)

# Select dates

subIndex = ibm.index[(ibm['close'] > 95) & (ibm['close'] < 100)]
msftOnSameDates = msft.reindex(subIndex)

# Insert columns

msft['hi-lo spread'] = msft['high'] - msft['low']
ibm['hi-lo spread'] = ibm['high'] - ibm['low']

# Aggregate monthly

def toMonthly(frame, how):
    offset = BMonthEnd()

    return frame.groupby(offset.rollforward).aggregate(how)

msftMonthly = toMonthly(msft, np.mean)
ibmMonthly = toMonthly(ibm, np.mean)

# Statistics

stdev = DataMatrix({
    'MSFT' : msft.std(),
    'IBM'  : ibm.std()
})

# Arithmetic

ratios = ibm / msft

# Works with different indices

ratio = ibm / ibmMonthly
monthlyRatio = ratio.reindex(ibmMonthly.index)

# Ratio relative to past month average

filledRatio = ibm / ibmMonthly.reindex(ibm.index, fillMethod='pad')
