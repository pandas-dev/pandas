"""A collection of random tools for dealing with dates in Python"""

from pandas.tseries.tools import *
from pandas.tseries.offsets import *
from pandas.tseries.frequencies import *
from dateutil import parser

day = DateOffset()
bday = BDay()
businessDay = bday
try:
    cday = CDay()
    customBusinessDay = CustomBusinessDay()
except NotImplementedError:
    cday = None
    customBusinessDay = None
monthEnd = MonthEnd()
yearEnd = YearEnd()
yearBegin = YearBegin()
bmonthEnd = BMonthEnd()
businessMonthEnd = bmonthEnd
bquarterEnd = BQuarterEnd()
quarterEnd = QuarterEnd()
byearEnd = BYearEnd()
week = Week()

# Functions/offsets to roll dates forward
thisMonthEnd = MonthEnd(0)
thisBMonthEnd = BMonthEnd(0)
thisYearEnd = YearEnd(0)
thisYearBegin = YearBegin(0)
thisBQuarterEnd = BQuarterEnd(0)
thisQuarterEnd = QuarterEnd(0)

# Functions to check where a date lies
isBusinessDay = BDay().onOffset
isMonthEnd = MonthEnd().onOffset
isBMonthEnd = BMonthEnd().onOffset
