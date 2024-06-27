from time import localtime
from datetime import date, datetime, time, timedelta


Date = date
Time = time
TimeDelta = timedelta
Timestamp = datetime


def DateFromTicks(ticks):
    return date(*localtime(ticks)[:3])


def TimeFromTicks(ticks):
    return time(*localtime(ticks)[3:6])


def TimestampFromTicks(ticks):
    return datetime(*localtime(ticks)[:6])
