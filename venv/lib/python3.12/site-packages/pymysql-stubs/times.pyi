from datetime import date, datetime, time, timedelta

Date: type[date]
Time: type[time]
TimeDelta: type[timedelta]
Timestamp: type[datetime]

def DateFromTicks(ticks: float | None) -> date: ...
def TimeFromTicks(ticks: float | None) -> time: ...
def TimestampFromTicks(ticks: float | None) -> datetime: ...
