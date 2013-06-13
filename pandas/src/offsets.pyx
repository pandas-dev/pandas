
ctypedef enum time_res:
    r_min = 0
    r_microsecond
    r_second
    r_minute
    r_hour
    r_day
    r_month
    r_year
    r_max = 98
    r_invalid = 99


cdef conversion_factor(time_res res1, time_res res2):
    cdef:
        time_res min_res, max_res
        int64_t factor

    min_res = min(res1, res2)
    max_res = max(res1, res2)
    factor = 1

    if min_res == max_res:
        return factor

    while min_res < max_res:
        if min_res < r_microsecond:
            raise "Cannot convert from less than us"
        elif min_res == r_microsecond:
            factor *= 1000000
            min_res = r_second
        elif min_res == r_second:
            factor *= 60
            min_res = r_minute
        elif min_res == r_minute:
            factor *= 60
            min_res = r_hour
        elif min_res == r_hour:
            factor *= 24
            min_res = r_day
        else:
            raise "Cannot convert to month or year"

    return factor

# Logic to generate ranges
# -----------------------------------------------------------------------------

cdef inline int64_t weekend_adjustment(int64_t dow, int bkwd):
    if dow > 4:                         # sat or sun?
        if bkwd:                        # roll back 1 or 2 days
            return (4 - dow)
        else:                           # roll forward 2 or 1 days
            return (7 - dow)
    return 0

cdef int64_t us_in_day = conversion_factor(r_microsecond, r_day)

cdef class _Offset:
    """
    Base class to generate timestamps. Set the anchor, and then move offsets
    with next & prev. Retrieve timestamp with ts attribute.
    """
    cdef:
        int64_t t, dow, biz, dayoffset
        object start
        _TSObject ts

    def __cinit__(self):
        self.t=0
        self.dow=0
        self.biz=0
        self.dayoffset=0

    cpdef anchor(self, object start=None):
        if start is not None:
            self.start = start
        self.ts = convert_to_tsobject(self.start, None, None)
        self._setup()

    cdef _setup(self):
        pass

    cpdef next(self):
        pass

    cpdef prev(self):
        pass

    cdef int64_t _ts(self):
        """
        Access the current timestamp value, with a possible weekday
        adjustment.
        """
        cdef int64_t adj

        if self.biz != 0:
            adj = weekend_adjustment(self.dow, self.biz < 0)
            return self.t + us_in_day * adj
        else:
            return self.t

    cdef int64_t _get_anchor(self):
        """
        Retrieve an anchor relating to current offset we're on.
        """
        return self.t - self.dayoffset * us_in_day

    property ts:
        def __get__(self):
            return self._ts()

cdef class YearOffset(_Offset):
    """
    Generate annual timestamps from provided start time; apply dayoffset to
    each timestamp. If biz > 0, we choose the next business day at each time;
    previous if < 0.

    Parameters
    ----------
    dayoffset : int
    biz : int
    """
    cdef:
        int64_t y, ly

    def __init__(self, int64_t dayoffset=0, int64_t biz=0, object anchor=None):
        self.dayoffset = dayoffset
        self.biz = biz

        if anchor is not None:
            self.anchor(anchor)

    cdef _setup(self):
        cdef _TSObject ts = self.ts

        self.t = ts.value + self.dayoffset * us_in_day
        self.y = ts.dts.year

        self.ly = (ts.dts.month > 2 or
                   ts.dts.month == 2 and ts.dts.day == 29)

        if self.biz != 0:
            self.dow = (ts_dayofweek(ts) + self.dayoffset) % 7

    cpdef next(self):
        cdef int64_t days

        days = 365 + is_leapyear(self.y + self.ly)

        self.t += days * us_in_day
        self.y += 1

        if self.biz != 0:
            self.dow = (self.dow + days) % 7

    cpdef prev(self):
        cdef int64_t days

        days = 365 + is_leapyear(self.y - (1-self.ly))

        self.t -= days * us_in_day
        self.y -= 1

        if self.biz != 0:
            self.dow = (self.dow - days) % 7

cdef class MonthOffset(_Offset):
    """
    Generate monthly timestamps from provided start time, and apply dayoffset
    to each timestamp.  Stride to construct strided timestamps (eg quarterly).
    If biz > 0, we choose the next business day at each time; previous if < 0.

    Parameters
    ----------
    dayoffset : int
    stride : int, > 0
    biz : int
    """
    cdef:
        Py_ssize_t stride, ly, m
        int64_t y

    def __init__(self, int64_t dayoffset=0, Py_ssize_t stride=1,
                 int64_t biz=0, object anchor=None):
        self.dayoffset = dayoffset
        self.stride = stride
        self.biz = biz

        if stride <= 0:
            raise ValueError("Stride must be positive")

        if anchor is not None:
            self.anchor(anchor)

    cdef _setup(self):
        cdef _TSObject ts = self.ts

        self.t = ts.value + (self.dayoffset * us_in_day)

        # for day counting
        self.m  = ts.dts.month - 1
        self.y  = ts.dts.year
        self.ly = is_leapyear(self.y)

        if self.biz != 0:
            self.dow = (ts_dayofweek(ts) + self.dayoffset) % 7

    cpdef next(self):
        cdef:
            int64_t tmp, days
            Py_ssize_t j

        days = 0
        for j in range(0, self.stride):
            if self.m >= 12:
                self.m -= 12
                self.y += 1
                self.ly = is_leapyear(self.y)
            days += days_per_month_table[self.ly][self.m]
            self.m += 1

        self.t += days * us_in_day

        if self.biz != 0:
            self.dow = (self.dow + days) % 7

    cpdef prev(self):
        cdef:
            int64_t tmp, days
            Py_ssize_t j

        days = 0
        for j in range(0, self.stride):
            self.m -= 1
            if self.m < 0:
                self.m += 12
                self.y -= 1
                self.ly = is_leapyear(self.y)
            days += days_per_month_table[self.ly][self.m]

        self.t -= days * us_in_day

        if self.biz != 0:
            self.dow = (self.dow - days) % 7

cdef class DayOfMonthOffset(_Offset):
    """
    Generate relative monthly timestamps from month & year of provided start
    time. For example, fridays of the third week of each month (week=3, day=4);
    or, thursdays of the last week of each month (week=-1, day=3).

    Parameters
    ----------
    week : int
    day : int, 0 to 6
    """
    cdef:
        Py_ssize_t ly, m
        int64_t y, day, week

    def __init__(self, int64_t week=0, int64_t day=0, object anchor=None):
        self.week = week
        self.day = day

        if self.day < 0 or self.day > 6:
            raise ValueError("Day offset must be 0 to 6")

        if anchor is not None:
            self.anchor(anchor)

    cdef _setup(self):
        cdef _TSObject ts = self.ts

        # rewind to beginning of month
        self.t = ts.value - (ts.dts.day - 1) * us_in_day
        self.dow = dayofweek(ts.dts.year, ts.dts.month, 1)

        # for day counting
        self.m = ts.dts.month - 1
        self.y = ts.dts.year
        self.ly = is_leapyear(self.y)

    cpdef next(self):
        cdef:
            int64_t tmp, days

        days = days_per_month_table[self.ly][self.m]
        self.t += days * us_in_day
        self.dow = (self.dow + days) % 7

        self.m += 1
        if self.m >= 12:
            self.m -= 12
            self.y += 1
            self.ly = is_leapyear(self.y)

    cpdef prev(self):
        cdef:
            int64_t tmp, days

        days = days_per_month_table[self.ly][(self.m - 1) % 12]
        self.t -= days * us_in_day
        self.dow = (self.dow - days) % 7

        self.m -= 1
        if self.m < 0:
            self.m += 12
            self.y -= 1
            self.ly = is_leapyear(self.y)

    cdef int64_t _ts(self):
        """
        Overwrite default adjustment
        """
        cdef int64_t adj = (self.week * 7) + (self.day - self.dow) % 7
        return self.t + us_in_day * adj

cdef class DayOffset(_Offset):
    """
    Generate daily timestamps beginning with first valid time >= start time. If
    biz != 0, we skip weekends. Stride, to construct weekly timestamps.

    Parameters
    ----------
    stride : int, > 0
    biz : boolean
    """
    cdef:
        Py_ssize_t stride

    def __init__(self, int64_t stride=1, int64_t biz=0, object anchor=None):
        self.stride = stride
        self.biz = biz

        if self.stride <= 0:
            raise ValueError("Stride must be positive")

        if anchor is not None:
            self.anchor(anchor)

    cdef _setup(self):
        cdef _TSObject ts = self.ts
        self.t = ts.value
        if self.biz != 0:
            self.dow = ts_dayofweek(ts)

    cpdef next(self):
        self.t += (self.stride * us_in_day)
        if self.biz != 0:
            self.dow = (self.dow + self.stride) % 7
            if self.dow >= 5:
                self.t += (7 - self.dow) * us_in_day
                self.dow = 0

    cpdef prev(self):
        self.t -= (self.stride * us_in_day)
        if self.biz != 0:
            self.dow = (self.dow - self.stride) % 7
            if self.dow >= 5:
                self.t += (4 - self.dow) * us_in_day
                self.dow = 4
