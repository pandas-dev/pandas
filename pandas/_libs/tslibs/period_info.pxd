# -*- coding: utf-8 -*-

import numpy as np
from numpy cimport int64_t

from util cimport INT32_MIN


ctypedef struct date_info:
    int64_t absdate
    double abstime
    double second
    int minute
    int hour
    int day
    int month
    int quarter
    int year
    int day_of_week
    int day_of_year
    int calendar


cdef enum OFFSETS:
    ORD_OFFSET = 719163LL   # days until 1970-01-01
    BDAY_OFFSET = 513689LL  # days until 1970-01-01
    WEEK_OFFSET = 102737LL

cdef enum CALENDARS:
    GREGORIAN_CALENDAR = 1
    JULIAN_CALENDAR = 2

cdef enum FREQS:
    FR_ANN = 1000       # Annual
    FR_ANNDEC = FR_ANN  # Annual - December year end
    FR_ANNJAN = 1001    # Annual - January year end
    FR_ANNFEB = 1002    # Annual - February year end
    FR_ANNMAR = 1003    # Annual - March year end
    FR_ANNAPR = 1004    # Annual - April year end
    FR_ANNMAY = 1005    # Annual - May year end
    FR_ANNJUN = 1006    # Annual - June year end
    FR_ANNJUL = 1007    # Annual - July year end
    FR_ANNAUG = 1008    # Annual - August year end
    FR_ANNSEP = 1009    # Annual - September year end
    FR_ANNOCT = 1010    # Annual - October year end
    FR_ANNNOV = 1011    # Annual - November year end

    # The standard quarterly frequencies with various fiscal year ends
    #   eg, Q42005 for Q@OCT runs Aug 1, 2005 to Oct 31, 2005
    FR_QTR = 2000       # Quarterly - December year end (default quarterly)
    FR_QTRDEC = FR_QTR  # Quarterly - December year end
    FR_QTRJAN = 2001    # Quarterly - January year end
    FR_QTRFEB = 2002    # Quarterly - February year end
    FR_QTRMAR = 2003    # Quarterly - March year end
    FR_QTRAPR = 2004    # Quarterly - April year end
    FR_QTRMAY = 2005    # Quarterly - May year end
    FR_QTRJUN = 2006    # Quarterly - June year end
    FR_QTRJUL = 2007    # Quarterly - July year end
    FR_QTRAUG = 2008    # Quarterly - August year end
    FR_QTRSEP = 2009    # Quarterly - September year end
    FR_QTROCT = 2010    # Quarterly - October year end
    FR_QTRNOV = 2011    # Quarterly - November year end

    FR_MTH = 3000     # Monthly

    FR_WK = 4000      # Weekly
    FR_WKSUN = FR_WK  # Weekly - Sunday end of week
    FR_WKMON = 4001   # Weekly - Monday end of week
    FR_WKTUE = 4002   # Weekly - Tuesday end of week
    FR_WKWED = 4003   # Weekly - Wednesday end of week
    FR_WKTHU = 4004   # Weekly - Thursday end of week
    FR_WKFRI = 4005   # Weekly - Friday end of week
    FR_WKSAT = 4006   # Weekly - Saturday end of week

    FR_BUS = 5000     # Business days
    FR_DAY = 6000     # Daily
    FR_HR = 7000      # Hourly
    FR_MIN = 8000     # Minutely
    FR_SEC = 9000     # Secondly
    FR_MS = 10000     # Millisecondly
    FR_US = 11000     # Microsecondly
    FR_NS = 12000     # Nanosecondly

    FR_UND = -10000   # Undefined


cdef int dInfoCalc_SetFromAbsDateTime(date_info *dinfo,
                                      int64_t absdate, double abstime,
                                      int calendar) nogil except -1
cdef int dInfoCalc_SetFromAbsDate(date_info *dinfo,
                                  int64_t absdate, int calendar) nogil
cdef int dInfoCalc_SetFromAbsTime(date_info *dinfo, double abstime) nogil

cdef int64_t absdate_from_ymd(int y, int m, int d) nogil
cdef int monthToQuarter(int month) nogil

cdef int dInfoCalc_YearOffset(int64_t year, int calendar) nogil except? -1
cdef int dInfoCalc_DayOfWeek(int64_t absdate) nogil
cdef bint dInfoCalc_Leapyear(int64_t year, int calendar) nogil
cdef int _ISOWeek(date_info *dinfo)

cdef int64_t get_period_ordinal(int year, int month, int day,
                                int hour, int minute, int second,
                                int microseconds, int picoseconds,
                                int freq) nogil except INT32_MIN
