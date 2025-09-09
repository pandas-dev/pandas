from numpy cimport int64_t

from pandas._libs.tslibs.np_datetime cimport NPY_DATETIMEUNIT


cdef str npy_unit_to_abbrev(NPY_DATETIMEUNIT unit)
cpdef NPY_DATETIMEUNIT abbrev_to_npy_unit(str abbrev)
cdef NPY_DATETIMEUNIT freq_group_code_to_npy_unit(int freq) noexcept nogil
cpdef int64_t periods_per_day(NPY_DATETIMEUNIT reso=*) except? -1
cpdef int64_t periods_per_second(NPY_DATETIMEUNIT reso) except? -1
cdef NPY_DATETIMEUNIT get_supported_reso(NPY_DATETIMEUNIT reso)
cdef bint is_supported_unit(NPY_DATETIMEUNIT reso)

cdef dict c_OFFSET_TO_PERIOD_FREQSTR
cdef dict c_PERIOD_TO_OFFSET_FREQSTR
cdef dict c_OFFSET_RENAMED_FREQSTR
cdef dict c_DEPR_UNITS
cdef dict c_PERIOD_AND_OFFSET_DEPR_FREQSTR
cdef dict attrname_to_abbrevs
cdef dict npy_unit_to_attrname
cdef dict attrname_to_npy_unit
cdef str INVALID_FREQ_ERR_MSG

cdef enum c_FreqGroup:
    # Mirrors FreqGroup in the .pyx file
    FR_ANN = 1000
    FR_QTR = 2000
    FR_MTH = 3000
    FR_WK = 4000
    FR_BUS = 5000
    FR_DAY = 6000
    FR_HR = 7000
    FR_MIN = 8000
    FR_SEC = 9000
    FR_MS = 10000
    FR_US = 11000
    FR_NS = 12000
    FR_UND = -10000  # undefined


cdef enum c_Resolution:
    # Mirrors Resolution in the .pyx file
    RESO_NS = 0
    RESO_US = 1
    RESO_MS = 2
    RESO_SEC = 3
    RESO_MIN = 4
    RESO_HR = 5
    RESO_DAY = 6
    RESO_MTH = 7
    RESO_QTR = 8
    RESO_YR = 9


cdef enum PeriodDtypeCode:
    # Annual freqs with various fiscal year ends.
    # eg, 2005 for A_FEB runs Mar 1, 2004 to Feb 28, 2005
    A = 1000      # Default alias
    A_DEC = 1000  # Annual - December year end
    A_JAN = 1001  # Annual - January year end
    A_FEB = 1002  # Annual - February year end
    A_MAR = 1003  # Annual - March year end
    A_APR = 1004  # Annual - April year end
    A_MAY = 1005  # Annual - May year end
    A_JUN = 1006  # Annual - June year end
    A_JUL = 1007  # Annual - July year end
    A_AUG = 1008  # Annual - August year end
    A_SEP = 1009  # Annual - September year end
    A_OCT = 1010  # Annual - October year end
    A_NOV = 1011  # Annual - November year end

    # Quarterly frequencies with various fiscal year ends.
    # eg, Q42005 for Q_OCT runs Aug 1, 2005 to Oct 31, 2005
    Q_DEC = 2000    # Quarterly - December year end
    Q_JAN = 2001    # Quarterly - January year end
    Q_FEB = 2002    # Quarterly - February year end
    Q_MAR = 2003    # Quarterly - March year end
    Q_APR = 2004    # Quarterly - April year end
    Q_MAY = 2005    # Quarterly - May year end
    Q_JUN = 2006    # Quarterly - June year end
    Q_JUL = 2007    # Quarterly - July year end
    Q_AUG = 2008    # Quarterly - August year end
    Q_SEP = 2009    # Quarterly - September year end
    Q_OCT = 2010    # Quarterly - October year end
    Q_NOV = 2011    # Quarterly - November year end

    M = 3000        # Monthly

    W_SUN = 4000    # Weekly - Sunday end of week
    W_MON = 4001    # Weekly - Monday end of week
    W_TUE = 4002    # Weekly - Tuesday end of week
    W_WED = 4003    # Weekly - Wednesday end of week
    W_THU = 4004    # Weekly - Thursday end of week
    W_FRI = 4005    # Weekly - Friday end of week
    W_SAT = 4006    # Weekly - Saturday end of week

    B = 5000        # Business days
    D = 6000        # Daily
    H = 7000        # Hourly
    T = 8000        # Minutely
    S = 9000        # Secondly
    L = 10000       # Millisecondly
    U = 11000       # Microsecondly
    N = 12000       # Nanosecondly

    UNDEFINED = -10_000


cdef class PeriodDtypeBase:
    cdef readonly:
        PeriodDtypeCode _dtype_code
        int64_t _n

    cpdef int _get_to_timestamp_base(self)
    cpdef bint _is_tick_like(self)
