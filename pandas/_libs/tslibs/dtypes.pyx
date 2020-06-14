# period frequency constants corresponding to scikits timeseries
# originals


cdef class PeriodDtypeBase:
    """
    Similar to an actual dtype, this contains all of the information
    describing a PeriodDtype in an integer code.
    """
    # cdef readonly:
    #    PeriodDtypeCode dtype_code

    def __cinit__(self, PeriodDtypeCode code):
        self.dtype_code = code

    def __eq__(self, other):
        if not isinstance(other, PeriodDtypeBase):
            return False
        if not isinstance(self, PeriodDtypeBase):
            # cython semantics, this is a reversed op
            return False
        return self.dtype_code == other.dtype_code

    @property
    def freq_group(self) -> int:
        # See also: libperiod.get_freq_group
        return (self.dtype_code // 1000) * 1000

    @property
    def date_offset(self):
        """
        Corresponding DateOffset object.

        This mapping is mainly for backward-compatibility.
        """
        from .offsets import to_offset

        freqstr = _reverse_period_code_map.get(self.dtype_code)
        # equiv: freqstr = libfrequencies.get_freq_str(self.dtype_code)

        return to_offset(freqstr)

    @classmethod
    def from_date_offset(cls, offset):
        code = offset._period_dtype_code
        return cls(code)


_period_code_map = {
    # Annual freqs with various fiscal year ends.
    # eg, 2005 for A-FEB runs Mar 1, 2004 to Feb 28, 2005
    "A-DEC": 1000,  # Annual - December year end
    "A-JAN": 1001,  # Annual - January year end
    "A-FEB": 1002,  # Annual - February year end
    "A-MAR": 1003,  # Annual - March year end
    "A-APR": 1004,  # Annual - April year end
    "A-MAY": 1005,  # Annual - May year end
    "A-JUN": 1006,  # Annual - June year end
    "A-JUL": 1007,  # Annual - July year end
    "A-AUG": 1008,  # Annual - August year end
    "A-SEP": 1009,  # Annual - September year end
    "A-OCT": 1010,  # Annual - October year end
    "A-NOV": 1011,  # Annual - November year end

    # Quarterly frequencies with various fiscal year ends.
    # eg, Q42005 for Q-OCT runs Aug 1, 2005 to Oct 31, 2005
    "Q-DEC": 2000,    # Quarterly - December year end
    "Q-JAN": 2001,    # Quarterly - January year end
    "Q-FEB": 2002,    # Quarterly - February year end
    "Q-MAR": 2003,    # Quarterly - March year end
    "Q-APR": 2004,    # Quarterly - April year end
    "Q-MAY": 2005,    # Quarterly - May year end
    "Q-JUN": 2006,    # Quarterly - June year end
    "Q-JUL": 2007,    # Quarterly - July year end
    "Q-AUG": 2008,    # Quarterly - August year end
    "Q-SEP": 2009,    # Quarterly - September year end
    "Q-OCT": 2010,    # Quarterly - October year end
    "Q-NOV": 2011,    # Quarterly - November year end

    "M": 3000,        # Monthly

    "W-SUN": 4000,    # Weekly - Sunday end of week
    "W-MON": 4001,    # Weekly - Monday end of week
    "W-TUE": 4002,    # Weekly - Tuesday end of week
    "W-WED": 4003,    # Weekly - Wednesday end of week
    "W-THU": 4004,    # Weekly - Thursday end of week
    "W-FRI": 4005,    # Weekly - Friday end of week
    "W-SAT": 4006,    # Weekly - Saturday end of week

    "B": 5000,        # Business days
    "D": 6000,        # Daily
    "H": 7000,        # Hourly
    "T": 8000,        # Minutely
    "S": 9000,        # Secondly
    "L": 10000,       # Millisecondly
    "U": 11000,       # Microsecondly
    "N": 12000,       # Nanosecondly
}

_reverse_period_code_map = {
    _period_code_map[key]: key for key in _period_code_map}

# Yearly aliases; careful not to put these in _reverse_period_code_map
_period_code_map.update({"Y" + key[1:]: _period_code_map[key]
                         for key in _period_code_map
                         if key.startswith("A-")})

_period_code_map.update({
    "Q": 2000,   # Quarterly - December year end (default quarterly)
    "A": 1000,   # Annual
    "W": 4000,   # Weekly
    "C": 5000,   # Custom Business Day
})


# Map attribute-name resolutions to resolution abbreviations
_attrname_to_abbrevs = {
    "year": "A",
    "quarter": "Q",
    "month": "M",
    "day": "D",
    "hour": "H",
    "minute": "T",
    "second": "S",
    "millisecond": "L",
    "microsecond": "U",
    "nanosecond": "N",
}
cdef dict attrname_to_abbrevs = _attrname_to_abbrevs


class FreqGroup:
    # Mirrors c_FreqGroup in the .pxd file
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

    @staticmethod
    def get_freq_group(code: int) -> int:
        # See also: PeriodDtypeBase.freq_group
        return (code // 1000) * 1000
