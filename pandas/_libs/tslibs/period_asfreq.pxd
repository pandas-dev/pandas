# -*- coding: utf-8 -*-

from numpy cimport int64_t

from period_info cimport date_info


ctypedef struct asfreq_info:
    int from_week_end
    int to_week_end

    int from_a_year_end
    int to_a_year_end

    int from_q_year_end
    int to_q_year_end

    int64_t intraday_conversion_factor


ctypedef int64_t (*freq_conv_func)(int64_t, char, asfreq_info*) nogil

cdef freq_conv_func get_asfreq_func(int fromFreq, int toFreq) nogil
cdef int64_t DtoQ_yq(int64_t ordinal, asfreq_info *af_info, int *year,
                     int *quarter) nogil
cdef void get_asfreq_info(int fromFreq, int toFreq, asfreq_info *af_info) nogil
cdef int64_t get_python_ordinal(int64_t period_ordinal, int freq) nogil
cdef int64_t asfreq(int64_t ordinal, int freq1, int freq2, char relation)

cdef int get_date_info(int64_t ordinal, int freq,
                       date_info *dinfo) nogil except -1

cdef int pqyear(int64_t ordinal, int freq)
cdef int pquarter(int64_t ordinal, int freq)
cdef int pday_of_year(int64_t ordinal, int freq)
cdef int pweek(int64_t ordinal, int freq)
cdef int pweekday(int64_t ordinal, int freq)
cdef int pyear(int64_t ordinal, int freq)
cdef int pmonth(int64_t ordinal, int freq)
cdef int pday(int64_t ordinal, int freq)
cdef int phour(int64_t ordinal, int freq)
cdef int pminute(int64_t ordinal, int freq)
cdef int psecond(int64_t ordinal, int freq)
cdef int pdays_in_month(int64_t ordinal, int freq)
