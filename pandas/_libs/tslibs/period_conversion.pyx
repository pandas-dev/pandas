# -*- coding: utf-8 -*-
# flake8: noqa:E241,E501
cimport cython

import numpy as np
cimport numpy as cnp
from numpy cimport ndarray, int64_t
cnp.import_array()


cdef enum FREQS:
    FR_ANN = 1000      # Annual
    FR_QTR = 2000      # Quarterly - December year end (default quarterly)
    FR_MTH = 3000      # Monthly
    FR_WK = 4000       # Weekly
    FR_BUS = 5000      # Business days
    FR_DAY = 6000      # Daily
    FR_HR = 7000       # Hourly
    FR_MIN = 8000      # Minutely
    FR_SEC = 9000      # Secondly
    FR_MS = 10000      # Millisecondly
    FR_US = 11000      # Microsecondly
    FR_NS = 12000      # Nanosecondly
    FR_UND = -10000    # Undefined


daytime_conversion_factors = [[FR_DAY, 1],   [FR_HR, 24],
                              [FR_MIN, 60],  [FR_SEC, 60],
                              [FR_MS, 1000], [FR_US, 1000],
                              [FR_NS, 1000], [0, 0]]


cdef int64_t[:, :] daytime_conversion_factor_matrix = np.array(
    [[<int64_t>val for val in row] for row in
    # Building on Windows seems to require super-explicit casting
    [
    [0, 0, 0, 0, 0, 0, 0,  0,    0,     0,        0,           0,              0],
    [0, 0, 0, 0, 0, 0, 0,  0,    0,     0,        0,           0,              0],
    [0, 0, 0, 0, 0, 0, 0,  0,    0,     0,        0,           0,              0],
    [0, 0, 0, 0, 0, 0, 0,  0,    0,     0,        0,           0,              0],
    [0, 0, 0, 0, 0, 0, 0,  0,    0,     0,        0,           0,              0],
    [0, 0, 0, 0, 0, 0, 0,  0,    0,     0,        0,           0,              0],
    [0, 0, 0, 0, 0, 0, 1, 24, 1440, 86400, 86400000, 86400000000, 86400000000000],
    [0, 0, 0, 0, 0, 0, 0,  1,   60,  3600,  3600000,  3600000000,  3600000000000],
    [0, 0, 0, 0, 0, 0, 0,  0,    1,    60,    60000,    60000000,    60000000000],
    [0, 0, 0, 0, 0, 0, 0,  0,    0,     1,     1000,     1000000,     1000000000],
    [0, 0, 0, 0, 0, 0, 0,  0,    0,     0,        1,        1000,        1000000],
    [0, 0, 0, 0, 0, 0, 0,  0,    0,     0,        0,           1,           1000],
    [0, 0, 0, 0, 0, 0, 0,  0,    0,     0,        0,           0,              1]]],
    dtype=np.int64)  # noqa




@cython.boundscheck(False)
cdef int64_t get_daytime_conversion_factor(int from_index, int to_index) nogil:
    cdef:
        Py_ssize_t row, col

    row = min_value(from_index, to_index)
    col = max_value(from_index, to_index)
    return daytime_conversion_factor_matrix[row][col]


cdef inline int max_value(int a, int b) nogil:
    if a > b:
        return a
    return b


cdef inline int min_value(int a, int b) nogil:
    if a < b:
        return a
    return b


@cython.cdivision
cdef double get_abs_time(int freq, int64_t date_ordinal,
                         int64_t ordinal) nogil:
    cdef:
        int freq_index, day_index, base_index
        int64_t per_day, start_ord
        double unit, result

    if freq <= FR_DAY:
        return 0

    freq_index = get_freq_group_index(freq)
    day_index = get_freq_group_index(FR_DAY)
    base_index = get_freq_group_index(FR_SEC)

    per_day = get_daytime_conversion_factor(day_index, freq_index)
    unit = get_daytime_conversion_factor(freq_index, base_index)

    if base_index < freq_index:
        unit = 1 / unit

    start_ord = date_ordinal * per_day
    result = <double>(unit * (ordinal - start_ord))
    return result


# ----------------------------------------------------------------------
# Conventions

@cython.cdivision
cdef inline int get_freq_group_index(int freq) nogil:
    return freq / 1000

@cython.cdivision
cdef inline int get_freq_group(int freq) nogil:
    return (freq / 1000) * 1000
