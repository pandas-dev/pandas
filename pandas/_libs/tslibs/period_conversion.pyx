# -*- coding: utf-8 -*-
cimport cython

import numpy as np
from numpy cimport int64_t


cdef enum FREQS:
    FR_ANN = 1000      # Annual
    FR_QTR = 2000      # Quarterly - December year end (default quarterly)
    FR_MTH = 3000 # Monthly
    FR_WK = 4000     # Weekly
    FR_BUS = 5000 # Business days
    FR_DAY = 6000 # Daily
    FR_HR = 7000  # Hourly
    FR_MIN = 8000 # Minutely
    FR_SEC = 9000 # Secondly
    FR_MS = 10000 # Millisecondly
    FR_US = 11000 # Microsecondly
    FR_NS = 12000 # Nanosecondly

    FR_UND = -10000 # Undefined


cdef int64_t[:, :] daytime_conversion_factor_matrix = np.array([
    [0, 0, 0, 0, 0, 0, 0,  0,    0,     0,        0,           0,              0],
    [0, 0, 0, 0, 0, 0, 0,  0,    0,     0,        0,           0,              0],
    [0, 0, 0, 0, 0, 0, 0,  0,    0,     0,        0,           0,              0],
    [0, 0, 0, 0, 0, 0, 0,  0,    0,     0,        0,           0,              0],
    [0, 0, 0, 0, 0, 0, 0,  0,    0,     0,        0,           0,              0],
    [0, 0, 0, 0, 0, 0, 0,  0,    0,     0,        0,           0,              0],
    [0, 0, 0, 0, 0, 0, 1, 24, 1440, 86400, 86400000, 86400000000, 86400000000000],
    [0, 0, 0, 0, 0, 0, 0,  1,   60,  3600,  3600000,  3600000000,  3600000000000],
    [0, 0, 0, 0, 0, 0, 0,  0,    1,   60,     60000,    60000000,    60000000000],
    [0, 0, 0, 0, 0, 0, 0,  0,    0,    1,      1000,     1000000,     1000000000],
    [0, 0, 0, 0, 0, 0, 0,  0,    0,    0,         1,        1000,        1000000],
    [0, 0, 0, 0, 0, 0, 0,  0,    0,    0,         0,           1,           1000],
    [0, 0, 0, 0, 0, 0, 0,  0,    0,    0,         0,           0,              1]],
    dtype=np.int64)


@cython.boundscheck(False)
cdef int64_t get_daytime_conversion_factor(int from_index, int to_index) nogil:
    return daytime_conversion_factor_matrix[min_value(from_index, to_index)][max_value(from_index, to_index)]


cdef inline int max_value(int a, int b) nogil:
    if a > b:
        return a
    return b


cdef inline int min_value(int a, int b) nogil:
    if a < b:
        return a
    return b


@cython.cdivision
cdef double get_abs_time(int freq, int64_t date_ordinal, int64_t ordinal) nogil:
    # // printf("get_abs_time %d %lld %lld\n", freq, date_ordinal, ordinal);
    cdef:
        int freq_index, day_index, base_index
        int64_t per_day, start_ord
        double unit, result

    if freq <= FR_DAY:
        return 0

    freq_index = get_freq_group_index(freq)
    day_index = get_freq_group_index(FR_DAY)
    base_index = get_freq_group_index(FR_SEC)

    #// printf("  indices: day %d, freq %d, base %d\n", day_index, freq_index,
    #// base_index);

    per_day = get_daytime_conversion_factor(day_index, freq_index)
    unit = get_daytime_conversion_factor(freq_index, base_index)

    #// printf("  per_day: %lld, unit: %f\n", per_day, unit);

    if base_index < freq_index:
        unit = 1 / unit
        #// printf("  corrected unit: %f\n", unit);

    start_ord = date_ordinal * per_day
    #// printf("start_ord: %lld\n", start_ord);
    result = <double>(unit * (ordinal - start_ord))
    #// printf("  result: %f\n", result);
    return result


# TODO: De-duplicate this with versions elsewhere
@cython.cdivision
cdef inline int get_freq_group_index(int freq) nogil:
    return freq / 1000


'''

daytime_conversion_factors = [[FR_DAY, 1],   [FR_HR, 24],
                              [FR_MIN, 60],  [FR_SEC, 60],
                              [FR_MS, 1000], [FR_US, 1000],
                              [FR_NS, 1000], [0, 0]]

def calculate_conversion_factor(int start_value, int end_value):
    conversion_factor = 0

    index = 0
    while True:
        freq_group = daytime_conversion_factors[index][0]

        if freq_group == 0:
            conversion_factor = 0
            break

        if freq_group == start_value:
            conversion_factor = 1
        else:
            conversion_factor *= daytime_conversion_factors[index][1]

        if freq_group == end_value:
            break

        index += 1
    return conversion_factor


def calc_conversion_factors_matrix_size():
    matrix_size = 0

    index = 0
    while True:
        period_value = get_freq_group_index(daytime_conversion_factors[index][0])
        if period_value == 0:
            break

        matrix_size = max_value(matrix_size, period_value)
        index += 1

    return matrix_size + 1


def populate_conversion_factors_matrix():
    N = len(daytime_conversion_factors)
    mat = np.zeros((13, 13), dtype=np.int64)

    row_index_index = 0
    while True:
        row_value = daytime_conversion_factors[row_index_index][0]
        if row_value == 0:
            break

        row_index = get_freq_group_index(row_value)

        column_index_index = row_index_index
        while True:
            column_value = daytime_conversion_factors[column_index_index][0]
            if column_value == 0:
                break

            column_index = get_freq_group_index(column_value)
            mat[row_index][column_index] = calculate_conversion_factor(row_value, column_value)
            column_index_index += 1

        row_index_index += 1
    return mat
'''
