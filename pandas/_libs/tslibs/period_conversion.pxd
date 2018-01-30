# -*- coding: utf-8 -*-

from numpy cimport int64_t

cdef int64_t[:, :] daytime_conversion_factor_matrix
cdef int64_t get_daytime_conversion_factor(int from_index, int to_index) nogil
cdef int max_value(int a, int b) nogil
cdef int min_value(int a, int b) nogil
cdef double get_abs_time(int freq, int64_t date_ordinal, int64_t ordinal) nogil
