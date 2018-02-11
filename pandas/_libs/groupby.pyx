# -*- coding: utf-8 -*-
# cython: profile=False

cimport numpy as cnp
import numpy as np

cimport cython

cnp.import_array()

from numpy cimport (ndarray,
                    double_t,
                    int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t,
                    uint32_t, uint64_t, float32_t, float64_t)

from libc.stdlib cimport malloc, free

from util cimport numeric, get_nat
from algos cimport (swap, TiebreakEnumType, TIEBREAK_AVERAGE, TIEBREAK_MIN,
                    TIEBREAK_MAX, TIEBREAK_FIRST, TIEBREAK_DENSE)
from algos import take_2d_axis1_float64_float64, groupsort_indexer, tiebreakers

cdef int64_t iNaT = get_nat()

cdef double NaN = <double> np.NaN
cdef double nan = NaN


cdef inline float64_t median_linear(float64_t* a, int n) nogil:
    cdef int i, j, na_count = 0
    cdef float64_t result
    cdef float64_t* tmp

    if n == 0:
        return NaN

    # count NAs
    for i in range(n):
        if a[i] != a[i]:
            na_count += 1

    if na_count:
        if na_count == n:
            return NaN

        tmp = <float64_t*> malloc((n - na_count) * sizeof(float64_t))

        j = 0
        for i in range(n):
            if a[i] == a[i]:
                tmp[j] = a[i]
                j += 1

        a = tmp
        n -= na_count

    if n % 2:
        result = kth_smallest_c( a, n / 2, n)
    else:
        result = (kth_smallest_c(a, n / 2, n) +
                  kth_smallest_c(a, n / 2 - 1, n)) / 2

    if na_count:
        free(a)

    return result


cdef inline float64_t kth_smallest_c(float64_t* a,
                                     Py_ssize_t k,
                                     Py_ssize_t n) nogil:
    cdef:
        Py_ssize_t i, j, l, m
        double_t x, t

    l = 0
    m = n -1
    while (l<m):
        x = a[k]
        i = l
        j = m

        while 1:
            while a[i] < x: i += 1
            while x < a[j]: j -= 1
            if i <= j:
                swap(&a[i], &a[j])
                i += 1; j -= 1

            if i > j: break

        if j < k: l = i
        if k < i: m = j
    return a[k]


# generated from template
include "groupby_helper.pxi"
