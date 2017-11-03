# -*- coding: utf-8 -*-
# cython: profile=False

from cpython cimport PyFloat_Check, PyComplex_Check

cimport numpy as cnp
cnp.import_array()

from util cimport is_datetime64_object, is_timedelta64_object

from datetime cimport get_timedelta64_value, get_datetime64_value

from pandas._libs.tslibs.nattype cimport NPY_NAT
from pandas._libs.tslibs.nattype import NaT


cdef inline bint checknull_with_nat(object val):
    """ utility to check if a value is a nat or not """
    return val is None or (
        PyFloat_Check(val) and val != val) or val is NaT


cdef inline bint check_all_nulls(object val):
    """ utility to check if a value is any type of null """
    cdef:
        bint res
    if PyFloat_Check(val) or PyComplex_Check(val):
        res = val != val
    elif val is NaT:
        res = 1
    elif val is None:
        res = 1
    elif is_datetime64_object(val):
        res = get_datetime64_value(val) == NPY_NAT
    elif is_timedelta64_object(val):
        res = get_timedelta64_value(val) == NPY_NAT
    else:
        res = 0
    return res
