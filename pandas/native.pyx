# cython: profile=False
# distutils: language = c++
# cython: embedsignature = True

from numpy cimport ndarray
cimport numpy as cnp

import numpy as np

from pandas.native cimport shared_ptr
cimport pandas.native as lp

cnp.import_array()


def wrap_numpy_array(ndarray arr):
    pass


cdef class PandasType:
    cdef:
        lp.TypePtr type

    cdef init(self, const lp.TypePtr& type):
        self.type = type

    def __repr__(self):
        return 'PandasType({0})'.format(self.type.get().ToString())


cdef class Category(PandasType):

    property categories:

        def __get__(self):
            pass


cdef class Array:
    cdef:
        lp.ArrayPtr arr

    cdef init(self, const lp.ArrayPtr& arr):
        self.arr = arr


cdef class NumericArray(Array):
    pass


cdef class FloatingArray(NumericArray):
    pass


cdef class IntegerArray(NumericArray):
    pass


cdef class Float32Array(FloatingArray):
    pass


cdef class BooleanArray(Array):
    cdef:
        lp.cBooleanArray* inst

    cdef init(self, const ArrayPtr& arr):
        Array.init(self, arr)


cdef class CategoryArray(Array):
    pass


cdef Array wrap_array(const lp.ArrayPtr& type):
    cdef:
        Array result

    if type.get().type_enum() == lp.CATEGORY:
        result = CategoryArray()
    else:
        result = Array()

    return result


cdef PandasType wrap_type(const lp.TypePtr& sp_type):
    cdef:
        lp.DataType* type = sp_type.get()
        PandasType result

    if type.type == lp.CATEGORY:
        result = Category()
    else:
        result = PandasType()

    result.init(sp_type)

    return result


cpdef PandasType convert_numpy_dtype(cnp.dtype dtype):
    cdef lp.TypePtr sp_type

    if dtype.type_num == cnp.NPY_INT8:
        sp_type.reset(new lp.Int8Type())
    else:
        raise NotImplementedError(dtype)

    return wrap_type(sp_type)
