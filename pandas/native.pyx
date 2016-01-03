# cython: profile=False
# distutils: language = c++
# cython: embedsignature = True

from numpy cimport ndarray
cimport numpy as cnp

import numpy as np

from pandas.native cimport shared_ptr, string
cimport pandas.native as lp

from cython.operator cimport dereference as deref

cnp.import_array()


UINT8 = lp.TypeEnum_UINT8
UINT16 = lp.TypeEnum_UINT16
UINT32 = lp.TypeEnum_UINT32
UINT64 = lp.TypeEnum_UINT64
INT8 = lp.TypeEnum_INT8
INT16 = lp.TypeEnum_INT16
INT32 = lp.TypeEnum_INT32
INT64 = lp.TypeEnum_INT64
BOOL = lp.TypeEnum_BOOL
FLOAT = lp.TypeEnum_FLOAT
DOUBLE = lp.TypeEnum_DOUBLE
PYOBJECT = lp.TypeEnum_PYOBJECT
CATEGORY = lp.TypeEnum_CATEGORY
TIMESTAMP = lp.TypeEnum_TIMESTAMP
TIMESTAMP_TZ = lp.TypeEnum_TIMESTAMP_TZ


class CPandasException(Exception):
    pass


class CPandasBadStatus(CPandasException):
    """
    libpandas operation returned an error Status
    """
    pass


class CPandasNotImplemented(CPandasBadStatus):
    pass


cdef check_status(const lp.Status& status):
    if status.ok():
        return

    message = status.ToString()

    if status.IsNotImplemented():
        raise CPandasNotImplemented(message)
    else:
        raise CPandasBadStatus(message)


cdef class Scalar:
    cdef readonly:
        lp.TypeEnum type


cdef class NAType(Scalar):

    def __cinit__(self):
        self.type = lp.TypeEnum_NA

NA = NAType()


cdef dict _primitive_type_aliases = {
    'u1': UINT8,
    'u2': UINT16,
    'u4': UINT32,
    'u8': UINT64,
    'uint8': UINT8,
    'uint16': UINT16,
    'uint32': UINT32,
    'uint64': UINT64,

    'i1': INT8,
    'i2': INT16,
    'i4': INT32,
    'i8': INT64,
    'int8': INT8,
    'int16': INT16,
    'int32': INT32,
    'int64': INT64,

    'f4': FLOAT,
    'f8': DOUBLE,
    'float32': FLOAT,
    'float64': DOUBLE,

    'b1': BOOL,
    'bool': BOOL,

    'O8': PYOBJECT,
    'object': PYOBJECT,
}


def wrap_numpy_array(ndarray arr):
    pass


cdef class PandasType:
    cdef:
        lp.TypePtr type

    cdef init(self, const lp.TypePtr& type):
        self.type = type

    def __repr__(self):
        return 'PandasType({0})'.format(self.type.get().ToString())

    def equals(PandasType self, PandasType other):
        return self.type.get().Equals(deref(other.type.get()))


def primitive_type(TypeEnum tp_enum):
    cdef:
        lp.TypePtr sp_type
        lp.DataType* type

    check_status(lp.primitive_type_from_enum(tp_enum, &type))
    sp_type.reset(type)
    return wrap_type(sp_type)


cdef class Category(PandasType):

    property categories:

        def __get__(self):
            pass


def category_type(categories):
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

    if type.get().type_enum() == lp.TypeEnum_CATEGORY:
        result = CategoryArray()
    else:
        result = Array()

    return result


cdef PandasType wrap_type(const lp.TypePtr& sp_type):
    cdef:
        lp.DataType* type = sp_type.get()
        PandasType result

    if type.type == lp.TypeEnum_CATEGORY:
        result = Category()
    else:
        result = PandasType()

    result.init(sp_type)

    return result


cpdef PandasType convert_numpy_dtype(cnp.dtype dtype):
    cdef lp.TypeEnum pandas_typenum

    check_status(lp.numpy_type_num_to_pandas(dtype.type_num,
                                             &pandas_typenum))

    return primitive_type(pandas_typenum)


def to_pandas_array(values):
    pass


def numpy_to_pandas_array(ndarray arr):
    pass
