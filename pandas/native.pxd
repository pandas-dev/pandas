# Copyright 2014 Cloudera, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# distutils: language = c++

from libc.stdint cimport *
from libcpp cimport bool as c_bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from cpython cimport PyObject

# This must be included for cerr and other things to work
cdef extern from "<iostream>":
    pass


cdef extern from "<memory>" namespace "std" nogil:

    cdef cppclass shared_ptr[T]:
        T* get()
        void reset()
        void reset(T* p)


cdef extern from "pandas/status.h" namespace "pandas" nogil:

    # We can later add more of the common status factory methods as needed
    cdef Status Status_OK "Status::OK"()

    cdef cppclass Status:
        Status()

        string ToString()

        c_bool ok()
        c_bool IsKeyError()
        c_bool IsOutOfMemory()
        c_bool IsInvalid()
        c_bool IsNotImplemented()


cdef extern from "pandas/api.h" namespace "pandas":

    enum TypeEnum:
        TypeEnum_NA " pandas::TypeEnum::NA"
        TypeEnum_UINT8 " pandas::TypeEnum::UINT8"
        TypeEnum_UINT16 " pandas::TypeEnum::UINT16"
        TypeEnum_UINT32 " pandas::TypeEnum::UINT32"
        TypeEnum_UINT64 " pandas::TypeEnum::UINT64"
        TypeEnum_INT8 " pandas::TypeEnum::INT8"
        TypeEnum_INT16 " pandas::TypeEnum::INT16"
        TypeEnum_INT32 " pandas::TypeEnum::INT32"
        TypeEnum_INT64 " pandas::TypeEnum::INT64"
        TypeEnum_BOOL " pandas::TypeEnum::BOOL"
        TypeEnum_FLOAT " pandas::TypeEnum::FLOAT"
        TypeEnum_DOUBLE " pandas::TypeEnum::DOUBLE"
        TypeEnum_PYOBJECT " pandas::TypeEnum::PYOBJECT"
        TypeEnum_CATEGORY " pandas::TypeEnum::CATEGORY"
        TypeEnum_TIMESTAMP " pandas::TypeEnum::TIMESTAMP"
        TypeEnum_TIMESTAMP_TZ " pandas::TypeEnum::TIMESTAMP_TZ"

    cdef cppclass DataType:
        TypeEnum type

        DataType()

        c_bool Equals(const DataType& other)
        string ToString()

    ctypedef shared_ptr[DataType] TypePtr

    cdef cppclass Int8Type(DataType):
        pass

    cdef cppclass Int16Type(DataType):
        pass

    cdef cppclass Int32Type(DataType):
        pass

    cdef cppclass Int64Type(DataType):
        pass

    cdef cppclass UInt8Type(DataType):
        pass

    cdef cppclass UInt16Type(DataType):
        pass

    cdef cppclass UInt32Type(DataType):
        pass

    cdef cppclass UInt64Type(DataType):
        pass

    cdef cppclass FloatType(DataType):
        pass

    cdef cppclass DoubleType(DataType):
        pass

    cdef cppclass PyObjectType(DataType):
        pass

    cdef cppclass CategoryType(DataType):
        pass

    cdef cppclass cArray" pandas::Array":

        const TypePtr& type()
        TypeEnum type_enum()
        size_t length()

    cdef cppclass cCategoryArray" pandas::CategoryArray"(cArray):
        pass

    cdef cppclass cBooleanArray" pandas::BooleanArray"(cArray):
        pass

    ctypedef shared_ptr[cArray] ArrayPtr

    Status numpy_type_num_to_pandas(int type_num, TypeEnum* pandas_type)
    Status primitive_type_from_enum(TypeEnum tp_enum, DataType** out)

    Status array_from_numpy(PyObject* arr, cArray** out)
    Status array_from_masked_numpy(PyObject* arr, cArray** out)

    void init_numpy()


cdef extern from "pandas/pytypes.h" namespace "pandas::py":
    void init_natype(object type_obj)
    c_bool is_na(object type_obj)
