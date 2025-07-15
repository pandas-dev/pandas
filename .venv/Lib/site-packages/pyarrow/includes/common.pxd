# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.

# distutils: language = c++

from libc.stdint cimport *

from libcpp cimport bool as c_bool, nullptr
from libcpp.functional cimport function
from libcpp.memory cimport (shared_ptr, unique_ptr, make_shared,
                            static_pointer_cast, dynamic_pointer_cast)
from libcpp.optional cimport nullopt, optional
from libcpp.string cimport string as c_string
from libcpp.utility cimport move, pair
from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map
from libcpp.unordered_set cimport unordered_set

from cpython cimport PyObject
from cpython.datetime cimport PyDateTime_DateTime
cimport cpython


cdef extern from "<string_view>" namespace "std" nogil:
    # Needed until https://github.com/cython/cython/issues/6651 is fixed
    cdef cppclass cpp_string_view "std::string_view":
        string_view()
        string_view(const char*)
        string_view(c_string&)
        size_t size()
        bint empty()
        const char* data()


cdef extern from * namespace "arrow::py" nogil:
    """
    #include <memory>
    #include <string>
    #include <string_view>
    #include <utility>

    namespace arrow {
    namespace py {

    template <typename T>
    std::shared_ptr<T> to_shared(std::unique_ptr<T>& t) {
        return std::move(t);
    }
    template <typename T>
    std::shared_ptr<T> to_shared(std::unique_ptr<T>&& t) {
        return std::move(t);
    }

    // Needed until https://github.com/cython/cython/issues/6651 is fixed
    inline std::string to_string(std::string_view s) {
        return std::string(s);
    }

    }  // namespace py
    }  // namespace arrow
    """
    cdef shared_ptr[T] to_shared" arrow::py::to_shared"[T](unique_ptr[T])
    cdef c_string to_string(cpp_string_view s)

cdef extern from "arrow/python/platform.h":
    pass

cdef extern from "<Python.h>":
    void Py_XDECREF(PyObject* o)
    Py_ssize_t Py_REFCNT(PyObject* o)

cdef extern from "numpy/halffloat.h":
    ctypedef uint16_t npy_half

cdef extern from "arrow/api.h" namespace "arrow" nogil:
    # We can later add more of the common status factory methods as needed
    cdef CStatus CStatus_OK "arrow::Status::OK"()

    cdef CStatus CStatus_Invalid "arrow::Status::Invalid"()
    cdef CStatus CStatus_NotImplemented \
        "arrow::Status::NotImplemented"(const c_string& msg)
    cdef CStatus CStatus_UnknownError \
        "arrow::Status::UnknownError"(const c_string& msg)

    cdef cppclass CStatus "arrow::Status":
        CStatus()

        c_string ToString()
        c_string message()
        shared_ptr[CStatusDetail] detail()

        c_bool ok()
        c_bool IsIOError()
        c_bool IsOutOfMemory()
        c_bool IsInvalid()
        c_bool IsKeyError()
        c_bool IsNotImplemented()
        c_bool IsTypeError()
        c_bool IsCapacityError()
        c_bool IsIndexError()
        c_bool IsSerializationError()
        c_bool IsCancelled()

        void Warn()

    cdef cppclass CStatusDetail "arrow::StatusDetail":
        c_string ToString()


cdef extern from "arrow/result.h" namespace "arrow" nogil:
    cdef cppclass CResult "arrow::Result"[T]:
        CResult()
        CResult(CStatus)
        CResult(T)
        c_bool ok()
        CStatus status()
        CStatus Value(T*)
        T operator*()


cdef extern from "arrow/util/future.h" namespace "arrow" nogil:
    cdef cppclass CFuture "arrow::Future"[T]:
        CFuture()


cdef extern from "arrow/python/async.h" namespace "arrow::py" nogil:
    # BindFuture's third argument is really a C++ callable with
    # the signature `object(T*)`, but Cython does not allow declaring that.
    # We use an ellipsis as a workaround.
    # Another possibility is to type-erase the argument by making it
    # `object(void*)`, but it would lose compile-time C++ type safety.
    void BindFuture[T](CFuture[T], object cb, ...)


cdef extern from "arrow/python/common.h" namespace "arrow::py" nogil:
    T GetResultValue[T](CResult[T]) except *
    cdef function[F] BindFunction[F](void* unbound, object bound, ...)


cdef inline object PyObject_to_object(PyObject* o):
    # Cast to "object" increments reference count
    cdef object result = <object> o
    cpython.Py_DECREF(result)
    return result
