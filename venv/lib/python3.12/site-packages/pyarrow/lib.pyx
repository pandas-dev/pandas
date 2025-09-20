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

# cython: profile = False
# cython: nonecheck = True
# distutils: language = c++

import datetime
import decimal as _pydecimal
try:
    import numpy as np
except ImportError:
    np = None
import os
import sys

from cython.operator cimport dereference as deref
from pyarrow.includes.libarrow cimport *
from pyarrow.includes.libarrow_python cimport *
from pyarrow.includes.common cimport PyObject_to_object
cimport pyarrow.includes.libarrow_python as libarrow_python
cimport cpython as cp


# Initialize NumPy C API only if numpy was able to be imported
if np is not None:
    arrow_init_numpy()

# Initialize PyArrow C++ API
# (used from some of our C++ code, see e.g. ARROW-5260)
import_pyarrow()


MonthDayNano = NewMonthDayNanoTupleType()


def cpu_count():
    """
    Return the number of threads to use in parallel operations.

    The number of threads is determined at startup by inspecting the
    ``OMP_NUM_THREADS`` and ``OMP_THREAD_LIMIT`` environment variables.
    If neither is present, it will default to the number of hardware threads
    on the system. It can be modified at runtime by calling
    :func:`set_cpu_count()`.

    See Also
    --------
    set_cpu_count : Modify the size of this pool.
    io_thread_count : The analogous function for the I/O thread pool.
    """
    return GetCpuThreadPoolCapacity()


def set_cpu_count(int count):
    """
    Set the number of threads to use in parallel operations.

    Parameters
    ----------
    count : int
        The number of concurrent threads that should be used.

    See Also
    --------
    cpu_count : Get the size of this pool.
    set_io_thread_count : The analogous function for the I/O thread pool.
    """
    if count < 1:
        raise ValueError("CPU count must be strictly positive")
    check_status(SetCpuThreadPoolCapacity(count))


def is_threading_enabled() -> bool:
    """
    Returns True if threading is enabled in libarrow.

    If it isn't enabled, then python shouldn't create any
    threads either, because we're probably on a system where
    threading doesn't work (e.g. Emscripten).
    """
    return libarrow_python.IsThreadingEnabled()


Type_NA = _Type_NA
Type_BOOL = _Type_BOOL
Type_UINT8 = _Type_UINT8
Type_INT8 = _Type_INT8
Type_UINT16 = _Type_UINT16
Type_INT16 = _Type_INT16
Type_UINT32 = _Type_UINT32
Type_INT32 = _Type_INT32
Type_UINT64 = _Type_UINT64
Type_INT64 = _Type_INT64
Type_HALF_FLOAT = _Type_HALF_FLOAT
Type_FLOAT = _Type_FLOAT
Type_DOUBLE = _Type_DOUBLE
Type_DECIMAL32 = _Type_DECIMAL32
Type_DECIMAL64 = _Type_DECIMAL64
Type_DECIMAL128 = _Type_DECIMAL128
Type_DECIMAL256 = _Type_DECIMAL256
Type_DATE32 = _Type_DATE32
Type_DATE64 = _Type_DATE64
Type_TIMESTAMP = _Type_TIMESTAMP
Type_TIME32 = _Type_TIME32
Type_TIME64 = _Type_TIME64
Type_DURATION = _Type_DURATION
Type_INTERVAL_MONTH_DAY_NANO = _Type_INTERVAL_MONTH_DAY_NANO
Type_BINARY = _Type_BINARY
Type_STRING = _Type_STRING
Type_LARGE_BINARY = _Type_LARGE_BINARY
Type_LARGE_STRING = _Type_LARGE_STRING
Type_FIXED_SIZE_BINARY = _Type_FIXED_SIZE_BINARY
Type_BINARY_VIEW = _Type_BINARY_VIEW
Type_STRING_VIEW = _Type_STRING_VIEW
Type_LIST = _Type_LIST
Type_LARGE_LIST = _Type_LARGE_LIST
Type_LIST_VIEW = _Type_LIST_VIEW
Type_LARGE_LIST_VIEW = _Type_LARGE_LIST_VIEW
Type_MAP = _Type_MAP
Type_FIXED_SIZE_LIST = _Type_FIXED_SIZE_LIST
Type_STRUCT = _Type_STRUCT
Type_SPARSE_UNION = _Type_SPARSE_UNION
Type_DENSE_UNION = _Type_DENSE_UNION
Type_DICTIONARY = _Type_DICTIONARY
Type_RUN_END_ENCODED = _Type_RUN_END_ENCODED

UnionMode_SPARSE = _UnionMode_SPARSE
UnionMode_DENSE = _UnionMode_DENSE

__pc = None
__pac = None
__cuda_loaded = None


def _pc():
    global __pc
    if __pc is None:
        import pyarrow.compute as pc
        __pc = pc
    return __pc


def _pac():
    global __pac
    if __pac is None:
        import pyarrow.acero as pac
        __pac = pac
    return __pac


def _ensure_cuda_loaded():
    # Try importing the cuda module to ensure libarrow_cuda gets loaded
    # to register the CUDA device for the C Data Interface import
    global __cuda_loaded
    if __cuda_loaded is None:
        try:
            import pyarrow.cuda  # no-cython-lint
            __cuda_loaded = True
        except ImportError as exc:
            __cuda_loaded = str(exc)

    if __cuda_loaded is not True:
        raise ImportError(
            "Trying to import data on a CUDA device, but PyArrow is not built with "
            f"CUDA support.\n(importing 'pyarrow.cuda' resulted in \"{__cuda_loaded}\")."
        )


def _gdb_test_session():
    GdbTestSession()


# Assorted compatibility helpers
include "compat.pxi"

# Exception types and Status handling
include "error.pxi"

# Configuration information
include "config.pxi"

# pandas API shim
include "pandas-shim.pxi"

# Memory pools and allocation
include "memory.pxi"

# Device type and memory manager
include "device.pxi"

# DataType, Field, Schema
include "types.pxi"

# Array scalar values
include "scalar.pxi"

# Array types
include "array.pxi"

# Builders
include "builder.pxi"

# Column, Table, Record Batch
include "table.pxi"

# Tensors
include "tensor.pxi"

# DLPack
include "_dlpack.pxi"

# File IO
include "io.pxi"

# IPC / Messaging
include "ipc.pxi"

# Micro-benchmark routines
include "benchmark.pxi"

# Public API
include "public-api.pxi"
