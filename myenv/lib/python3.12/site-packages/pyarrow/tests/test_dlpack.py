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

import ctypes
from functools import wraps
import pytest

import numpy as np

import pyarrow as pa
from pyarrow.vendored.version import Version


def PyCapsule_IsValid(capsule, name):
    return ctypes.pythonapi.PyCapsule_IsValid(ctypes.py_object(capsule), name) == 1


def check_dlpack_export(arr, expected_arr):
    DLTensor = arr.__dlpack__()
    assert PyCapsule_IsValid(DLTensor, b"dltensor") is True

    result = np.from_dlpack(arr)
    np.testing.assert_array_equal(result, expected_arr, strict=True)

    assert arr.__dlpack_device__() == (1, 0)


def check_bytes_allocated(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        allocated_bytes = pa.total_allocated_bytes()
        try:
            return f(*args, **kwargs)
        finally:
            assert pa.total_allocated_bytes() == allocated_bytes
    return wrapper


@check_bytes_allocated
@pytest.mark.parametrize(
    ('value_type', 'np_type'),
    [
        (pa.uint8(), np.uint8),
        (pa.uint16(), np.uint16),
        (pa.uint32(), np.uint32),
        (pa.uint64(), np.uint64),
        (pa.int8(), np.int8),
        (pa.int16(), np.int16),
        (pa.int32(), np.int32),
        (pa.int64(), np.int64),
        (pa.float16(), np.float16),
        (pa.float32(), np.float32),
        (pa.float64(), np.float64),
    ]
)
def test_dlpack(value_type, np_type):
    if Version(np.__version__) < Version("1.24.0"):
        pytest.skip("No dlpack support in numpy versions older than 1.22.0, "
                    "strict keyword in assert_array_equal added in numpy version "
                    "1.24.0")

    expected = np.array([1, 2, 3], dtype=np_type)
    arr = pa.array(expected, type=value_type)
    check_dlpack_export(arr, expected)

    arr_sliced = arr.slice(1, 1)
    expected = np.array([2], dtype=np_type)
    check_dlpack_export(arr_sliced, expected)

    arr_sliced = arr.slice(0, 1)
    expected = np.array([1], dtype=np_type)
    check_dlpack_export(arr_sliced, expected)

    arr_sliced = arr.slice(1)
    expected = np.array([2, 3], dtype=np_type)
    check_dlpack_export(arr_sliced, expected)

    arr_zero = pa.array([], type=value_type)
    expected = np.array([], dtype=np_type)
    check_dlpack_export(arr_zero, expected)


def test_dlpack_not_supported():
    if Version(np.__version__) < Version("1.22.0"):
        pytest.skip("No dlpack support in numpy versions older than 1.22.0.")

    arr = pa.array([1, None, 3])
    with pytest.raises(TypeError, match="Can only use DLPack "
                       "on arrays with no nulls."):
        np.from_dlpack(arr)

    arr = pa.array(
        [[0, 1], [3, 4]],
        type=pa.list_(pa.int32())
    )
    with pytest.raises(TypeError, match="DataType is not compatible with DLPack spec"):
        np.from_dlpack(arr)

    arr = pa.array([])
    with pytest.raises(TypeError, match="DataType is not compatible with DLPack spec"):
        np.from_dlpack(arr)

    # DLPack doesn't support bit-packed boolean values
    arr = pa.array([True, False, True])
    with pytest.raises(TypeError, match="Bit-packed boolean data type "
                       "not supported by DLPack."):
        np.from_dlpack(arr)


def test_dlpack_cuda_not_supported():
    cuda = pytest.importorskip("pyarrow.cuda")

    schema = pa.schema([pa.field('f0', pa.int16())])
    a0 = pa.array([1, 2, 3], type=pa.int16())
    batch = pa.record_batch([a0], schema=schema)

    cbuf = cuda.serialize_record_batch(batch, cuda.Context(0))
    cbatch = cuda.read_record_batch(cbuf, batch.schema)
    carr = cbatch["f0"]

    # CudaBuffers not yet supported
    with pytest.raises(NotImplementedError, match="DLPack support is implemented "
                       "only for buffers on CPU device."):
        np.from_dlpack(carr)

    with pytest.raises(NotImplementedError, match="DLPack support is implemented "
                       "only for buffers on CPU device."):
        carr.__dlpack_device__()
