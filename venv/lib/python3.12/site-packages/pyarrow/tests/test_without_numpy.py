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

import pytest

import pyarrow as pa

# Marks all of the tests in this module
# Ignore these with pytest ... -m 'not nonumpy'
pytestmark = pytest.mark.nonumpy


def test_array_to_np():
    arr = pa.array(range(10))

    msg = "Cannot return a numpy.ndarray if NumPy is not present"

    with pytest.raises(ImportError, match=msg):
        arr.to_numpy()


def test_chunked_array_to_np():
    data = pa.chunked_array([
        [1, 2, 3],
        [4, 5, 6],
        []
    ])
    msg = "Cannot return a numpy.ndarray if NumPy is not present"

    with pytest.raises(ImportError, match=msg):
        data.to_numpy()


def test_tensor_to_np():
    tensor_type = pa.fixed_shape_tensor(pa.int32(), [2, 2])
    arr = [[1, 2, 3, 4], [10, 20, 30, 40], [100, 200, 300, 400]]
    storage = pa.array(arr, pa.list_(pa.int32(), 4))
    tensor_array = pa.ExtensionArray.from_storage(tensor_type, storage)

    tensor = tensor_array.to_tensor()
    msg = "Cannot return a numpy.ndarray if NumPy is not present"

    with pytest.raises(ImportError, match=msg):
        tensor.to_numpy()
