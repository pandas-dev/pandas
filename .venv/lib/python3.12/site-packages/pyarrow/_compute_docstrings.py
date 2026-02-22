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

"""
Custom documentation additions for compute functions.
"""

function_doc_additions = {}

function_doc_additions["filter"] = """
    Examples
    --------
    >>> import pyarrow as pa
    >>> arr = pa.array(["a", "b", "c", None, "e"])
    >>> mask = pa.array([True, False, None, False, True])
    >>> arr.filter(mask)
    <pyarrow.lib.StringArray object at ...>
    [
      "a",
      "e"
    ]
    >>> arr.filter(mask, null_selection_behavior='emit_null')
    <pyarrow.lib.StringArray object at ...>
    [
      "a",
      null,
      "e"
    ]
    """

function_doc_additions["mode"] = """
    Examples
    --------
    >>> import pyarrow as pa
    >>> import pyarrow.compute as pc
    >>> arr = pa.array([1, 1, 2, 2, 3, 2, 2, 2])
    >>> modes = pc.mode(arr, 2)
    >>> modes[0]
    <pyarrow.StructScalar: [('mode', 2), ('count', 5)]>
    >>> modes[1]
    <pyarrow.StructScalar: [('mode', 1), ('count', 2)]>
    """

function_doc_additions["min"] = """
    Examples
    --------
    >>> import pyarrow as pa
    >>> import pyarrow.compute as pc
    >>> arr1 = pa.array([1, 1, 2, 2, 3, 2, 2, 2])
    >>> pc.min(arr1)
    <pyarrow.Int64Scalar: 1>

    Using ``skip_nulls`` to handle null values.

    >>> arr2 = pa.array([1.0, None, 2.0, 3.0])
    >>> pc.min(arr2)
    <pyarrow.DoubleScalar: 1.0>
    >>> pc.min(arr2, skip_nulls=False)
    <pyarrow.DoubleScalar: None>

    Using ``ScalarAggregateOptions`` to control minimum number of non-null values.

    >>> arr3 = pa.array([1.0, None, float("nan"), 3.0])
    >>> pc.min(arr3)
    <pyarrow.DoubleScalar: 1.0>
    >>> pc.min(arr3, options=pc.ScalarAggregateOptions(min_count=3))
    <pyarrow.DoubleScalar: 1.0>
    >>> pc.min(arr3, options=pc.ScalarAggregateOptions(min_count=4))
    <pyarrow.DoubleScalar: None>

    This function also works with string values.

    >>> arr4 = pa.array(["z", None, "y", "x"])
    >>> pc.min(arr4)
    <pyarrow.StringScalar: 'x'>
    """

function_doc_additions["max"] = """
    Examples
    --------
    >>> import pyarrow as pa
    >>> import pyarrow.compute as pc
    >>> arr1 = pa.array([1, 1, 2, 2, 3, 2, 2, 2])
    >>> pc.max(arr1)
    <pyarrow.Int64Scalar: 3>

    Using ``skip_nulls`` to handle null values.

    >>> arr2 = pa.array([1.0, None, 2.0, 3.0])
    >>> pc.max(arr2)
    <pyarrow.DoubleScalar: 3.0>
    >>> pc.max(arr2, skip_nulls=False)
    <pyarrow.DoubleScalar: None>

    Using ``ScalarAggregateOptions`` to control minimum number of non-null values.

    >>> arr3 = pa.array([1.0, None, float("nan"), 3.0])
    >>> pc.max(arr3)
    <pyarrow.DoubleScalar: 3.0>
    >>> pc.max(arr3, options=pc.ScalarAggregateOptions(min_count=3))
    <pyarrow.DoubleScalar: 3.0>
    >>> pc.max(arr3, options=pc.ScalarAggregateOptions(min_count=4))
    <pyarrow.DoubleScalar: None>

    This function also works with string values.

    >>> arr4 = pa.array(["z", None, "y", "x"])
    >>> pc.max(arr4)
    <pyarrow.StringScalar: 'z'>
    """

function_doc_additions["min_max"] = """
    Examples
    --------
    >>> import pyarrow as pa
    >>> import pyarrow.compute as pc
    >>> arr1 = pa.array([1, 1, 2, 2, 3, 2, 2, 2])
    >>> pc.min_max(arr1)
    <pyarrow.StructScalar: [('min', 1), ('max', 3)]>

    Using ``skip_nulls`` to handle null values.

    >>> arr2 = pa.array([1.0, None, 2.0, 3.0])
    >>> pc.min_max(arr2)
    <pyarrow.StructScalar: [('min', 1.0), ('max', 3.0)]>
    >>> pc.min_max(arr2, skip_nulls=False)
    <pyarrow.StructScalar: [('min', None), ('max', None)]>

    Using ``ScalarAggregateOptions`` to control minimum number of non-null values.

    >>> arr3 = pa.array([1.0, None, float("nan"), 3.0])
    >>> pc.min_max(arr3)
    <pyarrow.StructScalar: [('min', 1.0), ('max', 3.0)]>
    >>> pc.min_max(arr3, options=pc.ScalarAggregateOptions(min_count=3))
    <pyarrow.StructScalar: [('min', 1.0), ('max', 3.0)]>
    >>> pc.min_max(arr3, options=pc.ScalarAggregateOptions(min_count=4))
    <pyarrow.StructScalar: [('min', None), ('max', None)]>

    This function also works with string values.

    >>> arr4 = pa.array(["z", None, "y", "x"])
    >>> pc.min_max(arr4)
    <pyarrow.StructScalar: [('min', 'x'), ('max', 'z')]>
    """
