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

import numpy as np

import pyarrow as pa
from pyarrow import compute as pc

# UDFs are all tested with a dataset scan
pytestmark = pytest.mark.dataset

# For convenience, most of the test here doesn't care about udf func docs
empty_udf_doc = {"summary": "", "description": ""}

try:
    import pyarrow.dataset as ds
except ImportError:
    ds = None


def mock_udf_context(batch_length=10):
    from pyarrow._compute import _get_udf_context
    return _get_udf_context(pa.default_memory_pool(), batch_length)


class MyError(RuntimeError):
    pass


@pytest.fixture(scope="session")
def sum_agg_func_fixture():
    """
    Register a unary aggregate function (mean)
    """
    def func(ctx, x, *args):
        return pa.scalar(np.nansum(x))

    func_name = "sum_udf"
    func_doc = empty_udf_doc

    pc.register_aggregate_function(func,
                                   func_name,
                                   func_doc,
                                   {
                                       "x": pa.float64(),
                                   },
                                   pa.float64()
                                   )
    return func, func_name


@pytest.fixture(scope="session")
def exception_agg_func_fixture():
    def func(ctx, x):
        raise RuntimeError("Oops")
        return pa.scalar(len(x))

    func_name = "y=exception_len(x)"
    func_doc = empty_udf_doc

    pc.register_aggregate_function(func,
                                   func_name,
                                   func_doc,
                                   {
                                       "x": pa.int64(),
                                   },
                                   pa.int64()
                                   )
    return func, func_name


@pytest.fixture(scope="session")
def wrong_output_dtype_agg_func_fixture(scope="session"):
    def func(ctx, x):
        return pa.scalar(len(x), pa.int32())

    func_name = "y=wrong_output_dtype(x)"
    func_doc = empty_udf_doc

    pc.register_aggregate_function(func,
                                   func_name,
                                   func_doc,
                                   {
                                       "x": pa.int64(),
                                   },
                                   pa.int64()
                                   )
    return func, func_name


@pytest.fixture(scope="session")
def wrong_output_type_agg_func_fixture(scope="session"):
    def func(ctx, x):
        return len(x)

    func_name = "y=wrong_output_type(x)"
    func_doc = empty_udf_doc

    pc.register_aggregate_function(func,
                                   func_name,
                                   func_doc,
                                   {
                                       "x": pa.int64(),
                                   },
                                   pa.int64()
                                   )
    return func, func_name


@pytest.fixture(scope="session")
def binary_func_fixture():
    """
    Register a binary scalar function.
    """
    def binary_function(ctx, m, x):
        return pc.call_function("multiply", [m, x],
                                memory_pool=ctx.memory_pool)
    func_name = "y=mx"
    binary_doc = {"summary": "y=mx",
                  "description": "find y from y = mx"}
    pc.register_scalar_function(binary_function,
                                func_name,
                                binary_doc,
                                {"m": pa.int64(),
                                 "x": pa.int64(),
                                 },
                                pa.int64())
    return binary_function, func_name


@pytest.fixture(scope="session")
def ternary_func_fixture():
    """
    Register a ternary scalar function.
    """
    def ternary_function(ctx, m, x, c):
        mx = pc.call_function("multiply", [m, x],
                              memory_pool=ctx.memory_pool)
        return pc.call_function("add", [mx, c],
                                memory_pool=ctx.memory_pool)
    ternary_doc = {"summary": "y=mx+c",
                   "description": "find y from y = mx + c"}
    func_name = "y=mx+c"
    pc.register_scalar_function(ternary_function,
                                func_name,
                                ternary_doc,
                                {
                                    "array1": pa.int64(),
                                    "array2": pa.int64(),
                                    "array3": pa.int64(),
                                },
                                pa.int64())
    return ternary_function, func_name


@pytest.fixture(scope="session")
def varargs_func_fixture():
    """
    Register a varargs scalar function with at least two arguments.
    """
    def varargs_function(ctx, first, *values):
        acc = first
        for val in values:
            acc = pc.call_function("add", [acc, val],
                                   memory_pool=ctx.memory_pool)
        return acc
    func_name = "z=ax+by+c"
    varargs_doc = {"summary": "z=ax+by+c",
                   "description": "find z from z = ax + by + c"
                   }
    pc.register_scalar_function(varargs_function,
                                func_name,
                                varargs_doc,
                                {
                                    "array1": pa.int64(),
                                    "array2": pa.int64(),
                                },
                                pa.int64())
    return varargs_function, func_name


@pytest.fixture(scope="session")
def nullary_func_fixture():
    """
    Register a nullary scalar function.
    """
    def nullary_func(context):
        return pa.array([42] * context.batch_length, type=pa.int64(),
                        memory_pool=context.memory_pool)

    func_doc = {
        "summary": "random function",
        "description": "generates a random value"
    }
    func_name = "test_nullary_func"
    pc.register_scalar_function(nullary_func,
                                func_name,
                                func_doc,
                                {},
                                pa.int64())

    return nullary_func, func_name


@pytest.fixture(scope="session")
def wrong_output_type_func_fixture():
    """
    Register a scalar function which returns something that is neither
    a Arrow scalar or array.
    """
    def wrong_output_type(ctx):
        return 42

    func_name = "test_wrong_output_type"
    in_types = {}
    out_type = pa.int64()
    doc = {
        "summary": "return wrong output type",
        "description": ""
    }
    pc.register_scalar_function(wrong_output_type, func_name, doc,
                                in_types, out_type)
    return wrong_output_type, func_name


@pytest.fixture(scope="session")
def wrong_output_datatype_func_fixture():
    """
    Register a scalar function whose actual output DataType doesn't
    match the declared output DataType.
    """
    def wrong_output_datatype(ctx, array):
        return pc.call_function("add", [array, 1])
    func_name = "test_wrong_output_datatype"
    in_types = {"array": pa.int64()}
    # The actual output DataType will be int64.
    out_type = pa.int16()
    doc = {
        "summary": "return wrong output datatype",
        "description": ""
    }
    pc.register_scalar_function(wrong_output_datatype, func_name, doc,
                                in_types, out_type)
    return wrong_output_datatype, func_name


@pytest.fixture(scope="session")
def wrong_signature_func_fixture():
    """
    Register a scalar function with the wrong signature.
    """
    # Missing the context argument
    def wrong_signature():
        return pa.scalar(1, type=pa.int64())

    func_name = "test_wrong_signature"
    in_types = {}
    out_type = pa.int64()
    doc = {
        "summary": "UDF with wrong signature",
        "description": ""
    }
    pc.register_scalar_function(wrong_signature, func_name, doc,
                                in_types, out_type)
    return wrong_signature, func_name


@pytest.fixture(scope="session")
def raising_func_fixture():
    """
    Register a scalar function which raises a custom exception.
    """
    def raising_func(ctx):
        raise MyError("error raised by scalar UDF")
    func_name = "test_raise"
    doc = {
        "summary": "raising function",
        "description": ""
    }
    pc.register_scalar_function(raising_func, func_name, doc,
                                {}, pa.int64())
    return raising_func, func_name


@pytest.fixture(scope="session")
def unary_vector_func_fixture():
    """
    Register a vector function
    """
    def pct_rank(ctx, x):
        # copy here to get around pandas 1.0 issue
        return pa.array(x.to_pandas().copy().rank(pct=True))

    func_name = "y=pct_rank(x)"
    doc = empty_udf_doc
    pc.register_vector_function(pct_rank, func_name, doc, {
                                'x': pa.float64()}, pa.float64())

    return pct_rank, func_name


@pytest.fixture(scope="session")
def struct_vector_func_fixture():
    """
    Register a vector function that returns a struct array
    """
    def pivot(ctx, k, v, c):
        df = pa.RecordBatch.from_arrays([k, v, c], names=['k', 'v', 'c']).to_pandas()
        df_pivot = df.pivot(columns='c', values='v', index='k').reset_index()
        return pa.RecordBatch.from_pandas(df_pivot).to_struct_array()

    func_name = "y=pivot(x)"
    doc = empty_udf_doc
    pc.register_vector_function(
        pivot, func_name, doc,
        {'k': pa.int64(), 'v': pa.float64(), 'c': pa.utf8()},
        pa.struct([('k', pa.int64()), ('v1', pa.float64()), ('v2', pa.float64())])
    )

    return pivot, func_name


def check_scalar_function(func_fixture,
                          inputs, *,
                          run_in_dataset=True,
                          batch_length=None):
    function, name = func_fixture
    if batch_length is None:
        all_scalar = True
        for arg in inputs:
            if isinstance(arg, pa.Array):
                all_scalar = False
                batch_length = len(arg)
        if all_scalar:
            batch_length = 1

    func = pc.get_function(name)
    assert func.name == name

    result = pc.call_function(name, inputs, length=batch_length)
    expected_output = function(mock_udf_context(batch_length), *inputs)
    assert result == expected_output
    # At the moment there is an issue when handling nullary functions.
    # See: ARROW-15286 and ARROW-16290.
    if run_in_dataset:
        field_names = [f'field{index}' for index, in_arr in inputs]
        table = pa.Table.from_arrays(inputs, field_names)
        dataset = ds.dataset(table)
        func_args = [ds.field(field_name) for field_name in field_names]
        result_table = dataset.to_table(
            columns={'result': ds.field('')._call(name, func_args)})
        assert result_table.column(0).chunks[0] == expected_output


def test_udf_array_unary(unary_func_fixture):
    check_scalar_function(unary_func_fixture,
                          [
                              pa.array([10, 20], pa.int64())
                          ]
                          )


def test_udf_array_binary(binary_func_fixture):
    check_scalar_function(binary_func_fixture,
                          [
                              pa.array([10, 20], pa.int64()),
                              pa.array([2, 4], pa.int64())
                          ]
                          )


def test_udf_array_ternary(ternary_func_fixture):
    check_scalar_function(ternary_func_fixture,
                          [
                              pa.array([10, 20], pa.int64()),
                              pa.array([2, 4], pa.int64()),
                              pa.array([5, 10], pa.int64())
                          ]
                          )


def test_udf_array_varargs(varargs_func_fixture):
    check_scalar_function(varargs_func_fixture,
                          [
                              pa.array([2, 3], pa.int64()),
                              pa.array([10, 20], pa.int64()),
                              pa.array([3, 7], pa.int64()),
                              pa.array([20, 30], pa.int64()),
                              pa.array([5, 10], pa.int64())
                          ]
                          )


def test_registration_errors():
    # validate function name
    doc = {
        "summary": "test udf input",
        "description": "parameters are validated"
    }
    in_types = {"scalar": pa.int64()}
    out_type = pa.int64()

    def test_reg_function(context):
        return pa.array([10])

    with pytest.raises(TypeError):
        pc.register_scalar_function(test_reg_function,
                                    None, doc, in_types,
                                    out_type)

    # validate function
    with pytest.raises(TypeError, match="func must be a callable"):
        pc.register_scalar_function(None, "test_none_function", doc, in_types,
                                    out_type)

    # validate output type
    expected_expr = "DataType expected, got <class 'NoneType'>"
    with pytest.raises(TypeError, match=expected_expr):
        pc.register_scalar_function(test_reg_function,
                                    "test_output_function", doc, in_types,
                                    None)

    # validate input type
    expected_expr = "in_types must be a dictionary of DataType"
    with pytest.raises(TypeError, match=expected_expr):
        pc.register_scalar_function(test_reg_function,
                                    "test_input_function", doc, None,
                                    out_type)

    # register an already registered function
    # first registration
    pc.register_scalar_function(test_reg_function,
                                "test_reg_function", doc, {},
                                out_type)
    # second registration
    expected_expr = "Already have a function registered with name:" \
        + " test_reg_function"
    with pytest.raises(KeyError, match=expected_expr):
        pc.register_scalar_function(test_reg_function,
                                    "test_reg_function", doc, {},
                                    out_type)


def test_varargs_function_validation(varargs_func_fixture):
    _, func_name = varargs_func_fixture

    error_msg = r"VarArgs function 'z=ax\+by\+c' needs at least 2 arguments"

    with pytest.raises(ValueError, match=error_msg):
        pc.call_function(func_name, [42])


def test_function_doc_validation():
    # validate arity
    in_types = {"scalar": pa.int64()}
    out_type = pa.int64()

    # doc with no summary
    func_doc = {
        "description": "desc"
    }

    def add_const(ctx, scalar):
        return pc.call_function("add", [scalar, 1])

    with pytest.raises(ValueError,
                       match="Function doc must contain a summary"):
        pc.register_scalar_function(add_const, "test_no_summary",
                                    func_doc, in_types,
                                    out_type)

    # doc with no description
    func_doc = {
        "summary": "test summary"
    }

    with pytest.raises(ValueError,
                       match="Function doc must contain a description"):
        pc.register_scalar_function(add_const, "test_no_desc",
                                    func_doc, in_types,
                                    out_type)


def test_nullary_function(nullary_func_fixture):
    # XXX the Python compute layer API doesn't let us override batch_length,
    # so only test with the default value of 1.
    check_scalar_function(nullary_func_fixture, [], run_in_dataset=False,
                          batch_length=1)


def test_wrong_output_type(wrong_output_type_func_fixture):
    _, func_name = wrong_output_type_func_fixture

    with pytest.raises(TypeError,
                       match="Unexpected output type: int"):
        pc.call_function(func_name, [], length=1)


def test_wrong_output_datatype(wrong_output_datatype_func_fixture):
    _, func_name = wrong_output_datatype_func_fixture

    expected_expr = ("Expected output datatype int16, "
                     "but function returned datatype int64")

    with pytest.raises(TypeError, match=expected_expr):
        pc.call_function(func_name, [pa.array([20, 30])])


def test_wrong_signature(wrong_signature_func_fixture):
    _, func_name = wrong_signature_func_fixture

    expected_expr = (r"wrong_signature\(\) takes 0 positional arguments "
                     "but 1 was given")

    with pytest.raises(TypeError, match=expected_expr):
        pc.call_function(func_name, [], length=1)


def test_wrong_datatype_declaration():
    def identity(ctx, val):
        return val

    func_name = "test_wrong_datatype_declaration"
    in_types = {"array": pa.int64()}
    out_type = {}
    doc = {
        "summary": "test output value",
        "description": "test output"
    }
    with pytest.raises(TypeError,
                       match="DataType expected, got <class 'dict'>"):
        pc.register_scalar_function(identity, func_name,
                                    doc, in_types, out_type)


def test_wrong_input_type_declaration():
    def identity(ctx, val):
        return val

    func_name = "test_wrong_input_type_declaration"
    in_types = {"array": None}
    out_type = pa.int64()
    doc = {
        "summary": "test invalid input type",
        "description": "invalid input function"
    }
    with pytest.raises(TypeError,
                       match="DataType expected, got <class 'NoneType'>"):
        pc.register_scalar_function(identity, func_name, doc,
                                    in_types, out_type)


def test_scalar_udf_context(unary_func_fixture):
    # Check the memory_pool argument is properly propagated
    proxy_pool = pa.proxy_memory_pool(pa.default_memory_pool())
    _, func_name = unary_func_fixture

    res = pc.call_function(func_name,
                           [pa.array([1] * 1000, type=pa.int64())],
                           memory_pool=proxy_pool)
    assert res == pa.array([2] * 1000, type=pa.int64())
    assert proxy_pool.bytes_allocated() == 1000 * 8
    # Destroying Python array should destroy underlying C++ memory
    res = None
    assert proxy_pool.bytes_allocated() == 0


def test_raising_func(raising_func_fixture):
    _, func_name = raising_func_fixture
    with pytest.raises(MyError, match="error raised by scalar UDF"):
        pc.call_function(func_name, [], length=1)


def test_scalar_input(unary_func_fixture):
    function, func_name = unary_func_fixture
    res = pc.call_function(func_name, [pa.scalar(10)])
    assert res == pa.scalar(11)


def test_input_lifetime(unary_func_fixture):
    function, func_name = unary_func_fixture

    proxy_pool = pa.proxy_memory_pool(pa.default_memory_pool())
    assert proxy_pool.bytes_allocated() == 0

    v = pa.array([1] * 1000, type=pa.int64(), memory_pool=proxy_pool)
    assert proxy_pool.bytes_allocated() == 1000 * 8
    pc.call_function(func_name, [v])
    assert proxy_pool.bytes_allocated() == 1000 * 8
    # Calling a UDF should not have kept `v` alive longer than required
    v = None
    assert proxy_pool.bytes_allocated() == 0


def _record_batch_from_iters(schema, *iters):
    arrays = [pa.array(list(v), type=schema[i].type)
              for i, v in enumerate(iters)]
    return pa.RecordBatch.from_arrays(arrays=arrays, schema=schema)


def _record_batch_for_range(schema, n):
    return _record_batch_from_iters(schema,
                                    range(n, n + 10),
                                    range(n + 1, n + 11))


def make_udt_func(schema, batch_gen):
    def udf_func(ctx):
        class UDT:
            def __init__(self):
                self.caller = None

            def __call__(self, ctx):
                try:
                    if self.caller is None:
                        self.caller, ctx = batch_gen(ctx).send, None
                    batch = self.caller(ctx)
                except StopIteration:
                    arrays = [pa.array([], type=field.type)
                              for field in schema]
                    batch = pa.RecordBatch.from_arrays(
                        arrays=arrays, schema=schema)
                return batch.to_struct_array()
        return UDT()
    return udf_func


def datasource1_direct():
    """A short dataset"""
    schema = datasource1_schema()

    class Generator:
        def __init__(self):
            self.n = 3

        def __call__(self, ctx):
            if self.n == 0:
                batch = _record_batch_from_iters(schema, [], [])
            else:
                self.n -= 1
                batch = _record_batch_for_range(schema, self.n)
            return batch.to_struct_array()
    return lambda ctx: Generator()


def datasource1_generator():
    schema = datasource1_schema()

    def batch_gen(ctx):
        for n in range(3, 0, -1):
            # ctx =
            yield _record_batch_for_range(schema, n - 1)
    return make_udt_func(schema, batch_gen)


def datasource1_exception():
    schema = datasource1_schema()

    def batch_gen(ctx):
        for n in range(3, 0, -1):
            # ctx =
            yield _record_batch_for_range(schema, n - 1)
        raise RuntimeError("datasource1_exception")
    return make_udt_func(schema, batch_gen)


def datasource1_schema():
    return pa.schema([('', pa.int32()), ('', pa.int32())])


def datasource1_args(func, func_name):
    func_doc = {"summary": f"{func_name} UDT",
                "description": "test {func_name} UDT"}
    in_types = {}
    out_type = pa.struct([("", pa.int32()), ("", pa.int32())])
    return func, func_name, func_doc, in_types, out_type


def _test_datasource1_udt(func_maker):
    schema = datasource1_schema()
    func = func_maker()
    func_name = func_maker.__name__
    func_args = datasource1_args(func, func_name)
    pc.register_tabular_function(*func_args)
    n = 3
    for item in pc.call_tabular_function(func_name):
        n -= 1
        assert item == _record_batch_for_range(schema, n)


def test_udt_datasource1_direct():
    _test_datasource1_udt(datasource1_direct)


def test_udt_datasource1_generator():
    _test_datasource1_udt(datasource1_generator)


def test_udt_datasource1_exception():
    with pytest.raises(RuntimeError, match='datasource1_exception'):
        _test_datasource1_udt(datasource1_exception)


def test_scalar_agg_basic(unary_agg_func_fixture):
    arr = pa.array([10.0, 20.0, 30.0, 40.0, 50.0], pa.float64())
    result = pc.call_function("mean_udf", [arr])
    expected = pa.scalar(30.0)
    assert result == expected


def test_scalar_agg_empty(unary_agg_func_fixture):
    empty = pa.array([], pa.float64())

    with pytest.raises(pa.ArrowInvalid, match='empty inputs'):
        pc.call_function("mean_udf", [empty])


def test_scalar_agg_wrong_output_dtype(wrong_output_dtype_agg_func_fixture):
    arr = pa.array([10, 20, 30, 40, 50], pa.int64())
    with pytest.raises(pa.ArrowTypeError, match="output datatype"):
        pc.call_function("y=wrong_output_dtype(x)", [arr])


def test_scalar_agg_wrong_output_type(wrong_output_type_agg_func_fixture):
    arr = pa.array([10, 20, 30, 40, 50], pa.int64())
    with pytest.raises(pa.ArrowTypeError, match="output type"):
        pc.call_function("y=wrong_output_type(x)", [arr])


def test_scalar_agg_varargs(varargs_agg_func_fixture):
    arr1 = pa.array([10, 20, 30, 40, 50], pa.int64())
    arr2 = pa.array([1.0, 2.0, 3.0, 4.0, 5.0], pa.float64())

    result = pc.call_function(
        "sum_mean", [arr1, arr2]
    )
    expected = pa.scalar(33.0)
    assert result == expected


def test_scalar_agg_exception(exception_agg_func_fixture):
    arr = pa.array([10, 20, 30, 40, 50, 60], pa.int64())

    with pytest.raises(RuntimeError, match='Oops'):
        pc.call_function("y=exception_len(x)", [arr])


def test_hash_agg_basic(unary_agg_func_fixture):
    arr1 = pa.array([10.0, 20.0, 30.0, 40.0, 50.0], pa.float64())
    arr2 = pa.array([4, 2, 1, 2, 1], pa.int32())

    arr3 = pa.array([60.0, 70.0, 80.0, 90.0, 100.0], pa.float64())
    arr4 = pa.array([5, 1, 1, 4, 1], pa.int32())

    table1 = pa.table([arr2, arr1], names=["id", "value"])
    table2 = pa.table([arr4, arr3], names=["id", "value"])
    table = pa.concat_tables([table1, table2])

    result = table.group_by("id").aggregate([("value", "mean_udf")])
    expected = table.group_by("id").aggregate(
        [("value", "mean")]).rename_columns(['id', 'value_mean_udf'])

    assert result.sort_by('id') == expected.sort_by('id')


def test_hash_agg_empty(unary_agg_func_fixture):
    arr1 = pa.array([], pa.float64())
    arr2 = pa.array([], pa.int32())
    table = pa.table([arr2, arr1], names=["id", "value"])

    result = table.group_by("id").aggregate([("value", "mean_udf")])
    expected = pa.table([pa.array([], pa.int32()), pa.array(
        [], pa.float64())], names=['id', 'value_mean_udf'])

    assert result == expected


def test_hash_agg_wrong_output_dtype(wrong_output_dtype_agg_func_fixture):
    arr1 = pa.array([10, 20, 30, 40, 50], pa.int64())
    arr2 = pa.array([4, 2, 1, 2, 1], pa.int32())

    table = pa.table([arr2, arr1], names=["id", "value"])
    with pytest.raises(pa.ArrowTypeError, match="output datatype"):
        table.group_by("id").aggregate([("value", "y=wrong_output_dtype(x)")])


def test_hash_agg_wrong_output_type(wrong_output_type_agg_func_fixture):
    arr1 = pa.array([10, 20, 30, 40, 50], pa.int64())
    arr2 = pa.array([4, 2, 1, 2, 1], pa.int32())
    table = pa.table([arr2, arr1], names=["id", "value"])

    with pytest.raises(pa.ArrowTypeError, match="output type"):
        table.group_by("id").aggregate([("value", "y=wrong_output_type(x)")])


def test_hash_agg_exception(exception_agg_func_fixture):
    arr1 = pa.array([10, 20, 30, 40, 50], pa.int64())
    arr2 = pa.array([4, 2, 1, 2, 1], pa.int32())
    table = pa.table([arr2, arr1], names=["id", "value"])

    with pytest.raises(RuntimeError, match='Oops'):
        table.group_by("id").aggregate([("value", "y=exception_len(x)")])


def test_hash_agg_random(sum_agg_func_fixture):
    """Test hash aggregate udf with randomly sampled data"""

    value_num = 1000000
    group_num = 1000

    arr1 = pa.array(np.repeat(1, value_num), pa.float64())
    arr2 = pa.array(np.random.choice(group_num, value_num), pa.int32())

    table = pa.table([arr2, arr1], names=['id', 'value'])

    result = table.group_by("id").aggregate([("value", "sum_udf")])
    expected = table.group_by("id").aggregate(
        [("value", "sum")]).rename_columns(['id', 'value_sum_udf'])

    assert result.sort_by('id') == expected.sort_by('id')


@pytest.mark.pandas
def test_vector_basic(unary_vector_func_fixture):
    arr = pa.array([10.0, 20.0, 30.0, 40.0, 50.0], pa.float64())
    result = pc.call_function("y=pct_rank(x)", [arr])
    expected = unary_vector_func_fixture[0](None, arr)
    assert result == expected


@pytest.mark.pandas
def test_vector_empty(unary_vector_func_fixture):
    arr = pa.array([1], pa.float64())
    result = pc.call_function("y=pct_rank(x)", [arr])
    expected = unary_vector_func_fixture[0](None, arr)
    assert result == expected


@pytest.mark.pandas
def test_vector_struct(struct_vector_func_fixture):
    k = pa.array(
        [1, 1, 2, 2], pa.int64()
    )
    v = pa.array(
        [1.0, 2.0, 3.0, 4.0], pa.float64()
    )
    c = pa.array(
        ['v1', 'v2', 'v1', 'v2']
    )
    result = pc.call_function("y=pivot(x)", [k, v, c])
    expected = struct_vector_func_fixture[0](None, k, v, c)
    assert result == expected
