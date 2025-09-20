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

# cython: language_level = 3
from cython.operator cimport dereference as deref
from libcpp.vector cimport vector as std_vector

from pyarrow import Buffer, py_buffer
from pyarrow._compute cimport Expression
from pyarrow.lib import frombytes, tobytes
from pyarrow.lib cimport *
from pyarrow.includes.libarrow cimport *
from pyarrow.includes.libarrow_substrait cimport *

try:
    import substrait as py_substrait
except ImportError:
    py_substrait = None
else:
    import substrait.proto  # no-cython-lint


# TODO GH-37235: Fix exception handling
cdef CDeclaration _create_named_table_provider(
    dict named_args, const std_vector[c_string]& names, const CSchema& schema
) noexcept:
    cdef:
        c_string c_name
        shared_ptr[CTable] c_in_table
        shared_ptr[CTableSourceNodeOptions] c_tablesourceopts
        shared_ptr[CExecNodeOptions] c_input_node_opts
        vector[CDeclaration.Input] no_c_inputs

    py_names = []
    for i in range(names.size()):
        c_name = names[i]
        py_names.append(frombytes(c_name))
    py_schema = pyarrow_wrap_schema(make_shared[CSchema](schema))

    py_table = named_args["provider"](py_names, py_schema)
    c_in_table = pyarrow_unwrap_table(py_table)
    c_tablesourceopts = make_shared[CTableSourceNodeOptions](c_in_table)
    c_input_node_opts = static_pointer_cast[CExecNodeOptions, CTableSourceNodeOptions](
        c_tablesourceopts)
    return CDeclaration(tobytes("table_source"),
                        no_c_inputs, c_input_node_opts)


def run_query(plan, *, table_provider=None, use_threads=True):
    """
    Execute a Substrait plan and read the results as a RecordBatchReader.

    Parameters
    ----------
    plan : Union[Buffer, bytes]
        The serialized Substrait plan to execute.
    table_provider : object (optional)
        A function to resolve any NamedTable relation to a table.
        The function will receive two arguments which will be a list
        of strings representing the table name and a pyarrow.Schema representing
        the expected schema and should return a pyarrow.Table.
    use_threads : bool, default True
        If True then multiple threads will be used to run the query.  If False then
        all CPU intensive work will be done on the calling thread.

    Returns
    -------
    RecordBatchReader
        A reader containing the result of the executed query

    Examples
    --------
    >>> import pyarrow as pa
    >>> from pyarrow.lib import tobytes
    >>> import pyarrow.substrait as substrait
    >>> test_table_1 = pa.Table.from_pydict({"x": [1, 2, 3]})
    >>> test_table_2 = pa.Table.from_pydict({"x": [4, 5, 6]})
    >>> def table_provider(names, schema):
    ...     if not names:
    ...        raise Exception("No names provided")
    ...     elif names[0] == "t1":
    ...        return test_table_1
    ...     elif names[1] == "t2":
    ...        return test_table_2
    ...     else:
    ...        raise Exception("Unrecognized table name")
    ...
    >>> substrait_query = '''
    ...         {
    ...             "relations": [
    ...             {"rel": {
    ...                 "read": {
    ...                 "base_schema": {
    ...                     "struct": {
    ...                     "types": [
    ...                                 {"i64": {}}
    ...                             ]
    ...                     },
    ...                     "names": [
    ...                             "x"
    ...                             ]
    ...                 },
    ...                 "namedTable": {
    ...                         "names": ["t1"]
    ...                 }
    ...                 }
    ...             }}
    ...             ]
    ...         }
    ... '''
    >>> buf = pa._substrait._parse_json_plan(tobytes(substrait_query))
    >>> reader = pa.substrait.run_query(buf, table_provider=table_provider)
    >>> reader.read_all()
    pyarrow.Table
    x: int64
    ----
    x: [[1,2,3]]
    """

    cdef:
        CResult[shared_ptr[CRecordBatchReader]] c_res_reader
        shared_ptr[CRecordBatchReader] c_reader
        RecordBatchReader reader
        shared_ptr[CBuffer] c_buf_plan
        CConversionOptions c_conversion_options
        c_bool c_use_threads

    c_use_threads = use_threads
    if isinstance(plan, (bytes, memoryview)):
        c_buf_plan = pyarrow_unwrap_buffer(py_buffer(plan))
    elif isinstance(plan, Buffer):
        c_buf_plan = pyarrow_unwrap_buffer(plan)
    else:
        raise TypeError(
            f"Expected 'pyarrow.Buffer' or bytes, got '{type(plan)}'")

    if table_provider is not None:
        named_table_args = {
            "provider": table_provider
        }
        c_conversion_options.named_table_provider = BindFunction[CNamedTableProvider](
            &_create_named_table_provider, named_table_args)

    with nogil:
        c_res_reader = ExecuteSerializedPlan(
            deref(c_buf_plan), default_extension_id_registry(),
            GetFunctionRegistry(), c_conversion_options, c_use_threads)

    c_reader = GetResultValue(c_res_reader)

    reader = RecordBatchReader.__new__(RecordBatchReader)
    reader.reader = c_reader
    return reader


def _parse_json_plan(plan):
    """
    Parse a JSON plan into equivalent serialized Protobuf.

    Parameters
    ----------
    plan : bytes
        Substrait plan in JSON.

    Returns
    -------
    Buffer
        A buffer containing the serialized Protobuf plan.
    """

    cdef:
        CResult[shared_ptr[CBuffer]] c_res_buffer
        c_string c_str_plan
        shared_ptr[CBuffer] c_buf_plan

    c_str_plan = plan
    c_res_buffer = SerializeJsonPlan(c_str_plan)
    with nogil:
        c_buf_plan = GetResultValue(c_res_buffer)
    return pyarrow_wrap_buffer(c_buf_plan)


class SubstraitSchema:
    """A Schema encoded for Substrait usage.

    The SubstraitSchema contains a schema represented
    both as a substrait ``NamedStruct`` and as an
    ``ExtendedExpression``.

    The ``ExtendedExpression`` is available for cases where types
    used by the schema require extensions to decode them.
    In such case the schema will be the ``base_schema`` of the
    ``ExtendedExpression`` and all extensions will be provided.
    """

    def __init__(self, schema, expression):
        self.schema = schema
        self.expression = expression

    def to_pysubstrait(self):
        """Convert the schema to a substrait-python ExtendedExpression object."""
        if py_substrait is None:
            raise ImportError("The 'substrait' package is required.")
        return py_substrait.proto.ExtendedExpression.FromString(self.expression)


def serialize_schema(schema):
    """
    Serialize a schema into a SubstraitSchema object.

    Parameters
    ----------
    schema : Schema
        The schema to serialize

    Returns
    -------
    SubstraitSchema
        The schema stored in a SubstraitSchema object.
    """
    return SubstraitSchema(
        schema=_serialize_namedstruct_schema(schema),
        expression=serialize_expressions([], [], schema, allow_arrow_extensions=True)
    )


def _serialize_namedstruct_schema(schema):
    cdef:
        CResult[shared_ptr[CBuffer]] c_res_buffer
        shared_ptr[CBuffer] c_buffer
        CConversionOptions c_conversion_options
        CExtensionSet c_extensions

    with nogil:
        c_res_buffer = SerializeSchema(deref((<Schema> schema).sp_schema), &c_extensions, c_conversion_options)
        c_buffer = GetResultValue(c_res_buffer)

    return memoryview(pyarrow_wrap_buffer(c_buffer))


def deserialize_schema(buf):
    """
    Deserialize a ``NamedStruct`` Substrait message
    or a SubstraitSchema object into an Arrow Schema object

    Parameters
    ----------
    buf : Buffer or bytes or SubstraitSchema
        The message to deserialize

    Returns
    -------
    Schema
        The deserialized schema
    """
    cdef:
        shared_ptr[CBuffer] c_buffer
        CResult[shared_ptr[CSchema]] c_res_schema
        shared_ptr[CSchema] c_schema
        CConversionOptions c_conversion_options
        CExtensionSet c_extensions

    if isinstance(buf, SubstraitSchema):
        return deserialize_expressions(buf.expression).schema

    if isinstance(buf, (bytes, memoryview)):
        c_buffer = pyarrow_unwrap_buffer(py_buffer(buf))
    elif isinstance(buf, Buffer):
        c_buffer = pyarrow_unwrap_buffer(buf)
    else:
        raise TypeError(
            f"Expected 'pyarrow.Buffer' or bytes, got '{type(buf)}'")

    with nogil:
        c_res_schema = DeserializeSchema(
            deref(c_buffer), c_extensions, c_conversion_options)
        c_schema = GetResultValue(c_res_schema)

    return pyarrow_wrap_schema(c_schema)


def serialize_expressions(exprs, names, schema, *, allow_arrow_extensions=False):
    """
    Serialize a collection of expressions into Substrait

    Substrait expressions must be bound to a schema.  For example,
    the Substrait expression ``a:i32 + b:i32`` is different from the
    Substrait expression ``a:i64 + b:i64``.  Pyarrow expressions are
    typically unbound.  For example, both of the above expressions
    would be represented as ``a + b`` in pyarrow.

    This means a schema must be provided when serializing an expression.
    It also means that the serialization may fail if a matching function
    call cannot be found for the expression.

    Parameters
    ----------
    exprs : list of Expression
        The expressions to serialize
    names : list of str
        Names for the expressions
    schema : Schema
        The schema the expressions will be bound to
    allow_arrow_extensions : bool, default False
        If False then only functions that are part of the core Substrait function
        definitions will be allowed.  Set this to True to allow pyarrow-specific functions
        and user defined functions but the result may not be accepted by other
        compute libraries.

    Returns
    -------
    Buffer
        An ExtendedExpression message containing the serialized expressions
    """
    cdef:
        CResult[shared_ptr[CBuffer]] c_res_buffer
        shared_ptr[CBuffer] c_buffer
        CNamedExpression c_named_expr
        CBoundExpressions c_bound_exprs
        CConversionOptions c_conversion_options

    if len(exprs) != len(names):
        raise ValueError("exprs and names need to have the same length")
    for expr, name in zip(exprs, names):
        if not isinstance(expr, Expression):
            raise TypeError(f"Expected Expression, got '{type(expr)}' in exprs")
        if not isinstance(name, str):
            raise TypeError(f"Expected str, got '{type(name)}' in names")
        c_named_expr.expression = (<Expression> expr).unwrap()
        c_named_expr.name = tobytes(<str> name)
        c_bound_exprs.named_expressions.push_back(c_named_expr)

    c_bound_exprs.schema = (<Schema> schema).sp_schema

    c_conversion_options.allow_arrow_extensions = allow_arrow_extensions

    with nogil:
        c_res_buffer = SerializeExpressions(c_bound_exprs, c_conversion_options)
        c_buffer = GetResultValue(c_res_buffer)
    return memoryview(pyarrow_wrap_buffer(c_buffer))


cdef class BoundExpressions(_Weakrefable):
    """
    A collection of named expressions and the schema they are bound to

    This is equivalent to the Substrait ExtendedExpression message
    """

    cdef:
        CBoundExpressions c_bound_exprs

    def __init__(self):
        msg = 'BoundExpressions is an abstract class thus cannot be initialized.'
        raise TypeError(msg)

    cdef void init(self, CBoundExpressions bound_expressions):
        self.c_bound_exprs = bound_expressions

    @property
    def schema(self):
        """
        The common schema that all expressions are bound to
        """
        return pyarrow_wrap_schema(self.c_bound_exprs.schema)

    @property
    def expressions(self):
        """
        A dict from expression name to expression
        """
        expr_dict = {}
        for named_expr in self.c_bound_exprs.named_expressions:
            name = frombytes(named_expr.name)
            expr = Expression.wrap(named_expr.expression)
            expr_dict[name] = expr
        return expr_dict

    @staticmethod
    cdef wrap(const CBoundExpressions& bound_expressions):
        cdef BoundExpressions self = BoundExpressions.__new__(BoundExpressions)
        self.init(bound_expressions)
        return self

    @classmethod
    def from_substrait(cls, message):
        """
        Convert a Substrait message into a BoundExpressions object

        Parameters
        ----------
        message : Buffer or bytes or protobuf Message
            The message to convert to a BoundExpressions object

        Returns
        -------
        BoundExpressions
            The converted expressions, their names, and the bound schema
        """
        if isinstance(message, (bytes, memoryview)):
            return deserialize_expressions(message)
        elif isinstance(message, Buffer):
            return deserialize_expressions(message)
        else:
            try:
                return deserialize_expressions(message.SerializeToString())
            except AttributeError:
                raise TypeError(
                    f"Expected 'pyarrow.Buffer' or bytes or protobuf Message, got '{type(message)}'")


def deserialize_expressions(buf):
    """
    Deserialize an ExtendedExpression Substrait message into a BoundExpressions object

    Parameters
    ----------
    buf : Buffer or bytes
        The message to deserialize

    Returns
    -------
    BoundExpressions
        The deserialized expressions, their names, and the bound schema
    """
    cdef:
        shared_ptr[CBuffer] c_buffer
        CResult[CBoundExpressions] c_res_bound_exprs
        CBoundExpressions c_bound_exprs

    if isinstance(buf, (bytes, memoryview)):
        c_buffer = pyarrow_unwrap_buffer(py_buffer(buf))
    elif isinstance(buf, Buffer):
        c_buffer = pyarrow_unwrap_buffer(buf)
    else:
        raise TypeError(
            f"Expected 'pyarrow.Buffer' or bytes, got '{type(buf)}'")

    with nogil:
        c_res_bound_exprs = DeserializeExpressions(deref(c_buffer))
        c_bound_exprs = GetResultValue(c_res_bound_exprs)

    return BoundExpressions.wrap(c_bound_exprs)


def get_supported_functions():
    """
    Get a list of Substrait functions that the underlying
    engine currently supports.

    Returns
    -------
    list[str]
        A list of function ids encoded as '{uri}#{name}'
    """

    cdef:
        ExtensionIdRegistry* c_id_registry
        std_vector[c_string] c_ids

    c_id_registry = default_extension_id_registry()
    c_ids = c_id_registry.GetSupportedSubstraitFunctions()

    functions_list = []
    for c_id in c_ids:
        functions_list.append(frombytes(c_id))
    return functions_list
