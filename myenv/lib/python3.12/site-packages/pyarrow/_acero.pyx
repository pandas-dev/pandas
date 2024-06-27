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

# ---------------------------------------------------------------------
# Low-level Acero bindings

# cython: profile=False
# distutils: language = c++
# cython: language_level = 3

from pyarrow.includes.common cimport *
from pyarrow.includes.libarrow cimport *
from pyarrow.includes.libarrow_acero cimport *
from pyarrow.lib cimport (Table, pyarrow_unwrap_table, pyarrow_wrap_table,
                          RecordBatchReader)
from pyarrow.lib import frombytes, tobytes
from pyarrow._compute cimport (
    Expression, FunctionOptions, _ensure_field_ref, _true,
    unwrap_null_placement, unwrap_sort_order
)


cdef class ExecNodeOptions(_Weakrefable):
    """
    Base class for the node options.

    Use one of the subclasses to construct an options object.
    """
    __slots__ = ()  # avoid mistakingly creating attributes

    cdef void init(self, const shared_ptr[CExecNodeOptions]& sp):
        self.wrapped = sp

    cdef inline shared_ptr[CExecNodeOptions] unwrap(self) nogil:
        return self.wrapped


cdef class _TableSourceNodeOptions(ExecNodeOptions):

    def _set_options(self, Table table):
        cdef:
            shared_ptr[CTable] c_table

        c_table = pyarrow_unwrap_table(table)
        self.wrapped.reset(
            new CTableSourceNodeOptions(c_table)
        )


class TableSourceNodeOptions(_TableSourceNodeOptions):
    """
    A Source node which accepts a table.

    This is the option class for the "table_source" node factory.

    Parameters
    ----------
    table : pyarrow.Table
        The table which acts as the data source.
    """

    def __init__(self, Table table):
        self._set_options(table)


cdef class _FilterNodeOptions(ExecNodeOptions):

    def _set_options(self, Expression filter_expression not None):
        self.wrapped.reset(
            new CFilterNodeOptions(<CExpression>filter_expression.unwrap())
        )


class FilterNodeOptions(_FilterNodeOptions):
    """
    Make a node which excludes some rows from batches passed through it.

    This is the option class for the "filter" node factory.

    The "filter" operation provides an option to define data filtering
    criteria. It selects rows where the given expression evaluates to true.
    Filters can be written using pyarrow.compute.Expression, and the
    expression must have a return type of boolean.

    Parameters
    ----------
    filter_expression : pyarrow.compute.Expression
    """

    def __init__(self, Expression filter_expression):
        self._set_options(filter_expression)


cdef class _ProjectNodeOptions(ExecNodeOptions):

    def _set_options(self, expressions, names=None):
        cdef:
            Expression expr
            vector[CExpression] c_expressions
            vector[c_string] c_names

        for expr in expressions:
            c_expressions.push_back(expr.unwrap())

        if names is not None:
            if len(names) != len(expressions):
                raise ValueError(
                    "The number of names should be equal to the number of expressions"
                )

            for name in names:
                c_names.push_back(<c_string>tobytes(name))

            self.wrapped.reset(
                new CProjectNodeOptions(c_expressions, c_names)
            )
        else:
            self.wrapped.reset(
                new CProjectNodeOptions(c_expressions)
            )


class ProjectNodeOptions(_ProjectNodeOptions):
    """
    Make a node which executes expressions on input batches,
    producing batches of the same length with new columns.

    This is the option class for the "project" node factory.

    The "project" operation rearranges, deletes, transforms, and
    creates columns. Each output column is computed by evaluating
    an expression against the source record batch. These must be
    scalar expressions (expressions consisting of scalar literals,
    field references and scalar functions, i.e. elementwise functions
    that return one value for each input row independent of the value
    of all other rows).

    Parameters
    ----------
    expressions : list of pyarrow.compute.Expression
        List of expressions to evaluate against the source batch. This must
        be scalar expressions.
    names : list of str, optional
        List of names for each of the output columns (same length as
        `expressions`). If `names` is not provided, the string
        representations of exprs will be used.
    """

    def __init__(self, expressions, names=None):
        self._set_options(expressions, names)


cdef class _AggregateNodeOptions(ExecNodeOptions):

    def _set_options(self, aggregates, keys=None):
        cdef:
            CAggregate c_aggr
            vector[CAggregate] c_aggregations
            vector[CFieldRef] c_keys

        for arg_names, func_name, opts, name in aggregates:
            c_aggr.function = tobytes(func_name)
            if opts is not None:
                c_aggr.options = (<FunctionOptions?>opts).wrapped
            else:
                c_aggr.options = <shared_ptr[CFunctionOptions]>nullptr
            if not isinstance(arg_names, (list, tuple)):
                arg_names = [arg_names]
            for arg in arg_names:
                c_aggr.target.push_back(_ensure_field_ref(arg))
            c_aggr.name = tobytes(name)

            c_aggregations.push_back(move(c_aggr))

        if keys is None:
            keys = []
        for name in keys:
            c_keys.push_back(_ensure_field_ref(name))

        self.wrapped.reset(
            new CAggregateNodeOptions(c_aggregations, c_keys)
        )


class AggregateNodeOptions(_AggregateNodeOptions):
    """
    Make a node which aggregates input batches, optionally grouped by keys.

    This is the option class for the "aggregate" node factory.

    Acero supports two types of aggregates: "scalar" aggregates,
    and "hash" aggregates. Scalar aggregates reduce an array or scalar
    input to a single scalar output (e.g. computing the mean of a column).
    Hash aggregates act like GROUP BY in SQL and first partition data
    based on one or more key columns, then reduce the data in each partition.
    The aggregate node supports both types of computation, and can compute
    any number of aggregations at once.

    Parameters
    ----------
    aggregates : list of tuples
        Aggregations which will be applied to the targeted fields.
        Specified as a list of tuples, where each tuple is one aggregation
        specification and consists of: aggregation target column(s) followed
        by function name, aggregation function options object and the
        output field name.
        The target column(s) specification can be a single field reference,
        an empty list or a list of fields unary, nullary and n-ary aggregation
        functions respectively. Each field reference can be a string
        column name or expression.
    keys : list of field references, optional
        Keys by which aggregations will be grouped. Each key can reference
        a field using a string name or expression.
    """

    def __init__(self, aggregates, keys=None):
        self._set_options(aggregates, keys)


cdef class _OrderByNodeOptions(ExecNodeOptions):

    def _set_options(self, sort_keys, null_placement):
        cdef:
            vector[CSortKey] c_sort_keys

        for name, order in sort_keys:
            c_sort_keys.push_back(
                CSortKey(_ensure_field_ref(name), unwrap_sort_order(order))
            )

        self.wrapped.reset(
            new COrderByNodeOptions(
                COrdering(c_sort_keys, unwrap_null_placement(null_placement))
            )
        )


class OrderByNodeOptions(_OrderByNodeOptions):
    """
    Make a node which applies a new ordering to the data.

    Currently this node works by accumulating all data, sorting, and then
    emitting the new data with an updated batch index.
    Larger-than-memory sort is not currently supported.

    This is the option class for the "order_by" node factory.

    Parameters
    ----------
    sort_keys : sequence of (name, order) tuples
        Names of field/column keys to sort the input on,
        along with the order each field/column is sorted in.
        Accepted values for `order` are "ascending", "descending".
        Each field reference can be a string column name or expression.
    null_placement : str, default "at_end"
        Where nulls in input should be sorted, only applying to
        columns/fields mentioned in `sort_keys`.
        Accepted values are "at_start", "at_end".
    """

    def __init__(self, sort_keys=(), *, null_placement="at_end"):
        self._set_options(sort_keys, null_placement)


cdef class _HashJoinNodeOptions(ExecNodeOptions):

    def _set_options(
        self, join_type, left_keys, right_keys, left_output=None, right_output=None,
        output_suffix_for_left="", output_suffix_for_right="",
    ):
        cdef:
            CJoinType c_join_type
            vector[CFieldRef] c_left_keys
            vector[CFieldRef] c_right_keys
            vector[CFieldRef] c_left_output
            vector[CFieldRef] c_right_output

        # join type
        if join_type == "left semi":
            c_join_type = CJoinType_LEFT_SEMI
        elif join_type == "right semi":
            c_join_type = CJoinType_RIGHT_SEMI
        elif join_type == "left anti":
            c_join_type = CJoinType_LEFT_ANTI
        elif join_type == "right anti":
            c_join_type = CJoinType_RIGHT_ANTI
        elif join_type == "inner":
            c_join_type = CJoinType_INNER
        elif join_type == "left outer":
            c_join_type = CJoinType_LEFT_OUTER
        elif join_type == "right outer":
            c_join_type = CJoinType_RIGHT_OUTER
        elif join_type == "full outer":
            c_join_type = CJoinType_FULL_OUTER
        else:
            raise ValueError("Unsupported join type")

        # left/right keys
        if not isinstance(left_keys, (list, tuple)):
            left_keys = [left_keys]
        for key in left_keys:
            c_left_keys.push_back(_ensure_field_ref(key))
        if not isinstance(right_keys, (list, tuple)):
            right_keys = [right_keys]
        for key in right_keys:
            c_right_keys.push_back(_ensure_field_ref(key))

        # left/right output fields
        if left_output is not None and right_output is not None:
            for colname in left_output:
                c_left_output.push_back(_ensure_field_ref(colname))
            for colname in right_output:
                c_right_output.push_back(_ensure_field_ref(colname))

            self.wrapped.reset(
                new CHashJoinNodeOptions(
                    c_join_type, c_left_keys, c_right_keys,
                    c_left_output, c_right_output,
                    _true,
                    <c_string>tobytes(output_suffix_for_left),
                    <c_string>tobytes(output_suffix_for_right)
                )
            )
        else:
            self.wrapped.reset(
                new CHashJoinNodeOptions(
                    c_join_type, c_left_keys, c_right_keys,
                    _true,
                    <c_string>tobytes(output_suffix_for_left),
                    <c_string>tobytes(output_suffix_for_right)
                )
            )


class HashJoinNodeOptions(_HashJoinNodeOptions):
    """
    Make a node which implements join operation using hash join strategy.

    This is the option class for the "hashjoin" node factory.

    Parameters
    ----------
    join_type : str
        Type of join. One of "left semi", "right semi", "left anti",
        "right anti", "inner", "left outer", "right outer", "full outer".
    left_keys : str, Expression or list
        Key fields from left input. Each key can be a string column name
        or a field expression, or a list of such field references.
    right_keys : str, Expression or list
        Key fields from right input. See `left_keys` for details.
    left_output : list, optional
        List of output fields passed from left input. If left and right
        output fields are not specified, all valid fields from both left and
        right input will be output. Each field can be a string column name
        or a field expression.
    right_output : list, optional
        List of output fields passed from right input. If left and right
        output fields are not specified, all valid fields from both left and
        right input will be output. Each field can be a string column name
        or a field expression.
    output_suffix_for_left : str
        Suffix added to names of output fields coming from left input
        (used to distinguish, if necessary, between fields of the same
        name in left and right input and can be left empty if there are
        no name collisions).
    output_suffix_for_right : str
        Suffix added to names of output fields coming from right input,
        see `output_suffix_for_left` for details.
    """

    def __init__(
        self, join_type, left_keys, right_keys, left_output=None, right_output=None,
        output_suffix_for_left="", output_suffix_for_right=""
    ):
        self._set_options(
            join_type, left_keys, right_keys, left_output, right_output,
            output_suffix_for_left, output_suffix_for_right
        )


cdef class _AsofJoinNodeOptions(ExecNodeOptions):

    def _set_options(self, left_on, left_by, right_on, right_by, tolerance):
        cdef:
            vector[CFieldRef] c_left_by
            vector[CFieldRef] c_right_by
            CAsofJoinKeys c_left_keys
            CAsofJoinKeys c_right_keys
            vector[CAsofJoinKeys] c_input_keys

        # Prepare left AsofJoinNodeOption::Keys
        if not isinstance(left_by, (list, tuple)):
            left_by = [left_by]
        for key in left_by:
            c_left_by.push_back(_ensure_field_ref(key))

        c_left_keys.on_key = _ensure_field_ref(left_on)
        c_left_keys.by_key = c_left_by

        c_input_keys.push_back(c_left_keys)

        # Prepare right AsofJoinNodeOption::Keys
        if not isinstance(right_by, (list, tuple)):
            right_by = [right_by]
        for key in right_by:
            c_right_by.push_back(_ensure_field_ref(key))

        c_right_keys.on_key = _ensure_field_ref(right_on)
        c_right_keys.by_key = c_right_by

        c_input_keys.push_back(c_right_keys)

        self.wrapped.reset(
            new CAsofJoinNodeOptions(
                c_input_keys,
                tolerance,
            )
        )


class AsofJoinNodeOptions(_AsofJoinNodeOptions):
    """
    Make a node which implements 'as of join' operation.

    This is the option class for the "asofjoin" node factory.

    Parameters
    ----------
    left_on : str, Expression
        The left key on which the join operation should be performed.
        Can be a string column name or a field expression.

        An inexact match is used on the "on" key, i.e. a row is considered a
        match if and only if left_on - tolerance <= right_on <= left_on.

        The input dataset must be sorted by the "on" key. Must be a single
        field of a common type.

        Currently, the "on" key must be an integer, date, or timestamp type.
    left_by: str, Expression or list
        The left keys on which the join operation should be performed.
        Exact equality is used for each field of the "by" keys.
        Each key can be a string column name or a field expression,
        or a list of such field references.
    right_on : str, Expression
        The right key on which the join operation should be performed.
        See `left_on` for details.
    right_by: str, Expression or list
        The right keys on which the join operation should be performed.
        See `left_by` for details.
    tolerance : int
        The tolerance to use for the asof join. The tolerance is interpreted in
        the same units as the "on" key.
    """

    def __init__(self, left_on, left_by, right_on, right_by, tolerance):
        self._set_options(left_on, left_by, right_on, right_by, tolerance)


cdef class Declaration(_Weakrefable):
    """
    Helper class for declaring the nodes of an ExecPlan.

    A Declaration represents an unconstructed ExecNode, and potentially
    more since its inputs may also be Declarations or when constructed
    with ``from_sequence``.

    The possible ExecNodes to use are registered with a name,
    the "factory name", and need to be specified using this name, together
    with its corresponding ExecNodeOptions subclass.

    Parameters
    ----------
    factory_name : str
        The ExecNode factory name, such as "table_source", "filter",
        "project" etc. See the ExecNodeOptions subclasses for the exact
        factory names to use.
    options : ExecNodeOptions
        Corresponding ExecNodeOptions subclass (matching the factory name).
    inputs : list of Declaration, optional
        Input nodes for this declaration. Optional if the node is a source
        node, or when the declaration gets combined later with
        ``from_sequence``.

    Returns
    -------
    Declaration
    """
    cdef void init(self, const CDeclaration& c_decl):
        self.decl = c_decl

    @staticmethod
    cdef wrap(const CDeclaration& c_decl):
        cdef Declaration self = Declaration.__new__(Declaration)
        self.init(c_decl)
        return self

    cdef inline CDeclaration unwrap(self) nogil:
        return self.decl

    def __init__(self, factory_name, ExecNodeOptions options, inputs=None):
        cdef:
            c_string c_factory_name
            CDeclaration c_decl
            vector[CDeclaration.Input] c_inputs

        c_factory_name = tobytes(factory_name)

        if inputs is not None:
            for ipt in inputs:
                c_inputs.push_back(
                    CDeclaration.Input((<Declaration>ipt).unwrap())
                )

        c_decl = CDeclaration(c_factory_name, c_inputs, options.unwrap())
        self.init(c_decl)

    @staticmethod
    def from_sequence(decls):
        """
        Convenience factory for the common case of a simple sequence of nodes.

        Each of the declarations will be appended to the inputs of the
        subsequent declaration, and the final modified declaration will
        be returned.

        Parameters
        ----------
        decls : list of Declaration

        Returns
        -------
        Declaration
        """
        cdef:
            vector[CDeclaration] c_decls
            CDeclaration c_decl

        for decl in decls:
            c_decls.push_back((<Declaration> decl).unwrap())

        c_decl = CDeclaration.Sequence(c_decls)
        return Declaration.wrap(c_decl)

    def __str__(self):
        return frombytes(GetResultValue(DeclarationToString(self.decl)))

    def __repr__(self):
        return "<pyarrow.acero.Declaration>\n{0}".format(str(self))

    def to_table(self, bint use_threads=True):
        """
        Run the declaration and collect the results into a table.

        This method will implicitly add a sink node to the declaration
        to collect results into a table. It will then create an ExecPlan
        from the declaration, start the exec plan, block until the plan
        has finished, and return the created table.

        Parameters
        ----------
        use_threads : bool, default True
            If set to False, then all CPU work will be done on the calling
            thread. I/O tasks will still happen on the I/O executor
            and may be multi-threaded (but should not use significant CPU
            resources).

        Returns
        -------
        pyarrow.Table
        """
        cdef:
            shared_ptr[CTable] c_table

        with nogil:
            c_table = GetResultValue(DeclarationToTable(self.unwrap(), use_threads))
        return pyarrow_wrap_table(c_table)

    def to_reader(self, bint use_threads=True):
        """Run the declaration and return results as a RecordBatchReader.

        For details about the parameters, see `to_table`.

        Returns
        -------
        pyarrow.RecordBatchReader
        """
        cdef:
            RecordBatchReader reader
        reader = RecordBatchReader.__new__(RecordBatchReader)
        reader.reader.reset(
            GetResultValue(DeclarationToReader(self.unwrap(), use_threads)).release()
        )
        return reader
