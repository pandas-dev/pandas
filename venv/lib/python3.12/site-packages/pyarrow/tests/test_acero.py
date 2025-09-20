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
import pyarrow.compute as pc
from pyarrow.compute import field

try:
    from pyarrow.acero import (
        Declaration,
        TableSourceNodeOptions,
        FilterNodeOptions,
        ProjectNodeOptions,
        AggregateNodeOptions,
        OrderByNodeOptions,
        HashJoinNodeOptions,
        AsofJoinNodeOptions,
    )
except ImportError:
    pass

try:
    import pyarrow.dataset as ds
    from pyarrow.acero import ScanNodeOptions
except ImportError:
    ds = None

pytestmark = pytest.mark.acero


@pytest.fixture
def table_source():
    table = pa.table({'a': [1, 2, 3], 'b': [4, 5, 6]})
    table_opts = TableSourceNodeOptions(table)
    table_source = Declaration("table_source", options=table_opts)
    return table_source


def test_declaration():

    table = pa.table({'a': [1, 2, 3], 'b': [4, 5, 6]})
    table_opts = TableSourceNodeOptions(table)
    filter_opts = FilterNodeOptions(field('a') > 1)

    # using sequence
    decl = Declaration.from_sequence([
        Declaration("table_source", options=table_opts),
        Declaration("filter", options=filter_opts)
    ])
    result = decl.to_table()
    assert result.equals(table.slice(1, 2))

    # using explicit inputs
    table_source = Declaration("table_source", options=table_opts)
    filtered = Declaration("filter", options=filter_opts, inputs=[table_source])
    result = filtered.to_table()
    assert result.equals(table.slice(1, 2))


def test_declaration_repr(table_source):

    assert "TableSourceNode" in str(table_source)
    assert "TableSourceNode" in repr(table_source)


def test_declaration_to_reader(table_source):
    with table_source.to_reader() as reader:
        assert reader.schema == pa.schema([("a", pa.int64()), ("b", pa.int64())])
        result = reader.read_all()
    expected = pa.table({'a': [1, 2, 3], 'b': [4, 5, 6]})
    assert result.equals(expected)


def test_table_source():
    with pytest.raises(TypeError):
        TableSourceNodeOptions(pa.record_batch([pa.array([1, 2, 3])], ["a"]))

    table_source = TableSourceNodeOptions(None)
    decl = Declaration("table_source", table_source)
    with pytest.raises(
        ValueError, match="TableSourceNode requires table which is not null"
    ):
        _ = decl.to_table()


def test_filter(table_source):
    # referencing unknown field
    decl = Declaration.from_sequence([
        table_source,
        Declaration("filter", options=FilterNodeOptions(field("c") > 1))
    ])
    with pytest.raises(ValueError, match=r"No match for FieldRef.Name\(c\)"):
        _ = decl.to_table()

    # requires a pyarrow Expression
    with pytest.raises(TypeError):
        FilterNodeOptions(pa.array([True, False, True]))
    with pytest.raises(TypeError):
        FilterNodeOptions(None)


@pytest.mark.parametrize('source', [
    pa.record_batch({"number": [1, 2, 3]}),
    pa.table({"number": [1, 2, 3]})
])
def test_filter_all_rows(source):
    # GH-46057: filtering all rows should return empty RecordBatch with same schema
    result_expr = source.filter(pc.field("number") < 0)

    assert result_expr.num_rows == 0
    assert type(result_expr) is type(source)
    assert result_expr.schema.equals(source.schema)

    result_mask = source.filter(pa.array([False, False, False]))

    assert result_mask.num_rows == 0
    assert type(result_mask) is type(source)
    assert result_mask.schema.equals(source.schema)


def test_project(table_source):
    # default name from expression
    decl = Declaration.from_sequence([
        table_source,
        Declaration("project", ProjectNodeOptions([pc.multiply(field("a"), 2)]))
    ])
    result = decl.to_table()
    assert result.schema.names == ["multiply(a, 2)"]
    assert result[0].to_pylist() == [2, 4, 6]

    # provide name
    decl = Declaration.from_sequence([
        table_source,
        Declaration("project", ProjectNodeOptions([pc.multiply(field("a"), 2)], ["a2"]))
    ])
    result = decl.to_table()
    assert result.schema.names == ["a2"]
    assert result["a2"].to_pylist() == [2, 4, 6]

    # input validation
    with pytest.raises(ValueError):
        ProjectNodeOptions([pc.multiply(field("a"), 2)], ["a2", "b2"])

    # no scalar expression
    decl = Declaration.from_sequence([
        table_source,
        Declaration("project", ProjectNodeOptions([pc.sum(field("a"))]))
    ])
    with pytest.raises(ValueError, match="cannot Execute non-scalar expression"):
        _ = decl.to_table()


def test_aggregate_scalar(table_source):
    decl = Declaration.from_sequence([
        table_source,
        Declaration("aggregate", AggregateNodeOptions([("a", "sum", None, "a_sum")]))
    ])
    result = decl.to_table()
    assert result.schema.names == ["a_sum"]
    assert result["a_sum"].to_pylist() == [6]

    # with options class
    table = pa.table({'a': [1, 2, None]})
    aggr_opts = AggregateNodeOptions(
        [("a", "sum", pc.ScalarAggregateOptions(skip_nulls=False), "a_sum")]
    )
    decl = Declaration.from_sequence([
        Declaration("table_source", TableSourceNodeOptions(table)),
        Declaration("aggregate", aggr_opts),
    ])
    result = decl.to_table()
    assert result.schema.names == ["a_sum"]
    assert result["a_sum"].to_pylist() == [None]

    # test various ways of specifying the target column
    for target in ["a", field("a"), 0, field(0), ["a"], [field("a")], [0]]:
        aggr_opts = AggregateNodeOptions([(target, "sum", None, "a_sum")])
        decl = Declaration.from_sequence(
            [table_source, Declaration("aggregate", aggr_opts)]
        )
        result = decl.to_table()
        assert result.schema.names == ["a_sum"]
        assert result["a_sum"].to_pylist() == [6]

    # proper error when specifying the wrong number of target columns
    aggr_opts = AggregateNodeOptions([(["a", "b"], "sum", None, "a_sum")])
    decl = Declaration.from_sequence(
        [table_source, Declaration("aggregate", aggr_opts)]
    )
    with pytest.raises(
        ValueError, match="Function 'sum' accepts 1 arguments but 2 passed"
    ):
        _ = decl.to_table()

    # proper error when using hash aggregation without keys
    aggr_opts = AggregateNodeOptions([("a", "hash_sum", None, "a_sum")])
    decl = Declaration.from_sequence(
        [table_source, Declaration("aggregate", aggr_opts)]
    )
    with pytest.raises(ValueError, match="is a hash aggregate function"):
        _ = decl.to_table()


def test_aggregate_hash():
    table = pa.table({'a': [1, 2, None], 'b': ["foo", "bar", "foo"]})
    table_opts = TableSourceNodeOptions(table)
    table_source = Declaration("table_source", options=table_opts)

    # default options
    aggr_opts = AggregateNodeOptions(
        [("a", "hash_count", None, "count(a)")], keys=["b"])
    decl = Declaration.from_sequence([
        table_source, Declaration("aggregate", aggr_opts)
    ])
    result = decl.to_table()
    expected = pa.table({"b": ["foo", "bar"], "count(a)": [1, 1]})
    assert result.equals(expected)

    # specify function options
    aggr_opts = AggregateNodeOptions(
        [("a", "hash_count", pc.CountOptions("all"), "count(a)")], keys=["b"]
    )
    decl = Declaration.from_sequence([
        table_source, Declaration("aggregate", aggr_opts)
    ])
    result = decl.to_table()
    expected_all = pa.table({"b": ["foo", "bar"], "count(a)": [2, 1]})
    assert result.equals(expected_all)

    # specify keys as field references
    aggr_opts = AggregateNodeOptions(
        [("a", "hash_count", None, "count(a)")], keys=[field("b")]
    )
    decl = Declaration.from_sequence([
        table_source, Declaration("aggregate", aggr_opts)
    ])
    result = decl.to_table()
    assert result.equals(expected)

    # wrong type of (aggregation) function
    # TODO test with kernel that matches number of arguments (arity) -> avoid segfault
    aggr_opts = AggregateNodeOptions([("a", "sum", None, "a_sum")], keys=["b"])
    decl = Declaration.from_sequence([
        table_source, Declaration("aggregate", aggr_opts)
    ])
    with pytest.raises(ValueError):
        _ = decl.to_table()


def test_order_by():
    table = pa.table({'a': [1, 2, 3, 4], 'b': [1, 3, None, 2]})
    table_source = Declaration("table_source", TableSourceNodeOptions(table))

    ord_opts = OrderByNodeOptions([("b", "ascending")])
    decl = Declaration.from_sequence([table_source, Declaration("order_by", ord_opts)])
    result = decl.to_table()
    expected = pa.table({"a": [1, 4, 2, 3], "b": [1, 2, 3, None]})
    assert result.equals(expected)

    ord_opts = OrderByNodeOptions([(field("b"), "descending")])
    decl = Declaration.from_sequence([table_source, Declaration("order_by", ord_opts)])
    result = decl.to_table()
    expected = pa.table({"a": [2, 4, 1, 3], "b": [3, 2, 1, None]})
    assert result.equals(expected)

    ord_opts = OrderByNodeOptions([(1, "descending")], null_placement="at_start")
    decl = Declaration.from_sequence([table_source, Declaration("order_by", ord_opts)])
    result = decl.to_table()
    expected = pa.table({"a": [3, 2, 4, 1], "b": [None, 3, 2, 1]})
    assert result.equals(expected)

    # empty ordering
    ord_opts = OrderByNodeOptions([])
    decl = Declaration.from_sequence([table_source, Declaration("order_by", ord_opts)])
    with pytest.raises(
        ValueError, match="`ordering` must be an explicit non-empty ordering"
    ):
        _ = decl.to_table()

    with pytest.raises(ValueError, match="\"decreasing\" is not a valid sort order"):
        _ = OrderByNodeOptions([("b", "decreasing")])

    with pytest.raises(ValueError, match="\"start\" is not a valid null placement"):
        _ = OrderByNodeOptions([("b", "ascending")], null_placement="start")


def test_hash_join():
    left = pa.table({'key': [1, 2, 3], 'a': [4, 5, 6]})
    left_source = Declaration("table_source", options=TableSourceNodeOptions(left))
    right = pa.table({'key': [2, 3, 4], 'b': [4, 5, 6]})
    right_source = Declaration("table_source", options=TableSourceNodeOptions(right))

    # inner join
    join_opts = HashJoinNodeOptions("inner", left_keys="key", right_keys="key")
    joined = Declaration(
        "hashjoin", options=join_opts, inputs=[left_source, right_source])
    result = joined.to_table()
    expected = pa.table(
        [[2, 3], [5, 6], [2, 3], [4, 5]],
        names=["key", "a", "key", "b"])
    assert result.equals(expected)

    for keys in [field("key"), ["key"], [field("key")]]:
        join_opts = HashJoinNodeOptions("inner", left_keys=keys, right_keys=keys)
        joined = Declaration(
            "hashjoin", options=join_opts, inputs=[left_source, right_source])
        result = joined.to_table()
        assert result.equals(expected)

    # left join
    join_opts = HashJoinNodeOptions(
        "left outer", left_keys="key", right_keys="key")
    joined = Declaration(
        "hashjoin", options=join_opts, inputs=[left_source, right_source])
    result = joined.to_table()
    expected = pa.table(
        [[1, 2, 3], [4, 5, 6], [None, 2, 3], [None, 4, 5]],
        names=["key", "a", "key", "b"]
    )
    assert result.sort_by("a").equals(expected)

    # suffixes
    join_opts = HashJoinNodeOptions(
        "left outer", left_keys="key", right_keys="key",
        output_suffix_for_left="_left", output_suffix_for_right="_right")
    joined = Declaration(
        "hashjoin", options=join_opts, inputs=[left_source, right_source])
    result = joined.to_table()
    expected = pa.table(
        [[1, 2, 3], [4, 5, 6], [None, 2, 3], [None, 4, 5]],
        names=["key_left", "a", "key_right", "b"]
    )
    assert result.sort_by("a").equals(expected)

    # manually specifying output columns
    join_opts = HashJoinNodeOptions(
        "left outer", left_keys="key", right_keys="key",
        left_output=["key", "a"], right_output=[field("b")])
    joined = Declaration(
        "hashjoin", options=join_opts, inputs=[left_source, right_source])
    result = joined.to_table()
    expected = pa.table(
        [[1, 2, 3], [4, 5, 6], [None, 4, 5]],
        names=["key", "a", "b"]
    )
    assert result.sort_by("a").equals(expected)


def test_hash_join_with_residual_filter():
    left = pa.table({'key': [1, 2, 3], 'a': [4, 5, 6]})
    left_source = Declaration("table_source", options=TableSourceNodeOptions(left))
    right = pa.table({'key': [2, 3, 4], 'b': [4, 5, 6]})
    right_source = Declaration("table_source", options=TableSourceNodeOptions(right))

    join_opts = HashJoinNodeOptions(
        "inner", left_keys="key", right_keys="key",
        filter_expression=pc.equal(pc.field('a'), 5))
    joined = Declaration(
        "hashjoin", options=join_opts, inputs=[left_source, right_source])
    result = joined.to_table()
    expected = pa.table(
        [[2], [5], [2], [4]],
        names=["key", "a", "key", "b"])
    assert result.equals(expected)

    # test filter expression referencing columns from both side
    join_opts = HashJoinNodeOptions(
        "left outer", left_keys="key", right_keys="key",
        filter_expression=pc.equal(pc.field("a"), 5) | pc.equal(pc.field("b"), 10)
    )
    joined = Declaration(
        "hashjoin", options=join_opts, inputs=[left_source, right_source])
    result = joined.to_table()
    expected = pa.table(
        [[2, 1, 3], [5, 4, 6], [2, None, None], [4, None, None]],
        names=["key", "a", "key", "b"])
    assert result.equals(expected)

    # test with always true
    always_true = pc.scalar(True)
    join_opts = HashJoinNodeOptions(
        "inner", left_keys="key", right_keys="key",
        filter_expression=always_true)
    joined = Declaration(
        "hashjoin", options=join_opts, inputs=[left_source, right_source])
    result = joined.to_table()
    expected = pa.table(
        [[2, 3], [5, 6], [2, 3], [4, 5]],
        names=["key", "a", "key", "b"]
    )
    assert result.equals(expected)

    # test with always false
    always_false = pc.scalar(False)
    join_opts = HashJoinNodeOptions(
        "inner", left_keys="key", right_keys="key",
        filter_expression=always_false)
    joined = Declaration(
        "hashjoin", options=join_opts, inputs=[left_source, right_source])
    result = joined.to_table()
    expected = pa.table(
        [
            pa.array([], type=pa.int64()),
            pa.array([], type=pa.int64()),
            pa.array([], type=pa.int64()),
            pa.array([], type=pa.int64())
        ],
        names=["key", "a", "key", "b"]
    )
    assert result.equals(expected)


def test_asof_join():
    left = pa.table({'key': [1, 2, 3], 'ts': [1, 1, 1], 'a': [4, 5, 6]})
    left_source = Declaration("table_source", options=TableSourceNodeOptions(left))
    right = pa.table({'key': [2, 3, 4], 'ts': [2, 5, 2], 'b': [4, 5, 6]})
    right_source = Declaration("table_source", options=TableSourceNodeOptions(right))

    # asof join
    join_opts = AsofJoinNodeOptions(
        left_on="ts", left_by=["key"],
        right_on="ts", right_by=["key"],
        tolerance=1,
    )
    joined = Declaration(
        "asofjoin", options=join_opts, inputs=[left_source, right_source]
    )
    result = joined.to_table()
    expected = pa.table(
        [[1, 2, 3], [1, 1, 1], [4, 5, 6], [None, 4, None]],
        names=["key", "ts", "a", "b"])
    assert result == expected

    for by in [field("key"), ["key"], [field("key")]]:
        for on in [field("ts"), "ts"]:
            join_opts = AsofJoinNodeOptions(
                left_on=on, left_by=by,
                right_on=on, right_by=by,
                tolerance=1,
            )
            joined = Declaration(
                "asofjoin", options=join_opts, inputs=[left_source, right_source])
            result = joined.to_table()
            assert result == expected


@pytest.mark.dataset
def test_scan(tempdir):
    table = pa.table({'a': [1, 2, 3], 'b': [4, 5, 6]})
    ds.write_dataset(table, tempdir / "dataset", format="parquet")
    dataset = ds.dataset(tempdir / "dataset", format="parquet")
    decl = Declaration("scan", ScanNodeOptions(dataset))
    result = decl.to_table()
    assert result.schema.names == [
        "a", "b", "__fragment_index", "__batch_index",
        "__last_in_fragment", "__filename"
    ]
    assert result.select(["a", "b"]).equals(table)

    # using a filter only does pushdown (depending on file format), not actual filter

    scan_opts = ScanNodeOptions(dataset, filter=field('a') > 1)
    decl = Declaration("scan", scan_opts)
    # fragment not filtered based on min/max statistics
    assert decl.to_table().num_rows == 3

    scan_opts = ScanNodeOptions(dataset, filter=field('a') > 4)
    decl = Declaration("scan", scan_opts)
    # full fragment filtered based on min/max statistics
    assert decl.to_table().num_rows == 0

    # projection scan option

    scan_opts = ScanNodeOptions(dataset, columns={"a2": pc.multiply(field("a"), 2)})
    decl = Declaration("scan", scan_opts)
    result = decl.to_table()
    # "a" is included in the result (needed later on for the actual projection)
    assert result["a"].to_pylist() == [1, 2, 3]
    # "b" is still included, but without data as it will be removed by the projection
    assert pc.all(result["b"].is_null()).as_py()
