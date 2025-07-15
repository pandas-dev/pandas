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
# Implement Internal ExecPlan bindings

# cython: profile=False
# distutils: language = c++
# cython: language_level = 3

from pyarrow.lib import Table, RecordBatch
from pyarrow.compute import Expression, field

try:
    from pyarrow._acero import (  # noqa
        Declaration,
        ExecNodeOptions,
        TableSourceNodeOptions,
        FilterNodeOptions,
        ProjectNodeOptions,
        AggregateNodeOptions,
        OrderByNodeOptions,
        HashJoinNodeOptions,
        AsofJoinNodeOptions,
    )
except ImportError as exc:
    raise ImportError(
        f"The pyarrow installation is not built with support for 'acero' ({str(exc)})"
    ) from None


try:
    import pyarrow.dataset as ds
    from pyarrow._dataset import ScanNodeOptions
except ImportError:
    class DatasetModuleStub:
        class Dataset:
            pass

        class InMemoryDataset:
            pass
    ds = DatasetModuleStub


def _dataset_to_decl(dataset, use_threads=True, implicit_ordering=False):
    decl = Declaration("scan", ScanNodeOptions(
        dataset, use_threads=use_threads,
        implicit_ordering=implicit_ordering))

    # Get rid of special dataset columns
    # "__fragment_index", "__batch_index", "__last_in_fragment", "__filename"
    projections = [field(f) for f in dataset.schema.names]
    decl = Declaration.from_sequence(
        [decl, Declaration("project", ProjectNodeOptions(projections))]
    )

    filter_expr = dataset._scan_options.get("filter")
    if filter_expr is not None:
        # Filters applied in CScanNodeOptions are "best effort" for the scan node itself
        # so we always need to inject an additional Filter node to apply them for real.
        decl = Declaration.from_sequence(
            [decl, Declaration("filter", FilterNodeOptions(filter_expr))]
        )

    return decl


def _perform_join(join_type, left_operand, left_keys,
                  right_operand, right_keys,
                  left_suffix=None, right_suffix=None,
                  use_threads=True, coalesce_keys=False,
                  output_type=Table):
    """
    Perform join of two tables or datasets.

    The result will be an output table with the result of the join operation

    Parameters
    ----------
    join_type : str
        One of supported join types.
    left_operand : Table or Dataset
        The left operand for the join operation.
    left_keys : str or list[str]
        The left key (or keys) on which the join operation should be performed.
    right_operand : Table or Dataset
        The right operand for the join operation.
    right_keys : str or list[str]
        The right key (or keys) on which the join operation should be performed.
    left_suffix : str, default None
        Which suffix to add to left column names. This prevents confusion
        when the columns in left and right operands have colliding names.
    right_suffix : str, default None
        Which suffix to add to the right column names. This prevents confusion
        when the columns in left and right operands have colliding names.
    use_threads : bool, default True
        Whether to use multithreading or not.
    coalesce_keys : bool, default False
        If the duplicated keys should be omitted from one of the sides
        in the join result.
    output_type: Table or InMemoryDataset
        The output type for the exec plan result.

    Returns
    -------
    result_table : Table or InMemoryDataset
    """
    if not isinstance(left_operand, (Table, ds.Dataset)):
        raise TypeError(f"Expected Table or Dataset, got {type(left_operand)}")
    if not isinstance(right_operand, (Table, ds.Dataset)):
        raise TypeError(f"Expected Table or Dataset, got {type(right_operand)}")

    # Prepare left and right tables Keys to send them to the C++ function
    left_keys_order = {}
    if not isinstance(left_keys, (tuple, list)):
        left_keys = [left_keys]
    for idx, key in enumerate(left_keys):
        left_keys_order[key] = idx

    right_keys_order = {}
    if not isinstance(right_keys, (list, tuple)):
        right_keys = [right_keys]
    for idx, key in enumerate(right_keys):
        right_keys_order[key] = idx

    # By default expose all columns on both left and right table
    left_columns = left_operand.schema.names
    right_columns = right_operand.schema.names

    # Pick the join type
    if join_type == "left semi" or join_type == "left anti":
        right_columns = []
    elif join_type == "right semi" or join_type == "right anti":
        left_columns = []
    elif join_type == "inner" or join_type == "left outer":
        right_columns = [
            col for col in right_columns if col not in right_keys_order
        ]
    elif join_type == "right outer":
        left_columns = [
            col for col in left_columns if col not in left_keys_order
        ]

    # Turn the columns to vectors of FieldRefs
    # and set aside indices of keys.
    left_column_keys_indices = {}
    for idx, colname in enumerate(left_columns):
        if colname in left_keys:
            left_column_keys_indices[colname] = idx
    right_column_keys_indices = {}
    for idx, colname in enumerate(right_columns):
        if colname in right_keys:
            right_column_keys_indices[colname] = idx

    # Add the join node to the execplan
    if isinstance(left_operand, ds.Dataset):
        left_source = _dataset_to_decl(left_operand, use_threads=use_threads)
    else:
        left_source = Declaration("table_source", TableSourceNodeOptions(left_operand))
    if isinstance(right_operand, ds.Dataset):
        right_source = _dataset_to_decl(right_operand, use_threads=use_threads)
    else:
        right_source = Declaration(
            "table_source", TableSourceNodeOptions(right_operand)
        )

    if coalesce_keys:
        join_opts = HashJoinNodeOptions(
            join_type, left_keys, right_keys, left_columns, right_columns,
            output_suffix_for_left=left_suffix or "",
            output_suffix_for_right=right_suffix or "",
        )
    else:
        join_opts = HashJoinNodeOptions(
            join_type, left_keys, right_keys,
            output_suffix_for_left=left_suffix or "",
            output_suffix_for_right=right_suffix or "",
        )
    decl = Declaration(
        "hashjoin", options=join_opts, inputs=[left_source, right_source]
    )

    if coalesce_keys and join_type == "full outer":
        # In case of full outer joins, the join operation will output all columns
        # so that we can coalesce the keys and exclude duplicates in a subsequent
        # projection.
        left_columns_set = set(left_columns)
        right_columns_set = set(right_columns)
        # Where the right table columns start.
        right_operand_index = len(left_columns)
        projected_col_names = []
        projections = []
        for idx, col in enumerate(left_columns + right_columns):
            if idx < len(left_columns) and col in left_column_keys_indices:
                # Include keys only once and coalesce left+right table keys.
                projected_col_names.append(col)
                # Get the index of the right key that is being paired
                # with this left key. We do so by retrieving the name
                # of the right key that is in the same position in the provided keys
                # and then looking up the index for that name in the right table.
                right_key_index = right_column_keys_indices[
                    right_keys[left_keys_order[col]]]
                projections.append(
                    Expression._call("coalesce", [
                        Expression._field(idx), Expression._field(
                            right_operand_index+right_key_index)
                    ])
                )
            elif idx >= right_operand_index and col in right_column_keys_indices:
                # Do not include right table keys. As they would lead to duplicated keys
                continue
            else:
                # For all the other columns include them as they are.
                # Just recompute the suffixes that the join produced as the projection
                # would lose them otherwise.
                if (
                    left_suffix and idx < right_operand_index
                    and col in right_columns_set
                ):
                    col += left_suffix
                if (
                    right_suffix and idx >= right_operand_index
                    and col in left_columns_set
                ):
                    col += right_suffix
                projected_col_names.append(col)
                projections.append(
                    Expression._field(idx)
                )
        projection = Declaration(
            "project", ProjectNodeOptions(projections, projected_col_names)
        )
        decl = Declaration.from_sequence([decl, projection])

    result_table = decl.to_table(use_threads=use_threads)

    if output_type == Table:
        return result_table
    elif output_type == ds.InMemoryDataset:
        return ds.InMemoryDataset(result_table)
    else:
        raise TypeError("Unsupported output type")


def _perform_join_asof(left_operand, left_on, left_by,
                       right_operand, right_on, right_by,
                       tolerance, use_threads=True,
                       output_type=Table):
    """
    Perform asof join of two tables or datasets.

    The result will be an output table with the result of the join operation

    Parameters
    ----------
    left_operand : Table or Dataset
        The left operand for the join operation.
    left_on : str
        The left key (or keys) on which the join operation should be performed.
    left_by: str or list[str]
        The left key (or keys) on which the join operation should be performed.
    right_operand : Table or Dataset
        The right operand for the join operation.
    right_on : str or list[str]
        The right key (or keys) on which the join operation should be performed.
    right_by: str or list[str]
        The right key (or keys) on which the join operation should be performed.
    tolerance : int
        The tolerance to use for the asof join. The tolerance is interpreted in
        the same units as the "on" key.
    output_type: Table or InMemoryDataset
        The output type for the exec plan result.

    Returns
    -------
    result_table : Table or InMemoryDataset
    """
    if not isinstance(left_operand, (Table, ds.Dataset)):
        raise TypeError(f"Expected Table or Dataset, got {type(left_operand)}")
    if not isinstance(right_operand, (Table, ds.Dataset)):
        raise TypeError(f"Expected Table or Dataset, got {type(right_operand)}")

    if not isinstance(left_by, (tuple, list)):
        left_by = [left_by]
    if not isinstance(right_by, (tuple, list)):
        right_by = [right_by]

    # AsofJoin does not return on or by columns for right_operand.
    right_columns = [
        col for col in right_operand.schema.names
        if col not in [right_on] + right_by
    ]
    columns_collisions = set(left_operand.schema.names) & set(right_columns)
    if columns_collisions:
        raise ValueError(
            "Columns {} present in both tables. AsofJoin does not support "
            "column collisions.".format(columns_collisions),
        )

    # Add the join node to the execplan
    if isinstance(left_operand, ds.Dataset):
        left_source = _dataset_to_decl(
            left_operand,
            use_threads=use_threads,
            implicit_ordering=True)
    else:
        left_source = Declaration(
            "table_source", TableSourceNodeOptions(left_operand),
        )
    if isinstance(right_operand, ds.Dataset):
        right_source = _dataset_to_decl(
            right_operand, use_threads=use_threads,
            implicit_ordering=True)
    else:
        right_source = Declaration(
            "table_source", TableSourceNodeOptions(right_operand)
        )

    join_opts = AsofJoinNodeOptions(
        left_on, left_by, right_on, right_by, tolerance
    )
    decl = Declaration(
        "asofjoin", options=join_opts, inputs=[left_source, right_source]
    )

    result_table = decl.to_table(use_threads=use_threads)

    if output_type == Table:
        return result_table
    elif output_type == ds.InMemoryDataset:
        return ds.InMemoryDataset(result_table)
    else:
        raise TypeError("Unsupported output type")


def _filter_table(table, expression):
    """Filter rows of a table based on the provided expression.

    The result will be an output table with only the rows matching
    the provided expression.

    Parameters
    ----------
    table : Table or RecordBatch
        Table that should be filtered.
    expression : Expression
        The expression on which rows should be filtered.

    Returns
    -------
    Table
    """
    is_batch = False
    if isinstance(table, RecordBatch):
        table = Table.from_batches([table])
        is_batch = True

    decl = Declaration.from_sequence([
        Declaration("table_source", options=TableSourceNodeOptions(table)),
        Declaration("filter", options=FilterNodeOptions(expression))
    ])
    result = decl.to_table(use_threads=True)
    if is_batch:
        result = result.combine_chunks().to_batches()[0]
    return result


def _sort_source(table_or_dataset, sort_keys, output_type=Table, **kwargs):

    if isinstance(table_or_dataset, ds.Dataset):
        data_source = _dataset_to_decl(table_or_dataset, use_threads=True)
    else:
        data_source = Declaration(
            "table_source", TableSourceNodeOptions(table_or_dataset)
        )

    order_by = Declaration("order_by", OrderByNodeOptions(sort_keys, **kwargs))

    decl = Declaration.from_sequence([data_source, order_by])
    result_table = decl.to_table(use_threads=True)

    if output_type == Table:
        return result_table
    elif output_type == ds.InMemoryDataset:
        return ds.InMemoryDataset(result_table)
    else:
        raise TypeError("Unsupported output type")


def _group_by(table, aggregates, keys, use_threads=True):

    decl = Declaration.from_sequence([
        Declaration("table_source", TableSourceNodeOptions(table)),
        Declaration("aggregate", AggregateNodeOptions(aggregates, keys=keys))
    ])
    return decl.to_table(use_threads=use_threads)
