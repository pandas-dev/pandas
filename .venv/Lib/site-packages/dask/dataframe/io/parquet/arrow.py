from __future__ import annotations

import itertools
import json
import operator
import textwrap
from collections import defaultdict
from datetime import datetime
from functools import reduce

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# Check PyArrow version for feature support
from fsspec.core import expand_paths_if_needed, stringify_path
from fsspec.implementations.arrow import ArrowFSWrapper
from pyarrow import dataset as pa_ds
from pyarrow import fs as pa_fs

from dask.core import flatten
from dask.dataframe._compat import PANDAS_GE_220
from dask.dataframe.backends import pyarrow_schema_dispatch
from dask.dataframe.io.parquet.utils import (
    Engine,
    _get_aggregation_depth,
    _infer_split_row_groups,
    _normalize_index_columns,
    _process_open_file_options,
    _row_groups_to_parts,
    _set_gather_statistics,
    _set_metadata_task_size,
    _sort_and_analyze_paths,
)
from dask.dataframe.io.utils import _get_pyarrow_dtypes, _is_local_fs, _open_input_files
from dask.dataframe.utils import clear_known_categories, pyarrow_strings_enabled
from dask.delayed import Delayed
from dask.tokenize import normalize_token, tokenize
from dask.utils import getargspec, natural_sort_key

PYARROW_NULLABLE_DTYPE_MAPPING = {
    pa.int8(): pd.Int8Dtype(),
    pa.int16(): pd.Int16Dtype(),
    pa.int32(): pd.Int32Dtype(),
    pa.int64(): pd.Int64Dtype(),
    pa.uint8(): pd.UInt8Dtype(),
    pa.uint16(): pd.UInt16Dtype(),
    pa.uint32(): pd.UInt32Dtype(),
    pa.uint64(): pd.UInt64Dtype(),
    pa.bool_(): pd.BooleanDtype(),
    pa.string(): pd.StringDtype(),
    pa.float32(): pd.Float32Dtype(),
    pa.float64(): pd.Float64Dtype(),
}


@normalize_token.register(pa_fs.FileSystem)
def tokenize_arrowfs(obj):
    return obj.__reduce__()


#
#  Helper Utilities
#


def _wrapped_fs(fs):
    """Return the wrapped filesystem if fs is ArrowFSWrapper"""
    return fs.fs if isinstance(fs, ArrowFSWrapper) else fs


def _append_row_groups(metadata, md):
    """Append row-group metadata and include a helpful
    error message if an inconsistent schema is detected.
    """
    try:
        metadata.append_row_groups(md)
    except RuntimeError as err:
        if "requires equal schemas" in str(err):
            raise RuntimeError(
                "Schemas are inconsistent, try using "
                '`to_parquet(..., schema="infer")`, or pass an explicit '
                "pyarrow schema. Such as "
                '`to_parquet(..., schema={"column1": pa.string()})`'
            ) from err
        else:
            raise err


def _write_partitioned(
    table,
    df,
    root_path,
    filename,
    partition_cols,
    fs,
    pandas_to_arrow_table,
    preserve_index,
    index_cols=(),
    return_metadata=True,
    **kwargs,
):
    """Write table to a partitioned dataset with pyarrow.

    Logic copied from pyarrow.parquet.
    (arrow/python/pyarrow/parquet.py::write_to_dataset)

    TODO: Remove this in favor of pyarrow's `write_to_dataset`
          once ARROW-8244 is addressed.
    """
    fs.mkdirs(root_path, exist_ok=True)

    if preserve_index:
        df.reset_index(inplace=True)
    df = df[table.schema.names]

    index_cols = list(index_cols) if index_cols else []
    preserve_index = False
    if index_cols:
        df.set_index(index_cols, inplace=True)
        preserve_index = True

    partition_keys = [df[col] for col in partition_cols]
    data_df = df.drop(partition_cols, axis="columns")
    data_cols = df.columns.drop(partition_cols)
    if len(data_cols) == 0 and not index_cols:
        raise ValueError("No data left to save outside partition columns")

    subschema = table.schema
    for col in table.schema.names:
        if col in partition_cols:
            subschema = subschema.remove(subschema.get_field_index(col))

    md_list = []
    partition_keys = partition_keys[0] if len(partition_keys) == 1 else partition_keys
    gb = data_df.groupby(partition_keys, dropna=False, observed=False)
    for keys, subgroup in gb:
        if not isinstance(keys, tuple):
            keys = (keys,)
        subdir = fs.sep.join(
            [_hive_dirname(name, val) for name, val in zip(partition_cols, keys)]
        )
        subtable = pandas_to_arrow_table(
            subgroup, preserve_index=preserve_index, schema=subschema
        )
        prefix = fs.sep.join([root_path, subdir])
        fs.mkdirs(prefix, exist_ok=True)
        full_path = fs.sep.join([prefix, filename])
        with fs.open(full_path, "wb") as f:
            pq.write_table(
                subtable,
                f,
                metadata_collector=md_list if return_metadata else None,
                **kwargs,
            )
        if return_metadata:
            md_list[-1].set_file_path(fs.sep.join([subdir, filename]))

    return md_list


def _index_in_schema(index, schema):
    """Simple utility to check if all `index` columns are included
    in the known `schema`.
    """
    if index and schema is not None:
        # Make sure all index columns are in user-defined schema
        return len(set(index).intersection(schema.names)) == len(index)
    elif index:
        return True  # Schema is not user-specified, all good
    else:
        return False  # No index to check


class PartitionObj:
    """Simple class providing a `name` and `keys` attribute
    for a single partition column.

    This class was originally designed as a mechanism to build a
    duck-typed version of pyarrow's deprecated `ParquetPartitions`
    class. Now that `ArrowLegacyEngine` is deprecated, this class
    can be modified/removed, but it is still used as a convenience.
    """

    def __init__(self, name, keys):
        self.name = name
        self.keys = pd.Index(keys.sort_values(), copy=False)

    def __dask_tokenize__(self):
        return tokenize(self.name, self.keys)


def _frag_subset(old_frag, row_groups):
    """Create new fragment with row-group subset."""
    return old_frag.format.make_fragment(
        old_frag.path,
        old_frag.filesystem,
        old_frag.partition_expression,
        row_groups=row_groups,
    )


def _get_pandas_metadata(schema):
    """Get pandas-specific metadata from schema."""

    has_pandas_metadata = schema.metadata is not None and b"pandas" in schema.metadata
    if has_pandas_metadata:
        return json.loads(schema.metadata[b"pandas"].decode("utf8"))
    else:
        return {}


def _read_table_from_path(
    path,
    fs,
    row_groups,
    columns,
    schema,
    filters,
    **kwargs,
):
    """Read arrow table from file path.

    Used by `ArrowDatasetEngine._read_table` when no filters
    are specified (otherwise fragments are converted directly
    into tables).
    """

    # Define file-opening options
    read_kwargs = kwargs.get("read", {}).copy()
    precache_options, open_file_options = _process_open_file_options(
        read_kwargs.pop("open_file_options", {}),
        **(
            {
                "allow_precache": False,
                "default_cache": "none",
            }
            if _is_local_fs(fs)
            else {
                "columns": columns,
                "row_groups": row_groups if row_groups == [None] else [row_groups],
                "default_engine": "pyarrow",
                "default_cache": "none",
            }
        ),
    )

    # Use `pre_buffer=True` if the option is supported and an optimized
    # "pre-caching" method isn't already specified in `precache_options`
    # (The distinct fsspec and pyarrow optimizations will conflict)
    pre_buffer_default = precache_options.get("method", None) is None
    pre_buffer = {"pre_buffer": read_kwargs.pop("pre_buffer", pre_buffer_default)}

    with _open_input_files(
        [path],
        fs=fs,
        precache_options=precache_options,
        **open_file_options,
    )[0] as fil:
        if row_groups == [None]:
            return pq.ParquetFile(fil, **pre_buffer).read(
                columns=columns,
                use_threads=False,
                use_pandas_metadata=True,
                **read_kwargs,
            )
        else:
            return pq.ParquetFile(fil, **pre_buffer).read_row_groups(
                row_groups,
                columns=columns,
                use_threads=False,
                use_pandas_metadata=True,
                **read_kwargs,
            )


def _get_rg_statistics(row_group, col_names):
    """Custom version of pyarrow's RowGroupInfo.statistics method
    (https://github.com/apache/arrow/blob/master/python/pyarrow/_dataset.pyx)

    We use column names to specify the specific subset of columns
    that we need statistics for.  This is more optimal than the
    upstream `RowGroupInfo.statistics` method, which will return
    statistics for all columns.
    """

    row_group_schema = dict(
        zip(
            row_group.schema.names,
            itertools.accumulate(
                [
                    # Need to account for multi-field struct columns
                    max(row_group.schema.types[i].num_fields, 1)
                    for i in range(len(row_group.schema.names) - 1)
                ],
                initial=0,
            ),
        )
    )

    def name_stats(column_name):
        col = row_group.metadata.column(row_group_schema[column_name])

        stats = col.statistics
        if stats is None or not stats.has_min_max:
            return None, None

        name = col.path_in_schema
        field_index = row_group.schema.get_field_index(name)
        if field_index < 0:
            return None, None

        return col.path_in_schema, {
            "min": stats.min,
            "max": stats.max,
            "null_count": stats.null_count,
        }

    return {
        name: stats for name, stats in map(name_stats, col_names) if stats is not None
    }


def _need_filtering(filters, partition_keys):
    # Check if we need to generate a fragment for filtering.
    # We only need to do this if we are applying filters to
    # columns that were not already filtered by "partition".

    partition_cols = (
        {v[0] for v in flatten(partition_keys, container=list) if len(v)}
        if partition_keys
        else set()
    )
    filtered_cols = (
        {v[0] for v in flatten(filters, container=list) if len(v)} if filters else set()
    )

    return bool(filtered_cols - partition_cols)


def _hive_dirname(name, val):
    # Simple utility to produce hive directory name.
    # Note that "__HIVE_DEFAULT_PARTITION__" is the
    # conventional "null" label in other platforms
    val = "__HIVE_DEFAULT_PARTITION__" if pd.isna(val) else val
    return f"{name}={val}"


def _process_kwargs(partitioning=None, **kwargs):
    # Pre-process a dict of `pyarrow.dataset.dataset`` key-word
    # arguments. Primary purpose is to convert a dictionary-based
    # "partitioning" option into a proper `Partitioning` object
    return {
        "partitioning": (
            pa_ds.partitioning(**partitioning)
            if isinstance(partitioning, dict)
            else partitioning
        ),
        **kwargs,
    }


def _filters_to_expression(filters, propagate_null=False, nan_is_null=True):
    # Mostly copied from: pq.filters_to_expression
    # TODO: Use pq.filters_to_expression if/when null-value
    # handling is resolved.
    # See: https://github.com/dask/dask/issues/9845

    if isinstance(filters, pa_ds.Expression):
        return filters

    if filters is not None:
        if len(filters) == 0 or any(len(f) == 0 for f in filters):
            raise ValueError("Malformed filters")
        if isinstance(filters[0][0], str):
            # We have encountered the situation where we have one nesting level
            # too few:
            #   We have [(,,), ..] instead of [[(,,), ..]]
            filters = [filters]

    def convert_single_predicate(col, op, val):
        field = pa_ds.field(col)

        # Handle null-value comparison
        if val is None or (nan_is_null and val is np.nan):
            if op == "is":
                return field.is_null(nan_is_null=nan_is_null)
            elif op == "is not":
                return ~field.is_null(nan_is_null=nan_is_null)
            else:
                raise ValueError(
                    f'"{(col, op, val)}" is not a supported predicate '
                    f'Please use "is" or "is not" for null comparison.'
                )

        if op == "=" or op == "==":
            expr = field == val
        elif op == "!=":
            expr = field != val
        elif op == "<":
            expr = field < val
        elif op == ">":
            expr = field > val
        elif op == "<=":
            expr = field <= val
        elif op == ">=":
            expr = field >= val
        elif op == "in":
            expr = field.isin(val)
        elif op == "not in":
            expr = ~field.isin(val)
        else:
            raise ValueError(
                f'"{(col, op, val)}" is not a valid operator in predicates.'
            )

        # (Optionally) Avoid null-value propagation
        if not propagate_null and op in ("!=", "not in"):
            return field.is_null(nan_is_null=nan_is_null) | expr
        return expr

    disjunction_members = []

    for conjunction in filters:
        conjunction_members = [
            convert_single_predicate(col, op, val) for col, op, val in conjunction
        ]

        disjunction_members.append(reduce(operator.and_, conjunction_members))

    return reduce(operator.or_, disjunction_members)


#
#  ArrowDatasetEngine
#


class ArrowDatasetEngine(Engine):
    #
    # Public Class Methods
    #

    @classmethod
    def extract_filesystem(
        cls,
        urlpath,
        filesystem,
        dataset_options,
        open_file_options,
        storage_options,
    ):
        # Check if filesystem was specified as a dataset option
        if filesystem is None:
            fs = dataset_options.pop("filesystem", "fsspec")
        else:
            if "filesystem" in dataset_options:
                raise ValueError(
                    "Cannot specify a filesystem argument if the "
                    "'filesystem' dataset option is also defined."
                )
            fs = filesystem

        # Handle pyarrow-based filesystem
        if isinstance(fs, pa_fs.FileSystem) or fs in ("arrow", "pyarrow"):
            if isinstance(urlpath, (list, tuple, set)):
                if not urlpath:
                    raise ValueError("empty urlpath sequence")
                urlpath = [stringify_path(u) for u in urlpath]
            else:
                urlpath = [stringify_path(urlpath)]

            if fs in ("arrow", "pyarrow"):
                fs = pa_fs.FileSystem.from_uri(urlpath[0])[0]
                if storage_options:
                    # Use inferred region as the default
                    region = (
                        {} if "region" in storage_options else {"region": fs.region}
                    )
                    fs = type(fs)(**region, **storage_options)

            fsspec_fs = ArrowFSWrapper(fs)
            if urlpath[0].startswith("C:") and isinstance(fs, pa_fs.LocalFileSystem):
                # ArrowFSWrapper._strip_protocol not reliable on windows
                # See: https://github.com/fsspec/filesystem_spec/issues/1137
                from fsspec.implementations.local import LocalFileSystem

                fs_strip = LocalFileSystem()
            else:
                fs_strip = fsspec_fs
            paths = expand_paths_if_needed(urlpath, "rb", 1, fsspec_fs, None)
            return (
                fsspec_fs,
                [fs_strip._strip_protocol(u) for u in paths],
                dataset_options,
                {"open_file_func": fs.open_input_file},
            )

        # Use default file-system initialization
        return Engine.extract_filesystem(
            urlpath,
            fs,
            dataset_options,
            open_file_options,
            storage_options,
        )

    @classmethod
    def multi_support(cls):
        return cls == ArrowDatasetEngine

    @classmethod
    def read_partition(
        cls,
        fs,
        pieces,
        columns,
        index,
        dtype_backend=None,
        categories=(),
        partitions=(),
        filters=None,
        schema=None,
        **kwargs,
    ):
        """Read in a single output partition"""

        if isinstance(index, list):
            for level in index:
                # unclear if we can use set ops here. I think the order matters.
                # Need the membership test to avoid duplicating index when
                # we slice with `columns` later on.
                if level not in columns:
                    columns.append(level)

        # Ensure `columns` and `partitions` do not overlap
        columns_and_parts = columns.copy()
        if not isinstance(partitions, (list, tuple)):
            if columns_and_parts and partitions:
                for part_name in partitions.partition_names:
                    if part_name in columns:
                        columns.remove(part_name)
                    else:
                        columns_and_parts.append(part_name)
                columns = columns or None

        # Always convert pieces to list
        if not isinstance(pieces, list):
            pieces = [pieces]

        tables = []
        multi_read = len(pieces) > 1
        for piece in pieces:
            if isinstance(piece, str):
                # `piece` is a file-path string
                path_or_frag = piece
                row_group = None
                partition_keys = None
            else:
                # `piece` contains (path, row_group, partition_keys)
                (path_or_frag, row_group, partition_keys) = piece

            # Convert row_group to a list and be sure to
            # check if msgpack converted it to a tuple
            if isinstance(row_group, tuple):
                row_group = list(row_group)
            if not isinstance(row_group, list):
                row_group = [row_group]

            # Read in arrow table and convert to pandas
            arrow_table = cls._read_table(
                path_or_frag,
                fs,
                row_group,
                columns,
                schema,
                filters,
                partitions,
                partition_keys,
                **kwargs,
            )
            if multi_read:
                tables.append(arrow_table)

        if multi_read:
            arrow_table = pa.concat_tables(tables)

        # Convert to pandas
        df = cls._arrow_table_to_pandas(
            arrow_table, categories, dtype_backend=dtype_backend, **kwargs
        )

        # For pyarrow.dataset api, need to convert partition columns
        # to categorigal manually for integer types.
        if partitions and isinstance(partitions, list):
            for partition in partitions:
                if len(partition.keys) and df[partition.name].dtype.name != "category":
                    # We read directly from fragments, so the partition
                    # columns are already in our dataframe.  We just
                    # need to convert non-categorical types.
                    df[partition.name] = pd.Series(
                        pd.Categorical(
                            categories=partition.keys,
                            values=df[partition.name].values,
                        ),
                        index=df.index,
                    )

        # Note that `to_pandas(ignore_metadata=False)` means
        # pyarrow will use the pandas metadata to set the index.
        index_in_columns_and_parts = set(df.index.names).issubset(
            set(columns_and_parts)
        )
        if not index:
            if index_in_columns_and_parts:
                # User does not want to set index and a desired
                # column/partition has been set to the index
                df.reset_index(drop=False, inplace=True)
            else:
                # User does not want to set index and an
                # "unwanted" column has been set to the index
                df.reset_index(drop=True, inplace=True)
        else:
            if set(df.index.names) != set(index) and index_in_columns_and_parts:
                # The wrong index has been set and it contains
                # one or more desired columns/partitions
                df.reset_index(drop=False, inplace=True)
            elif index_in_columns_and_parts:
                # The correct index has already been set
                index = False
                columns_and_parts = list(set(columns_and_parts) - set(df.index.names))
        df = df[list(columns_and_parts)]

        if index:
            df = df.set_index(index)
        return df

    @classmethod
    def initialize_write(
        cls,
        df,
        fs,
        path,
        append=False,
        partition_on=None,
        ignore_divisions=False,
        division_info=None,
        schema="infer",
        index_cols=None,
        **kwargs,
    ):
        if schema == "infer" or isinstance(schema, dict):
            # Start with schema from _meta_nonempty
            inferred_schema = pyarrow_schema_dispatch(
                df._meta_nonempty.set_index(index_cols)
                if index_cols
                else df._meta_nonempty
            ).remove_metadata()

            # Use dict to update our inferred schema
            if isinstance(schema, dict):
                schema = pa.schema(schema)
                for name in schema.names:
                    i = inferred_schema.get_field_index(name)
                    j = schema.get_field_index(name)
                    inferred_schema = inferred_schema.set(i, schema.field(j))
            schema = inferred_schema

        # Check that target directory exists
        fs.mkdirs(path, exist_ok=True)
        if append and division_info is None:
            ignore_divisions = True

        full_metadata = None  # metadata for the full dataset, from _metadata
        tail_metadata = None  # metadata for at least the last file in the dataset
        i_offset = 0
        metadata_file_exists = False
        if append:
            # Extract metadata and get file offset if appending
            ds = pa_ds.dataset(path, filesystem=_wrapped_fs(fs), format="parquet")
            i_offset = len(ds.files)
            if i_offset > 0:
                try:
                    with fs.open(fs.sep.join([path, "_metadata"]), mode="rb") as fil:
                        full_metadata = pq.read_metadata(fil)
                    tail_metadata = full_metadata
                    metadata_file_exists = True
                except OSError:
                    try:
                        with fs.open(
                            sorted(ds.files, key=natural_sort_key)[-1], mode="rb"
                        ) as fil:
                            tail_metadata = pq.read_metadata(fil)
                    except OSError:
                        pass
            else:
                append = False  # No existing files, can skip the append logic

        # If appending, validate against the initial metadata file (if present)
        if append and tail_metadata is not None:
            arrow_schema = tail_metadata.schema.to_arrow_schema()
            names = arrow_schema.names
            has_pandas_metadata = (
                arrow_schema.metadata is not None and b"pandas" in arrow_schema.metadata
            )
            if has_pandas_metadata:
                pandas_metadata = json.loads(
                    arrow_schema.metadata[b"pandas"].decode("utf8")
                )
                categories = [
                    c["name"]
                    for c in pandas_metadata["columns"]
                    if c["pandas_type"] == "categorical"
                ]
            else:
                categories = None
            dtypes = _get_pyarrow_dtypes(arrow_schema, categories)
            if set(names) != set(df.columns) - set(partition_on):
                raise ValueError(
                    "Appended columns not the same.\n"
                    "Previous: {} | New: {}".format(names, list(df.columns))
                )
            elif (pd.Series(dtypes).loc[names] != df[names].dtypes).any():
                # TODO Coerce values for compatible but different dtypes
                raise ValueError(
                    "Appended dtypes differ.\n{}".format(
                        set(dtypes.items()) ^ set(df.dtypes.items())
                    )
                )

            # Check divisions if necessary
            if division_info["name"] not in names:
                ignore_divisions = True
            if not ignore_divisions:
                old_end = None
                row_groups = (
                    tail_metadata.row_group(i)
                    for i in range(tail_metadata.num_row_groups)
                )
                index_col_i = names.index(division_info["name"])
                for row_group in row_groups:
                    column = row_group.column(index_col_i)
                    if column.statistics:
                        if old_end is None:
                            old_end = column.statistics.max
                        elif column.statistics.max > old_end:
                            old_end = column.statistics.max
                        else:
                            # Existing column on disk isn't sorted, set
                            # `old_end = None` to skip check below
                            old_end = None
                            break

                divisions = division_info["divisions"]
                if old_end is not None and divisions[0] <= old_end:
                    raise ValueError(
                        "The divisions of the appended dataframe overlap with "
                        "previously written divisions. If this is desired, set "
                        "``ignore_divisions=True`` to append anyway.\n"
                        "- End of last written partition: {old_end}\n"
                        "- Start of first new partition: {divisions[0]}"
                    )

        extra_write_kwargs = {"schema": schema, "index_cols": index_cols}
        return i_offset, full_metadata, metadata_file_exists, extra_write_kwargs

    @classmethod
    def _pandas_to_arrow_table(
        cls, df: pd.DataFrame, preserve_index=False, schema=None
    ) -> pa.Table:
        try:
            return pa.Table.from_pandas(
                df, nthreads=1, preserve_index=preserve_index, schema=schema
            )
        except pa.ArrowException as exc:
            if schema is None:
                raise
            df_schema = pa.Schema.from_pandas(df)
            expected = textwrap.indent(
                schema.to_string(show_schema_metadata=False), "    "
            )
            actual = textwrap.indent(
                df_schema.to_string(show_schema_metadata=False), "    "
            )
            raise ValueError(
                f"Failed to convert partition to expected pyarrow schema:\n"
                f"    `{exc!r}`\n"
                f"\n"
                f"Expected partition schema:\n"
                f"{expected}\n"
                f"\n"
                f"Received partition schema:\n"
                f"{actual}\n"
                f"\n"
                f"This error *may* be resolved by passing in schema information for\n"
                f"the mismatched column(s) using the `schema` keyword in `to_parquet`."
            ) from None

    @classmethod
    def write_partition(
        cls,
        df,
        path,
        fs,
        filename,
        partition_on,
        return_metadata,
        fmd=None,
        compression=None,
        index_cols=None,
        schema=None,
        head=False,
        custom_metadata=None,
        **kwargs,
    ):
        _meta = None
        preserve_index = False
        if _index_in_schema(index_cols, schema):
            df.set_index(index_cols, inplace=True)
            preserve_index = True
        else:
            index_cols = []

        t = cls._pandas_to_arrow_table(df, preserve_index=preserve_index, schema=schema)
        if custom_metadata:
            _md = t.schema.metadata
            _md.update(custom_metadata)
            t = t.replace_schema_metadata(metadata=_md)

        if partition_on:
            md_list = _write_partitioned(
                t,
                df,
                path,
                filename,
                partition_on,
                fs,
                cls._pandas_to_arrow_table,
                preserve_index,
                index_cols=index_cols,
                compression=compression,
                return_metadata=return_metadata,
                **kwargs,
            )
            if md_list:
                _meta = md_list[0]
                for i in range(1, len(md_list)):
                    _append_row_groups(_meta, md_list[i])
        else:
            md_list = []
            with fs.open(fs.sep.join([path, filename]), "wb") as fil:
                pq.write_table(
                    t,
                    fil,
                    compression=compression,
                    metadata_collector=md_list if return_metadata else None,
                    **kwargs,
                )
            if md_list:
                _meta = md_list[0]
                _meta.set_file_path(filename)
        # Return the schema needed to write the metadata
        if return_metadata:
            d = {"meta": _meta}
            if head:
                # Only return schema if this is the "head" partition
                d["schema"] = t.schema
            return [d]
        else:
            return []

    @classmethod
    def write_metadata(cls, parts, meta, fs, path, append=False, **kwargs):
        schema = parts[0][0].get("schema", None)
        parts = [p for p in parts if p[0]["meta"] is not None]
        if parts:
            if not append:
                # Get only arguments specified in the function
                common_metadata_path = fs.sep.join([path, "_common_metadata"])
                keywords = getargspec(pq.write_metadata).args
                kwargs_meta = {k: v for k, v in kwargs.items() if k in keywords}
                with fs.open(common_metadata_path, "wb") as fil:
                    pq.write_metadata(schema, fil, **kwargs_meta)

            # Aggregate metadata and write to _metadata file
            metadata_path = fs.sep.join([path, "_metadata"])
            if append and meta is not None:
                _meta = meta
                i_start = 0
            else:
                _meta = parts[0][0]["meta"]
                i_start = 1
            for i in range(i_start, len(parts)):
                _append_row_groups(_meta, parts[i][0]["meta"])
            with fs.open(metadata_path, "wb") as fil:
                _meta.write_metadata_file(fil)

    #
    # Private Class Methods
    #

    @classmethod
    def _collect_dataset_info(
        cls,
        paths,
        fs,
        categories,
        index,
        gather_statistics,
        filters,
        split_row_groups,
        blocksize,
        aggregate_files,
        ignore_metadata_file,
        metadata_task_size,
        parquet_file_extension,
        kwargs,
    ):
        """pyarrow.dataset version of _collect_dataset_info
        Use pyarrow.dataset API to construct a dictionary of all
        general information needed to read the dataset.
        """

        # Use pyarrow.dataset API
        ds = None
        valid_paths = None  # Only used if `paths` is a list containing _metadata

        # Extract dataset-specific options
        _dataset_kwargs = kwargs.pop("dataset", {})

        if "partitioning" not in _dataset_kwargs:
            _dataset_kwargs["partitioning"] = "hive"

        if "format" not in _dataset_kwargs:
            _dataset_kwargs["format"] = pa_ds.ParquetFileFormat()
        _processed_dataset_kwargs = _process_kwargs(**_dataset_kwargs)

        # Case-dependent pyarrow.dataset creation
        has_metadata_file = False
        if len(paths) == 1 and fs.isdir(paths[0]):
            # Use _analyze_paths to avoid relative-path
            # problems (see GH#5608)
            paths, base, fns = _sort_and_analyze_paths(paths, fs)
            paths = fs.sep.join([base, fns[0]])

            meta_path = fs.sep.join([paths, "_metadata"])
            if not ignore_metadata_file and fs.exists(meta_path):
                # Use _metadata file
                ds = pa_ds.parquet_dataset(
                    meta_path,
                    filesystem=_wrapped_fs(fs),
                    **_processed_dataset_kwargs,
                )
                has_metadata_file = True
            elif parquet_file_extension:
                # Need to materialize all paths if we are missing the _metadata file
                # Raise error if all files have been filtered by extension
                len0 = len(paths)
                paths = [
                    path
                    for path in fs.find(paths)
                    if path.endswith(parquet_file_extension)
                ]
                if len0 and paths == []:
                    raise ValueError(
                        "No files satisfy the `parquet_file_extension` criteria "
                        f"(files must end with {parquet_file_extension})."
                    )

        elif len(paths) > 1:
            paths, base, fns = _sort_and_analyze_paths(paths, fs)
            meta_path = fs.sep.join([base, "_metadata"])
            if "_metadata" in fns:
                # Pyarrow cannot handle "_metadata" when `paths` is a list
                # Use _metadata file
                if not ignore_metadata_file:
                    ds = pa_ds.parquet_dataset(
                        meta_path,
                        filesystem=_wrapped_fs(fs),
                        **_processed_dataset_kwargs,
                    )
                    has_metadata_file = True

                # Populate valid_paths, since the original path list
                # must be used to filter the _metadata-based dataset
                fns.remove("_metadata")
                valid_paths = fns

        # Final "catch-all" pyarrow.dataset call
        if ds is None:
            ds = pa_ds.dataset(
                paths,
                filesystem=_wrapped_fs(fs),
                **_processed_dataset_kwargs,
            )

        # Get file_frag sample and extract physical_schema
        try:
            file_frag = next(ds.get_fragments())
            physical_schema = file_frag.physical_schema
        except StopIteration:
            file_frag = None
            physical_schema = ds.schema

        # Set split_row_groups for desired partitioning behavior
        #
        # Expected behavior for split_row_groups + blocksize combinations:
        # +======+==================+===========+=============================+
        # | Case | split_row_groups | blocksize | Behavior                    |
        # +======+==================+===========+=============================+
        # |  A   |  "infer"         |  not None | Go to E or G (using md)     |
        # +------+------------------+-----------+-----------------------------+
        # |  B   |  "infer"         |  None     | Go to H                     |
        # +------+------------------+-----------+-----------------------------+
        # |  C   |  "adaptive"      |  not None | Go to E                     |
        # +------+------------------+-----------+-----------------------------+
        # |  D   |  "adaptive"      |  None     | Go to H                     |
        # +======+==================+===========+=============================+
        # |  E*  |  True            |  not None | Adaptive partitioning       |
        # +------+------------------+-----------+-----------------------------+
        # |  F   |  True            |  None     | 1 row-group per partition   |
        # +------+------------------+-----------+-----------------------------+
        # |  G*  |  False           |  not None | 1+ full files per partition |
        # +------+------------------+-----------+-----------------------------+
        # |  H   |  False           |  None     | 1 full file per partition   |
        # +------+------------------+-----------+-----------------------------+
        # |  I   |  n               |  N/A      | n row-groups per partition  |
        # +======+==================+===========+=============================+
        # NOTES:
        # - Adaptive partitioning (E) means that the individual size of each
        #   row-group will be accounted for when deciding how many row-groups
        #   to map to each output partition.
        # - E, G and I will only aggregate data from multiple files into the
        #   same output partition if `bool(aggregate_files) == True`.
        # - Default partitioning will correspond to either E or G. All other
        #   behavior requires user input.

        if split_row_groups == "infer":
            if blocksize:
                # Sample row-group sizes in first file
                if file_frag is None:
                    # Empty dataset
                    split_row_groups = False
                else:
                    split_row_groups = _infer_split_row_groups(
                        [rg.total_byte_size for rg in file_frag.row_groups],
                        blocksize,
                        bool(aggregate_files),
                    )
            else:
                split_row_groups = False

        if split_row_groups == "adaptive":
            if blocksize:
                split_row_groups = True
            else:
                split_row_groups = False

        # Note on (hive) partitioning information:
        #
        #    - "partition_obj" : (list of PartitionObj) This is a list of
        #          simple objects providing `name` and `keys` attributes
        #          for each partition column.
        #    - "partition_names" : (list)  This is a list containing the
        #          names of partitioned columns.
        #
        partition_obj, partition_names = [], []
        if ds.partitioning and ds.partitioning.schema:
            partition_names = list(ds.partitioning.schema.names)
            for i, name in enumerate(partition_names):
                dictionary = (
                    ds.partitioning.dictionaries[i]
                    if ds.partitioning.dictionaries
                    else None
                )
                partition_obj.append(
                    PartitionObj(
                        name,
                        (
                            pd.Series([], dtype="object")
                            if dictionary is None
                            else dictionary.to_pandas()
                        ),
                    )
                )

        # Check the `aggregate_files` setting
        aggregation_depth = _get_aggregation_depth(aggregate_files, partition_names)

        return {
            "ds": ds,
            "physical_schema": physical_schema,
            "has_metadata_file": has_metadata_file,
            "schema": ds.schema,
            "fs": fs,
            "valid_paths": valid_paths,
            "gather_statistics": gather_statistics,
            "categories": categories,
            "index": index,
            "filters": filters,
            "split_row_groups": split_row_groups,
            "blocksize": blocksize,
            "aggregate_files": aggregate_files,
            "aggregation_depth": aggregation_depth,
            "partitions": partition_obj,
            "partition_names": partition_names,
            "metadata_task_size": metadata_task_size,
            "kwargs": {
                "dataset": _dataset_kwargs,
                "convert_string": pyarrow_strings_enabled(),
                **kwargs,
            },
        }

    @classmethod
    def _create_dd_meta(cls, dataset_info):
        """Use parquet schema and hive-partition information
        (stored in dataset_info) to construct DataFrame metadata.
        """

        # Collect necessary information from dataset_info
        schema = dataset_info["schema"]
        index = dataset_info["index"]
        categories = dataset_info["categories"]
        partition_obj = dataset_info["partitions"]
        partitions = dataset_info["partition_names"]
        physical_column_names = dataset_info.get("physical_schema", schema).names
        columns = None

        # Use pandas metadata to update categories
        pandas_metadata = _get_pandas_metadata(schema) or {}
        if pandas_metadata:
            if categories is None:
                categories = []
                for col in pandas_metadata["columns"]:
                    if (col["pandas_type"] == "categorical") and (
                        col["name"] not in categories
                    ):
                        categories.append(col["name"])

        # Use _arrow_table_to_pandas to generate meta
        arrow_to_pandas = dataset_info["kwargs"].get("arrow_to_pandas", {}).copy()
        convert_string = dataset_info["kwargs"]["convert_string"]
        dtype_backend = dataset_info["kwargs"]["dtype_backend"]
        meta = cls._arrow_table_to_pandas(
            schema.empty_table(),
            categories,
            arrow_to_pandas=arrow_to_pandas,
            dtype_backend=dtype_backend,
            convert_string=convert_string,
        )
        index_names = list(meta.index.names)
        column_names = list(meta.columns)
        if index_names and index_names != [None]:
            # Reset the index if non-null index name
            meta.reset_index(inplace=True)

        # Use index specified in the pandas metadata if
        # the index column was not specified by the user
        if (
            index is None
            and index_names
            and (
                # Only set to `[None]` if pandas metadata includes an index
                index_names != [None]
                or pandas_metadata.get("index_columns", None)
            )
        ):
            index = index_names

        # Set proper index for meta
        index_cols = index or ()
        if index_cols and index_cols != [None]:
            meta.set_index(index_cols, inplace=True)

        # Ensure that there is no overlap between partition columns
        # and explicit column storage
        if partitions:
            _partitions = [p for p in partitions if p not in physical_column_names]
            if not _partitions:
                partitions = []
                dataset_info["partitions"] = []
                dataset_info["partition_keys"] = {}
                dataset_info["partition_names"] = partitions
            elif len(_partitions) != len(partitions):
                raise ValueError(
                    "No partition-columns should be written in the \n"
                    "file unless they are ALL written in the file.\n"
                    "physical columns: {} | partitions: {}".format(
                        physical_column_names, partitions
                    )
                )

        # Get all available column names
        column_names, index_names = _normalize_index_columns(
            columns, column_names + partitions, index, index_names
        )
        all_columns = index_names + column_names

        if categories:
            # Check that categories are included in columns
            if not set(categories).intersection(all_columns):
                raise ValueError(
                    "categories not in available columns.\n"
                    "categories: {} | columns: {}".format(categories, list(all_columns))
                )

            # Make sure all categories are set to "unknown".
            # Cannot include index names in the `cols` argument.
            meta = clear_known_categories(
                meta,
                cols=[c for c in categories if c not in meta.index.names],
                dtype_backend=dtype_backend,
            )

        if partition_obj:
            # Update meta dtypes for partitioned columns
            for partition in partition_obj:
                if not len(partition.keys):
                    continue
                if isinstance(index, list) and partition.name == index[0]:
                    # Index from directory structure
                    meta.index = pd.CategoricalIndex(
                        [], categories=partition.keys, name=index[0]
                    )
                elif partition.name == meta.index.name:
                    # Index created from a categorical column
                    meta.index = pd.CategoricalIndex(
                        [], categories=partition.keys, name=meta.index.name
                    )
                elif partition.name in meta.columns:
                    meta[partition.name] = pd.Series(
                        pd.Categorical(categories=partition.keys, values=[]),
                        index=meta.index,
                    )

        # Update `dataset_info` and return `meta`
        dataset_info["index"] = index
        dataset_info["index_cols"] = index_cols
        dataset_info["categories"] = categories

        return meta

    @classmethod
    def _construct_collection_plan(cls, dataset_info):
        """pyarrow.dataset version of _construct_collection_plan
        Use dataset_info to construct the general plan for
        generating the output DataFrame collection.

        The "plan" is essentially a list (called `parts`) of
        information that is needed to produce each output partition.
        After this function is returned, the information in each
        element of `parts` will be used to produce a single Dask-
        DataFrame partition (unless some elements of `parts`
        are aggregated together in a follow-up step).

        This method also returns ``stats`` (which is a list of
        parquet-metadata statistics for each element of parts),
        and ``common_metadata`` (which is a dictionary of kwargs
        that should be passed to the ``read_partition`` call for
        every output partition).
        """

        # Collect necessary dataset information from dataset_info
        ds = dataset_info["ds"]
        fs = dataset_info["fs"]
        filters = dataset_info["filters"]
        split_row_groups = dataset_info["split_row_groups"]
        gather_statistics = dataset_info["gather_statistics"]
        blocksize = dataset_info["blocksize"]
        aggregation_depth = dataset_info["aggregation_depth"]
        index_cols = dataset_info["index_cols"]
        schema = dataset_info["schema"]
        partition_names = dataset_info["partition_names"]
        partitions = dataset_info["partitions"]
        categories = dataset_info["categories"]
        has_metadata_file = dataset_info["has_metadata_file"]
        valid_paths = dataset_info["valid_paths"]
        kwargs = dataset_info["kwargs"]

        # Ensure metadata_task_size is set
        # (Using config file or defaults)
        metadata_task_size = _set_metadata_task_size(
            dataset_info["metadata_task_size"], fs
        )

        # Make sure that any `in`-predicate filters have iterable values
        filter_columns = set()
        if filters is not None:
            for filter in flatten(filters, container=list):
                col, op, val = filter
                if op == "in" and not isinstance(val, (set, list, tuple)):
                    raise TypeError(
                        "Value of 'in' filter must be a list, set or tuple."
                    )
                filter_columns.add(col)

        # Determine which columns need statistics.
        # At this point, gather_statistics is only True if
        # the user specified calculate_divisions=True
        stat_col_indices = {}
        _index_cols = index_cols if (gather_statistics and len(index_cols) == 1) else []
        for i, name in enumerate(schema.names):
            if name in _index_cols or name in filter_columns:
                if name in partition_names:
                    # Partition columns won't have statistics
                    continue
                stat_col_indices[name] = i

        # Decide final `gather_statistics` setting
        gather_statistics = _set_gather_statistics(
            gather_statistics,
            blocksize,
            split_row_groups,
            aggregation_depth,
            filter_columns,
            set(stat_col_indices),
        )

        # Add common kwargs
        common_kwargs = {
            "partitions": partitions,
            "categories": categories,
            "filters": filters,
            "schema": schema,
            **kwargs,
        }

        # Check if this is a very simple case where we can just return
        # the path names
        if gather_statistics is False and not (split_row_groups or filters):
            return (
                [
                    {"piece": (full_path, None, None)}
                    for full_path in sorted(ds.files, key=natural_sort_key)
                ],
                [],
                common_kwargs,
            )

        # Get/translate filters
        ds_filters = None
        if filters is not None:
            ds_filters = _filters_to_expression(filters)

        # Define subset of `dataset_info` required by _collect_file_parts
        dataset_info_kwargs = {
            "fs": fs,
            "split_row_groups": split_row_groups,
            "gather_statistics": gather_statistics,
            "filters": filters,
            "ds_filters": ds_filters,
            "schema": schema,
            "stat_col_indices": stat_col_indices,
            "aggregation_depth": aggregation_depth,
            "blocksize": blocksize,
            "partitions": partitions,
            "dataset_options": kwargs["dataset"],
        }

        # Main parts/stats-construction
        if (
            has_metadata_file
            or metadata_task_size == 0
            or metadata_task_size > len(ds.files)
        ):
            # We have a global _metadata file to work with.
            # Therefore, we can just loop over fragments on the client.

            # Start with sorted (by path) list of file-based fragments
            file_frags = sorted(
                (frag for frag in ds.get_fragments(ds_filters)),
                key=lambda x: natural_sort_key(x.path),
            )
            parts, stats = cls._collect_file_parts(file_frags, dataset_info_kwargs)
        else:
            # We DON'T have a global _metadata file to work with.
            # We should loop over files in parallel

            if filters and partitions:
                # Start with sorted (by path) list of file-based fragments
                all_files = sorted(
                    (frag for frag in ds.get_fragments(ds_filters)),
                    key=lambda x: natural_sort_key(x.path),
                )
            else:
                # Collect list of file paths.
                # If valid_paths is not None, the user passed in a list
                # of files containing a _metadata file.  Since we used
                # the _metadata file to generate our dataset object , we need
                # to ignore any file fragments that are not in the list.
                all_files = sorted(ds.files, key=natural_sort_key)
                if valid_paths:
                    all_files = [
                        filef
                        for filef in all_files
                        if filef.split(fs.sep)[-1] in valid_paths
                    ]

            parts, stats = [], []
            if all_files:
                # Build and compute a task graph to construct stats/parts
                gather_parts_dsk = {}
                name = "gather-pq-parts-" + tokenize(all_files, dataset_info_kwargs)
                finalize_list = []
                for task_i, file_i in enumerate(
                    range(0, len(all_files), metadata_task_size)
                ):
                    finalize_list.append((name, task_i))
                    gather_parts_dsk[finalize_list[-1]] = (
                        cls._collect_file_parts,
                        all_files[file_i : file_i + metadata_task_size],
                        dataset_info_kwargs,
                    )

                def _combine_parts(parts_and_stats):
                    parts, stats = [], []
                    for part, stat in parts_and_stats:
                        parts += part
                        if stat:
                            stats += stat
                    return parts, stats

                gather_parts_dsk["final-" + name] = (_combine_parts, finalize_list)
                parts, stats = Delayed("final-" + name, gather_parts_dsk).compute()

        return parts, stats, common_kwargs

    @classmethod
    def _collect_file_parts(
        cls,
        files_or_frags,
        dataset_info_kwargs,
    ):
        # Collect necessary information from dataset_info
        fs = dataset_info_kwargs["fs"]
        split_row_groups = dataset_info_kwargs["split_row_groups"]
        gather_statistics = dataset_info_kwargs["gather_statistics"]
        partitions = dataset_info_kwargs["partitions"]
        dataset_options = dataset_info_kwargs["dataset_options"]

        # Make sure we are processing a non-empty list
        if not isinstance(files_or_frags, list):
            files_or_frags = [files_or_frags]
        elif not files_or_frags:
            return [], []

        # Make sure we are starting with file fragments
        if isinstance(files_or_frags[0], str):
            # Check if we are using a simple file-partition map
            # without requiring any file or row-group statistics
            if not (split_row_groups or partitions) and gather_statistics is False:
                # Cool - We can return immediately
                return [
                    {"piece": (file_or_frag, None, None)}
                    for file_or_frag in files_or_frags
                ], None

            # Need more information - convert the path to a fragment
            file_frags = list(
                pa_ds.dataset(
                    files_or_frags,
                    filesystem=_wrapped_fs(fs),
                    **_process_kwargs(**dataset_options),
                ).get_fragments()
            )
        else:
            file_frags = files_or_frags

        # Collect settings from dataset_info
        filters = dataset_info_kwargs["filters"]
        ds_filters = dataset_info_kwargs["ds_filters"]
        schema = dataset_info_kwargs["schema"]
        stat_col_indices = dataset_info_kwargs["stat_col_indices"]
        aggregation_depth = dataset_info_kwargs["aggregation_depth"]
        blocksize = dataset_info_kwargs["blocksize"]

        # Initialize row-group and statistics data structures
        file_row_groups = defaultdict(list)
        file_row_group_stats = defaultdict(list)
        file_row_group_column_stats = defaultdict(list)
        single_rg_parts = int(split_row_groups) == 1
        hive_partition_keys = {}
        cmax_last = {}
        for file_frag in file_frags:
            fpath = file_frag.path

            # Extract hive-partition keys, and make sure they
            # are orederd the same as they are in `partitions`
            raw_keys = pa_ds._get_partition_keys(file_frag.partition_expression)
            hive_partition_keys[fpath] = [
                (hive_part.name, raw_keys[hive_part.name]) for hive_part in partitions
            ]

            for frag in file_frag.split_by_row_group(ds_filters, schema=schema):
                row_group_info = frag.row_groups
                if gather_statistics or split_row_groups:
                    # If we are gathering statistics or splitting by
                    # row-group, we may need to ensure our fragment
                    # metadata is complete.
                    if row_group_info is None:
                        frag.ensure_complete_metadata()
                        row_group_info = frag.row_groups
                    if not len(row_group_info):
                        continue
                else:
                    file_row_groups[fpath] = [None]
                    continue
                for row_group in row_group_info:
                    file_row_groups[fpath].append(row_group.id)
                    if gather_statistics:
                        statistics = _get_rg_statistics(
                            row_group, list(stat_col_indices)
                        )
                        if single_rg_parts:
                            s = {
                                "file_path_0": fpath,
                                "num-rows": row_group.num_rows,
                                "total_byte_size": row_group.total_byte_size,
                                "columns": [],
                            }
                        else:
                            s = {
                                "num-rows": row_group.num_rows,
                                "total_byte_size": row_group.total_byte_size,
                            }
                        cstats = []
                        for name in stat_col_indices.keys():
                            if name in statistics:
                                cmin = statistics[name]["min"]
                                cmax = statistics[name]["max"]
                                null_count = statistics[name]["null_count"]
                                cmin = (
                                    pd.Timestamp(cmin)
                                    if isinstance(cmin, datetime)
                                    else cmin
                                )
                                cmax = (
                                    pd.Timestamp(cmax)
                                    if isinstance(cmax, datetime)
                                    else cmax
                                )
                                last = cmax_last.get(name)
                                if not (
                                    filters
                                    or (blocksize and split_row_groups is True)
                                    or aggregation_depth
                                ):
                                    # Only think about bailing if we don't need
                                    # stats for filtering
                                    if cmin is None or (last and cmin < last):
                                        # We are collecting statistics for divisions
                                        # only (no filters) - Column isn't sorted, or
                                        # we have an all-null partition, so lets bail.
                                        #
                                        # Note: This assumes ascending order.
                                        #
                                        gather_statistics = False
                                        file_row_group_stats = {}
                                        file_row_group_column_stats = {}
                                        break

                                if single_rg_parts:
                                    s["columns"].append(
                                        {
                                            "name": name,
                                            "min": cmin,
                                            "max": cmax,
                                            "null_count": null_count,
                                        }
                                    )
                                else:
                                    cstats += [cmin, cmax, null_count]
                                cmax_last[name] = cmax
                            else:
                                if single_rg_parts:
                                    s["columns"].append({"name": name})
                                else:
                                    cstats += [None, None, None]
                        if gather_statistics:
                            file_row_group_stats[fpath].append(s)
                            if not single_rg_parts:
                                file_row_group_column_stats[fpath].append(tuple(cstats))

        # Check if we have empty parts to return
        if not file_row_groups:
            return [], []

        # Convert organized row-groups to parts
        return _row_groups_to_parts(
            gather_statistics,
            split_row_groups,
            aggregation_depth,
            file_row_groups,
            file_row_group_stats,
            file_row_group_column_stats,
            stat_col_indices,
            cls._make_part,
            make_part_kwargs={
                "fs": fs,
                "partition_keys": hive_partition_keys,
                "partition_obj": partitions,
                "data_path": "",
            },
        )

    @classmethod
    def _make_part(
        cls,
        filename,
        rg_list,
        fs=None,
        partition_keys=None,
        partition_obj=None,
        data_path=None,
    ):
        """Generate a partition-specific element of `parts`."""

        # Get full path (empty strings should be ignored)
        full_path = fs.sep.join([p for p in [data_path, filename] if p != ""])

        pkeys = partition_keys.get(full_path, None)
        if partition_obj and pkeys is None:
            return None  # This partition was filtered
        return {"piece": (full_path, rg_list, pkeys)}

    @classmethod
    def _read_table(
        cls,
        path_or_frag,
        fs,
        row_groups,
        columns,
        schema,
        filters,
        partitions,
        partition_keys,
        **kwargs,
    ):
        """Read in a pyarrow table"""

        if isinstance(path_or_frag, pa_ds.ParquetFileFragment):
            frag = path_or_frag

        else:
            frag = None

            # Check if we have partitioning information.
            partitioning = kwargs.get("dataset", {}).get("partitioning", None)

            # Check if we need to generate a fragment.
            # NOTE: We only need a fragment if we are doing row-wise
            # filtering, or if we are missing necessary information
            # about the hive/directory partitioning. For the case
            # of filtering, we only need a fragment if we are applying
            # filters to "un-partitioned" columns. Partitioned-column
            # filters should have been applied earlier.
            missing_partitioning_info = (
                # Need to discover partition_keys
                (partitions and partition_keys is None)
                # Need to apply custom partitioning schema
                or (partitioning and not isinstance(partitioning, (str, list)))
            )
            if missing_partitioning_info or _need_filtering(filters, partition_keys):
                # Convert the path and row-group IDs to a single fragment
                ds = pa_ds.dataset(
                    path_or_frag,
                    filesystem=_wrapped_fs(fs),
                    **_process_kwargs(**kwargs.get("dataset", {})),
                )
                frags = list(ds.get_fragments())
                assert len(frags) == 1
                frag = (
                    _frag_subset(frags[0], row_groups)
                    if row_groups != [None]
                    else frags[0]
                )

                # Extract hive-partition keys, and make sure they
                # are ordered the same as they are in `partitions`
                raw_keys = pa_ds._get_partition_keys(frag.partition_expression)
                partition_keys = [
                    (hive_part.name, raw_keys[hive_part.name])
                    for hive_part in partitions
                ]

        if frag:
            cols = []
            for name in columns:
                if name is None:
                    if "__index_level_0__" in schema.names:
                        columns.append("__index_level_0__")
                else:
                    cols.append(name)

            arrow_table = frag.to_table(
                use_threads=False,
                schema=schema,
                columns=cols,
                filter=_filters_to_expression(filters) if filters else None,
            )
        else:
            arrow_table = _read_table_from_path(
                path_or_frag,
                fs,
                row_groups,
                columns,
                schema,
                filters,
                **kwargs,
            )

        # For pyarrow.dataset api, if we did not read directly from
        # fragments, we need to add the partitioned columns here.
        if partitions and isinstance(partitions, list):
            keys_dict = {k: v for (k, v) in partition_keys}
            for partition in partitions:
                if partition.name not in arrow_table.schema.names:
                    # We read from file paths, so the partition
                    # columns may NOT be in our table yet.
                    cat = keys_dict.get(partition.name, None)
                    if not len(partition.keys):
                        arr = pa.array(np.full(len(arrow_table), cat))
                    else:
                        cat_ind = np.full(
                            len(arrow_table), partition.keys.get_loc(cat), dtype="i4"
                        )
                        arr = pa.DictionaryArray.from_arrays(
                            cat_ind, pa.array(partition.keys)
                        )
                    arrow_table = arrow_table.append_column(partition.name, arr)

        return arrow_table

    @classmethod
    def _determine_type_mapper(
        cls, *, dtype_backend=None, convert_string=False, **kwargs
    ):
        user_mapper = kwargs.get("arrow_to_pandas", {}).get("types_mapper")
        type_mappers = []

        def pyarrow_type_mapper(pyarrow_dtype):
            # Special case pyarrow strings to use more feature complete dtype
            # See https://github.com/pandas-dev/pandas/issues/50074
            if PANDAS_GE_220 and pyarrow_dtype == pa.large_string():
                return pd.StringDtype("pyarrow")
            if pyarrow_dtype == pa.string():
                return pd.StringDtype("pyarrow")
            else:
                return pd.ArrowDtype(pyarrow_dtype)

        # always use the user-defined mapper first, if available
        if user_mapper is not None:
            type_mappers.append(user_mapper)

        # next in priority is converting strings
        if convert_string:
            type_mappers.append({pa.string(): pd.StringDtype("pyarrow")}.get)
            if PANDAS_GE_220:
                type_mappers.append({pa.large_string(): pd.StringDtype("pyarrow")}.get)
            type_mappers.append({pa.date32(): pd.ArrowDtype(pa.date32())}.get)
            type_mappers.append({pa.date64(): pd.ArrowDtype(pa.date64())}.get)

            def _convert_decimal_type(type):
                if pa.types.is_decimal(type):
                    return pd.ArrowDtype(type)
                return None

            type_mappers.append(_convert_decimal_type)

        # and then nullable types
        if dtype_backend == "numpy_nullable":
            type_mappers.append(PYARROW_NULLABLE_DTYPE_MAPPING.get)
        elif dtype_backend == "pyarrow":
            type_mappers.append(pyarrow_type_mapper)

        def default_types_mapper(pyarrow_dtype):
            """Try all type mappers in order, starting from the user type mapper."""
            for type_converter in type_mappers:
                converted_type = type_converter(pyarrow_dtype)
                if converted_type is not None:
                    return converted_type

        if len(type_mappers) > 0:
            return default_types_mapper

    @classmethod
    def _arrow_table_to_pandas(
        cls,
        arrow_table: pa.Table,
        categories,
        dtype_backend=None,
        convert_string=False,
        **kwargs,
    ) -> pd.DataFrame:
        _kwargs = kwargs.get("arrow_to_pandas", {})
        _kwargs.update({"use_threads": False, "ignore_metadata": False})

        types_mapper = cls._determine_type_mapper(
            dtype_backend=dtype_backend,
            convert_string=convert_string,
            **kwargs,
        )
        if types_mapper is not None:
            _kwargs["types_mapper"] = types_mapper

        res = arrow_table.to_pandas(categories=categories, **_kwargs)
        # TODO: remove this when fixed in pyarrow: https://github.com/apache/arrow/issues/34283
        if (
            convert_string
            and isinstance(res.index, pd.Index)
            and not isinstance(res.index, pd.MultiIndex)
            and pd.api.types.is_string_dtype(res.index.dtype)
            and res.index.dtype
            not in (pd.StringDtype("pyarrow"), pd.ArrowDtype(pa.string()))
        ):
            res.index = res.index.astype(pd.StringDtype("pyarrow"))
        return res

    @classmethod
    def collect_file_metadata(cls, path, fs, file_path):
        with fs.open(path, "rb") as f:
            meta = pq.ParquetFile(f).metadata
        if file_path:
            meta.set_file_path(file_path)
        return meta

    @classmethod
    def aggregate_metadata(cls, meta_list, fs, out_path):
        meta = None
        for _meta in meta_list:
            if meta:
                _append_row_groups(meta, _meta)
            else:
                meta = _meta
        if out_path:
            metadata_path = fs.sep.join([out_path, "_metadata"])
            with fs.open(metadata_path, "wb") as fil:
                if not meta:
                    raise ValueError("Cannot write empty metadata!")
                meta.write_metadata_file(fil)
            return None
        else:
            return meta
