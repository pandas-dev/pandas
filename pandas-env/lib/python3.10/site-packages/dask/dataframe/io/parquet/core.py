from __future__ import annotations

import math
from typing import Literal

import pandas as pd
import tlz as toolz
from fsspec.core import get_fs_token_paths

from dask.base import tokenize
from dask.dataframe.io.parquet.utils import Engine, _sort_and_analyze_paths
from dask.dataframe.io.utils import DataFrameIOFunction
from dask.dataframe.methods import concat
from dask.delayed import Delayed
from dask.highlevelgraph import HighLevelGraph
from dask.utils import import_required, parse_bytes

NONE_LABEL = "__null_dask_index__"

# ----------------------------------------------------------------------
# User API


class ParquetFunctionWrapper(DataFrameIOFunction):
    """
    Parquet Function-Wrapper Class
    Reads parquet data from disk to produce a partition
    (given a `part` argument).
    """

    def __init__(
        self,
        engine,
        fs,
        meta,
        columns,
        index,
        dtype_backend,
        kwargs,
        common_kwargs,
    ):
        self.engine = engine
        self.fs = fs
        self.meta = meta
        self._columns = columns
        self.index = index
        self.dtype_backend = dtype_backend

        # `kwargs` = user-defined kwargs to be passed
        #            identically for all partitions.
        #
        # `common_kwargs` = kwargs set by engine to be
        #                   passed identically for all
        #                   partitions.
        self.common_kwargs = toolz.merge(common_kwargs, kwargs or {})

    @property
    def columns(self):
        return self._columns

    def project_columns(self, columns):
        """Return a new ParquetFunctionWrapper object
        with a sub-column projection.
        """
        if columns == self.columns:
            return self
        return ParquetFunctionWrapper(
            self.engine,
            self.fs,
            self.meta,
            columns,
            self.index,
            self.dtype_backend,
            None,  # Already merged into common_kwargs
            self.common_kwargs,
        )

    def __call__(self, part):
        if not isinstance(part, list):
            part = [part]

        return read_parquet_part(
            self.fs,
            self.engine,
            self.meta,
            [
                # Temporary workaround for HLG serialization bug
                # (see: https://github.com/dask/dask/issues/8581)
                (
                    (p.data["piece"], p.data.get("kwargs", {}))
                    if hasattr(p, "data")
                    else (p["piece"], p.get("kwargs", {}))
                )
                for p in part
            ],
            self.columns,
            self.index,
            self.common_kwargs,
        )


class ToParquetFunctionWrapper:
    """
    Parquet Function-Wrapper Class

    Writes a DataFrame partition into a distinct parquet
    file. When called, the function also requires the
    current block index (via ``blockwise.BlockIndex``).
    """

    def __init__(
        self,
        engine,
        path,
        fs,
        partition_on,
        write_metadata_file,
        i_offset,
        name_function,
        kwargs_pass,
    ):
        self.engine = engine
        self.path = path
        self.fs = fs
        self.partition_on = partition_on
        self.write_metadata_file = write_metadata_file
        self.i_offset = i_offset
        self.name_function = name_function
        self.kwargs_pass = kwargs_pass

        # NOTE: __name__ must be with "to-parquet"
        # for the name of the resulting `Blockwise`
        # layer to begin with "to-parquet"
        self.__name__ = "to-parquet"

    def __dask_tokenize__(self):
        return (
            self.engine,
            self.path,
            self.fs,
            self.partition_on,
            self.write_metadata_file,
            self.i_offset,
            self.name_function,
            self.kwargs_pass,
        )

    def __call__(self, df, block_index: tuple[int]):
        # Get partition index from block index tuple
        part_i = block_index[0]
        filename = (
            f"part.{part_i + self.i_offset}.parquet"
            if self.name_function is None
            else self.name_function(part_i + self.i_offset)
        )

        # Write out data
        return self.engine.write_partition(
            df,
            self.path,
            self.fs,
            filename,
            self.partition_on,
            self.write_metadata_file,
            **(dict(self.kwargs_pass, head=True) if part_i == 0 else self.kwargs_pass),
        )


def check_multi_support(engine):
    # Helper function to check that the engine
    # supports a multi-partition read
    return hasattr(engine, "multi_support") and engine.multi_support()


def read_parquet_part(fs, engine, meta, part, columns, index, kwargs):
    """Read a part of a parquet dataset

    This function is used by `read_parquet`."""
    if isinstance(part, list):
        if len(part) == 1 or part[0][1] or not check_multi_support(engine):
            # Part kwargs expected
            func = engine.read_partition
            dfs = [
                func(
                    fs,
                    rg,
                    columns.copy(),
                    index,
                    **toolz.merge(kwargs, kw),
                )
                for (rg, kw) in part
            ]
            df = concat(dfs, axis=0) if len(dfs) > 1 else dfs[0]
        else:
            # No part specific kwargs, let engine read
            # list of parts at once
            df = engine.read_partition(
                fs,
                [p[0] for p in part],
                columns.copy(),
                index,
                **kwargs,
            )
    else:
        # NOTE: `kwargs` are the same for all parts, while `part_kwargs` may
        #       be different for each part.
        rg, part_kwargs = part
        df = engine.read_partition(
            fs,
            rg,
            columns,
            index,
            **toolz.merge(kwargs, part_kwargs),
        )

    if meta.columns.name:
        df.columns.name = meta.columns.name
    columns = columns or []
    index = index or []
    df = df[[c for c in columns if c not in index]]
    if index == [NONE_LABEL]:
        df.index.name = None
    return df


def create_metadata_file(
    paths,
    root_dir=None,
    out_dir=None,
    engine="pyarrow",
    storage_options=None,
    split_every=32,
    compute=True,
    compute_kwargs=None,
    fs=None,
):
    """Construct a global _metadata file from a list of parquet files.

    Dask's read_parquet function is designed to leverage a global
    _metadata file whenever one is available.  The to_parquet
    function will generate this file automatically by default, but it
    may not exist if the dataset was generated outside of Dask.  This
    utility provides a mechanism to generate a _metadata file from a
    list of existing parquet files.

    Parameters
    ----------
    paths : list(string)
        List of files to collect footer metadata from.
    root_dir : string, optional
        Root directory of dataset.  The `file_path` fields in the new
        _metadata file will relative to this directory.  If None, a common
        root directory will be inferred.
    out_dir : string or False, optional
        Directory location to write the final _metadata file.  By default,
        this will be set to `root_dir`.  If False is specified, the global
        metadata will be returned as an in-memory object (and will not be
        written to disk).
    engine : str or Engine, default 'pyarrow'
        Parquet Engine to use. Only 'pyarrow' is supported if a string
        is passed.
    storage_options : dict, optional
        Key/value pairs to be passed on to the file-system backend, if any.
    split_every : int, optional
        The final metadata object that is written to _metadata can be much
        smaller than the list of footer metadata. In order to avoid the
        aggregation of all metadata within a single task, a tree reduction
        is used.  This argument specifies the maximum number of metadata
        inputs to be handled by any one task in the tree. Defaults to 32.
    compute : bool, optional
        If True (default) then the result is computed immediately. If False
        then a ``dask.delayed`` object is returned for future computation.
    compute_kwargs : dict, optional
        Options to be passed in to the compute method
    fs : fsspec object, optional
        File-system instance to use for file handling. If prefixes have
        been removed from the elements of ``paths`` before calling this
        function, an ``fs`` argument must be provided to ensure correct
        behavior on remote file systems ("naked" paths cannot be used
        to infer file-system information).
    """
    if isinstance(engine, str):
        engine = get_engine(engine)

    # Process input path list
    if fs is None:
        # Only do this if an fsspec file-system object is not
        # already defined. The prefixes may already be stripped.
        fs, _, paths = get_fs_token_paths(
            paths, mode="rb", storage_options=storage_options
        )
    ap_kwargs = {"root": root_dir} if root_dir else {}
    paths, root_dir, fns = _sort_and_analyze_paths(paths, fs, **ap_kwargs)
    out_dir = root_dir if out_dir is None else out_dir

    # Start constructing a raw graph
    dsk = {}
    name = "gen-metadata-" + tokenize(paths, fs)
    collect_name = "collect-" + name
    agg_name = "agg-" + name

    # Define a "collect" task for each file in the input list.
    # Each tasks will:
    #   1. Extract the footer metadata from a distinct file
    #   2. Populate the `file_path` field in the metadata
    #   3. Return the extracted/modified metadata
    for p, (fn, path) in enumerate(zip(fns, paths)):
        key = (collect_name, p, 0)
        dsk[key] = (engine.collect_file_metadata, path, fs, fn)

    # Build a reduction tree to aggregate all footer metadata
    # into a single metadata object.  Each task in the tree
    # will take in a list of metadata objects as input, and will
    # usually output a single (aggregated) metadata object.
    # The final task in the tree will write the result to disk
    # instead of returning it (this behavior is triggered by
    # passing a file path to `engine.aggregate_metadata`).
    parts = len(paths)
    widths = [parts]
    while parts > 1:
        parts = math.ceil(parts / split_every)
        widths.append(parts)
    height = len(widths)
    for depth in range(1, height):
        for group in range(widths[depth]):
            p_max = widths[depth - 1]
            lstart = split_every * group
            lstop = min(lstart + split_every, p_max)
            dep_task_name = collect_name if depth == 1 else agg_name
            node_list = [(dep_task_name, p, depth - 1) for p in range(lstart, lstop)]
            if depth == height - 1:
                assert group == 0
                dsk[name] = (engine.aggregate_metadata, node_list, fs, out_dir)
            else:
                dsk[(agg_name, group, depth)] = (
                    engine.aggregate_metadata,
                    node_list,
                    None,
                    None,
                )

    # There will be no aggregation tasks if there is only one file
    if len(paths) == 1:
        dsk[name] = (engine.aggregate_metadata, [(collect_name, 0, 0)], fs, out_dir)

    # Convert the raw graph to a `Delayed` object
    graph = HighLevelGraph.from_collections(name, dsk, dependencies=[])
    out = Delayed(name, graph)

    # Optionally compute the result
    if compute:
        if compute_kwargs is None:
            compute_kwargs = dict()
        out = out.compute(**compute_kwargs)
    return out


_ENGINES: dict[str, type[Engine]] = {}


def get_engine(engine: Literal["auto", "pyarrow"] | type[Engine]) -> type[Engine]:
    """Get the parquet engine backend implementation.

    Parameters
    ----------
    engine : {'auto', 'pyarrow'} or Engine subclass
        Parquet library to use. Defaults to 'pyarrow'.

    This can be used to inject third-party engine; e.g. from dask_cudf.
    """
    if isinstance(engine, type) and issubclass(engine, Engine):
        return engine

    if engine in ("auto", "pyarrow"):
        if "pyarrow" in _ENGINES:
            return _ENGINES["pyarrow"]

        try:
            import_required("pyarrow", "`pyarrow` not installed")
        except RuntimeError:
            if engine != "auto":
                raise
        else:
            from dask.dataframe.io.parquet.arrow import ArrowDatasetEngine

            _ENGINES["pyarrow"] = eng = ArrowDatasetEngine
            return eng

    if engine == "auto":
        raise RuntimeError("`pyarrow` not installed")

    raise ValueError(
        f'Unsupported engine: "{engine}". Valid choices are "pyarrow" or "auto".'
    )


#####################
# Utility Functions #
#####################


def sorted_columns(statistics, columns=None):
    """Find sorted columns given row-group statistics

    This finds all columns that are sorted, along with the
    appropriate ``divisions`` for those columns. If the (optional)
    ``columns`` argument is used, the search will be restricted
    to the specified column set.

    Returns
    -------
    out: List of {'name': str, 'divisions': List[str]} dictionaries
    """
    if not statistics:
        return []

    out = []
    for i, c in enumerate(statistics[0]["columns"]):
        if columns and c["name"] not in columns:
            continue
        if not all(
            "min" in s["columns"][i] and "max" in s["columns"][i] for s in statistics
        ):
            continue
        divisions = [c["min"]]
        max = c["max"]
        success = c["min"] is not None
        for stats in statistics[1:]:
            c = stats["columns"][i]
            if c["min"] is None:
                success = False
                break
            if c["min"] >= max:
                divisions.append(c["min"])
                max = c["max"]
            else:
                success = False
                break

        if success:
            divisions.append(max)
            assert divisions == sorted(divisions)
            out.append({"name": c["name"], "divisions": divisions})

    return out


def apply_filters(parts, statistics, filters):
    """Apply filters onto parts/statistics pairs

    Parameters
    ----------
    parts: list
        Tokens corresponding to row groups to read in the future
    statistics: List[dict]
        List of statistics for each part, including min and max values
    filters: Union[List[Tuple[str, str, Any]], List[List[Tuple[str, str, Any]]]]
        List of filters to apply, like ``[[('x', '=', 0), ...], ...]``. This
        implements partition-level (hive) filtering only, i.e., to prevent the
        loading of some row-groups and/or files.

        Predicates can be expressed in disjunctive normal form (DNF). This means
        that the innermost tuple describes a single column predicate. These
        inner predicates are combined with an AND conjunction into a larger
        predicate. The outer-most list then combines all of the combined
        filters with an OR disjunction.

        Predicates can also be expressed as a List[Tuple]. These are evaluated
        as an AND conjunction. To express OR in predicates, one must use the
        (preferred) List[List[Tuple]] notation.

    Returns
    -------
    parts, statistics: the same as the input, but possibly a subset
    """

    # Supported predicate operators
    _supported_operators = {
        "=",
        "==",
        "!=",
        "<",
        "<=",
        ">",
        ">=",
        "is",
        "is not",
        "in",
        "not in",
    }

    def apply_conjunction(parts, statistics, conjunction):
        for column, operator, value in conjunction:
            if operator not in _supported_operators:
                # Use same error message as `_filters_to_expression`
                raise ValueError(
                    f'"{(column, operator, value)}" is not a valid operator in predicates.'
                )
            elif operator in ("in", "not in") and not isinstance(
                value, (list, set, tuple)
            ):
                raise TypeError("Value of 'in' filter must be a list, set, or tuple.")
            out_parts = []
            out_statistics = []
            for part, stats in zip(parts, statistics):
                if "filter" in stats and stats["filter"]:
                    continue  # Filtered by engine
                try:
                    c = toolz.groupby("name", stats["columns"])[column][0]
                    min = c["min"]
                    max = c["max"]
                    null_count = c.get("null_count", None)
                except KeyError:
                    out_parts.append(part)
                    out_statistics.append(stats)
                else:
                    if (
                        # Must allow row-groups with "missing" stats
                        (min is None and max is None and not null_count)
                        # Check "is" and "is not" filters first
                        or operator == "is"
                        and null_count
                        or operator == "is not"
                        and (not pd.isna(min) or not pd.isna(max))
                        # Allow all-null row-groups if not filtering out nulls
                        or operator != "is not"
                        and min is None
                        and max is None
                        and null_count
                        # Start conventional (non-null) filtering
                        # (main/max cannot be None for remaining checks)
                        or operator in ("==", "=")
                        and min <= value <= max
                        or operator == "!="
                        and (null_count or min != value or max != value)
                        or operator == "<"
                        and min < value
                        or operator == "<="
                        and min <= value
                        or operator == ">"
                        and max > value
                        or operator == ">="
                        and max >= value
                        or operator == "in"
                        and any(min <= item <= max for item in value)
                        or operator == "not in"
                        and not any(min == max == item for item in value)
                    ):
                        out_parts.append(part)
                        out_statistics.append(stats)

            parts, statistics = out_parts, out_statistics

        return parts, statistics

    conjunction, *disjunction = filters if isinstance(filters[0], list) else [filters]

    out_parts, out_statistics = apply_conjunction(parts, statistics, conjunction)
    for conjunction in disjunction:
        for part, stats in zip(*apply_conjunction(parts, statistics, conjunction)):
            if part not in out_parts:
                out_parts.append(part)
                out_statistics.append(stats)

    return out_parts, out_statistics


def set_index_columns(meta, index, columns, auto_index_allowed):
    """Handle index/column arguments, and modify `meta`
    Used in read_parquet.
    """
    ignore_index_column_intersection = False
    if columns is None:
        # User didn't specify columns, so ignore any intersection
        # of auto-detected values with the index (if necessary)
        ignore_index_column_intersection = True
        # Do not allow "un-named" fields to be read in as columns.
        # These were intended to be un-named indices at write time.
        _index = index or []
        columns = [
            c for c in meta.columns if c not in (None, NONE_LABEL) or c in _index
        ]

    if not set(columns).issubset(set(meta.columns)):
        raise ValueError(
            "The following columns were not found in the dataset %s\n"
            "The following columns were found %s"
            % (set(columns) - set(meta.columns), meta.columns)
        )

    if index:
        if isinstance(index, str):
            index = [index]
        if isinstance(columns, str):
            columns = [columns]

        if ignore_index_column_intersection:
            columns = [col for col in columns if col not in index]
        if set(index).intersection(columns):
            if auto_index_allowed:
                raise ValueError(
                    "Specified index and column arguments must not intersect"
                    " (set index=False or remove the detected index from columns).\n"
                    "index: {} | column: {}".format(index, columns)
                )
            else:
                raise ValueError(
                    "Specified index and column arguments must not intersect.\n"
                    "index: {} | column: {}".format(index, columns)
                )

    return meta[list(columns)], index, columns


def aggregate_row_groups(
    parts, stats, blocksize, split_row_groups, fs, aggregation_depth
):
    if not stats or not stats[0].get("file_path_0", None):
        return parts, stats

    parts_agg = []
    stats_agg = []

    use_row_group_criteria = split_row_groups and int(split_row_groups) > 1
    use_blocksize_criteria = bool(blocksize)
    if use_blocksize_criteria:
        blocksize = parse_bytes(blocksize)
    next_part, next_stat = [parts[0].copy()], stats[0].copy()
    for i in range(1, len(parts)):
        stat, part = stats[i], parts[i]

        # Criteria #1 for aggregating parts: parts are within the same file
        same_path = stat["file_path_0"] == next_stat["file_path_0"]
        multi_path_allowed = False

        if aggregation_depth:
            # Criteria #2 for aggregating parts: The part does not include
            # row-group information, or both parts include the same kind
            # of row_group aggregation (all None, or all indices)
            multi_path_allowed = len(part["piece"]) == 1
            if not (same_path or multi_path_allowed):
                rgs = set(list(part["piece"][1]) + list(next_part[-1]["piece"][1]))
                multi_path_allowed = (rgs == {None}) or (None not in rgs)

            # Criteria #3 for aggregating parts: The parts share a
            # directory at the "depth" allowed by `aggregation_depth`
            if not same_path and multi_path_allowed:
                if aggregation_depth is True:
                    multi_path_allowed = True
                elif isinstance(aggregation_depth, int):
                    # Make sure files share the same directory
                    root = stat["file_path_0"].split(fs.sep)[:-aggregation_depth]
                    next_root = next_stat["file_path_0"].split(fs.sep)[
                        :-aggregation_depth
                    ]
                    multi_path_allowed = root == next_root
                else:
                    raise ValueError(
                        f"{aggregation_depth} not supported for `aggregation_depth`"
                    )

        def _check_row_group_criteria(stat, next_stat):
            if use_row_group_criteria:
                return (next_stat["num-row-groups"] + stat["num-row-groups"]) <= int(
                    split_row_groups
                )
            else:
                return False

        def _check_blocksize_criteria(stat, next_stat):
            if use_blocksize_criteria:
                return (
                    next_stat["total_byte_size"] + stat["total_byte_size"]
                ) <= blocksize
            else:
                return False

        stat["num-row-groups"] = stat.get("num-row-groups", 1)
        next_stat["num-row-groups"] = next_stat.get("num-row-groups", 1)

        if (same_path or multi_path_allowed) and (
            _check_row_group_criteria(stat, next_stat)
            or _check_blocksize_criteria(stat, next_stat)
        ):
            # Update part list
            next_piece = next_part[-1]["piece"]
            this_piece = part["piece"]
            if (
                same_path
                and len(next_piece) > 1
                and next_piece[1] != [None]
                and this_piece[1] != [None]
            ):
                next_piece[1].extend(this_piece[1])
            else:
                next_part.append(part)

            # Update Statistics
            next_stat["total_byte_size"] += stat["total_byte_size"]
            next_stat["num-rows"] += stat["num-rows"]
            next_stat["num-row-groups"] += stat["num-row-groups"]
            for col, col_add in zip(next_stat["columns"], stat["columns"]):
                if col["name"] != col_add["name"]:
                    raise ValueError("Columns are different!!")
                if "min" in col:
                    col["min"] = min(col["min"], col_add["min"])
                if "max" in col:
                    col["max"] = max(col["max"], col_add["max"])
        else:
            parts_agg.append(next_part)
            stats_agg.append(next_stat)
            next_part, next_stat = [part.copy()], stat.copy()

    parts_agg.append(next_part)
    stats_agg.append(next_stat)

    return parts_agg, stats_agg
