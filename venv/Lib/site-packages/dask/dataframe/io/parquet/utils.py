from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from fsspec.core import expand_paths_if_needed, get_fs_token_paths, stringify_path
from fsspec.spec import AbstractFileSystem

from dask import config
from dask.dataframe.io.utils import _is_local_fs
from dask.utils import natural_sort_key, parse_bytes


class Engine:
    """The API necessary to provide a new Parquet reader/writer"""

    @classmethod
    def extract_filesystem(
        cls,
        urlpath,
        filesystem,
        dataset_options,
        open_file_options,
        storage_options,
    ):
        """Extract filesystem object from urlpath or user arguments

        This classmethod should only be overridden for engines that need
        to handle filesystem implementations other than ``fsspec``
        (e.g. ``pyarrow.fs.S3FileSystem``).

        Parameters
        ----------
        urlpath: str or List[str]
            Source directory for data, or path(s) to individual parquet files.
        filesystem: "fsspec" or fsspec.AbstractFileSystem
            Filesystem backend to use. Default is "fsspec"
        dataset_options: dict
            Engine-specific dataset options.
        open_file_options: dict
            Options to be used for file-opening at read time.
        storage_options: dict
            Options to be passed on to the file-system backend.

        Returns
        -------
        fs: Any
            A global filesystem object to be used for metadata
            processing and file-opening by the engine.
        paths: List[str]
            List of data-source paths.
        dataset_options: dict
            Engine-specific dataset options.
        open_file_options: dict
            Options to be used for file-opening at read time.
        """

        # Check if fs was specified as a dataset option
        if filesystem is None:
            fs = dataset_options.pop("fs", "fsspec")
        else:
            if "fs" in dataset_options:
                raise ValueError(
                    "Cannot specify a filesystem argument if the "
                    "'fs' dataset option is also defined."
                )
            fs = filesystem

        if fs in (None, "fsspec"):
            # Use fsspec to infer a filesystem by default
            fs, _, paths = get_fs_token_paths(
                urlpath, mode="rb", storage_options=storage_options
            )
            return fs, paths, dataset_options, open_file_options

        else:
            # Check that an initialized filesystem object was provided
            if not isinstance(fs, AbstractFileSystem):
                raise ValueError(
                    f"Expected fsspec.AbstractFileSystem or 'fsspec'. Got {fs}"
                )

            if storage_options:
                # The filesystem was already specified. Can't pass in
                # any storage options
                raise ValueError(
                    f"Cannot specify storage_options when an explicit "
                    f"filesystem object is specified. Got: {storage_options}"
                )

            if isinstance(urlpath, (list, tuple, set)):
                if not urlpath:
                    raise ValueError("empty urlpath sequence")
                urlpath = [stringify_path(u) for u in urlpath]
            else:
                urlpath = [stringify_path(urlpath)]

            paths = expand_paths_if_needed(urlpath, "rb", 1, fs, None)
            return (
                fs,
                [fs._strip_protocol(u) for u in paths],
                dataset_options,
                open_file_options,
            )

    @classmethod
    def default_blocksize(cls):
        return "256 MiB"

    @classmethod
    def read_partition(
        cls, fs, piece, columns, index, use_nullable_dtypes=False, **kwargs
    ):
        """Read a single piece of a Parquet dataset into a Pandas DataFrame

        This function is called many times in individual tasks

        Parameters
        ----------
        fs: FileSystem
        piece: object
            This is some token that is returned by Engine.read_metadata.
            Typically it represents a row group in a Parquet dataset
        columns: List[str]
            List of column names to pull out of that row group
        index: str, List[str], or False
            The index name(s).
        use_nullable_dtypes: boolean
            Whether to use pandas nullable dtypes (like "string" or "Int64")
            where appropriate when reading parquet files.
        dtype_backend: {"numpy_nullable", "pyarrow"}
            Whether to use pandas nullable dtypes (like "string" or "Int64")
            where appropriate when reading parquet files.
        convert_string: boolean
            Whether to use pyarrow strings when reading parquet files.
        **kwargs:
            Includes `"kwargs"` values stored within the `parts` output
            of `engine.read_metadata`. May also include arguments to be
            passed to the backend (if stored under a top-level `"read"` key).

        Returns
        -------
        A Pandas DataFrame
        """
        raise NotImplementedError()

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
        **kwargs,
    ):
        """Perform engine-specific initialization steps for this dataset

        Parameters
        ----------
        df: dask.dataframe.DataFrame
        fs: FileSystem
        path: str
            Destination directory for data.  Prepend with protocol like ``s3://``
            or ``hdfs://`` for remote data.
        append: bool
            If True, may use existing metadata (if any) and perform checks
            against the new data being stored.
        partition_on: List(str)
            Column(s) to use for dataset partitioning in parquet.
        ignore_divisions: bool
            Whether or not to ignore old divisions when appending.  Otherwise,
            overlapping divisions will lead to an error being raised.
        division_info: dict
            Dictionary containing the divisions and corresponding column name.
        **kwargs: dict
            Other keyword arguments (including `index_cols`)

        Returns
        -------
        tuple:
            engine-specific instance
            list of filenames, one per partition
        """
        raise NotImplementedError

    @classmethod
    def write_partition(
        cls, df, path, fs, filename, partition_on, return_metadata, **kwargs
    ):
        """
        Output a partition of a dask.DataFrame. This will correspond to
        one output file, unless partition_on is set, in which case, it will
        correspond to up to one file in each sub-directory.

        Parameters
        ----------
        df: dask.dataframe.DataFrame
        path: str
            Destination directory for data.  Prepend with protocol like ``s3://``
            or ``hdfs://`` for remote data.
        fs: FileSystem
        filename: str
        partition_on: List(str)
            Column(s) to use for dataset partitioning in parquet.
        return_metadata : bool
            Whether to return list of instances from this write, one for each
            output file. These will be passed to write_metadata if an output
            metadata file is requested.
        **kwargs: dict
            Other keyword arguments (including `fmd` and `index_cols`)

        Returns
        -------
        List of metadata-containing instances (if `return_metadata` is `True`)
        or empty list
        """
        raise NotImplementedError

    @classmethod
    def write_metadata(cls, parts, meta, fs, path, append=False, **kwargs):
        """
        Write the shared metadata file for a parquet dataset.

        Parameters
        ----------
        parts: List
            Contains metadata objects to write, of the type undrestood by the
            specific implementation
        meta: non-chunk metadata
            Details that do not depend on the specifics of each chunk write,
            typically the schema and pandas metadata, in a format the writer
            can use.
        fs: FileSystem
        path: str
            Output file to write to, usually ``"_metadata"`` in the root of
            the output dataset
        append: boolean
            Whether or not to consolidate new metadata with existing (True)
            or start from scratch (False)
        **kwargs: dict
            Other keyword arguments (including `compression`)
        """
        raise NotImplementedError()

    @classmethod
    def collect_file_metadata(cls, path, fs, file_path):
        """
        Collect parquet metadata from a file and set the file_path.

        Parameters
        ----------
        path: str
            Parquet-file path to extract metadata from.
        fs: FileSystem
        file_path: str
            Relative path to set as `file_path` in the metadata.

        Returns
        -------
        A metadata object.  The specific type should be recognized
        by the aggregate_metadata method.
        """
        raise NotImplementedError()

    @classmethod
    def aggregate_metadata(cls, meta_list, fs, out_path):
        """
        Aggregate a list of metadata objects and optionally
        write out the final result as a _metadata file.

        Parameters
        ----------
        meta_list: list
            List of metadata objects to be aggregated into a single
            metadata object, and optionally written to disk. The
            specific element type can be engine specific.
        fs: FileSystem
        out_path: str or None
            Directory to write the final _metadata file. If None
            is specified, the aggregated metadata will be returned,
            and nothing will be written to disk.

        Returns
        -------
        If out_path is None, an aggregate metadata object is returned.
        Otherwise, None is returned.
        """
        raise NotImplementedError()


def _normalize_index_columns(user_columns, data_columns, user_index, data_index):
    """Normalize user and file-provided column and index names

    Parameters
    ----------
    user_columns : None, str or list of str
    data_columns : list of str
    user_index : None, str, or list of str
    data_index : list of str

    Returns
    -------
    column_names : list of str
    index_names : list of str
    """
    specified_columns = user_columns is not None
    specified_index = user_index is not None

    if user_columns is None:
        user_columns = list(data_columns)
    elif isinstance(user_columns, str):
        user_columns = [user_columns]
    else:
        user_columns = list(user_columns)

    if user_index is None:
        user_index = data_index
    elif user_index is False:
        # When index is False, use no index and all fields should be treated as
        # columns (unless `columns` provided).
        user_index = []
        data_columns = data_index + data_columns
    elif isinstance(user_index, str):
        user_index = [user_index]
    else:
        user_index = list(user_index)

    if specified_index and not specified_columns:
        # Only `index` provided. Use specified index, and all column fields
        # that weren't specified as indices
        index_names = user_index
        column_names = [x for x in data_columns if x not in index_names]
    elif specified_columns and not specified_index:
        # Only `columns` provided. Use specified columns, and all index fields
        # that weren't specified as columns
        column_names = user_columns
        index_names = [x for x in data_index if x not in column_names]
    elif specified_index and specified_columns:
        # Both `index` and `columns` provided. Use as specified, but error if
        # they intersect.
        column_names = user_columns
        index_names = user_index
        if set(column_names).intersection(index_names):
            raise ValueError("Specified index and column names must not intersect")
    else:
        # Use default columns and index from the metadata
        column_names = data_columns
        index_names = data_index

    return column_names, index_names


def _sort_and_analyze_paths(file_list, fs, root=False):
    file_list = sorted(file_list, key=natural_sort_key)
    base, fns = _analyze_paths(file_list, fs, root=root)
    return file_list, base, fns


def _analyze_paths(file_list, fs, root=False):
    """Consolidate list of file-paths into parquet relative paths

    Note: This function was mostly copied from dask/fastparquet to
    use in ArrowEngine`."""

    def _join_path(*path):
        def _scrub(i, p):
            # Convert path to standard form
            # this means windows path separators are converted to linux
            p = p.replace(fs.sep, "/")
            if p == "":  # empty path is assumed to be a relative path
                return "."
            if p[-1] == "/":  # trailing slashes are not allowed
                p = p[:-1]
            if i > 0 and p[0] == "/":  # only the first path can start with /
                p = p[1:]
            return p

        abs_prefix = ""
        if path and path[0]:
            if path[0][0] == "/":
                abs_prefix = "/"
                path = list(path)
                path[0] = path[0][1:]
            elif fs.sep == "\\" and path[0][1:].startswith(":/"):
                # If windows, then look for the "c:/" prefix
                abs_prefix = path[0][0:3]
                path = list(path)
                path[0] = path[0][3:]

        _scrubbed = []
        for i, p in enumerate(path):
            _scrubbed.extend(_scrub(i, p).split("/"))
        simpler = []
        for s in _scrubbed:
            if s == ".":
                pass
            elif s == "..":
                if simpler:
                    if simpler[-1] == "..":
                        simpler.append(s)
                    else:
                        simpler.pop()
                elif abs_prefix:
                    raise Exception("can not get parent of root")
                else:
                    simpler.append(s)
            else:
                simpler.append(s)

        if not simpler:
            if abs_prefix:
                joined = abs_prefix
            else:
                joined = "."
        else:
            joined = abs_prefix + ("/".join(simpler))
        return joined

    path_parts_list = [_join_path(fn).split("/") for fn in file_list]
    if root is False:
        basepath = path_parts_list[0][:-1]
        for path_parts in path_parts_list:
            j = len(path_parts) - 1
            for k, (base_part, path_part) in enumerate(zip(basepath, path_parts)):
                if base_part != path_part:
                    j = k
                    break
            basepath = basepath[:j]
        l = len(basepath)
    else:
        basepath = _join_path(root).split("/")
        l = len(basepath)
        assert all(
            p[:l] == basepath for p in path_parts_list
        ), "All paths must begin with the given root"
    out_list = []
    for path_parts in path_parts_list:
        out_list.append(
            "/".join(path_parts[l:])
        )  # use '/'.join() instead of _join_path to be consistent with split('/')

    return (
        "/".join(basepath),
        out_list,
    )  # use '/'.join() instead of _join_path to be consistent with split('/')


def _aggregate_stats(
    file_path,
    file_row_group_stats,
    file_row_group_column_stats,
    stat_col_indices,
):
    """Utility to aggregate the statistics for N row-groups
    into a single dictionary.

    Used by `Engine._construct_parts`
    """
    if len(file_row_group_stats) < 1:
        # Empty statistics
        return {}
    elif len(file_row_group_column_stats) == 0:
        assert len(file_row_group_stats) == 1
        return file_row_group_stats[0]
    else:
        # Note: It would be better to avoid df_rgs and df_cols
        #       construction altogether. It makes it fast to aggregate
        #       the statistics for many row groups, but isn't
        #       worthwhile for a small number of row groups.
        if len(file_row_group_stats) > 1:
            df_rgs = pd.DataFrame(file_row_group_stats)
            s = {
                "file_path_0": file_path,
                "num-rows": df_rgs["num-rows"].sum(),
                "num-row-groups": df_rgs["num-rows"].count(),
                "total_byte_size": df_rgs["total_byte_size"].sum(),
                "columns": [],
            }
        else:
            s = {
                "file_path_0": file_path,
                "num-rows": file_row_group_stats[0]["num-rows"],
                "num-row-groups": 1,
                "total_byte_size": file_row_group_stats[0]["total_byte_size"],
                "columns": [],
            }

        df_cols = None
        if len(file_row_group_column_stats) > 1:
            df_cols = pd.DataFrame(file_row_group_column_stats)
        for ind, name in enumerate(stat_col_indices):
            i = ind * 3
            if df_cols is None:
                minval = file_row_group_column_stats[0][i]
                maxval = file_row_group_column_stats[0][i + 1]
                null_count = file_row_group_column_stats[0][i + 2]
                if minval == maxval and null_count:
                    # Remove "dangerous" stats (min == max, but null values exist)
                    s["columns"].append({"null_count": null_count})
                else:
                    s["columns"].append(
                        {
                            "name": name,
                            "min": minval,
                            "max": maxval,
                            "null_count": null_count,
                        }
                    )
            else:
                minval = df_cols.iloc[:, i].dropna().min()
                maxval = df_cols.iloc[:, i + 1].dropna().max()
                null_count = df_cols.iloc[:, i + 2].sum()
                if minval == maxval and null_count:
                    s["columns"].append({"null_count": null_count})
                else:
                    s["columns"].append(
                        {
                            "name": name,
                            "min": minval,
                            "max": maxval,
                            "null_count": null_count,
                        }
                    )
        return s


def _row_groups_to_parts(
    gather_statistics,
    split_row_groups,
    aggregation_depth,
    file_row_groups,
    file_row_group_stats,
    file_row_group_column_stats,
    stat_col_indices,
    make_part_func,
    make_part_kwargs,
):
    # Construct `parts` and `stats`
    parts = []
    stats = []
    if split_row_groups:
        # Create parts from each file,
        # limiting the number of row_groups in each piece
        split_row_groups = int(split_row_groups)
        residual = 0
        for filename, row_groups in file_row_groups.items():
            row_group_count = len(row_groups)
            if residual:
                _rgs = [0] + list(range(residual, row_group_count, split_row_groups))
            else:
                _rgs = list(range(residual, row_group_count, split_row_groups))

            for i in _rgs:
                i_end = i + split_row_groups
                if aggregation_depth is True:
                    if residual and i == 0:
                        i_end = residual
                        residual = 0
                    _residual = i_end - row_group_count
                    if _residual > 0:
                        residual = _residual

                rg_list = row_groups[i:i_end]

                part = make_part_func(
                    filename,
                    rg_list,
                    **make_part_kwargs,
                )
                if part is None:
                    continue

                parts.append(part)
                if gather_statistics:
                    stat = _aggregate_stats(
                        filename,
                        file_row_group_stats[filename][i:i_end],
                        file_row_group_column_stats[filename][i:i_end],
                        stat_col_indices,
                    )
                    stats.append(stat)
    else:
        for filename, row_groups in file_row_groups.items():
            part = make_part_func(
                filename,
                row_groups,
                **make_part_kwargs,
            )
            if part is None:
                continue

            parts.append(part)
            if gather_statistics:
                stat = _aggregate_stats(
                    filename,
                    file_row_group_stats[filename],
                    file_row_group_column_stats[filename],
                    stat_col_indices,
                )
                stats.append(stat)

    return parts, stats


def _get_aggregation_depth(aggregate_files, partition_names):
    # Use `aggregate_files` to set `aggregation_depth`
    #
    # Note that `partition_names` must be ordered. `True` means that we allow
    # aggregation of any two files. `False` means that we will never aggregate
    # files.  If a string is specified, it must be the name of a partition
    # column, and the "partition depth" of that column will be used for
    # aggregation.  Note that we always convert the string into the partition
    # "depth" to simplify the aggregation logic.

    # Summary of output `aggregation_depth` settings:
    #
    # True  : Free-for-all aggregation (any two files may be aggregated)
    # False : No file aggregation allowed
    # <int> : Allow aggregation within this partition-hierarchy depth

    aggregation_depth = aggregate_files
    if isinstance(aggregate_files, str):
        if aggregate_files in partition_names:
            # aggregate_files corresponds to a partition column. Reset the
            # value of this variable to reflect the partition "depth" (in the
            # range of 1 to the total number of partition levels)
            aggregation_depth = len(partition_names) - partition_names.index(
                aggregate_files
            )
        else:
            raise ValueError(
                f"{aggregate_files} is not a recognized directory partition."
            )

    return aggregation_depth


def _set_metadata_task_size(metadata_task_size, fs):
    # Set metadata_task_size using the config file
    # if the kwarg value was not specified
    if metadata_task_size is None:
        # If a default value is not specified in the config file,
        # otherwise we use "0"
        config_str = "dataframe.parquet.metadata-task-size-" + (
            "local" if _is_local_fs(fs) else "remote"
        )
        return config.get(config_str, 0)

    return metadata_task_size


def _process_open_file_options(
    open_file_options,
    metadata=None,
    columns=None,
    row_groups=None,
    default_engine=None,
    default_cache="readahead",
    allow_precache=True,
):
    # Process `open_file_options`.
    # Set default values and extract `precache_options`
    open_file_options = (open_file_options or {}).copy()
    precache_options = open_file_options.pop("precache_options", {}).copy()
    if not allow_precache:
        # Precaching not allowed
        # (probably because the file system is local)
        precache_options = {}
    if "open_file_func" not in open_file_options:
        if precache_options.get("method", None) == "parquet":
            open_file_options["cache_type"] = open_file_options.get(
                "cache_type", "parts"
            )
            precache_options.update(
                {
                    "metadata": metadata,
                    "columns": columns,
                    "row_groups": row_groups,
                    "engine": precache_options.get("engine", default_engine),
                }
            )
        else:
            open_file_options["cache_type"] = open_file_options.get(
                "cache_type", default_cache
            )
            open_file_options["mode"] = open_file_options.get("mode", "rb")
    return precache_options, open_file_options


def _split_user_options(**kwargs):
    # Check user-defined options.
    # Split into "dataset"-specific kwargs
    user_kwargs = kwargs.copy()

    if "file" in user_kwargs:
        # Deprecation warning to move toward a single `dataset` key
        warnings.warn(
            "Passing user options with the 'file' argument is now deprecated."
            " Please use 'dataset' instead.",
            FutureWarning,
        )

    dataset_options = {
        **user_kwargs.pop("file", {}).copy(),
        **user_kwargs.pop("dataset", {}).copy(),
    }
    read_options = user_kwargs.pop("read", {}).copy()
    open_file_options = user_kwargs.pop("open_file_options", {}).copy()
    return (
        dataset_options,
        read_options,
        open_file_options,
        user_kwargs,
    )


def _set_gather_statistics(
    gather_statistics,
    blocksize,
    split_row_groups,
    aggregation_depth,
    filter_columns,
    stat_columns,
):
    # Use available information about the current read options
    # and target dataset to decide if we need to gather metadata
    # statistics to construct the graph for a `read_parquet` op.

    # If the user has specified `calculate_divisions=True`, then
    # we will be starting with `gather_statistics=True` here.
    if (
        (blocksize and split_row_groups is True)
        or (int(split_row_groups) > 1 and aggregation_depth)
        or filter_columns.intersection(stat_columns)
    ):
        # Need to gather statistics if we are aggregating files
        # or filtering
        # NOTE: Should avoid gathering statistics when the agg
        # does not depend on a row-group statistic
        gather_statistics = True
    elif not stat_columns:
        # Not aggregating files/row-groups.
        # We only need to gather statistics if `stat_columns`
        # is populated
        gather_statistics = False

    return bool(gather_statistics)


def _infer_split_row_groups(row_group_sizes, blocksize, aggregate_files=False):
    # Use blocksize to choose an appropriate split_row_groups value
    if row_group_sizes:
        blocksize = parse_bytes(blocksize)
        if aggregate_files or np.sum(row_group_sizes) > 2 * blocksize:
            # If we are aggregating files, or the file is larger
            # than `blocksize`, set split_row_groups to "adaptive"
            return "adaptive"
    return False
