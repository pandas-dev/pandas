from __future__ import annotations

import copy
import pickle
import threading
import warnings
from collections import OrderedDict, defaultdict
from contextlib import ExitStack

import numpy as np
import pandas as pd
import tlz as toolz
from packaging.version import Version

from dask.core import flatten
from dask.dataframe._compat import PANDAS_GE_201

try:
    import fastparquet
    from fastparquet import ParquetFile
    from fastparquet.util import ex_from_sep, get_file_scheme, groupby_types, val_to_num
    from fastparquet.writer import make_part_file, partition_on_columns
except ImportError:
    pass

from dask.base import tokenize

#########################
# Fastparquet interface #
#########################
from dask.dataframe.io.parquet.utils import (
    Engine,
    _get_aggregation_depth,
    _infer_split_row_groups,
    _normalize_index_columns,
    _parse_pandas_metadata,
    _process_open_file_options,
    _row_groups_to_parts,
    _set_gather_statistics,
    _set_metadata_task_size,
    _sort_and_analyze_paths,
)
from dask.dataframe.io.utils import _is_local_fs, _meta_from_dtypes, _open_input_files
from dask.dataframe.utils import UNKNOWN_CATEGORIES
from dask.delayed import Delayed
from dask.utils import natural_sort_key

# Thread lock required to reset row-groups
_FP_FILE_LOCK = threading.RLock()


def _paths_to_cats(paths, file_scheme):
    """
    Extract categorical fields and labels from hive- or drill-style paths.
    FixMe: This has been pasted from https://github.com/dask/fastparquet/pull/471
    Use fastparquet.api.paths_to_cats from fastparquet>0.3.2 instead.

    Parameters
    ----------
    paths (Iterable[str]): file paths relative to root
    file_scheme (str):

    Returns
    -------
    cats (OrderedDict[str, List[Any]]): a dict of field names and their values
    """
    if file_scheme in ["simple", "flat", "other"]:
        cats = {}
        return cats

    cats = OrderedDict()
    raw_cats = OrderedDict()
    s = ex_from_sep("/")
    paths = toolz.unique(paths)
    if file_scheme == "hive":
        partitions = toolz.unique((k, v) for path in paths for k, v in s.findall(path))
        for key, val in partitions:
            cats.setdefault(key, set()).add(val_to_num(val))
            raw_cats.setdefault(key, set()).add(val)
    else:
        i_val = toolz.unique(
            (i, val) for path in paths for i, val in enumerate(path.split("/")[:-1])
        )
        for i, val in i_val:
            key = "dir%i" % i
            cats.setdefault(key, set()).add(val_to_num(val))
            raw_cats.setdefault(key, set()).add(val)

    for key, v in cats.items():
        # Check that no partition names map to the same value after transformation by val_to_num
        raw = raw_cats[key]
        if len(v) != len(raw):
            conflicts_by_value = OrderedDict()
            for raw_val in raw_cats[key]:
                conflicts_by_value.setdefault(val_to_num(raw_val), set()).add(raw_val)
            conflicts = [
                c for k in conflicts_by_value.values() if len(k) > 1 for c in k
            ]
            raise ValueError("Partition names map to the same value: %s" % conflicts)
        vals_by_type = groupby_types(v)

        # Check that all partition names map to the same type after transformation by val_to_num
        if len(vals_by_type) > 1:
            examples = [x[0] for x in vals_by_type.values()]
            warnings.warn(
                "Partition names coerce to values of different types, e.g. %s"
                % examples
            )

    cats = OrderedDict([(key, list(v)) for key, v in cats.items()])
    return cats


paths_to_cats = (
    _paths_to_cats  # FixMe: use fastparquet.api.paths_to_cats for fastparquet>0.3.2
)


class FastParquetEngine(Engine):
    @classmethod
    def _organize_row_groups(
        cls,
        pf,
        split_row_groups,
        gather_statistics,
        stat_col_indices,
        filters,
        dtypes,
        base_path,
        has_metadata_file,
        blocksize,
        aggregation_depth,
    ):
        """Organize row-groups by file."""

        # Get partitioning metadata
        pqpartitions = list(pf.cats)

        # Fastparquet does not use a natural sorting
        # order for partitioned data. Re-sort by path
        if (
            pqpartitions
            and aggregation_depth
            and pf.row_groups
            and pf.row_groups[0].columns[0].file_path
        ):
            pf.row_groups = sorted(
                pf.row_groups,
                key=lambda x: natural_sort_key(x.columns[0].file_path),
            )

        # Store types specified in pandas metadata
        pandas_type = {}
        if pf.row_groups and pf.pandas_metadata:
            for c in pf.pandas_metadata.get("columns", []):
                if "field_name" in c:
                    pandas_type[c["field_name"]] = c.get("pandas_type", None)

        # Get the number of row groups per file
        single_rg_parts = int(split_row_groups) == 1
        file_row_groups = defaultdict(list)
        file_row_group_stats = defaultdict(list)
        file_row_group_column_stats = defaultdict(list)
        cmax_last = {}
        for rg, row_group in enumerate(pf.row_groups):
            # We can filter partition columns here without dealing
            # with statistics
            if (
                pqpartitions
                and filters
                and fastparquet.api.filter_out_cats(row_group, filters)
            ):
                continue

            # NOTE: Here we assume that all column chunks are stored
            # in the same file. This is not strictly required by the
            # parquet spec.
            fp = row_group.columns[0].file_path
            fpath = fp.decode() if isinstance(fp, bytes) else fp
            if fpath is None:
                if not has_metadata_file:
                    # There doesn't need to be a file_path if the
                    # row group is in the same file as the metadata.
                    # Assume this is a single-file dataset.
                    fpath = pf.fn
                    base_path = base_path or ""
                else:
                    raise ValueError(
                        "Global metadata structure is missing a file_path string. "
                        "If the dataset includes a _metadata file, that file may "
                        "have one or more missing file_path fields."
                    )

            # Append a tuple to file_row_groups. This tuple will
            # be structured as: `(<local-row-group-id>, <global-row-group-id>)`
            if file_row_groups[fpath]:
                file_row_groups[fpath].append((file_row_groups[fpath][-1][0] + 1, rg))
            else:
                file_row_groups[fpath].append((0, rg))

            if gather_statistics:
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
                for name, i in stat_col_indices.items():
                    column = row_group.columns[i]
                    if column.meta_data.statistics:
                        cmin = None
                        cmax = None
                        null_count = None
                        # TODO: Avoid use of `pf.statistics`
                        if pf.statistics["min"][name][0] is not None:
                            cmin = pf.statistics["min"][name][rg]
                            cmax = pf.statistics["max"][name][rg]
                            null_count = pf.statistics["null_count"][name][rg]
                        elif dtypes[name] == "object":
                            cmin = column.meta_data.statistics.min_value
                            cmax = column.meta_data.statistics.max_value
                            null_count = column.meta_data.statistics.null_count
                            # Older versions may not have cmin/cmax_value
                            if cmin is None:
                                cmin = column.meta_data.statistics.min
                            if cmax is None:
                                cmax = column.meta_data.statistics.max
                            # Decode bytes as long as "bytes" is not the
                            # expected `pandas_type` for this column
                            if (
                                isinstance(cmin, (bytes, bytearray))
                                and pandas_type.get(name, None) != "bytes"
                            ):
                                cmin = cmin.decode("utf-8")
                                cmax = cmax.decode("utf-8")
                            if isinstance(null_count, (bytes, bytearray)):
                                null_count = null_count.decode("utf-8")
                        if isinstance(cmin, np.datetime64):
                            tz = getattr(dtypes[name], "tz", None)
                            cmin = pd.Timestamp(cmin, tz=tz)
                            cmax = pd.Timestamp(cmax, tz=tz)
                        last = cmax_last.get(name, None)

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
                        if (
                            not (
                                filters
                                or (blocksize and split_row_groups is True)
                                or aggregation_depth
                            )
                            and column.meta_data.num_values > 0
                        ):
                            # We are collecting statistics for divisions
                            # only (no filters) - Lets bail.
                            gather_statistics = False
                            file_row_group_stats = {}
                            file_row_group_column_stats = {}
                            break

                        if single_rg_parts:
                            s["columns"].append({"name": name})
                        else:
                            cstats += [None, None, None]
                if gather_statistics:
                    file_row_group_stats[fpath].append(s)
                    if not single_rg_parts:
                        file_row_group_column_stats[fpath].append(tuple(cstats))

        return (
            file_row_groups,
            file_row_group_stats,
            file_row_group_column_stats,
            gather_statistics,
            base_path,
        )

    @classmethod
    def _get_thrift_row_groups(
        cls,
        pf,
        filename,
        row_groups,
    ):
        """Turn a set of row-groups into bytes-serialized form
        using thrift via pickle.
        """

        real_row_groups = []
        for _, rg_global in row_groups:
            row_group = pf.row_groups[rg_global]
            columns = row_group.columns
            for c, col in enumerate(columns):
                if c:
                    col.file_path = None
                md = col.meta_data
                md.key_value_metadata = None
                # NOTE: Fastparquet may need the null count in the
                # statistics, so we cannot just set statistics
                # to none.  Set attributes separately:
                st = md.statistics
                if st:
                    st.distinct_count = None
                    st.max = None
                    st.min = None
                    st.max_value = None
                    st.min_value = None
                md.encodings = None
                md.total_uncompressed_size = None
                md.encoding_stats = None
            row_group.columns = columns
            real_row_groups.append(row_group)
        return real_row_groups

    @classmethod
    def _make_part(
        cls,
        filename,
        rg_list,
        fs=None,
        pf=None,
        base_path=None,
        partitions=None,
    ):
        """Generate a partition-specific element of `parts`."""

        if partitions:
            real_row_groups = cls._get_thrift_row_groups(
                pf,
                filename,
                rg_list,
            )
            part = {"piece": (real_row_groups,)}
        else:
            # Get full path (empty strings should be ignored)
            full_path = fs.sep.join([p for p in [base_path, filename] if p != ""])
            row_groups = [rg[0] for rg in rg_list]  # Don't need global IDs
            part = {"piece": (full_path, row_groups)}

        return part

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
        # Define the parquet-file (pf) object to use for metadata,
        # Also, initialize `parts`.  If `parts` is populated here,
        # then each part will correspond to a file.  Otherwise, each part will
        # correspond to a row group (populated later).

        # Extract dataset-specific options
        dataset_kwargs = kwargs.pop("dataset", {})

        parts = []
        _metadata_exists = False
        if len(paths) == 1 and fs.isdir(paths[0]):
            # This is a directory.
            # Check if _metadata and/or _common_metadata files exists
            base = paths[0]
            _metadata_exists = True
            if not ignore_metadata_file:
                _metadata_exists = fs.isfile(fs.sep.join([base, "_metadata"]))

            # Find all files if we are not using a _metadata file
            if ignore_metadata_file or not _metadata_exists:
                # For now, we need to discover every file under paths[0]
                paths, base, fns = _sort_and_analyze_paths(fs.find(base), fs, root=base)
                _update_paths = False
                for fn in ["_metadata", "_common_metadata"]:
                    try:
                        fns.remove(fn)
                        _update_paths = True
                    except ValueError:
                        pass
                if _update_paths:
                    paths = [fs.sep.join([base, fn]) for fn in fns]
                _metadata_exists = False
            if _metadata_exists:
                # Using _metadata file (best-case scenario)
                pf = ParquetFile(
                    fs.sep.join([base, "_metadata"]),
                    open_with=fs.open,
                    **dataset_kwargs,
                )
            else:
                # Use 0th file
                # Note that "_common_metadata" can cause issues for
                # partitioned datasets.
                if parquet_file_extension:
                    # Raise error if all files have been filtered by extension
                    len0 = len(paths)
                    paths = [
                        path for path in paths if path.endswith(parquet_file_extension)
                    ]
                    fns = [fn for fn in fns if fn.endswith(parquet_file_extension)]
                    if len0 and paths == []:
                        raise ValueError(
                            "No files satisfy the `parquet_file_extension` criteria "
                            f"(files must end with {parquet_file_extension})."
                        )
                pf = ParquetFile(
                    paths[:1], open_with=fs.open, root=base, **dataset_kwargs
                )
                scheme = get_file_scheme(fns)
                pf.file_scheme = scheme
                pf.cats = paths_to_cats(fns, scheme)
                if not gather_statistics:
                    parts = [fs.sep.join([base, fn]) for fn in fns]
        else:
            # This is a list of files
            paths, base, fns = _sort_and_analyze_paths(paths, fs)

            # Check if _metadata is in paths, and
            # remove it if ignore_metadata_file=True
            _metadata_exists = "_metadata" in fns
            if _metadata_exists and ignore_metadata_file:
                fns.remove("_metadata")
                _metadata_exists = False
            paths = [fs.sep.join([base, fn]) for fn in fns]

            if _metadata_exists:
                # We have a _metadata file, lets use it
                pf = ParquetFile(
                    fs.sep.join([base, "_metadata"]),
                    open_with=fs.open,
                    **dataset_kwargs,
                )
            else:
                # Rely on metadata for 0th file.
                # Will need to pass a list of paths to read_partition
                scheme = get_file_scheme(fns)
                pf = ParquetFile(
                    paths[:1], open_with=fs.open, root=base, **dataset_kwargs
                )
                pf.file_scheme = scheme
                pf.cats = paths_to_cats(fns, scheme)
                if not gather_statistics:
                    parts = paths.copy()

        # Check the `aggregate_files` setting
        aggregation_depth = _get_aggregation_depth(
            aggregate_files,
            list(pf.cats),
        )

        # Ensure that there is no overlap between partition columns
        # and explicit columns in `pf`
        if pf.cats:
            _partitions = [p for p in pf.cats if p not in pf.columns]
            if not _partitions:
                pf.cats = {}
            elif len(_partitions) != len(pf.cats):
                raise ValueError(
                    "No partition-columns should be written in the \n"
                    "file unless they are ALL written in the file.\n"
                    "This restriction is removed as of fastparquet 0.8.4\n"
                    "columns: {} | partitions: {}".format(pf.columns, pf.cats.keys())
                )

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
                pf_sample = ParquetFile(
                    paths[0],
                    open_with=fs.open,
                    **dataset_kwargs,
                )
                split_row_groups = _infer_split_row_groups(
                    [rg.total_byte_size for rg in pf_sample.row_groups],
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

        return {
            "pf": pf,
            "paths": paths,
            "has_metadata_file": _metadata_exists,
            "parts": parts,
            "base": base,
            "fs": fs,
            "gather_statistics": gather_statistics,
            "categories": categories,
            "index": index,
            "filters": filters,
            "split_row_groups": split_row_groups,
            "blocksize": blocksize,
            "aggregate_files": aggregate_files,
            "aggregation_depth": aggregation_depth,
            "metadata_task_size": metadata_task_size,
            "kwargs": {
                "dataset": dataset_kwargs,
                **kwargs,
            },
        }

    @classmethod
    def _create_dd_meta(cls, dataset_info):
        # Collect necessary information from dataset_info
        pf = dataset_info["pf"]
        index = dataset_info["index"]
        categories = dataset_info["categories"]

        columns = None
        pandas_md = pf.pandas_metadata

        if pandas_md:
            (
                index_names,
                column_names,
                storage_name_mapping,
                column_index_names,
            ) = _parse_pandas_metadata(pandas_md)
            #  auto-ranges should not be created by fastparquet
            column_names.extend(pf.cats)
        else:
            index_names = []
            column_names = pf.columns + list(pf.cats)
            storage_name_mapping = {k: k for k in column_names}
            column_index_names = [None]
        if index is None and len(index_names) > 0:
            if len(index_names) == 1 and index_names[0] is not None:
                index = index_names[0]
            else:
                index = index_names

        # Normalize user inputs
        column_names, index_names = _normalize_index_columns(
            columns, column_names, index, index_names
        )

        all_columns = index_names + column_names

        categories_dict = None
        if isinstance(categories, dict):
            categories_dict = categories

        if categories is None:
            categories = pf.categories
        elif isinstance(categories, str):
            categories = [categories]
        else:
            categories = list(categories)

        # Check that categories are included in columns
        if categories and not set(categories).intersection(all_columns):
            raise ValueError(
                "categories not in available columns.\n"
                "categories: {} | columns: {}".format(categories, list(all_columns))
            )

        dtypes = pf._dtypes(categories)
        dtypes = {storage_name_mapping.get(k, k): v for k, v in dtypes.items()}

        index_cols = index or ()
        if isinstance(index_cols, str):
            index_cols = [index_cols]
        for ind in index_cols:
            if getattr(dtypes.get(ind), "numpy_dtype", None):
                # index does not support masked types
                dtypes[ind] = dtypes[ind].numpy_dtype
        for cat in categories:
            if cat in all_columns:
                dtypes[cat] = pd.CategoricalDtype(categories=[UNKNOWN_CATEGORIES])

        for catcol in pf.cats:
            if catcol in all_columns:
                dtypes[catcol] = pd.CategoricalDtype(categories=pf.cats[catcol])

        meta = _meta_from_dtypes(all_columns, dtypes, index_cols, column_index_names)

        # Update `dataset_info` and return `meta`
        dataset_info["dtypes"] = dtypes
        dataset_info["index"] = index
        dataset_info["index_cols"] = index_cols
        dataset_info["categories"] = categories
        dataset_info["categories_dict"] = categories_dict

        return meta

    @classmethod
    def _construct_collection_plan(cls, dataset_info):
        # Collect necessary information from dataset_info
        fs = dataset_info["fs"]
        parts = dataset_info["parts"]
        paths = dataset_info["paths"]
        filters = dataset_info["filters"]
        pf = dataset_info["pf"]
        split_row_groups = dataset_info["split_row_groups"]
        blocksize = dataset_info["blocksize"]
        gather_statistics = dataset_info["gather_statistics"]
        base_path = dataset_info["base"]
        aggregation_depth = dataset_info["aggregation_depth"]
        index_cols = dataset_info["index_cols"]
        categories = dataset_info["categories"]
        dtypes = dataset_info["dtypes"]
        categories_dict = dataset_info["categories_dict"]
        has_metadata_file = dataset_info["has_metadata_file"]
        metadata_task_size = dataset_info["metadata_task_size"]
        kwargs = dataset_info["kwargs"]

        # Ensure metadata_task_size is set
        # (Using config file or defaults)
        metadata_task_size = _set_metadata_task_size(
            dataset_info["metadata_task_size"], fs
        )

        # Determine which columns need statistics.
        # At this point, gather_statistics is only True if
        # the user specified calculate_divisions=True
        filter_columns = {t[0] for t in flatten(filters or [], container=list)}
        stat_col_indices = {}
        _index_cols = index_cols if (gather_statistics and len(index_cols) == 1) else []
        for i, name in enumerate(pf.columns):
            if name in _index_cols or name in filter_columns:
                stat_col_indices[name] = i

        # Decide final `gather_statistics` setting.
        # NOTE: The "fastparquet" engine requires statistics for
        # filtering even if the filter is on a partitioned column
        gather_statistics = _set_gather_statistics(
            gather_statistics,
            blocksize,
            split_row_groups,
            aggregation_depth,
            filter_columns,
            set(stat_col_indices) | filter_columns,
        )

        # Define common_kwargs
        common_kwargs = {
            "categories": categories_dict or categories,
            "root_cats": pf.cats,
            "root_file_scheme": pf.file_scheme,
            "base_path": base_path,
            **kwargs,
        }

        # Check if this is a very simple case where we can just
        # return the path names. This requires that `parts`
        # already be a list of paths. Also, we cannot be splitting
        # by row-group or collecting statistics.
        if (
            gather_statistics is False
            and not split_row_groups
            and isinstance(parts, list)
            and len(parts)
            and isinstance(parts[0], str)
        ):
            return (
                [{"piece": (full_path, None)} for full_path in parts],
                [],
                common_kwargs,
            )

        dataset_info_kwargs = {
            "fs": fs,
            "split_row_groups": split_row_groups,
            "gather_statistics": gather_statistics,
            "filters": filters,
            "dtypes": dtypes,
            "stat_col_indices": stat_col_indices,
            "aggregation_depth": aggregation_depth,
            "blocksize": blocksize,
            "root_cats": pf.cats,
            "root_file_scheme": pf.file_scheme,
            "base_path": "" if base_path is None else base_path,
            "has_metadata_file": has_metadata_file,
        }

        if (
            has_metadata_file
            or metadata_task_size == 0
            or metadata_task_size > len(paths)
        ):
            # Construct the output-partitioning plan on the
            # client process (in serial).  This means we have
            # a global _metadata file, or that `metadata_task_size`
            # is zero or larger than the number of files.
            pf_or_paths = pf if has_metadata_file else paths
            parts, stats = cls._collect_file_parts(pf_or_paths, dataset_info_kwargs)

        else:
            # We DON'T have a global _metadata file to work with.
            # We should loop over files in parallel
            parts, stats = [], []
            if paths:
                # Build and compute a task graph to construct stats/parts
                gather_parts_dsk = {}
                name = "gather-pq-parts-" + tokenize(paths, dataset_info_kwargs)
                finalize_list = []
                for task_i, file_i in enumerate(
                    range(0, len(paths), metadata_task_size)
                ):
                    finalize_list.append((name, task_i))
                    gather_parts_dsk[finalize_list[-1]] = (
                        cls._collect_file_parts,
                        paths[file_i : file_i + metadata_task_size],
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
        pf_or_files,
        dataset_info_kwargs,
    ):
        # Collect necessary information from dataset_info
        fs = dataset_info_kwargs["fs"]
        split_row_groups = dataset_info_kwargs["split_row_groups"]
        gather_statistics = dataset_info_kwargs["gather_statistics"]
        stat_col_indices = dataset_info_kwargs["stat_col_indices"]
        filters = dataset_info_kwargs["filters"]
        dtypes = dataset_info_kwargs["dtypes"]
        blocksize = dataset_info_kwargs["blocksize"]
        aggregation_depth = dataset_info_kwargs["aggregation_depth"]
        base_path = dataset_info_kwargs.get("base_path", None)
        root_cats = dataset_info_kwargs.get("root_cats", None)
        root_file_scheme = dataset_info_kwargs.get("root_file_scheme", None)
        has_metadata_file = dataset_info_kwargs["has_metadata_file"]

        # Get ParquetFile
        if not isinstance(pf_or_files, fastparquet.api.ParquetFile):
            # Construct local `ParquetFile` object
            pf = ParquetFile(
                pf_or_files,
                open_with=fs.open,
                root=base_path,
            )
            # Update hive-partitioning to match global cats/scheme
            pf.cats = root_cats or {}
            if root_cats:
                pf.file_scheme = root_file_scheme
        else:
            # We already have a ParquetFile object to work with
            pf = pf_or_files

        # Organize row-groups by file
        (
            file_row_groups,
            file_row_group_stats,
            file_row_group_column_stats,
            gather_statistics,
            base_path,
        ) = cls._organize_row_groups(
            pf,
            split_row_groups,
            gather_statistics,
            stat_col_indices,
            filters,
            dtypes,
            base_path,
            has_metadata_file,
            blocksize,
            aggregation_depth,
        )

        # Convert organized row-groups to parts
        parts, stats = _row_groups_to_parts(
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
                "pf": pf,
                "base_path": base_path,
                "partitions": list(pf.cats),
            },
        )

        return parts, stats

    @classmethod
    def read_metadata(
        cls,
        fs,
        paths,
        categories=None,
        index=None,
        use_nullable_dtypes=None,
        dtype_backend=None,
        gather_statistics=None,
        filters=None,
        split_row_groups="adaptive",
        blocksize=None,
        aggregate_files=None,
        ignore_metadata_file=False,
        metadata_task_size=None,
        parquet_file_extension=None,
        **kwargs,
    ):
        if use_nullable_dtypes is not None:
            raise ValueError(
                "`use_nullable_dtypes` is not supported by the fastparquet engine"
            )
        if dtype_backend is not None:
            raise ValueError(
                "`dtype_backend` is not supported by the fastparquet engine"
            )

        # Stage 1: Collect general dataset information
        dataset_info = cls._collect_dataset_info(
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
        )

        # Stage 2: Generate output `meta`
        meta = cls._create_dd_meta(dataset_info)

        # Stage 3: Generate parts and stats
        parts, stats, common_kwargs = cls._construct_collection_plan(dataset_info)

        # Cannot allow `None` in columns if the user has specified index=False
        index = dataset_info["index"]
        if index is False and None in meta.columns:
            meta.drop(columns=[None], inplace=True)

        # Add `common_kwargs` to the first element of `parts`.
        # We can return as a separate element in the future, but
        # should avoid breaking the API for now.
        if len(parts):
            parts[0]["common_kwargs"] = common_kwargs
            parts[0]["aggregation_depth"] = dataset_info["aggregation_depth"]
            parts[0]["split_row_groups"] = dataset_info["split_row_groups"]

        if len(parts) and len(parts[0]["piece"]) == 1:
            # Strip all partition-dependent or unnecessary
            # data from the `ParquetFile` object
            pf = dataset_info["pf"]
            pf.row_groups = None
            pf.fmd.row_groups = None
            pf._statistics = None
            parts[0]["common_kwargs"]["parquet_file"] = pf

        return (meta, stats, parts, index)

    @classmethod
    def multi_support(cls):
        return cls == FastParquetEngine

    @classmethod
    def read_partition(
        cls,
        fs,
        pieces,
        columns,
        index,
        dtype_backend=None,
        categories=(),
        root_cats=None,
        root_file_scheme=None,
        base_path=None,
        **kwargs,
    ):
        null_index_name = False
        base_path = False if not root_cats else base_path
        if isinstance(index, list):
            if index == [None]:
                # Handling a None-labeled index...
                # The pandas metadata told us to read in an index
                # labeled `None`. If this corresponds to a `RangeIndex`,
                # fastparquet will need use the pandas metadata to
                # construct the index. Otherwise, the index will correspond
                # to a column named "__index_level_0__".  We will need to
                # check the `ParquetFile` object for this column below.
                index = []
                null_index_name = True
            columns += index

        # Use global `parquet_file` object.  Need to reattach
        # the desired row_group
        parquet_file = kwargs.pop("parquet_file", None)

        # Always convert pieces to list
        if not isinstance(pieces, list):
            pieces = [pieces]

        sample = pieces[0]
        if isinstance(sample, tuple):
            if isinstance(sample[0], str):
                # We have paths to read from
                assert parquet_file is None

                row_groups = []
                rg_offset = 0
                parquet_file = ParquetFile(
                    [p[0] for p in pieces],
                    open_with=fs.open,
                    root=base_path or False,
                    **kwargs.get("dataset", {}),
                )
                for piece in pieces:
                    _pf = (
                        parquet_file
                        if len(pieces) == 1
                        else ParquetFile(
                            piece[0],
                            open_with=fs.open,
                            root=base_path or False,
                            **kwargs.get("dataset", {}),
                        )
                    )
                    n_local_row_groups = len(_pf.row_groups)
                    local_rg_indices = piece[1] or list(range(n_local_row_groups))
                    row_groups += [
                        parquet_file.row_groups[rg + rg_offset]
                        for rg in local_rg_indices
                    ]
                    rg_offset += n_local_row_groups
                update_parquet_file = len(row_groups) < len(parquet_file.row_groups)

            elif parquet_file:
                row_groups = []
                for piece in pieces:
                    # `piece[1]` will contain actual row-group objects,
                    # but they may be pickled
                    rgs = piece[0]
                    if isinstance(rgs, bytes):
                        rgs = pickle.loads(rgs)
                    row_groups += rgs
                update_parquet_file = True

            else:
                raise ValueError("Neither path nor ParquetFile detected!")

            if update_parquet_file:
                with _FP_FILE_LOCK:
                    for rg in row_groups:
                        for chunk in rg.columns:
                            s = chunk.file_path
                            if s and isinstance(s, bytes):
                                chunk.file_path = s.decode()

                    parquet_file.fmd.row_groups = row_groups
                    # NOTE: May lose cats after `_set_attrs` call
                    save_cats = parquet_file.cats
                    parquet_file._set_attrs()
                    parquet_file.cats = save_cats

            if null_index_name:
                if "__index_level_0__" in parquet_file.columns:
                    # See "Handling a None-labeled index" comment above
                    index = ["__index_level_0__"]
                    columns += index

            # Update hive-partitioning information if necessary
            parquet_file.cats = root_cats or {}
            if root_cats:
                parquet_file.file_scheme = root_file_scheme

            parquet_file._dtypes = (
                lambda *args: parquet_file.dtypes
            )  # ugly patch, could be fixed

            # Convert ParquetFile to pandas
            return cls.pf_to_pandas(
                parquet_file,
                fs=fs,
                columns=columns,
                categories=categories,
                index=index,
                **kwargs.get("read", {}),
            )

        else:
            # `sample` is NOT a tuple
            raise ValueError(f"Expected tuple, got {type(sample)}")

    @classmethod
    def pf_to_pandas(
        cls,
        pf,
        fs=None,
        columns=None,
        categories=None,
        index=None,
        open_file_options=None,
        **kwargs,
    ):
        # This method was mostly copied from the fastparquet
        # `ParquetFile.to_pandas` definition. We maintain our
        # own implementation in Dask to enable better remote
        # file-handling control

        # Handle selected columns
        if columns is not None:
            columns = columns[:]
        else:
            columns = pf.columns + list(pf.cats)
        if index:
            columns += [i for i in index if i not in columns]

        # Extract row-groups and pre-allocate df
        rgs = pf.row_groups
        size = sum(rg.num_rows for rg in rgs)
        df, views = pf.pre_allocate(size, columns, categories, index)
        if (
            Version(fastparquet.__version__) <= Version("2023.02.0")
            and PANDAS_GE_201
            and df.columns.empty
        ):
            df.columns = pd.Index([], dtype=object)
        start = 0

        # Get a map of file names -> row-groups
        fn_rg_map = defaultdict(list)
        for rg in rgs:
            fn = pf.row_group_filename(rg)
            fn_rg_map[fn].append(rg)

        # Define file-opening options
        precache_options, open_file_options = _process_open_file_options(
            open_file_options,
            **(
                {
                    "allow_precache": False,
                    "default_cache": "readahead",
                }
                if _is_local_fs(fs)
                else {
                    "metadata": pf,
                    "columns": list(set(columns).intersection(pf.columns)),
                    "row_groups": [rgs for rgs in fn_rg_map.values()],
                    "default_engine": "fastparquet",
                    "default_cache": "readahead",
                }
            ),
        )

        with ExitStack() as stack:
            for fn, infile in zip(
                fn_rg_map.keys(),
                _open_input_files(
                    list(fn_rg_map.keys()),
                    fs=fs,
                    context_stack=stack,
                    precache_options=precache_options,
                    **open_file_options,
                ),
            ):
                for rg in fn_rg_map[fn]:
                    thislen = rg.num_rows
                    parts = {
                        name: (
                            v
                            if name.endswith("-catdef")
                            else v[start : start + thislen]
                        )
                        for (name, v) in views.items()
                    }

                    # Add row-group data to df
                    pf.read_row_group_file(
                        rg,
                        columns,
                        categories,
                        index,
                        assign=parts,
                        partition_meta=pf.partition_meta,
                        infile=infile,
                        **kwargs,
                    )
                    start += thislen
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
        object_encoding="utf8",
        index_cols=None,
        custom_metadata=None,
        **kwargs,
    ):
        if index_cols is None:
            index_cols = []
        if append and division_info is None:
            ignore_divisions = True
        fs.mkdirs(path, exist_ok=True)
        if object_encoding == "infer" or (
            isinstance(object_encoding, dict) and "infer" in object_encoding.values()
        ):
            raise ValueError(
                '"infer" not allowed as object encoding, '
                "because this required data in memory."
            )

        metadata_file_exists = False
        if append:
            try:
                # to append to a dataset without _metadata, need to load
                # _common_metadata or any data file here
                pf = fastparquet.api.ParquetFile(path, open_with=fs.open)
                metadata_file_exists = fs.exists(fs.sep.join([path, "_metadata"]))
            except (OSError, ValueError):
                # append for create
                append = False
        if append:
            from dask.dataframe._pyarrow import to_object_string

            if pf.file_scheme not in ["hive", "empty", "flat"]:
                raise ValueError(
                    "Requested file scheme is hive, but existing file scheme is not."
                )
            elif (set(pf.columns) != set(df.columns) - set(partition_on)) or (
                set(partition_on) != set(pf.cats)
            ):
                raise ValueError(
                    "Appended columns not the same.\n"
                    "Previous: {} | New: {}".format(pf.columns, list(df.columns))
                )
            elif (
                pd.Series(pf.dtypes).loc[pf.columns]
                != to_object_string(df[pf.columns]._meta).dtypes
            ).any():
                raise ValueError(
                    "Appended dtypes differ.\n{}".format(
                        set(pf.dtypes.items()) ^ set(df.dtypes.items())
                    )
                )
            else:
                df = df[pf.columns + partition_on]

            fmd = pf.fmd
            i_offset = fastparquet.writer.find_max_part(fmd.row_groups)
            if not ignore_divisions:
                if not set(index_cols).intersection([division_info["name"]]):
                    ignore_divisions = True
            if not ignore_divisions:
                minmax = fastparquet.api.sorted_partitioned_columns(pf)
                # If fastparquet detects that a partitioned column isn't sorted, it won't
                # appear in the resulting min/max dictionary
                old_end = (
                    minmax[index_cols[0]]["max"][-1]
                    if index_cols[0] in minmax
                    else None
                )
                divisions = division_info["divisions"]
                if old_end is not None and divisions[0] <= old_end:
                    raise ValueError(
                        "The divisions of the appended dataframe overlap with "
                        "previously written divisions. If this is desired, set "
                        "``ignore_divisions=True`` to append anyway.\n"
                        "- End of last written partition: {old_end}\n"
                        "- Start of first new partition: {divisions[0]}"
                    )
        else:
            fmd = fastparquet.writer.make_metadata(
                df._meta,
                object_encoding=object_encoding,
                index_cols=index_cols,
                ignore_columns=partition_on,
                **kwargs,
            )
            i_offset = 0
        if custom_metadata is not None:
            kvm = fmd.key_value_metadata or []
            kvm.extend(
                [
                    fastparquet.parquet_thrift.KeyValue(key=key, value=value)
                    for key, value in custom_metadata.items()
                ]
            )
            fmd.key_value_metadata = kvm

        extra_write_kwargs = {"fmd": fmd}
        return i_offset, fmd, metadata_file_exists, extra_write_kwargs

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
        custom_metadata=None,
        **kwargs,
    ):
        # Update key/value metadata if necessary
        fmd = copy.copy(fmd)
        for s in fmd.schema:
            try:
                # can be coerced to bytes on copy
                s.name = s.name.decode()
            except AttributeError:
                pass
        if custom_metadata and fmd is not None:
            fmd.key_value_metadata = fmd.key_value_metadata + (
                [
                    fastparquet.parquet_thrift.KeyValue(key=key, value=value)
                    for key, value in custom_metadata.items()
                ]
            )

        if not len(df):
            # Write nothing for empty partitions
            rgs = []
        elif partition_on:
            mkdirs = lambda x: fs.mkdirs(x, exist_ok=True)
            if Version(fastparquet.__version__) >= Version("0.1.4"):
                rgs = partition_on_columns(
                    df, partition_on, path, filename, fmd, compression, fs.open, mkdirs
                )
            else:
                rgs = partition_on_columns(
                    df,
                    partition_on,
                    path,
                    filename,
                    fmd,
                    fs.sep,
                    compression,
                    fs.open,
                    mkdirs,
                )
        else:
            with fs.open(fs.sep.join([path, filename]), "wb") as fil:
                fmd.num_rows = len(df)
                rg = make_part_file(
                    fil, df, fmd.schema, compression=compression, fmd=fmd
                )
            for chunk in rg.columns:
                chunk.file_path = filename
            rgs = [rg]
        if return_metadata:
            return rgs
        else:
            return []

    @classmethod
    def write_metadata(cls, parts, meta, fs, path, append=False, **kwargs):
        _meta = copy.copy(meta)
        rgs = meta.row_groups
        if parts:
            for rg in parts:
                if rg is not None:
                    if isinstance(rg, list):
                        for r in rg:
                            rgs.append(r)
                    else:
                        rgs.append(rg)
            _meta.row_groups = rgs
            fn = fs.sep.join([path, "_metadata"])
            fastparquet.writer.write_common_metadata(
                fn, _meta, open_with=fs.open, no_row_groups=False
            )

        # if appending, could skip this, but would need to check existence
        fn = fs.sep.join([path, "_common_metadata"])
        fastparquet.writer.write_common_metadata(fn, _meta, open_with=fs.open)
