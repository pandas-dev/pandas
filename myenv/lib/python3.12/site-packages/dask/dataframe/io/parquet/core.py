from __future__ import annotations

import contextlib
import math
import warnings
from typing import Literal

import pandas as pd
import tlz as toolz
from fsspec.core import get_fs_token_paths
from fsspec.utils import stringify_path

import dask
from dask.base import tokenize
from dask.blockwise import BlockIndex
from dask.dataframe.backends import dataframe_creation_dispatch
from dask.dataframe.core import DataFrame, Scalar
from dask.dataframe.io.io import from_map
from dask.dataframe.io.parquet.utils import (
    Engine,
    _sort_and_analyze_paths,
    _split_user_options,
)
from dask.dataframe.io.utils import DataFrameIOFunction, _is_local_fs
from dask.dataframe.methods import concat
from dask.delayed import Delayed
from dask.highlevelgraph import HighLevelGraph
from dask.layers import DataFrameIOLayer
from dask.utils import apply, import_required, natural_sort_key, parse_bytes

__all__ = ("read_parquet", "to_parquet")

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
                (p.data["piece"], p.data.get("kwargs", {}))
                if hasattr(p, "data")
                else (p["piece"], p.get("kwargs", {}))
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


@dataframe_creation_dispatch.register_inplace("pandas")
def read_parquet(
    path,
    columns=None,
    filters=None,
    categories=None,
    index=None,
    storage_options=None,
    engine="auto",
    use_nullable_dtypes: bool | None = None,
    dtype_backend=None,
    calculate_divisions=None,
    ignore_metadata_file=False,
    metadata_task_size=None,
    split_row_groups="infer",
    blocksize="default",
    aggregate_files=None,
    parquet_file_extension=(".parq", ".parquet", ".pq"),
    filesystem=None,
    **kwargs,
):
    """
    Read a Parquet file into a Dask DataFrame

    This reads a directory of Parquet data into a Dask.dataframe, one file per
    partition.  It selects the index among the sorted columns if any exist.

    Parameters
    ----------
    path : str or list
        Source directory for data, or path(s) to individual parquet files.
        Prefix with a protocol like ``s3://`` to read from alternative
        filesystems. To read from multiple files you can pass a globstring or a
        list of paths, with the caveat that they must all have the same
        protocol.
    columns : str or list, default None
        Field name(s) to read in as columns in the output. By default all
        non-index fields will be read (as determined by the pandas parquet
        metadata, if present). Provide a single field name instead of a list to
        read in the data as a Series.
    filters : Union[List[Tuple[str, str, Any]], List[List[Tuple[str, str, Any]]]], default None
        List of filters to apply, like ``[[('col1', '==', 0), ...], ...]``.
        Using this argument will result in row-wise filtering of the final partitions.

        Predicates can be expressed in disjunctive normal form (DNF). This means that
        the inner-most tuple describes a single column predicate. These inner predicates
        are combined with an AND conjunction into a larger predicate. The outer-most
        list then combines all of the combined filters with an OR disjunction.

        Predicates can also be expressed as a ``List[Tuple]``. These are evaluated
        as an AND conjunction. To express OR in predicates, one must use the
        (preferred for "pyarrow") ``List[List[Tuple]]`` notation.
    index : str, list or False, default None
        Field name(s) to use as the output frame index. By default will be
        inferred from the pandas parquet file metadata, if present. Use ``False``
        to read all fields as columns.
    categories : list or dict, default None
        For any fields listed here, if the parquet encoding is Dictionary,
        the column will be created with dtype category. Use only if it is
        guaranteed that the column is encoded as dictionary in all row-groups.
        If a list, assumes up to 2**16-1 labels; if a dict, specify the number
        of labels expected; if None, will load categories automatically for
        data written by dask, not otherwise.
    storage_options : dict, default None
        Key/value pairs to be passed on to the file-system backend, if any.
        Note that the default file-system backend can be configured with the
        ``filesystem`` argument, described below.
    open_file_options : dict, default None
        Key/value arguments to be passed along to ``AbstractFileSystem.open``
        when each parquet data file is open for reading. Experimental
        (optimized) "precaching" for remote file systems (e.g. S3, GCS) can
        be enabled by adding ``{"method": "parquet"}`` under the
        ``"precache_options"`` key. Also, a custom file-open function can be
        used (instead of ``AbstractFileSystem.open``), by specifying the
        desired function under the ``"open_file_func"`` key.
    engine : {'auto', 'pyarrow'}
        Parquet library to use. Defaults to 'auto', which uses ``pyarrow`` if
        it is installed, and falls back to the deprecated ``fastparquet`` otherwise.
        Note that ``fastparquet`` does not support all functionality offered by
        ``pyarrow``.
        This is also used by third-party packages (e.g. CuDF) to inject bespoke engines.
    use_nullable_dtypes : {False, True}
        Whether to use extension dtypes for the resulting ``DataFrame``.

        .. note::

            This option is deprecated. Use "dtype_backend" instead.

    dtype_backend : {'numpy_nullable', 'pyarrow'}, defaults to NumPy backed DataFrames
        Which dtype_backend to use, e.g. whether a DataFrame should have NumPy arrays,
        nullable dtypes are used for all dtypes that have a nullable implementation
        when 'numpy_nullable' is set, pyarrow is used for all dtypes if 'pyarrow'
        is set.
        ``dtype_backend="pyarrow"`` requires ``pandas`` 1.5+.
    calculate_divisions : bool, default False
        Whether to use min/max statistics from the footer metadata (or global
        ``_metadata`` file) to calculate divisions for the output DataFrame
        collection. Divisions will not be calculated if statistics are missing.
        This option will be ignored if ``index`` is not specified and there is
        no physical index column specified in the custom "pandas" Parquet
        metadata. Note that ``calculate_divisions=True`` may be extremely slow
        when no global ``_metadata`` file is present, especially when reading
        from remote storage. Set this to ``True`` only when known divisions
        are needed for your workload (see :ref:`dataframe-design-partitions`).
    ignore_metadata_file : bool, default False
        Whether to ignore the global ``_metadata`` file (when one is present).
        If ``True``, or if the global ``_metadata`` file is missing, the parquet
        metadata may be gathered and processed in parallel. Parallel metadata
        processing is currently supported for ``ArrowDatasetEngine`` only.
    metadata_task_size : int, default configurable
        If parquet metadata is processed in parallel (see ``ignore_metadata_file``
        description above), this argument can be used to specify the number of
        dataset files to be processed by each task in the Dask graph.  If this
        argument is set to ``0``, parallel metadata processing will be disabled.
        The default values for local and remote filesystems can be specified
        with the "metadata-task-size-local" and "metadata-task-size-remote"
        config fields, respectively (see "dataframe.parquet").
    split_row_groups : 'infer', 'adaptive', bool, or int, default 'infer'
        If True, then each output dataframe partition will correspond to a single
        parquet-file row-group. If False, each partition will correspond to a
        complete file.  If a positive integer value is given, each dataframe
        partition will correspond to that number of parquet row-groups (or fewer).
        If 'adaptive', the metadata of each file will be used to ensure that every
        partition satisfies ``blocksize``. If 'infer' (the default), the
        uncompressed storage-size metadata in the first file will be used to
        automatically set ``split_row_groups`` to either 'adaptive' or ``False``.
    blocksize : int or str, default 'default'
        The desired size of each output ``DataFrame`` partition in terms of total
        (uncompressed) parquet storage space. This argument is currently used to
        set the default value of ``split_row_groups`` (using row-group metadata
        from a single file), and will be ignored if ``split_row_groups`` is not
        set to 'infer' or 'adaptive'. Default is 256 MiB.
    aggregate_files : bool or str, default None
        WARNING: Passing a string argument to ``aggregate_files`` will result
        in experimental behavior. This behavior may change in the future.

        Whether distinct file paths may be aggregated into the same output
        partition. This parameter is only used when `split_row_groups` is set to
        'infer', 'adaptive' or to an integer >1. A setting of True means that any
        two file paths may be aggregated into the same output partition, while
        False means that inter-file aggregation is prohibited.

        For "hive-partitioned" datasets, a "partition"-column name can also be
        specified. In this case, we allow the aggregation of any two files
        sharing a file path up to, and including, the corresponding directory name.
        For example, if ``aggregate_files`` is set to ``"section"`` for the
        directory structure below, ``03.parquet`` and ``04.parquet`` may be
        aggregated together, but ``01.parquet`` and ``02.parquet`` cannot be.
        If, however, ``aggregate_files`` is set to ``"region"``, ``01.parquet``
        may be aggregated with ``02.parquet``, and ``03.parquet`` may be aggregated
        with ``04.parquet``::

            dataset-path/
            ├── region=1/
            │   ├── section=a/
            │   │   └── 01.parquet
            │   ├── section=b/
            │   └── └── 02.parquet
            └── region=2/
                ├── section=a/
                │   ├── 03.parquet
                └── └── 04.parquet

        Note that the default behavior of ``aggregate_files`` is ``False``.
    parquet_file_extension: str, tuple[str], or None, default (".parq", ".parquet", ".pq")
        A file extension or an iterable of extensions to use when discovering
        parquet files in a directory. Files that don't match these extensions
        will be ignored. This argument only applies when ``paths`` corresponds
        to a directory and no ``_metadata`` file is present (or
        ``ignore_metadata_file=True``). Passing in ``parquet_file_extension=None``
        will treat all files in the directory as parquet files.

        The purpose of this argument is to ensure that the engine will ignore
        unsupported metadata files (like Spark's '_SUCCESS' and 'crc' files).
        It may be necessary to change this argument if the data files in your
        parquet dataset do not end in ".parq", ".parquet", or ".pq".
    filesystem: "fsspec", "arrow", or fsspec.AbstractFileSystem backend to use.
    dataset: dict, default None
        Dictionary of options to use when creating a ``pyarrow.dataset.Dataset`` object.
        These options may include a "filesystem" key to configure the desired
        file-system backend. However, the top-level ``filesystem`` argument will always
        take precedence.

        **Note**: The ``dataset`` options may include a "partitioning" key.
        However, since ``pyarrow.dataset.Partitioning``
        objects cannot be serialized, the value can be a dict of key-word
        arguments for the ``pyarrow.dataset.partitioning`` API
        (e.g. ``dataset={"partitioning": {"flavor": "hive", "schema": ...}}``).
        Note that partitioned columns will not be converted to categorical
        dtypes when a custom partitioning schema is specified in this way.
    read: dict, default None
        Dictionary of options to pass through to ``engine.read_partitions``
        using the ``read`` key-word argument.
    arrow_to_pandas: dict, default None
        Dictionary of options to use when converting from ``pyarrow.Table`` to
        a pandas ``DataFrame`` object. Only used by the "arrow" engine.
    **kwargs: dict (of dicts)
        Options to pass through to ``engine.read_partitions`` as stand-alone
        key-word arguments. Note that these options will be ignored by the
        engines defined in ``dask.dataframe``, but may be used by other custom
        implementations.

    Examples
    --------
    >>> df = dd.read_parquet('s3://bucket/my-parquet-data')  # doctest: +SKIP

    See Also
    --------
    to_parquet
    pyarrow.parquet.ParquetDataset
    """
    # Handle `chunksize` deprecation
    if "chunksize" in kwargs:
        if blocksize != "default":
            raise ValueError(
                "The `chunksize` argument is now deprecated. "
                "Please use the `blocksize` argument instead."
            )
        blocksize = kwargs.pop("chunksize")
        warnings.warn(
            "The `chunksize` argument is deprecated, and will be "
            "removed in a future release. Setting the `blocksize` "
            "argument instead. Please see documentation on the "
            "`blocksize` argument for more information.",
            FutureWarning,
        )

    # FutureWarning for `aggregate_files`
    if aggregate_files and isinstance(aggregate_files, str):
        warnings.warn(
            "String support for `aggregate_files` is experimental. "
            "Behavior may change in the future. ",
            FutureWarning,
        )

    if "read_from_paths" in kwargs:
        kwargs.pop("read_from_paths")
        warnings.warn(
            "`read_from_paths` is no longer supported and will be ignored.",
            FutureWarning,
        )

    # Handle gather_statistics deprecation
    if "gather_statistics" in kwargs:
        if calculate_divisions is None:
            calculate_divisions = kwargs.pop("gather_statistics")
            warnings.warn(
                "``gather_statistics`` is deprecated and will be removed in a "
                "future release. Please use ``calculate_divisions`` instead.",
                FutureWarning,
            )
        else:
            warnings.warn(
                f"``gather_statistics`` is deprecated. Ignoring this option "
                f"in favor of ``calculate_divisions={calculate_divisions}``",
                FutureWarning,
            )
    calculate_divisions = bool(calculate_divisions)

    # We support a top-level `parquet_file_extension` kwarg, but
    # must check if the deprecated `require_extension` option is
    # being passed to the engine. If `parquet_file_extension` is
    # set to the default value, and `require_extension` was also
    # specified, we will use `require_extension` but warn the user.
    if (
        "dataset" in kwargs
        and "require_extension" in kwargs["dataset"]
        and parquet_file_extension == (".parq", ".parquet", ".pq")
    ):
        parquet_file_extension = kwargs["dataset"].pop("require_extension")
        warnings.warn(
            "require_extension is deprecated, and will be removed from "
            "read_parquet in a future release. Please use the top-level "
            "parquet_file_extension argument instead.",
            FutureWarning,
        )

    # Store initial function arguments
    input_kwargs = {
        "columns": columns,
        "filters": filters,
        "categories": categories,
        "index": index,
        "storage_options": storage_options,
        "engine": engine,
        "use_nullable_dtypes": use_nullable_dtypes,
        "dtype_backend": dtype_backend,
        "calculate_divisions": calculate_divisions,
        "ignore_metadata_file": ignore_metadata_file,
        "metadata_task_size": metadata_task_size,
        "split_row_groups": split_row_groups,
        "blocksize": blocksize,
        "aggregate_files": aggregate_files,
        "parquet_file_extension": parquet_file_extension,
        "filesystem": filesystem,
        **kwargs,
    }

    if isinstance(columns, str):
        input_kwargs["columns"] = [columns]
        df = read_parquet(path, **input_kwargs)
        return df[columns]

    if columns is not None:
        columns = list(columns)

    engine = get_engine(engine)

    if hasattr(path, "name"):
        path = stringify_path(path)

    # Update input_kwargs
    input_kwargs.update({"columns": columns, "engine": engine})

    # Process and split user options
    (
        dataset_options,
        read_options,
        open_file_options,
        other_options,
    ) = _split_user_options(**kwargs)

    # Extract global filesystem and paths
    fs, paths, dataset_options, open_file_options = engine.extract_filesystem(
        path,
        filesystem,
        dataset_options,
        open_file_options,
        storage_options,
    )
    read_options["open_file_options"] = open_file_options
    paths = sorted(paths, key=natural_sort_key)  # numeric rather than glob ordering

    auto_index_allowed = False
    if index is None:
        # User is allowing auto-detected index
        auto_index_allowed = True
    if index and isinstance(index, str):
        index = [index]

    if split_row_groups in ("infer", "adaptive"):
        # Using blocksize to plan partitioning
        blocksize = engine.default_blocksize() if blocksize == "default" else blocksize
    else:
        # Not using blocksize - Set to `None`
        blocksize = None

    read_metadata_result = engine.read_metadata(
        fs,
        paths,
        categories=categories,
        index=index,
        use_nullable_dtypes=use_nullable_dtypes,
        dtype_backend=dtype_backend,
        gather_statistics=calculate_divisions,
        filters=filters,
        split_row_groups=split_row_groups,
        blocksize=blocksize,
        aggregate_files=aggregate_files,
        ignore_metadata_file=ignore_metadata_file,
        metadata_task_size=metadata_task_size,
        parquet_file_extension=parquet_file_extension,
        dataset=dataset_options,
        read=read_options,
        **other_options,
    )

    # In the future, we may want to give the engine the
    # option to return a dedicated element for `common_kwargs`.
    # However, to avoid breaking the API, we just embed this
    # data in the first element of `parts` for now.
    # The logic below is intended to handle backward and forward
    # compatibility with a user-defined engine.
    meta, statistics, parts, index = read_metadata_result[:4]
    common_kwargs = {}
    aggregation_depth = False
    if len(parts):
        # For now, `common_kwargs` and `aggregation_depth`
        # may be stored in the first element of `parts`
        common_kwargs = parts[0].pop("common_kwargs", {})
        aggregation_depth = parts[0].pop("aggregation_depth", aggregation_depth)
        split_row_groups = parts[0].pop("split_row_groups", split_row_groups)

    # Parse dataset statistics from metadata (if available)
    parts, divisions, index = process_statistics(
        parts,
        statistics,
        filters,
        index,
        (blocksize if split_row_groups is True else None),
        split_row_groups,
        fs,
        aggregation_depth,
    )

    # Account for index and columns arguments.
    # Modify `meta` dataframe accordingly
    meta, index, columns = set_index_columns(meta, index, columns, auto_index_allowed)
    if meta.index.name == NONE_LABEL:
        meta.index.name = None

    if len(divisions) < 2:
        # empty dataframe - just use meta
        divisions = (None, None)
        io_func = lambda x: x
        parts = [meta]
    else:
        # Use IO function wrapper
        io_func = ParquetFunctionWrapper(
            engine,
            fs,
            meta,
            columns,
            index,
            dtype_backend,
            {},  # All kwargs should now be in `common_kwargs`
            common_kwargs,
        )

    # If we are using a remote filesystem and retries is not set, bump it
    # to be more fault tolerant, as transient transport errors can occur.
    # The specific number 5 isn't hugely motivated: it's less than ten and more
    # than two.
    annotations = dask.get_annotations()
    if "retries" not in annotations and not _is_local_fs(fs):
        ctx = dask.annotate(retries=5)
    else:
        ctx = contextlib.nullcontext()  # type: ignore

    with ctx:
        # Construct the output collection with from_map
        return from_map(
            io_func,
            parts,
            meta=meta,
            divisions=divisions,
            label="read-parquet",
            token=tokenize(
                path, dask.config.get("dataframe.convert-string"), **input_kwargs
            ),
            enforce_metadata=False,
            creation_info={
                "func": read_parquet,
                "args": (path,),
                "kwargs": input_kwargs,
            },
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


def to_parquet(
    df,
    path,
    engine="auto",
    compression="snappy",
    write_index=True,
    append=False,
    overwrite=False,
    ignore_divisions=False,
    partition_on=None,
    storage_options=None,
    custom_metadata=None,
    write_metadata_file=None,
    compute=True,
    compute_kwargs=None,
    schema="infer",
    name_function=None,
    filesystem=None,
    **kwargs,
):
    """Store Dask.dataframe to Parquet files

    Notes
    -----
    Each partition will be written to a separate file.

    Parameters
    ----------
    df : dask.dataframe.DataFrame
    path : string or pathlib.Path
        Destination directory for data.  Prepend with protocol like ``s3://``
        or ``hdfs://`` for remote data.
    engine : {'auto', 'pyarrow'}
        Parquet library to use. Defaults to 'auto', which uses ``pyarrow`` if
        it is installed, and falls back to the deprecated ``fastparquet`` otherwise.
        Note that ``fastparquet`` does not support all functionality offered by
        ``pyarrow``.
        This is also used by third-party packages (e.g. CuDF) to inject bespoke engines.
    compression : string or dict, default 'snappy'
        Either a string like ``"snappy"`` or a dictionary mapping column names
        to compressors like ``{"name": "gzip", "values": "snappy"}``. Defaults
        to ``"snappy"``.
    write_index : boolean, default True
        Whether or not to write the index. Defaults to True.
    append : bool, default False
        If False (default), construct data-set from scratch. If True, add new
        row-group(s) to an existing data-set. In the latter case, the data-set
        must exist, and the schema must match the input data.
    overwrite : bool, default False
        Whether or not to remove the contents of `path` before writing the dataset.
        The default is False.  If True, the specified path must correspond to
        a directory (but not the current working directory).  This option cannot
        be set to True if `append=True`.
        NOTE: `overwrite=True` will remove the original data even if the current
        write operation fails.  Use at your own risk.
    ignore_divisions : bool, default False
        If False (default) raises error when previous divisions overlap with
        the new appended divisions. Ignored if append=False.
    partition_on : list, default None
        Construct directory-based partitioning by splitting on these fields'
        values. Each dask partition will result in one or more datafiles,
        there will be no global groupby.
    storage_options : dict, default None
        Key/value pairs to be passed on to the file-system backend, if any.
    custom_metadata : dict, default None
        Custom key/value metadata to include in all footer metadata (and
        in the global "_metadata" file, if applicable).  Note that the custom
        metadata may not contain the reserved b"pandas" key.
    write_metadata_file : bool or None, default None
        Whether to write the special ``_metadata`` file. If ``None`` (the
        default), a ``_metadata`` file will only be written if ``append=True``
        and the dataset already has a ``_metadata`` file.
    compute : bool, default True
        If ``True`` (default) then the result is computed immediately. If
        ``False`` then a ``dask.dataframe.Scalar`` object is returned for
        future computation.
    compute_kwargs : dict, default True
        Options to be passed in to the compute method
    schema : pyarrow.Schema, dict, "infer", or None, default "infer"
        Global schema to use for the output dataset. Defaults to "infer", which
        will infer the schema from the dask dataframe metadata. This is usually
        sufficient for common schemas, but notably will fail for ``object``
        dtype columns that contain things other than strings. These columns
        will require an explicit schema be specified. The schema for a subset
        of columns can be overridden by passing in a dict of column names to
        pyarrow types (for example ``schema={"field": pa.string()}``); columns
        not present in this dict will still be automatically inferred.
        Alternatively, a full ``pyarrow.Schema`` may be passed, in which case
        no schema inference will be done. Passing in ``schema=None`` will
        disable the use of a global file schema - each written file may use a
        different schema dependent on the dtypes of the corresponding
        partition.
    name_function : callable, default None
        Function to generate the filename for each output partition.
        The function should accept an integer (partition index) as input and
        return a string which will be used as the filename for the corresponding
        partition. Should preserve the lexicographic order of partitions.
        If not specified, files will created using the convention
        ``part.0.parquet``, ``part.1.parquet``, ``part.2.parquet``, ...
        and so on for each partition in the DataFrame.
    filesystem: "fsspec", "arrow", or fsspec.AbstractFileSystem backend to use.
    **kwargs :
        Extra options to be passed on to the specific backend.

    Examples
    --------
    >>> df = dd.read_csv(...)  # doctest: +SKIP
    >>> df.to_parquet('/path/to/output/', ...)  # doctest: +SKIP

    By default, files will be created in the specified output directory using the
    convention ``part.0.parquet``, ``part.1.parquet``, ``part.2.parquet``, ... and so on for
    each partition in the DataFrame. To customize the names of each file, you can use the
    ``name_function=`` keyword argument. The function passed to ``name_function`` will be
    used to generate the filename for each partition and should expect a partition's index
    integer as input and return a string which will be used as the filename for the corresponding
    partition. Strings produced by ``name_function`` must preserve the order of their respective
    partition indices.

    For example:

    >>> name_function = lambda x: f"data-{x}.parquet"
    >>> df.to_parquet('/path/to/output/', name_function=name_function)  # doctest: +SKIP

    will result in the following files being created::

        /path/to/output/
            ├── data-0.parquet
            ├── data-1.parquet
            ├── data-2.parquet
            └── ...

    See Also
    --------
    read_parquet: Read parquet data to dask.dataframe
    """
    compute_kwargs = compute_kwargs or {}

    partition_on = partition_on or []
    if isinstance(partition_on, str):
        partition_on = [partition_on]

    if set(partition_on) - set(df.columns):
        raise ValueError(
            "Partitioning on non-existent column. "
            "partition_on=%s ."
            "columns=%s" % (str(partition_on), str(list(df.columns)))
        )

    if df.columns.inferred_type not in {"string", "empty"}:
        raise ValueError("parquet doesn't support non-string column names")

    engine = get_engine(engine)

    if hasattr(path, "name"):
        path = stringify_path(path)

    fs, _paths, _, _ = engine.extract_filesystem(
        path,
        filesystem=filesystem,
        dataset_options={},
        open_file_options={},
        storage_options=storage_options,
    )
    assert len(_paths) == 1, "only one path"
    path = _paths[0]

    if overwrite:
        if append:
            raise ValueError("Cannot use both `overwrite=True` and `append=True`!")

        if fs.exists(path) and fs.isdir(path):
            # Check for any previous parquet layers reading from a file in the
            # output directory, since deleting those files now would result in
            # errors or incorrect results.
            for layer_name, layer in df.dask.layers.items():
                if layer_name.startswith("read-parquet-") and isinstance(
                    layer, DataFrameIOLayer
                ):
                    path_with_slash = path.rstrip("/") + "/"  # ensure trailing slash
                    for input in layer.inputs:
                        # Note that `input` may be either `dict` or `List[dict]`
                        for piece_dict in input if isinstance(input, list) else [input]:
                            if piece_dict["piece"][0].startswith(path_with_slash):
                                raise ValueError(
                                    "Reading and writing to the same parquet file within "
                                    "the same task graph is not supported."
                                )

            # Don't remove the directory if it's the current working directory
            if _is_local_fs(fs):
                working_dir = fs.expand_path(".")[0]
                if path.rstrip("/") == working_dir.rstrip("/"):
                    raise ValueError(
                        "Cannot clear the contents of the current working directory!"
                    )

            # It's safe to clear the output directory
            fs.rm(path, recursive=True)

    # Always skip divisions checks if divisions are unknown
    if not df.known_divisions:
        ignore_divisions = True

    # Save divisions and corresponding index name. This is necessary,
    # because we may be resetting the index to write the file
    division_info = {"divisions": df.divisions, "name": df.index.name}
    if division_info["name"] is None:
        # As of 0.24.2, pandas will rename an index with name=None
        # when df.reset_index() is called.  The default name is "index",
        # but dask will always change the name to the NONE_LABEL constant
        if NONE_LABEL not in df.columns:
            division_info["name"] = NONE_LABEL
        elif write_index:
            raise ValueError(
                "Index must have a name if __null_dask_index__ is a column."
            )
        else:
            warnings.warn(
                "If read back by Dask, column named __null_dask_index__ "
                "will be set to the index (and renamed to None)."
            )

    # There are some "reserved" names that may be used as the default column
    # name after resetting the index. However, we don't want to treat it as
    # a "special" name if the string is already used as a "real" column name.
    reserved_names = []
    for name in ["index", "level_0"]:
        if name not in df.columns:
            reserved_names.append(name)

    # If write_index==True (default), reset the index and record the
    # name of the original index in `index_cols` (we will set the name
    # to the NONE_LABEL constant if it is originally `None`).
    # `fastparquet` will use `index_cols` to specify the index column(s)
    # in the metadata.  `pyarrow` will revert the `reset_index` call
    # below if `index_cols` is populated (because pyarrow will want to handle
    # index preservation itself).  For both engines, the column index
    # will be written to "pandas metadata" if write_index=True
    index_cols = []
    if write_index:
        real_cols = set(df.columns)
        none_index = list(df._meta.index.names) == [None]
        df = df.reset_index()
        if none_index:
            df.columns = [
                c if c not in reserved_names else NONE_LABEL for c in df.columns
            ]
        index_cols = [c for c in set(df.columns) - real_cols]
    else:
        # Not writing index - might as well drop it
        df = df.reset_index(drop=True)

    if custom_metadata and b"pandas" in custom_metadata.keys():
        raise ValueError(
            "User-defined key/value metadata (custom_metadata) can not "
            "contain a b'pandas' key.  This key is reserved by Pandas, "
            "and overwriting the corresponding value can render the "
            "entire dataset unreadable."
        )

    # Engine-specific initialization steps to write the dataset.
    # Possibly create parquet metadata, and load existing stuff if appending
    i_offset, fmd, metadata_file_exists, extra_write_kwargs = engine.initialize_write(
        df,
        fs,
        path,
        append=append,
        ignore_divisions=ignore_divisions,
        partition_on=partition_on,
        division_info=division_info,
        index_cols=index_cols,
        schema=schema,
        custom_metadata=custom_metadata,
        **kwargs,
    )

    # By default we only write a metadata file when appending if one already
    # exists
    if append and write_metadata_file is None:
        write_metadata_file = metadata_file_exists

    # Check that custom name_function is valid,
    # and that it will produce unique names
    if name_function is not None:
        if not callable(name_function):
            raise ValueError("``name_function`` must be a callable with one argument.")
        filenames = [name_function(i + i_offset) for i in range(df.npartitions)]
        if len(set(filenames)) < len(filenames):
            raise ValueError("``name_function`` must produce unique filenames.")

    # If we are using a remote filesystem and retries is not set, bump it
    # to be more fault tolerant, as transient transport errors can occur.
    # The specific number 5 isn't hugely motivated: it's less than ten and more
    # than two.
    annotations = dask.get_annotations()
    if "retries" not in annotations and not _is_local_fs(fs):
        ctx = dask.annotate(retries=5)
    else:
        ctx = contextlib.nullcontext()

    # Create Blockwise layer for parquet-data write
    with ctx:
        data_write = df.map_partitions(
            ToParquetFunctionWrapper(
                engine,
                path,
                fs,
                partition_on,
                write_metadata_file,
                i_offset,
                name_function,
                toolz.merge(
                    kwargs,
                    {"compression": compression, "custom_metadata": custom_metadata},
                    extra_write_kwargs,
                ),
            ),
            BlockIndex((df.npartitions,)),
            # Pass in the original metadata to avoid
            # metadata emulation in `map_partitions`.
            # This is necessary, because we are not
            # expecting a dataframe-like output.
            meta=df._meta,
            enforce_metadata=False,
            transform_divisions=False,
            align_dataframes=False,
        )

    # Collect metadata and write _metadata.
    # TODO: Use tree-reduction layer (when available)
    if write_metadata_file:
        final_name = "metadata-" + data_write._name
        dsk = {
            (final_name, 0): (
                apply,
                engine.write_metadata,
                [
                    data_write.__dask_keys__(),
                    fmd,
                    fs,
                    path,
                ],
                {"append": append, "compression": compression},
            )
        }
    else:
        # NOTE: We still define a single task to tie everything together
        # when we are not writing a _metadata file. We do not want to
        # return `data_write` (or a `data_write.to_bag()`), because calling
        # `compute()` on a multi-partition collection requires the overhead
        # of trying to concatenate results on the client.
        final_name = "store-" + data_write._name
        dsk = {(final_name, 0): (lambda x: None, data_write.__dask_keys__())}

    # Convert data_write + dsk to computable collection
    graph = HighLevelGraph.from_collections(final_name, dsk, dependencies=(data_write,))
    out = Scalar(graph, final_name, "")

    if compute:
        out = out.compute(**compute_kwargs)

    # Invalidate the filesystem listing cache for the output path after write.
    # We do this before returning, even if `compute=False`. This helps ensure
    # that reading files that were just written succeeds.
    fs.invalidate_cache(path)

    return out


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

    # Get engine.
    # Note that "fastparquet" is not supported.
    if isinstance(engine, str):
        engine = get_engine(engine)
        if engine is not _ENGINES.get("pyarrow"):
            raise ValueError(
                "fastparquet is not a supported engine for create_metadata_file."
                "Please install pyarrow."
            )

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
    engine : {'auto', 'pyarrow', 'fastparquet'} or Engine subclass
        Parquet library to use. Defaults to 'auto', which uses ``pyarrow`` if
        it is installed, and falls back to the deprecated ``fastparquet`` otherwise.

    This can be used to inject third-party engine; e.g. from dask_cudf.
    """
    if isinstance(engine, type) and issubclass(engine, Engine):
        return engine

    if engine in ("arrow", "pyarrow-dataset"):
        warnings.warn(
            f"engine='{engine}' has been deprecated. "
            "Please use engine='pyarrow' instead.",
            category=FutureWarning,
        )
        engine = "pyarrow"

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

    if engine in ("auto", "fastparquet"):

        def warn_fastparquet():
            if engine == "fastparquet":
                warnings.warn(
                    "The fastparquet engine is deprecated and will be removed "
                    "in a future release. Please install pyarrow.",
                    category=FutureWarning,
                )
            else:
                warnings.warn(
                    "Could not find pyarrow; falling back to fastparquet, which "
                    "is deprecated and will be removed in a future release.",
                    category=FutureWarning,
                )

        if "fastparquet" in _ENGINES:
            warn_fastparquet()
            return _ENGINES["fastparquet"]

        try:
            import_required("fastparquet", "`fastparquet` not installed")
        except RuntimeError:
            if engine == "auto":
                raise RuntimeError("`pyarrow` not installed")
            else:
                raise
        else:
            from dask.dataframe.io.parquet.fastparquet import FastParquetEngine

            _ENGINES["fastparquet"] = fpq_eng = FastParquetEngine
            warn_fastparquet()
            return fpq_eng

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

        Note that the "fastparquet" engine does not currently support DNF for
        the filtering of partitioned columns (List[Tuple] is required).
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


def process_statistics(
    parts,
    statistics,
    filters,
    index,
    blocksize,
    split_row_groups,
    fs,
    aggregation_depth,
):
    """Process row-group column statistics in metadata
    Used in read_parquet.
    """
    if statistics and len(parts) != len(statistics):
        # It is up to the Engine to guarantee that these
        # lists are the same length (if statistics are defined).
        # This misalignment may be indicative of a bug or
        # incorrect read_parquet usage, so throw a warning.
        warnings.warn(
            f"Length of partition statistics ({len(statistics)}) "
            f"does not match the partition count ({len(parts)}). "
            f"This may indicate a bug or incorrect read_parquet "
            f"usage. We must ignore the statistics and disable: "
            f"filtering, divisions, and/or file aggregation."
        )
        statistics = []

    divisions = None
    if statistics:
        result = list(
            zip(
                *[
                    (part, stats)
                    for part, stats in zip(parts, statistics)
                    if stats["num-rows"] > 0
                ]
            )
        )
        parts, statistics = result or [[], []]
        if filters:
            parts, statistics = apply_filters(parts, statistics, filters)

        # Aggregate parts/statistics if we are splitting by row-group
        if blocksize or (split_row_groups and int(split_row_groups) > 1):
            parts, statistics = aggregate_row_groups(
                parts, statistics, blocksize, split_row_groups, fs, aggregation_depth
            )

        # Convert str index to list
        index = [index] if isinstance(index, str) else index

        # TODO: Remove `filters` criteria below after deprecation cycle.
        # We can then remove the `sorted_col_names` logic and warning.
        # See: https://github.com/dask/dask/pull/9661
        process_columns = index if index and len(index) == 1 else None
        if filters:
            process_columns = None

        # Use statistics to define divisions
        if process_columns or filters:
            sorted_col_names = []
            for sorted_column_info in sorted_columns(
                statistics, columns=process_columns
            ):
                if index and sorted_column_info["name"] in index:
                    divisions = sorted_column_info["divisions"]
                    break
                else:
                    # Filtered columns may also be sorted
                    sorted_col_names.append(sorted_column_info["name"])

            if index is None and sorted_col_names:
                assert bool(filters)  # Should only get here when filtering
                warnings.warn(
                    f"Sorted columns detected: {sorted_col_names}\n"
                    f"Use the `index` argument to set a sorted column as your "
                    f"index to create a DataFrame collection with known `divisions`.",
                    UserWarning,
                )

    divisions = divisions or (None,) * (len(parts) + 1)
    return parts, divisions, index


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


DataFrame.to_parquet.__doc__ = to_parquet.__doc__
