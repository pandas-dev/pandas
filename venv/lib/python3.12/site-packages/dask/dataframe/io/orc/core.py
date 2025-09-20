from __future__ import annotations

from typing import TYPE_CHECKING, Literal

from fsspec.core import get_fs_token_paths
from fsspec.utils import stringify_path

import dask.dataframe as dd
from dask.base import compute_as_if_collection, tokenize
from dask.dataframe.backends import dataframe_creation_dispatch
from dask.dataframe.io.orc.utils import ORCEngine
from dask.highlevelgraph import HighLevelGraph
from dask.utils import apply

if TYPE_CHECKING:
    from dask.dataframe.io.orc.arrow import ArrowORCEngine


def _get_engine(
    engine: Literal["pyarrow"] | ORCEngine,
) -> type[ArrowORCEngine] | ORCEngine:
    if engine == "pyarrow":
        from dask.dataframe.io.orc.arrow import ArrowORCEngine

        return ArrowORCEngine
    elif not isinstance(engine, ORCEngine):
        raise TypeError("engine must be 'pyarrow' or an ORCEngine object")
    return engine


@dataframe_creation_dispatch.register_inplace("pandas")
def read_orc(
    path,
    engine="pyarrow",
    columns=None,
    index=None,
    split_stripes=1,
    aggregate_files=None,
    storage_options=None,
):
    """Read dataframe from ORC file(s)

    Parameters
    ----------
    path: str or list(str)
        Location of file(s), which can be a full URL with protocol
        specifier, and may include glob character if a single string.
    engine: 'pyarrow' or ORCEngine
        Backend ORC engine to use for I/O. Default is "pyarrow".
    columns: None or list(str)
        Columns to load. If None, loads all.
    index: str
        Column name to set as index.
    split_stripes: int or False
        Maximum number of ORC stripes to include in each output-DataFrame
        partition. Use False to specify a 1-to-1 mapping between files
        and partitions. Default is 1.
    aggregate_files : bool, default False
        Whether distinct file paths may be aggregated into the same output
        partition. A setting of True means that any two file paths may be
        aggregated into the same output partition, while False means that
        inter-file aggregation is prohibited.
    storage_options: None or dict
        Further parameters to pass to the bytes backend.

    Returns
    -------
    Dask.DataFrame (even if there is only one column)

    Examples
    --------
    >>> df = dd.read_orc('https://github.com/apache/orc/raw/'
    ...                  'master/examples/demo-11-zlib.orc')  # doctest: +SKIP
    """

    # Get engine
    engine = _get_engine(engine)

    # Process file path(s)
    storage_options = storage_options or {}
    fs, fs_token, paths = get_fs_token_paths(
        path, mode="rb", storage_options=storage_options
    )

    # Let backend engine generate a list of parts
    # from the ORC metadata.  The backend should also
    # return the schema and DataFrame-collection metadata
    parts, schema, meta = engine.read_metadata(
        fs,
        paths,
        columns,
        index,
        split_stripes,
        aggregate_files,
    )

    if columns is not None and index in columns:
        columns = [col for col in columns if col != index]
    return dd.from_map(
        _read_orc,
        parts,
        engine=engine,
        fs=fs,
        schema=schema,
        index=index,
        meta=meta,
        columns=columns,
    )


def _read_orc(parts, *, engine, fs, schema, index, columns=None):
    if index is not None and columns is not None:
        columns.append(index)
    _df = engine.read_partition(
        fs,
        parts,
        schema,
        columns,
    )
    if index:
        _df = _df.set_index(index)

    return _df


def to_orc(
    df,
    path,
    engine="pyarrow",
    write_index=True,
    storage_options=None,
    compute=True,
    compute_kwargs=None,
):
    """Store Dask.dataframe to ORC files

    Notes
    -----
    Each partition will be written to a separate file.

    Parameters
    ----------
    df : dask.dataframe.DataFrame
    path : string or pathlib.Path
        Destination directory for data.  Prepend with protocol like ``s3://``
        or ``hdfs://`` for remote data.
    engine: 'pyarrow' or ORCEngine
        Backend ORC engine to use for I/O. Default is "pyarrow".
    write_index : boolean, default True
        Whether or not to write the index. Defaults to True.
    storage_options : dict, default None
        Key/value pairs to be passed on to the file-system backend, if any.
    compute : bool, default True
        If True (default) then the result is computed immediately. If False
        then a ``dask.delayed`` object is returned for future computation.
    compute_kwargs : dict, default True
        Options to be passed in to the compute method

    Examples
    --------
    >>> df = dd.read_csv(...)  # doctest: +SKIP
    >>> df.to_orc('/path/to/output/', ...)  # doctest: +SKIP

    See Also
    --------
    read_orc: Read ORC data to dask.dataframe
    """
    df = df.optimize()

    # Get engine
    engine = _get_engine(engine)

    if hasattr(path, "name"):
        path = stringify_path(path)
    fs, _, _ = get_fs_token_paths(path, mode="wb", storage_options=storage_options)
    # Trim any protocol information from the path before forwarding
    path = fs._strip_protocol(path)

    if not write_index:
        # Not writing index - might as well drop it
        df = df.reset_index(drop=True)

    # Use df.npartitions to define file-name list
    fs.mkdirs(path, exist_ok=True)
    filenames = [f"part.{i}.orc" for i in range(df.npartitions)]

    # Construct IO graph
    dsk = {}
    name = "to-orc-" + tokenize(
        df,
        fs,
        path,
        engine,
        write_index,
        storage_options,
    )
    final_name = name + "-final"
    for d, filename in enumerate(filenames):
        dsk[(name, d)] = (
            apply,
            engine.write_partition,
            [
                (df._name, d),
                path,
                fs,
                filename,
            ],
        )
    part_tasks = list(dsk.keys())
    dsk[(final_name, 0)] = (lambda x: None, part_tasks)
    graph = HighLevelGraph.from_collections((final_name, 0), dsk, dependencies=[df])

    from dask.dataframe import DataFrame

    # Compute or return future
    if compute:
        if compute_kwargs is None:
            compute_kwargs = dict()
        return compute_as_if_collection(DataFrame, graph, part_tasks, **compute_kwargs)

    from dask.dataframe import from_graph

    dsk.update(df.optimize().__dask_graph__())
    return from_graph(dsk, None, (None, None), [(final_name, 0)], "to-orc-persister")
