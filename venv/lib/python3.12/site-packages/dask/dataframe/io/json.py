from __future__ import annotations

import io
import os
from functools import partial
from itertools import zip_longest

import pandas as pd
from fsspec.core import open_files

import dask.dataframe as dd
from dask.base import compute as dask_compute
from dask.bytes import read_bytes
from dask.core import flatten
from dask.dataframe.backends import dataframe_creation_dispatch
from dask.dataframe.utils import insert_meta_param_description, make_meta
from dask.delayed import delayed


def to_json(
    df,
    url_path,
    orient="records",
    lines=None,
    storage_options=None,
    compute=True,
    encoding="utf-8",
    errors="strict",
    compression=None,
    compute_kwargs=None,
    name_function=None,
    **kwargs,
):
    """Write dataframe into JSON text files

    This utilises ``pandas.DataFrame.to_json()``, and most parameters are
    passed through - see its docstring.

    Differences: orient is 'records' by default, with lines=True; this
    produces the kind of JSON output that is most common in big-data
    applications, and which can be chunked when reading (see ``read_json()``).

    Parameters
    ----------
    df: dask.DataFrame
        Data to save
    url_path: str, list of str
        Location to write to. If a string, and there are more than one
        partitions in df, should include a glob character to expand into a
        set of file names, or provide a ``name_function=`` parameter.
        Supports protocol specifications such as ``"s3://"``.
    encoding, errors:
        The text encoding to implement, e.g., "utf-8" and how to respond
        to errors in the conversion (see ``str.encode()``).
    orient, lines, kwargs
        passed to pandas; if not specified, lines=True when orient='records',
        False otherwise.
    storage_options: dict
        Passed to backend file-system implementation
    compute: bool
        If true, immediately executes. If False, returns a set of delayed
        objects, which can be computed at a later time.
    compute_kwargs : dict, optional
        Options to be passed in to the compute method
    compression : string or None
        String like 'gzip' or 'xz'.
    name_function : callable, default None
        Function accepting an integer (partition index) and producing a
        string to replace the asterisk in the given filename globstring.
        Should preserve the lexicographic order of partitions.
    """
    if lines is None:
        lines = orient == "records"
    if orient != "records" and lines:
        raise ValueError('Line-delimited JSON is only available with orient="records".')
    kwargs["orient"] = orient
    kwargs["lines"] = lines and orient == "records"
    outfiles = open_files(
        url_path,
        "wt",
        encoding=encoding,
        errors=errors,
        name_function=name_function,
        num=df.npartitions,
        compression=compression,
        **(storage_options or {}),
    )
    parts = [
        delayed(write_json_partition)(d, outfile, kwargs)
        for outfile, d in zip(outfiles, df.to_delayed())
    ]
    if compute:
        if compute_kwargs is None:
            compute_kwargs = dict()
        return list(dask_compute(*parts, **compute_kwargs))
    else:
        return parts


def write_json_partition(df, openfile, kwargs):
    with openfile as f:
        df.to_json(f, **kwargs)
    return os.path.normpath(openfile.path)


@dataframe_creation_dispatch.register_inplace("pandas")
@insert_meta_param_description
def read_json(
    url_path,
    orient="records",
    lines=None,
    storage_options=None,
    blocksize=None,
    sample=2**20,
    encoding="utf-8",
    errors="strict",
    compression="infer",
    meta=None,
    engine=pd.read_json,
    include_path_column=False,
    path_converter=None,
    **kwargs,
):
    """Create a dataframe from a set of JSON files

    This utilises ``pandas.read_json()``, and most parameters are
    passed through - see its docstring.

    Differences: orient is 'records' by default, with lines=True; this
    is appropriate for line-delimited "JSON-lines" data, the kind of JSON output
    that is most common in big-data scenarios, and which can be chunked when
    reading (see ``read_json()``). All other options require blocksize=None,
    i.e., one partition per input file.

    Parameters
    ----------
    url_path: str, list of str
        Location to read from. If a string, can include a glob character to
        find a set of file names.
        Supports protocol specifications such as ``"s3://"``.
    encoding, errors:
        The text encoding to implement, e.g., "utf-8" and how to respond
        to errors in the conversion (see ``str.encode()``).
    orient, lines, kwargs
        passed to pandas; if not specified, lines=True when orient='records',
        False otherwise.
    storage_options: dict
        Passed to backend file-system implementation
    blocksize: None or int
        If None, files are not blocked, and you get one partition per input
        file. If int, which can only be used for line-delimited JSON files,
        each partition will be approximately this size in bytes, to the nearest
        newline character.
    sample: int
        Number of bytes to pre-load, to provide an empty dataframe structure
        to any blocks without data. Only relevant when using blocksize.
    encoding, errors:
        Text conversion, ``see bytes.decode()``
    compression : string or None
        String like 'gzip' or 'xz'.
    engine : callable or str, default ``pd.read_json``
        The underlying function that dask will use to read JSON files. By
        default, this will be the pandas JSON reader (``pd.read_json``).
        If a string is specified, this value will be passed under the ``engine``
        key-word argument to ``pd.read_json`` (only supported for pandas>=2.0).
    include_path_column : bool or str, optional
        Include a column with the file path where each row in the dataframe
        originated. If ``True``, a new column is added to the dataframe called
        ``path``. If ``str``, sets new column name. Default is ``False``.
    path_converter : function or None, optional
        A function that takes one argument and returns a string. Used to convert
        paths in the ``path`` column, for instance, to strip a common prefix from
        all the paths.
    $META

    Returns
    -------
    dask.DataFrame

    Examples
    --------
    Load single file

    >>> dd.read_json('myfile.1.json')  # doctest: +SKIP

    Load multiple files

    >>> dd.read_json('myfile.*.json')  # doctest: +SKIP

    >>> dd.read_json(['myfile.1.json', 'myfile.2.json'])  # doctest: +SKIP

    Load large line-delimited JSON files using partitions of approx
    256MB size

    >> dd.read_json('data/file*.csv', blocksize=2**28)
    """
    if lines is None:
        lines = orient == "records"
    if orient != "records" and lines:
        raise ValueError('Line-delimited JSON is only available with orient="records".')
    if blocksize and (orient != "records" or not lines):
        raise ValueError(
            "JSON file chunking only allowed for JSON-lines"
            "input (orient='records', lines=True)."
        )
    storage_options = storage_options or {}
    if include_path_column is True:
        include_path_column = "path"

    if path_converter is None:
        path_converter = lambda x: x

    # Handle engine string
    if isinstance(engine, str):
        engine = partial(pd.read_json, engine=engine)

    if blocksize:
        b_out = read_bytes(
            url_path,
            b"\n",
            blocksize=blocksize,
            sample=sample,
            compression=compression,
            include_path=include_path_column,
            **storage_options,
        )
        if include_path_column:
            first, chunks, paths = b_out
            first_path = path_converter(paths[0])
            path_dtype = pd.CategoricalDtype(path_converter(p) for p in paths)
            flat_paths = flatten(
                [path_converter(p)] * len(chunk) for p, chunk in zip(paths, chunks)
            )
        else:
            first, chunks = b_out
            first_path = None
            flat_paths = (None,)
            path_dtype = None

        flat_chunks = flatten(chunks)
        if meta is None:
            meta = read_json_chunk(
                first,
                encoding,
                errors,
                engine,
                include_path_column,
                first_path,
                path_dtype,
                kwargs,
            )
        meta = make_meta(meta)
        parts = [
            delayed(read_json_chunk)(
                chunk,
                encoding,
                errors,
                engine,
                include_path_column,
                path,
                path_dtype,
                kwargs,
                meta=meta,
            )
            for chunk, path in zip_longest(flat_chunks, flat_paths)
        ]
    else:
        files = open_files(
            url_path,
            "rt",
            encoding=encoding,
            errors=errors,
            compression=compression,
            **storage_options,
        )
        path_dtype = pd.CategoricalDtype(path_converter(f.path) for f in files)
        parts = [
            delayed(read_json_file)(
                f,
                orient,
                lines,
                engine,
                include_path_column,
                path_converter(f.path),
                path_dtype,
                kwargs,
            )
            for f in files
        ]

    return dd.from_delayed(parts, meta=meta)


def read_json_chunk(
    chunk, encoding, errors, engine, column_name, path, path_dtype, kwargs, meta=None
):
    s = io.StringIO(chunk.decode(encoding, errors))
    s.seek(0)
    df = engine(s, orient="records", lines=True, **kwargs)
    if meta is not None and df.empty:
        return meta

    if column_name:
        df = add_path_column(df, column_name, path, path_dtype)

    return df


def read_json_file(f, orient, lines, engine, column_name, path, path_dtype, kwargs):
    with f as open_file:
        df = engine(open_file, orient=orient, lines=lines, **kwargs)
    if column_name:
        df = add_path_column(df, column_name, path, path_dtype)
    return df


def add_path_column(df, column_name, path, dtype):
    if column_name in df.columns:
        raise ValueError(
            f"Files already contain the column name: '{column_name}', so the path "
            "column cannot use this name. Please set `include_path_column` to a "
            "unique name."
        )
    return df.assign(**{column_name: pd.Series([path] * len(df), dtype=dtype)})
