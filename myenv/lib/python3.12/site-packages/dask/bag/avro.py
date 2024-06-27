from __future__ import annotations

import io
import uuid

from fsspec.core import OpenFile, get_fs_token_paths, open_files
from fsspec.utils import read_block
from fsspec.utils import tokenize as fs_tokenize

from dask.highlevelgraph import HighLevelGraph

MAGIC = b"Obj\x01"
SYNC_SIZE = 16


def read_long(fo):
    """variable-length, zig-zag encoding."""
    c = fo.read(1)
    b = ord(c)
    n = b & 0x7F
    shift = 7
    while (b & 0x80) != 0:
        b = ord(fo.read(1))
        n |= (b & 0x7F) << shift
        shift += 7
    return (n >> 1) ^ -(n & 1)


def read_bytes(fo):
    """a long followed by that many bytes of data."""
    size = read_long(fo)
    return fo.read(size)


def read_header(fo):
    """Extract an avro file's header

    fo: file-like
        This should be in bytes mode, e.g., io.BytesIO

    Returns dict representing the header

    Parameters
    ----------
    fo: file-like
    """
    assert fo.read(len(MAGIC)) == MAGIC, "Magic avro bytes missing"
    meta = {}
    out = {"meta": meta}
    while True:
        n_keys = read_long(fo)
        if n_keys == 0:
            break
        for _ in range(n_keys):
            # ignore dtype mapping for bag version
            read_bytes(fo)  # schema keys
            read_bytes(fo)  # schema values
    out["sync"] = fo.read(SYNC_SIZE)
    out["header_size"] = fo.tell()
    fo.seek(0)
    out["head_bytes"] = fo.read(out["header_size"])
    return out


def open_head(fs, path, compression):
    """Open a file just to read its head and size"""
    with OpenFile(fs, path, compression=compression) as f:
        head = read_header(f)
    size = fs.info(path)["size"]
    return head, size


def read_avro(urlpath, blocksize=100000000, storage_options=None, compression=None):
    """Read set of avro files

    Use this with arbitrary nested avro schemas. Please refer to the
    fastavro documentation for its capabilities:
    https://github.com/fastavro/fastavro

    Parameters
    ----------
    urlpath: string or list
        Absolute or relative filepath, URL (may include protocols like
        ``s3://``), or globstring pointing to data.
    blocksize: int or None
        Size of chunks in bytes. If None, there will be no chunking and each
        file will become one partition.
    storage_options: dict or None
        passed to backend file-system
    compression: str or None
        Compression format of the targe(s), like 'gzip'. Should only be used
        with blocksize=None.
    """
    from dask import compute, delayed
    from dask.bag import from_delayed
    from dask.utils import import_required

    import_required(
        "fastavro", "fastavro is a required dependency for using bag.read_avro()."
    )

    storage_options = storage_options or {}
    if blocksize is not None:
        fs, fs_token, paths = get_fs_token_paths(
            urlpath, mode="rb", storage_options=storage_options
        )
        dhead = delayed(open_head)
        out = compute(*[dhead(fs, path, compression) for path in paths])
        heads, sizes = zip(*out)
        dread = delayed(read_chunk)

        offsets = []
        lengths = []
        for size in sizes:
            off = list(range(0, size, blocksize))
            length = [blocksize] * len(off)
            offsets.append(off)
            lengths.append(length)

        out = []
        for path, offset, length, head in zip(paths, offsets, lengths, heads):
            delimiter = head["sync"]
            f = OpenFile(fs, path, compression=compression)
            token = fs_tokenize(
                fs_token, delimiter, path, fs.ukey(path), compression, offset
            )
            keys = [f"read-avro-{o}-{token}" for o in offset]
            values = [
                dread(f, o, l, head, dask_key_name=key)
                for o, key, l in zip(offset, keys, length)
            ]
            out.extend(values)

        return from_delayed(out)
    else:
        files = open_files(urlpath, compression=compression, **storage_options)
        dread = delayed(read_file)
        chunks = [dread(fo) for fo in files]
        return from_delayed(chunks)


def read_chunk(fobj, off, l, head):
    """Get rows from raw bytes block"""
    import fastavro

    if hasattr(fastavro, "iter_avro"):
        reader = fastavro.iter_avro
    else:
        reader = fastavro.reader

    with fobj as f:
        chunk = read_block(f, off, l, head["sync"])
    head_bytes = head["head_bytes"]
    if not chunk.startswith(MAGIC):
        chunk = head_bytes + chunk
    i = io.BytesIO(chunk)
    return list(reader(i))


def read_file(fo):
    """Get rows from file-like"""
    import fastavro

    if hasattr(fastavro, "iter_avro"):
        reader = fastavro.iter_avro
    else:
        reader = fastavro.reader

    with fo as f:
        return list(reader(f))


def to_avro(
    b,
    filename,
    schema,
    name_function=None,
    storage_options=None,
    codec="null",
    sync_interval=16000,
    metadata=None,
    compute=True,
    **kwargs,
):
    """Write bag to set of avro files

    The schema is a complex dictionary describing the data, see
    https://avro.apache.org/docs/1.8.2/gettingstartedpython.html#Defining+a+schema
    and https://fastavro.readthedocs.io/en/latest/writer.html .
    It's structure is as follows::

        {'name': 'Test',
         'namespace': 'Test',
         'doc': 'Descriptive text',
         'type': 'record',
         'fields': [
            {'name': 'a', 'type': 'int'},
         ]}

    where the "name" field is required, but "namespace" and "doc" are optional
    descriptors; "type" must always be "record". The list of fields should
    have an entry for every key of the input records, and the types are
    like the primitive, complex or logical types of the Avro spec
    ( https://avro.apache.org/docs/1.8.2/spec.html ).

    Results in one avro file per input partition.

    Parameters
    ----------
    b: dask.bag.Bag
    filename: list of str or str
        Filenames to write to. If a list, number must match the number of
        partitions. If a string, must include a glob character "*", which will
        be expanded using name_function
    schema: dict
        Avro schema dictionary, see above
    name_function: None or callable
        Expands integers into strings, see
        ``dask.bytes.utils.build_name_function``
    storage_options: None or dict
        Extra key/value options to pass to the backend file-system
    codec: 'null', 'deflate', or 'snappy'
        Compression algorithm
    sync_interval: int
        Number of records to include in each block within a file
    metadata: None or dict
        Included in the file header
    compute: bool
        If True, files are written immediately, and function blocks. If False,
        returns delayed objects, which can be computed by the user where
        convenient.
    kwargs: passed to compute(), if compute=True

    Examples
    --------
    >>> import dask.bag as db
    >>> b = db.from_sequence([{'name': 'Alice', 'value': 100},
    ...                       {'name': 'Bob', 'value': 200}])
    >>> schema = {'name': 'People', 'doc': "Set of people's scores",
    ...           'type': 'record',
    ...           'fields': [
    ...               {'name': 'name', 'type': 'string'},
    ...               {'name': 'value', 'type': 'int'}]}
    >>> b.to_avro('my-data.*.avro', schema)  # doctest: +SKIP
    ['my-data.0.avro', 'my-data.1.avro']
    """
    # TODO infer schema from first partition of data
    from dask.utils import import_required

    import_required(
        "fastavro", "fastavro is a required dependency for using bag.to_avro()."
    )
    _verify_schema(schema)

    storage_options = storage_options or {}
    files = open_files(
        filename,
        "wb",
        name_function=name_function,
        num=b.npartitions,
        **storage_options,
    )
    name = "to-avro-" + uuid.uuid4().hex
    dsk = {
        (name, i): (
            _write_avro_part,
            (b.name, i),
            f,
            schema,
            codec,
            sync_interval,
            metadata,
        )
        for i, f in enumerate(files)
    }
    graph = HighLevelGraph.from_collections(name, dsk, dependencies=[b])
    out = type(b)(graph, name, b.npartitions)
    if compute:
        out.compute(**kwargs)
        return [f.path for f in files]
    else:
        return out.to_delayed()


def _verify_schema(s):
    assert isinstance(s, dict), "Schema must be dictionary"
    for field in ["name", "type", "fields"]:
        assert field in s, "Schema missing '%s' field" % field
    assert s["type"] == "record", "Schema must be of type 'record'"
    assert isinstance(s["fields"], list), "Fields entry must be a list"
    for f in s["fields"]:
        assert "name" in f and "type" in f, "Field spec incomplete: %s" % f


def _write_avro_part(part, f, schema, codec, sync_interval, metadata):
    """Create single avro file from list of dictionaries"""
    import fastavro

    with f as f:
        fastavro.writer(f, schema, part, codec, sync_interval, metadata)
