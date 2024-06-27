from __future__ import annotations

import copy
import os

from fsspec.core import OpenFile, get_fs_token_paths
from fsspec.utils import infer_compression, read_block

from dask.base import tokenize
from dask.delayed import delayed
from dask.utils import is_integer, parse_bytes


def read_bytes(
    urlpath,
    delimiter=None,
    not_zero=False,
    blocksize="128 MiB",
    sample="10 kiB",
    compression=None,
    include_path=False,
    **kwargs,
):
    """Given a path or paths, return delayed objects that read from those paths.

    The path may be a filename like ``'2015-01-01.csv'`` or a globstring
    like ``'2015-*-*.csv'``.

    The path may be preceded by a protocol, like ``s3://`` or ``hdfs://`` if
    those libraries are installed.

    This cleanly breaks data by a delimiter if given, so that block boundaries
    start directly after a delimiter and end on the delimiter.

    Parameters
    ----------
    urlpath : string or list
        Absolute or relative filepath(s). Prefix with a protocol like ``s3://``
        to read from alternative filesystems. To read from multiple files you
        can pass a globstring or a list of paths, with the caveat that they
        must all have the same protocol.
    delimiter : bytes
        An optional delimiter, like ``b'\\n'`` on which to split blocks of
        bytes.
    not_zero : bool
        Force seek of start-of-file delimiter, discarding header.
    blocksize : int, str
        Chunk size in bytes, defaults to "128 MiB"
    compression : string or None
        String like 'gzip' or 'xz'.  Must support efficient random access.
    sample : int, string, or boolean
        Whether or not to return a header sample.
        Values can be ``False`` for "no sample requested"
        Or an integer or string value like ``2**20`` or ``"1 MiB"``
    include_path : bool
        Whether or not to include the path with the bytes representing a particular file.
        Default is False.
    **kwargs : dict
        Extra options that make sense to a particular storage connection, e.g.
        host, port, username, password, etc.

    Examples
    --------
    >>> sample, blocks = read_bytes('2015-*-*.csv', delimiter=b'\\n')  # doctest: +SKIP
    >>> sample, blocks = read_bytes('s3://bucket/2015-*-*.csv', delimiter=b'\\n')  # doctest: +SKIP
    >>> sample, paths, blocks = read_bytes('2015-*-*.csv', include_path=True)  # doctest: +SKIP

    Returns
    -------
    sample : bytes
        The sample header
    blocks : list of lists of ``dask.Delayed``
        Each list corresponds to a file, and each delayed object computes to a
        block of bytes from that file.
    paths : list of strings, only included if include_path is True
        List of same length as blocks, where each item is the path to the file
        represented in the corresponding block.

    """
    if not isinstance(urlpath, (str, list, tuple, os.PathLike)):
        raise TypeError("Path should be a string, os.PathLike, list or tuple")

    fs, fs_token, paths = get_fs_token_paths(urlpath, mode="rb", storage_options=kwargs)

    if len(paths) == 0:
        raise OSError("%s resolved to no files" % urlpath)

    if blocksize is not None:
        if isinstance(blocksize, str):
            blocksize = parse_bytes(blocksize)
        if not is_integer(blocksize):
            raise TypeError("blocksize must be an integer")
        blocksize = int(blocksize)

    if blocksize is None:
        offsets = [[0]] * len(paths)
        lengths = [[None]] * len(paths)
    else:
        offsets = []
        lengths = []
        for path in paths:
            if compression == "infer":
                comp = infer_compression(path)
            else:
                comp = compression
            if comp is not None:
                raise ValueError(
                    "Cannot do chunked reads on compressed files. "
                    "To read, set blocksize=None"
                )
            size = fs.info(path)["size"]
            if size is None:
                raise ValueError(
                    "Backing filesystem couldn't determine file size, cannot "
                    "do chunked reads. To read, set blocksize=None."
                )

            elif size == 0:
                # skip empty
                offsets.append([])
                lengths.append([])
            else:
                # shrink blocksize to give same number of parts
                if size % blocksize and size > blocksize:
                    blocksize1 = size / (size // blocksize)
                else:
                    blocksize1 = blocksize
                place = 0
                off = [0]
                length = []

                # figure out offsets, spreading around spare bytes
                while size - place > (blocksize1 * 2) - 1:
                    place += blocksize1
                    off.append(int(place))
                    length.append(off[-1] - off[-2])
                length.append(size - off[-1])

                if not_zero:
                    off[0] = 1
                    length[0] -= 1
                offsets.append(off)
                lengths.append(length)

    delayed_read = delayed(read_block_from_file)

    out = []
    for path, offset, length in zip(paths, offsets, lengths):
        token = tokenize(fs_token, delimiter, path, fs.ukey(path), compression, offset)
        keys = [f"read-block-{o}-{token}" for o in offset]
        values = [
            delayed_read(
                OpenFile(fs, path, compression=compression),
                o,
                l,
                delimiter,
                dask_key_name=key,
            )
            for o, key, l in zip(offset, keys, length)
        ]
        out.append(values)

    if sample:
        if sample is True:
            sample = "10 kiB"  # backwards compatibility
        if isinstance(sample, str):
            sample = parse_bytes(sample)
        with OpenFile(fs, paths[0], compression=compression) as f:
            # read block without seek (because we start at zero)
            if delimiter is None:
                sample = f.read(sample)
            else:
                sample_buff = f.read(sample)
                while True:
                    new = f.read(sample)
                    if not new:
                        break
                    if delimiter in new:
                        sample_buff = (
                            sample_buff + new.split(delimiter, 1)[0] + delimiter
                        )
                        break
                    sample_buff = sample_buff + new
                sample = sample_buff
    if include_path:
        return sample, out, paths
    return sample, out


def read_block_from_file(lazy_file, off, bs, delimiter):
    with copy.copy(lazy_file) as f:
        if off == 0 and bs is None:
            return f.read()
        return read_block(f, off, bs, delimiter)
