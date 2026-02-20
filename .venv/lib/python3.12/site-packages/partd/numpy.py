""" Store arrays

We put arrays on disk as raw bytes, extending along the first dimension.
Alongside each array x we ensure the value x.dtype which stores the string
description of the array's dtype.
"""
from contextlib import suppress
import pickle

import numpy as np
from toolz import valmap, identity, partial
from .core import Interface
from .file import File
from .utils import frame, framesplit, suffix


def serialize_dtype(dt):
    """ Serialize dtype to bytes

    >>> serialize_dtype(np.dtype('i4'))
    b'<i4'
    >>> serialize_dtype(np.dtype('M8[us]'))
    b'<M8[us]'
    """
    return dt.str.encode()


def parse_dtype(s):
    """ Parse text as numpy dtype

    >>> parse_dtype(b'i4')
    dtype('int32')

    >>> parse_dtype(b"[('a', 'i4')]")
    dtype([('a', '<i4')])
    """
    if s.startswith(b'['):
        return np.dtype(eval(s))  # Dangerous!
    else:
        return np.dtype(s)


class Numpy(Interface):
    def __init__(self, partd=None):
        if not partd or isinstance(partd, str):
            partd = File(partd)
        self.partd = partd
        Interface.__init__(self)

    def __getstate__(self):
        return {'partd': self.partd}

    def append(self, data, **kwargs):
        for k, v in data.items():
            self.partd.iset(suffix(k, '.dtype'), serialize_dtype(v.dtype))
        self.partd.append(valmap(serialize, data), **kwargs)

    def _get(self, keys, **kwargs):
        bytes = self.partd._get(keys, **kwargs)
        dtypes = self.partd._get([suffix(key, '.dtype') for key in keys],
                                 lock=False)
        dtypes = map(parse_dtype, dtypes)
        return list(map(deserialize, bytes, dtypes))

    def delete(self, keys, **kwargs):
        keys2 = [suffix(key, '.dtype') for key in keys]
        self.partd.delete(keys2, **kwargs)

    def _iset(self, key, value):
        return self.partd._iset(key, value)

    def drop(self):
        return self.partd.drop()

    def __del__(self):
        self.partd.__del__()

    @property
    def lock(self):
        return self.partd.lock

    def __exit__(self, *args):
        self.drop()
        self.partd.__exit__(self, *args)

try:
    from pandas import msgpack
except ImportError:
    try:
        import msgpack
    except ImportError:
        msgpack = False


def serialize(x):
    if x.dtype == 'O':
        l = x.flatten().tolist()
        with suppress(Exception):  # Try msgpack (faster on strings)
            return frame(msgpack.packb(l, use_bin_type=True))
        return frame(pickle.dumps(l, protocol=pickle.HIGHEST_PROTOCOL))
    else:
        return x.tobytes()


def deserialize(bytes, dtype, copy=False):
    if dtype == 'O':
        try:
            if msgpack.version >= (0, 5, 2):
                unpack_kwargs = {'raw': False}
            else:
                unpack_kwargs = {'encoding': 'utf-8'}

            blocks = [msgpack.unpackb(f, **unpack_kwargs)
                      for f in framesplit(bytes)]
        except Exception:
            blocks = [pickle.loads(f) for f in framesplit(bytes)]

        result = np.empty(sum(map(len, blocks)), dtype='O')
        i = 0
        for block in blocks:
            result[i:i + len(block)] = block
            i += len(block)
        return result
    else:
        result = np.frombuffer(bytes, dtype)
        if copy:
            result = result.copy()
        return result


compress_text = identity
decompress_text = identity
compress_bytes = lambda bytes, itemsize: bytes
decompress_bytes = identity

with suppress(ImportError):
    import blosc
    blosc.set_nthreads(1)

    compress_bytes = blosc.compress
    decompress_bytes = blosc.decompress

    compress_text = partial(blosc.compress, typesize=1)
    decompress_text = blosc.decompress

with suppress(ImportError):
    from snappy import compress as compress_text
    from snappy import decompress as decompress_text


def compress(bytes, dtype):
    if dtype == 'O':
        return compress_text(bytes)
    else:
        return compress_bytes(bytes, dtype.itemsize)


def decompress(bytes, dtype):
    if dtype == 'O':
        return decompress_text(bytes)
    else:
        return decompress_bytes(bytes)
