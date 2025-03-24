from contextlib import contextmanager
import os
import shutil
import tempfile
import struct


def raises(exc, lamda):
    try:
        lamda()
        return False
    except exc:
        return True


@contextmanager
def tmpfile(extension=''):
    extension = '.' + extension.lstrip('.')
    handle, filename = tempfile.mkstemp(extension)
    os.close(handle)
    os.remove(filename)

    try:
        yield filename
    finally:
        if os.path.exists(filename):
            if os.path.isdir(filename):
                shutil.rmtree(filename)
            else:
                os.remove(filename)


def frame(bytes):
    """ Pack the length of the bytes in front of the bytes

    TODO: This does a full copy.  This should maybe be inlined somehow
    wherever this gets used instead.  My laptop shows a data bandwidth of
    2GB/s
    """
    return struct.pack('Q', len(bytes)) + bytes


def framesplit(bytes):
    """ Split buffer into frames of concatenated chunks

    >>> data = frame(b'Hello') + frame(b'World')
    >>> list(framesplit(data))  # doctest: +SKIP
    [b'Hello', b'World']
    """
    i = 0; n = len(bytes)
    chunks = list()
    while i < n:
        nbytes = struct.unpack('Q', bytes[i:i+8])[0]
        i += 8
        yield bytes[i: i + nbytes]
        i += nbytes


def partition_all(n, bytes):
    """ Partition bytes into evenly sized blocks

    The final block holds the remainder and so may not be of equal size

    >>> list(partition_all(2, b'Hello'))
    [b'He', b'll', b'o']

    See Also:
        toolz.partition_all
    """
    if len(bytes) < n:  # zero copy fast common case
        yield bytes
    else:
        for i in range(0, len(bytes), n):
            yield bytes[i: i+n]


def nested_get(ind, coll, lazy=False):
    """ Get nested index from collection

    Examples
    --------

    >>> nested_get(1, 'abc')
    'b'
    >>> nested_get([1, 0], 'abc')
    ['b', 'a']
    >>> nested_get([[1, 0], [0, 1]], 'abc')
    [['b', 'a'], ['a', 'b']]
    """
    if isinstance(ind, list):
        if lazy:
            return (nested_get(i, coll, lazy=lazy) for i in ind)
        else:
            return [nested_get(i, coll, lazy=lazy) for i in ind]
    else:
        return coll[ind]


def flatten(seq):
    """

    >>> list(flatten([1]))
    [1]

    >>> list(flatten([[1, 2], [1, 2]]))
    [1, 2, 1, 2]

    >>> list(flatten([[[1], [2]], [[1], [2]]]))
    [1, 2, 1, 2]

    >>> list(flatten(((1, 2), (1, 2)))) # Don't flatten tuples
    [(1, 2), (1, 2)]

    >>> list(flatten((1, 2, [3, 4]))) # support heterogeneous
    [1, 2, 3, 4]
    """
    for item in seq:
        if isinstance(item, list):
            yield from flatten(item)
        else:
            yield item


def suffix(key, term):
    """ suffix a key with a suffix

    Works if they key is a string or a tuple

    >>> suffix('x', '.dtype')
    'x.dtype'
    >>> suffix(('a', 'b', 'c'), '.dtype')
    ('a', 'b', 'c.dtype')
    """
    if isinstance(key, str):
        return key + term
    elif isinstance(key, tuple):
        return key[:-1] + (suffix(key[-1], term),)
    else:
        return suffix(str(key), term)


def extend(key, term):
    """ extend a key with a another element in a tuple

    Works if they key is a string or a tuple

    >>> extend('x', '.dtype')
    ('x', '.dtype')
    >>> extend(('a', 'b', 'c'), '.dtype')
    ('a', 'b', 'c', '.dtype')
    """
    if isinstance(term, tuple):
        pass
    elif isinstance(term, str):
        term = (term,)
    else:
        term = (str(term),)

    if not isinstance(key, tuple):
        key = (key,)

    return key + term
