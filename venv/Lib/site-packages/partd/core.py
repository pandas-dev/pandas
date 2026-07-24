import os
import shutil
import locket
import string
from toolz import memoize
from contextlib import contextmanager
from .utils import nested_get, flatten



# http://stackoverflow.com/questions/295135/turn-a-string-into-a-valid-filename-in-python
valid_chars = "-_.() " + string.ascii_letters + string.digits + os.path.sep


def escape_filename(fn):
    """ Escape text so that it is a valid filename

    >>> escape_filename('Foo!bar?')
    'Foobar'

    """
    return ''.join(filter(valid_chars.__contains__, fn))


def filename(path, key):
    return os.path.join(path, escape_filename(token(key)))


def token(key):
    """

    >>> token('hello')
    'hello'
    >>> token(('hello', 'world'))  # doctest: +SKIP
    'hello/world'
    """
    if isinstance(key, str):
        return key
    elif isinstance(key, tuple):
        return os.path.join(*map(token, key))
    else:
        return str(key)


class Interface:
    def __init__(self):
        self._iset_seen = set()

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._iset_seen = set()

    def iset(self, key, value, **kwargs):
        if key in self._iset_seen:
            return
        else:
            self._iset(key, value, **kwargs)
            self._iset_seen.add(key)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.drop()

    def iget(self, key):
        return self._get([key], lock=False)[0]

    def get(self, keys, **kwargs):
        if not isinstance(keys, list):
            return self.get([keys], **kwargs)[0]
        elif any(isinstance(key, list) for key in keys):  # nested case
            flatkeys = list(flatten(keys))
            result = self.get(flatkeys, **kwargs)
            return nested_get(keys, dict(zip(flatkeys, result)))
        else:
            return self._get(keys, **kwargs)

    def delete(self, keys, **kwargs):
        if not isinstance(keys, list):
            return self._delete([keys], **kwargs)
        else:
            return self._delete(keys, **kwargs)

    def pop(self, keys, **kwargs):
        with self.partd.lock:
            result = self.partd.get(keys, lock=False)
            self.partd.delete(keys, lock=False)
        return result

