from contextlib import suppress
from functools import partial

from .encode import Encode

__all__ = []


def bytes_concat(L):
    return b''.join(L)


with suppress(ImportError, AttributeError):
    # In case snappy is not installed, or another package called snappy that does not implement compress / decompress.
    # For example, SnapPy (https://pypi.org/project/snappy/)
    import snappy
    Snappy = partial(Encode,
                     snappy.compress,
                     snappy.decompress,
                     bytes_concat)
    __all__.append('Snappy')


with suppress(ImportError):
    import zlib
    ZLib = partial(Encode,
                   zlib.compress,
                   zlib.decompress,
                   bytes_concat)
    __all__.append('ZLib')


with suppress(ImportError):
    import bz2
    BZ2 = partial(Encode,
                  bz2.compress,
                  bz2.decompress,
                  bytes_concat)
    __all__.append('BZ2')


with suppress(ImportError):
    import blosc
    Blosc = partial(Encode,
                    blosc.compress,
                    blosc.decompress,
                    bytes_concat)
    __all__.append('Blosc')
