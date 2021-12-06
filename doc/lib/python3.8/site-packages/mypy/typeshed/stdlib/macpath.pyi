import sys
from _typeshed import AnyPath, BytesPath, StrPath
from genericpath import (
    commonprefix as commonprefix,
    exists as exists,
    getatime as getatime,
    getctime as getctime,
    getmtime as getmtime,
    getsize as getsize,
    isdir as isdir,
    isfile as isfile,
)

if sys.version_info >= (3, 4):
    from genericpath import samefile as samefile, sameopenfile as sameopenfile, samestat as samestat

# Re-export common definitions from posixpath to reduce duplication
from posixpath import (
    abspath as abspath,
    curdir as curdir,
    defpath as defpath,
    devnull as devnull,
    expanduser as expanduser,
    expandvars as expandvars,
    extsep as extsep,
    isabs as isabs,
    lexists as lexists,
    pardir as pardir,
    pathsep as pathsep,
    sep as sep,
    splitdrive as splitdrive,
    splitext as splitext,
    supports_unicode_filenames as supports_unicode_filenames,
)
from typing import AnyStr, Optional, Text, Tuple, overload

altsep: Optional[str]

if sys.version_info >= (3, 6):
    from os import PathLike
    @overload
    def basename(s: PathLike[AnyStr]) -> AnyStr: ...
    @overload
    def basename(s: AnyStr) -> AnyStr: ...
    @overload
    def dirname(s: PathLike[AnyStr]) -> AnyStr: ...
    @overload
    def dirname(s: AnyStr) -> AnyStr: ...
    @overload
    def normcase(path: PathLike[AnyStr]) -> AnyStr: ...
    @overload
    def normcase(path: AnyStr) -> AnyStr: ...
    @overload
    def normpath(s: PathLike[AnyStr]) -> AnyStr: ...
    @overload
    def normpath(s: AnyStr) -> AnyStr: ...
    @overload
    def realpath(path: PathLike[AnyStr]) -> AnyStr: ...
    @overload
    def realpath(path: AnyStr) -> AnyStr: ...

else:
    def basename(s: AnyStr) -> AnyStr: ...
    def dirname(s: AnyStr) -> AnyStr: ...
    def normcase(path: AnyStr) -> AnyStr: ...
    def normpath(s: AnyStr) -> AnyStr: ...
    def realpath(path: AnyStr) -> AnyStr: ...

def islink(s: AnyPath) -> bool: ...

if sys.version_info >= (3, 6):
    # Mypy complains that the signatures overlap, but things seem to behave correctly anyway.
    @overload
    def join(s: StrPath, *paths: StrPath) -> Text: ...
    @overload
    def join(s: BytesPath, *paths: BytesPath) -> bytes: ...

elif sys.version_info >= (3, 0):
    def join(s: AnyStr, *paths: AnyStr) -> AnyStr: ...

else:
    # Make sure signatures are disjunct, and allow combinations of bytes and unicode.
    # (Since Python 2 allows that, too)
    # Note that e.g. os.path.join("a", "b", "c", "d", u"e") will still result in
    # a type error.
    @overload
    def join(__p1: bytes, *p: bytes) -> bytes: ...
    @overload
    def join(__p1: bytes, __p2: bytes, __p3: bytes, __p4: Text, *p: AnyPath) -> Text: ...
    @overload
    def join(__p1: bytes, __p2: bytes, __p3: Text, *p: AnyPath) -> Text: ...
    @overload
    def join(__p1: bytes, __p2: Text, *p: AnyPath) -> Text: ...
    @overload
    def join(__p1: Text, *p: AnyPath) -> Text: ...

if sys.version_info >= (3, 6):
    @overload
    def split(s: PathLike[AnyStr]) -> Tuple[AnyStr, AnyStr]: ...
    @overload
    def split(s: AnyStr) -> Tuple[AnyStr, AnyStr]: ...

else:
    def split(s: AnyStr) -> Tuple[AnyStr, AnyStr]: ...
