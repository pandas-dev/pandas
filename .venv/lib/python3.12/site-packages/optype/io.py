import sys
from typing import Literal, Protocol, TypeAlias as Alias

if sys.version_info >= (3, 13):
    from typing import TypeVar, runtime_checkable
else:
    from typing_extensions import TypeVar, runtime_checkable

__all__ = (
    "CanFSPath", "CanFileno",
    "CanFlush", "CanRead", "CanReadN", "CanReadline", "CanReadlineN", "CanWrite",
    "ModeAB", "ModeABU", "ModeAB_",
    "ModeAT", "ModeATU", "ModeAT_",
    "ModeOB", "ModeOBU", "ModeOB_",
    "ModeOT", "ModeOTU", "ModeOT_",
    "ModeRB", "ModeRBU", "ModeRB_",
    "ModeRT", "ModeRTU", "ModeRT_",
    "ModeWB", "ModeWBU", "ModeWB_",
    "ModeWT", "ModeWTU", "ModeWT_",
    "ModeXB", "ModeXBU", "ModeXB_",
    "ModeXT", "ModeXTU", "ModeXT_",
    "Mode_B", "Mode_BU", "Mode_B_",
    "Mode_T", "Mode_TU", "Mode_T_",
    "ToFileno", "ToPath",
)  # fmt: skip


def __dir__() -> tuple[str, ...]:
    return __all__


###

# not a type parameter
_StrOrBytes = TypeVar("_StrOrBytes", str, bytes, str | bytes, default=str | bytes)

_T_co = TypeVar("_T_co", covariant=True)
_T_contra = TypeVar("_T_contra", contravariant=True)
_RT_co = TypeVar("_RT_co", default=object, covariant=True)

_PathT_co = TypeVar("_PathT_co", bound=str | bytes, default=str | bytes, covariant=True)

###


@runtime_checkable
class CanFSPath(Protocol[_PathT_co]):
    """
    Similar to `os.PathLike`, but is is actually a protocol, doesn't incorrectly use a
    `TypeVar` with constraints", and therefore doesn't force type-unsafe usage of `Any`
    to express `str | bytes`.
    """

    def __fspath__(self, /) -> _PathT_co: ...


@runtime_checkable
class CanFileno(Protocol):
    """Runtime-checkable equivalent of `_typeshed.HasFileno`."""

    def fileno(self, /) -> int: ...


@runtime_checkable
class CanRead(Protocol[_T_co]):
    """
    Like `_typeshed.SupportsRead`, but without the required positional `int` argument,
    and is runtime-checkable.
    """

    def read(self, /) -> _T_co: ...


@runtime_checkable
class CanReadN(Protocol[_T_co]):
    """Runtime-checkable equivalent of `_typeshed.SupportsRead`."""

    def read(self, n: int = ..., /) -> _T_co: ...


@runtime_checkable
class CanReadline(Protocol[_T_co]):
    """
    Runtime-checkable equivalent of `_typeshed.SupportsNoArgReadline`, that
    additionally allows `self` to be positional-only.
    """

    def readline(self, /) -> _T_co: ...


@runtime_checkable
class CanReadlineN(Protocol[_T_co]):
    """Runtime-checkable equivalent of `_typeshed.SupportsReadline`."""

    def readline(self, n: int = ..., /) -> _T_co: ...


@runtime_checkable
class CanWrite(Protocol[_T_contra, _RT_co]):
    """
    Runtime-checkable equivalent of `_typeshed.SupportsWrite`, with an additional
    optional type parameter for the return type, that defaults to `object`.
    """

    def write(self, data: _T_contra, /) -> _RT_co: ...


@runtime_checkable
class CanFlush(Protocol[_RT_co]):
    """
    Runtime-checkable equivalent of `_typeshed.SupportsFlush`, with an additional
    optional type parameter for the return type, that defaults to `object`.
    """

    def flush(self, /) -> _RT_co: ...


###

# runtime-checkable `_typeshed.{Str,Bytes,StrOrBytes,Generic}Path` alternative
ToPath: Alias = _StrOrBytes | CanFSPath[_StrOrBytes]

# runtime-checkable `_typeshed.FileDescriptorLike` equivalent
ToFileno: Alias = int | CanFileno


###

# Literals for `builtins.open` file modes:
# https://docs.python.org/3/library/functions.html#open
#
# The names follow the pattern `Mode(R|W|X|A|O|_)(B|T|_)(U|_)?`, where each group
# matches a specific set of characters according to the following mapping:
#
# 1.
#   - R -> 'r'  (read)
#   - W -> 'w'  (write, truncate)
#   - X -> 'x'  (write, exclusive creation)
#   - A -> 'a'  (write, append)
#   - O -> 'w' | 'x' | 'a'  (write)
#   - _ -> 'r' | 'w' | 'x' | 'a'  (any)
# 2.
#   - B -> 'b'
#   - T -> 't' | ''
# 3.
#   - U -> '+'
#   - _ -> '+' | ''

ModeRB: Alias = Literal["rb", "br"]  # _typeshed.OpenBinaryModeReading (minus 'U')
ModeRBU: Alias = Literal["rb+", "r+b", "+rb", "br+", "b+r", "+br"]
ModeRB_: Alias = Literal[ModeRB, ModeRBU]
ModeRT: Alias = Literal["r", "rt", "tr"]  # _typeshed.OpenTextModeReading (minus 'U')
ModeRTU: Alias = Literal["r+", "+r", "rt+", "r+t", "+rt", "tr+", "t+r", "+tr"]
ModeRT_: Alias = Literal[ModeRT, ModeRTU]

ModeWB: Alias = Literal["wb", "bw"]
ModeWBU: Alias = Literal["wb+", "w+b", "+wb", "bw+", "b+w", "+bw"]
ModeWB_: Alias = Literal[ModeWB, ModeWBU]
ModeWT: Alias = Literal["w", "wt", "tw"]
ModeWTU: Alias = Literal["w+", "+w", "wt+", "w+t", "+wt", "tw+", "t+w", "+tw"]
ModeWT_: Alias = Literal[ModeWT, ModeWTU]

ModeXB: Alias = Literal["xb", "bx"]
ModeXBU: Alias = Literal["xb+", "x+b", "+xb", "bx+", "b+x", "+bx"]
ModeXB_: Alias = Literal[ModeXB, ModeXBU]
ModeXT: Alias = Literal["x", "xt", "tx"]
ModeXTU: Alias = Literal["x+", "+x", "xt+", "x+t", "+xt", "tx+", "t+x", "+tx"]
ModeXT_: Alias = Literal[ModeXT, ModeXTU]

ModeAB: Alias = Literal["ab", "ba"]
ModeABU: Alias = Literal["ab+", "a+b", "+ab", "ba+", "b+a", "+ba"]
ModeAB_: Alias = Literal[ModeAB, ModeABU]
ModeAT: Alias = Literal["a", "at", "ta"]
ModeATU: Alias = Literal["a+", "+a", "at+", "a+t", "+at", "ta+", "t+a", "+ta"]
ModeAT_: Alias = Literal[ModeAT, ModeATU]

ModeOB: Alias = Literal[ModeWB, ModeXB, ModeAB]  # typeshed.OpenBinaryModeWriting
ModeOBU: Alias = Literal[ModeWBU, ModeXBU, ModeABU]
ModeOB_: Alias = Literal[ModeWB_, ModeXB_, ModeAB_]
ModeOT: Alias = Literal[ModeWT, ModeXT, ModeAT]  # _typeshed.OpenTextModeWriting
ModeOTU: Alias = Literal[ModeWTU, ModeXTU, ModeATU]
ModeOT_: Alias = Literal[ModeWT_, ModeXT_, ModeAT_]

Mode_B: Alias = Literal[ModeRB, ModeOB]  # _typeshed.OpenBinaryModeUpdating
Mode_BU: Alias = Literal[ModeRBU, ModeOBU]  # _typeshed.OpenBinaryMode  (minus 'U')
Mode_B_: Alias = Literal[ModeRB_, ModeOB_]
Mode_T: Alias = Literal[ModeRT, ModeOT]
Mode_TU: Alias = Literal[ModeRTU, ModeOTU]  # _typeshed.OpenTextModeUpdating
Mode_T_: Alias = Literal[ModeRT_, ModeOT_]  # _typeshed.OpenTextMode  (minus 'U')
