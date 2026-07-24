import sys
from typing import Any, Literal, Protocol, TypeAliasType

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

_T_contra = TypeVar("_T_contra", contravariant=True)
_RT_co = TypeVar("_RT_co", default=Any, covariant=True)

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
class CanRead[T_co](Protocol):
    """
    Like `_typeshed.SupportsRead`, but without the required positional `int` argument,
    and is runtime-checkable.
    """

    def read(self, /) -> T_co: ...


@runtime_checkable
class CanReadN[T_co](Protocol):
    """Runtime-checkable equivalent of `_typeshed.SupportsRead`."""

    def read(self, n: int = ..., /) -> T_co: ...


@runtime_checkable
class CanReadline[T_co](Protocol):
    """
    Runtime-checkable equivalent of `_typeshed.SupportsNoArgReadline`, that
    additionally allows `self` to be positional-only.
    """

    def readline(self, /) -> T_co: ...


@runtime_checkable
class CanReadlineN[T_co](Protocol):
    """Runtime-checkable equivalent of `_typeshed.SupportsReadline`."""

    def readline(self, n: int = ..., /) -> T_co: ...


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
ToPath = TypeAliasType(  # noqa: UP040
    "ToPath",
    _StrOrBytes | CanFSPath[_StrOrBytes],
    type_params=(_StrOrBytes,),
)

# runtime-checkable `_typeshed.FileDescriptorLike` equivalent
type ToFileno = int | CanFileno


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

type ModeRB = Literal["rb", "br"]  # _typeshed.OpenBinaryModeReading (minus 'U')
type ModeRBU = Literal["rb+", "r+b", "+rb", "br+", "b+r", "+br"]
type ModeRB_ = Literal[ModeRB, ModeRBU]
type ModeRT = Literal["r", "rt", "tr"]  # _typeshed.OpenTextModeReading (minus 'U')
type ModeRTU = Literal["r+", "+r", "rt+", "r+t", "+rt", "tr+", "t+r", "+tr"]
type ModeRT_ = Literal[ModeRT, ModeRTU]

type ModeWB = Literal["wb", "bw"]
type ModeWBU = Literal["wb+", "w+b", "+wb", "bw+", "b+w", "+bw"]
type ModeWB_ = Literal[ModeWB, ModeWBU]
type ModeWT = Literal["w", "wt", "tw"]
type ModeWTU = Literal["w+", "+w", "wt+", "w+t", "+wt", "tw+", "t+w", "+tw"]
type ModeWT_ = Literal[ModeWT, ModeWTU]

type ModeXB = Literal["xb", "bx"]
type ModeXBU = Literal["xb+", "x+b", "+xb", "bx+", "b+x", "+bx"]
type ModeXB_ = Literal[ModeXB, ModeXBU]
type ModeXT = Literal["x", "xt", "tx"]
type ModeXTU = Literal["x+", "+x", "xt+", "x+t", "+xt", "tx+", "t+x", "+tx"]
type ModeXT_ = Literal[ModeXT, ModeXTU]

type ModeAB = Literal["ab", "ba"]
type ModeABU = Literal["ab+", "a+b", "+ab", "ba+", "b+a", "+ba"]
type ModeAB_ = Literal[ModeAB, ModeABU]
type ModeAT = Literal["a", "at", "ta"]
type ModeATU = Literal["a+", "+a", "at+", "a+t", "+at", "ta+", "t+a", "+ta"]
type ModeAT_ = Literal[ModeAT, ModeATU]

type ModeOB = Literal[ModeWB, ModeXB, ModeAB]  # typeshed.OpenBinaryModeWriting
type ModeOBU = Literal[ModeWBU, ModeXBU, ModeABU]
type ModeOB_ = Literal[ModeWB_, ModeXB_, ModeAB_]
type ModeOT = Literal[ModeWT, ModeXT, ModeAT]  # _typeshed.OpenTextModeWriting
type ModeOTU = Literal[ModeWTU, ModeXTU, ModeATU]
type ModeOT_ = Literal[ModeWT_, ModeXT_, ModeAT_]

type Mode_B = Literal[ModeRB, ModeOB]  # _typeshed.OpenBinaryModeUpdating
type Mode_BU = Literal[ModeRBU, ModeOBU]  # _typeshed.OpenBinaryMode  (minus 'U')
type Mode_B_ = Literal[ModeRB_, ModeOB_]
type Mode_T = Literal[ModeRT, ModeOT]
type Mode_TU = Literal[ModeRTU, ModeOTU]  # _typeshed.OpenTextModeUpdating
type Mode_T_ = Literal[ModeRT_, ModeOT_]  # _typeshed.OpenTextMode  (minus 'U')
