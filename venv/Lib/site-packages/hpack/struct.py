"""
Contains structures for representing header fields with associated metadata.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from typing_extensions import Self, TypeAlias  # noqa: UP035, pragma: no cover


class HeaderTuple(tuple[bytes, bytes]):
    """
    A data structure that stores a single header field.

    HTTP headers can be thought of as tuples of ``(field name, field value)``.
    A single header block is a sequence of such tuples.

    In HTTP/2, however, certain bits of additional information are required for
    compressing these headers: in particular, whether the header field can be
    safely added to the HPACK compression context.

    This class stores a header that can be added to the compression context. In
    all other ways it behaves exactly like a tuple.
    """

    __slots__ = ()

    indexable = True

    def __new__(cls, *args: Any) -> Self:
        return tuple.__new__(cls, args)


class NeverIndexedHeaderTuple(HeaderTuple):
    """
    A data structure that stores a single header field that cannot be added to
    a HTTP/2 header compression context.
    """

    __slots__ = ()

    indexable = False

    def __new__(cls, *args: Any) -> Self:
        return tuple.__new__(cls, args)


Header: TypeAlias = "HeaderTuple | NeverIndexedHeaderTuple | tuple[bytes, bytes]"
HeaderWeaklyTyped: TypeAlias = "HeaderTuple | NeverIndexedHeaderTuple | tuple[bytes | str, bytes | str]"
