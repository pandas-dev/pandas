from _typeshed import StrPath
from collections.abc import Iterable
from typing import overload
from typing_extensions import Self

from unidiff.constants import (
    DEFAULT_ENCODING as DEFAULT_ENCODING,
    LINE_TYPE_ADDED as LINE_TYPE_ADDED,
    LINE_TYPE_CONTEXT as LINE_TYPE_CONTEXT,
    LINE_TYPE_REMOVED as LINE_TYPE_REMOVED,
)
from unidiff.errors import UnidiffParseError as UnidiffParseError

class Line:
    source_line_no: int | None
    target_line_no: int | None
    diff_line_no: int | None
    line_type: str
    value: str
    def __init__(
        self,
        value: str,
        line_type: str,
        source_line_no: int | None = None,
        target_line_no: int | None = None,
        diff_line_no: int | None = None,
    ) -> None: ...
    def __eq__(self, other: Line) -> bool: ...  # type: ignore[override]
    @property
    def is_added(self) -> bool: ...
    @property
    def is_removed(self) -> bool: ...
    @property
    def is_context(self) -> bool: ...

class PatchInfo(list[str]): ...

class Hunk(list[Line]):
    source_start: int
    source_length: int
    target_start: int
    target_length: int
    section_header: str
    def __init__(
        self, src_start: int = 0, src_len: int | None = 0, tgt_start: int = 0, tgt_len: int | None = 0, section_header: str = ""
    ) -> None: ...
    def append(self, line: Line) -> None: ...
    @property
    def added(self) -> int: ...
    @property
    def removed(self) -> int: ...
    def is_valid(self) -> bool: ...
    def source_lines(self) -> Iterable[Line]: ...
    @property
    def source(self) -> Iterable[str]: ...
    def target_lines(self) -> Iterable[Line]: ...
    @property
    def target(self) -> Iterable[str]: ...

class PatchedFile(list[Hunk]):
    patch_info: PatchInfo | None
    source_file: str
    source_timestamp: str | None
    target_file: str
    target_timestamp: str | None
    is_binary_file: bool
    def __init__(
        self,
        patch_info: PatchInfo | None = None,
        source: str = "",
        target: str = "",
        source_timestamp: str | None = None,
        target_timestamp: str | None = None,
        is_binary_file: bool = False,
    ) -> None: ...
    @property
    def path(self) -> str: ...
    @property
    def added(self) -> int: ...
    @property
    def removed(self) -> int: ...
    @property
    def is_rename(self) -> bool: ...
    @property
    def is_added_file(self) -> bool: ...
    @property
    def is_removed_file(self) -> bool: ...
    @property
    def is_modified_file(self) -> bool: ...

class PatchSet(list[PatchedFile]):
    @overload
    def __init__(self, f: Iterable[str] | str, encoding: None = None, metadata_only: bool = False) -> None: ...
    @overload
    def __init__(self, f: Iterable[bytes] | bytes, encoding: str | None = None, metadata_only: bool = False) -> None: ...
    @classmethod
    def from_filename(
        cls, filename: StrPath, encoding: str = "UTF-8", errors: str | None = None, newline: str | None = None
    ) -> Self: ...
    @classmethod
    @overload
    def from_string(cls, data: str, encoding: None = None, errors: str | None = "strict") -> Self: ...
    @classmethod
    @overload
    def from_string(cls, data: bytes, encoding: str | None = None, errors: str | None = "strict") -> Self: ...
    @property
    def added_files(self) -> list[PatchedFile]: ...
    @property
    def removed_files(self) -> list[PatchedFile]: ...
    @property
    def modified_files(self) -> list[PatchedFile]: ...
    @property
    def added(self) -> int: ...
    @property
    def removed(self) -> int: ...
