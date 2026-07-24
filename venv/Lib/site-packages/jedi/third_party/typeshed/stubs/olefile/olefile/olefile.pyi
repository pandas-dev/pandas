import array
import datetime
import io
import logging
import traceback
from collections.abc import Sequence
from typing import IO, AnyStr, Generic
from typing_extensions import Self, TypeAlias

__date__: str
__version__: str
__author__: str

__all__ = [
    "isOleFile",
    "OleFileIO",
    "OleMetadata",
    "enable_logging",
    "MAGIC",
    "STGTY_EMPTY",
    "STGTY_STREAM",
    "STGTY_STORAGE",
    "STGTY_ROOT",
    "STGTY_PROPERTY",
    "STGTY_LOCKBYTES",
    "MINIMAL_OLEFILE_SIZE",
    "DEFECT_UNSURE",
    "DEFECT_POTENTIAL",
    "DEFECT_INCORRECT",
    "DEFECT_FATAL",
    "DEFAULT_PATH_ENCODING",
    "MAXREGSECT",
    "DIFSECT",
    "FATSECT",
    "ENDOFCHAIN",
    "FREESECT",
    "MAXREGSID",
    "NOSTREAM",
    "UNKNOWN_SIZE",
    "WORD_CLSID",
    "OleFileIONotClosed",
]

UINT32: str

DEFAULT_PATH_ENCODING: str | None

def get_logger(name: str, level: int = 51) -> logging.Logger: ...

log: logging.Logger

def enable_logging() -> None: ...

MAGIC: bytes

MAXREGSECT: int
DIFSECT: int
FATSECT: int
ENDOFCHAIN: int
FREESECT: int

MAXREGSID: int
NOSTREAM: int

STGTY_EMPTY: int
STGTY_STORAGE: int
STGTY_STREAM: int
STGTY_LOCKBYTES: int
STGTY_PROPERTY: int
STGTY_ROOT: int

UNKNOWN_SIZE: int

VT_EMPTY: int
VT_NULL: int
VT_I2: int
VT_I4: int
VT_R4: int
VT_R8: int
VT_CY: int
VT_DATE: int
VT_BSTR: int
VT_DISPATCH: int
VT_ERROR: int
VT_BOOL: int
VT_VARIANT: int
VT_UNKNOWN: int
VT_DECIMAL: int
VT_I1: int
VT_UI1: int
VT_UI2: int
VT_UI4: int
VT_I8: int
VT_UI8: int
VT_INT: int
VT_UINT: int
VT_VOID: int
VT_HRESULT: int
VT_PTR: int
VT_SAFEARRAY: int
VT_CARRAY: int
VT_USERDEFINED: int
VT_LPSTR: int
VT_LPWSTR: int
VT_FILETIME: int
VT_BLOB: int
VT_STREAM: int
VT_STORAGE: int
VT_STREAMED_OBJECT: int
VT_STORED_OBJECT: int
VT_BLOB_OBJECT: int
VT_CF: int
VT_CLSID: int
VT_VECTOR: int

VT: dict[int, str]

WORD_CLSID: str
DEFECT_UNSURE: int
DEFECT_POTENTIAL: int
DEFECT_INCORRECT: int
DEFECT_FATAL: int
MINIMAL_OLEFILE_SIZE: int

def isOleFile(filename: IO[bytes] | bytes | str | None = None, data: bytes | None = None) -> bool: ...
def i8(c: bytes | int) -> int: ...
def i16(c: bytes, o: int = 0) -> int: ...
def i32(c: bytes, o: int = 0) -> int: ...
def _clsid(clsid: bytes) -> str: ...
def filetime2datetime(filetime: int) -> datetime.datetime: ...

class OleFileError(IOError): ...
class NotOleFileError(OleFileError): ...

class OleMetadata:
    SUMMARY_ATTRIBS: list[str]
    DOCSUM_ATTRIBS: list[str]

    def __init__(self) -> None: ...
    def parse_properties(self, ole_file: OleFileIO[AnyStr]) -> None: ...
    def dump(self) -> None: ...

class OleFileIONotClosed(RuntimeWarning):
    def __init__(self, stack_of_open: traceback.FrameSummary | None = None) -> None: ...

class OleStream(io.BytesIO):
    def __init__(
        self,
        fp: IO[bytes],
        sect: int,
        size: int,
        offset: int,
        sectorsize: int,
        fat: list[int],
        filesize: int,
        olefileio: OleFileIO[AnyStr],
    ) -> None: ...

class OleDirectoryEntry(Generic[AnyStr]):
    STRUCT_DIRENTRY: str
    DIRENTRY_SIZE: int
    clsid: str

    def __init__(self, entry: bytes, sid: int, ole_file: OleFileIO[AnyStr]) -> None: ...
    def build_sect_chain(self, ole_file: OleFileIO[AnyStr]) -> None: ...
    def build_storage_tree(self) -> None: ...
    def append_kids(self, child_sid: int) -> None: ...
    def __eq__(self, other: OleDirectoryEntry[AnyStr]) -> bool: ...  # type: ignore[override]
    def __lt__(self, other: OleDirectoryEntry[AnyStr]) -> bool: ...
    def __ne__(self, other: OleDirectoryEntry[AnyStr]) -> bool: ...  # type: ignore[override]
    def __le__(self, other: OleDirectoryEntry[AnyStr]) -> bool: ...
    def dump(self, tab: int = 0) -> None: ...
    def getmtime(self) -> datetime.datetime | None: ...
    def getctime(self) -> datetime.datetime | None: ...

_Property: TypeAlias = int | str | bytes | bool | None

class OleFileIO(Generic[AnyStr]):
    root: OleDirectoryEntry[AnyStr] | None

    def __init__(
        self,
        filename: IO[bytes] | AnyStr | None = None,
        raise_defects: int = 40,
        write_mode: bool = False,
        debug: bool = False,
        path_encoding: str | None = DEFAULT_PATH_ENCODING,  # noqa: Y011
    ) -> None: ...
    def __del__(self) -> None: ...
    def __enter__(self) -> Self: ...
    def __exit__(self, *args: object) -> None: ...
    def _raise_defect(
        self, defect_level: int, message: str, exception_type: type[Exception] = OleFileError  # noqa: Y011
    ) -> None: ...
    def _decode_utf16_str(self, utf16_str: bytes, errors: str = "replace") -> str | bytes: ...
    def open(self, filename: IO[bytes] | AnyStr, write_mode: bool = False) -> None: ...
    def close(self) -> None: ...
    def _close(self, warn: bool = False) -> None: ...
    def _check_duplicate_stream(self, first_sect: int, minifat: bool = False) -> None: ...
    def dumpfat(self, fat: Sequence[int], firstindex: int = 0) -> None: ...
    def dumpsect(self, sector: bytes, firstindex: int = 0) -> None: ...
    def sect2array(self, sect: bytes) -> Sequence[int]: ...
    def loadfat_sect(self, sect: bytes | array.array[int]) -> int | None: ...
    def loadfat(self, header: bytes) -> None: ...
    def loadminifat(self) -> None: ...
    def getsect(self, sect: int) -> bytes: ...
    def write_sect(self, sect: int, data: bytes, padding: bytes = b"\x00") -> None: ...
    def _write_mini_sect(self, fp_pos: int, data: bytes, padding: bytes = b"\x00") -> None: ...
    def loaddirectory(self, sect: int) -> None: ...
    def _load_direntry(self, sid: int) -> OleDirectoryEntry[AnyStr]: ...
    def dumpdirectory(self) -> None: ...
    def _open(self, start: int, size: int = 0x7FFFFFFF, force_FAT: bool = False) -> OleStream: ...
    def _list(
        self,
        files: list[list[AnyStr]],
        prefix: list[AnyStr],
        node: OleDirectoryEntry[AnyStr],
        streams: bool = True,
        storages: bool = False,
    ) -> None: ...
    def listdir(self, streams: bool = True, storages: bool = False) -> list[list[AnyStr]]: ...
    def _find(self, filename: str | Sequence[str]) -> int: ...
    def openstream(self, filename: AnyStr | Sequence[AnyStr]) -> OleStream: ...
    def _write_mini_stream(self, entry: OleDirectoryEntry[AnyStr], data_to_write: bytes) -> None: ...
    def write_stream(self, stream_name: str | Sequence[str], data: bytes) -> None: ...
    def get_type(self, filename: AnyStr | Sequence[AnyStr]) -> bool | int: ...
    def getclsid(self, filename: AnyStr | Sequence[AnyStr]) -> str: ...
    def getmtime(self, filename: AnyStr | Sequence[AnyStr]) -> datetime.datetime | None: ...
    def getctime(self, filename: AnyStr | Sequence[AnyStr]) -> datetime.datetime | None: ...
    def exists(self, filename: AnyStr | Sequence[AnyStr]) -> bool: ...
    def get_size(self, filename: AnyStr | Sequence[AnyStr]) -> int: ...
    def get_rootentry_name(self) -> bytes: ...
    def getproperties(
        self, filename: AnyStr | Sequence[AnyStr], convert_time: bool = False, no_conversion: list[int] | None = None
    ) -> dict[int, list[_Property] | _Property]: ...
    def _parse_property(
        self, s: bytes, offset: int, property_id: int, property_type: int, convert_time: bool, no_conversion: list[int]
    ) -> list[_Property] | _Property: ...
    def _parse_property_basic(
        self, s: bytes, offset: int, property_id: int, property_type: int, convert_time: bool, no_conversion: list[int]
    ) -> tuple[_Property, int]: ...
    def get_metadata(self) -> OleMetadata: ...
    def get_userdefined_properties(
        self, filename: AnyStr | Sequence[AnyStr], convert_time: bool = False, no_conversion: list[int] | None = None
    ) -> list[dict[str, bytes | int | None]]: ...

def main() -> None: ...
