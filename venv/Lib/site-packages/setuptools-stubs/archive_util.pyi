from _typeshed import Incomplete
from collections.abc import Callable

from ._distutils.errors import DistutilsError

__all__ = [
    "unpack_archive",
    "unpack_zipfile",
    "unpack_tarfile",
    "default_filter",
    "UnrecognizedFormat",
    "extraction_drivers",
    "unpack_directory",
]

class UnrecognizedFormat(DistutilsError): ...

def default_filter(src, dst): ...
def unpack_archive(filename, extract_dir, progress_filter=..., drivers: Incomplete | None = None) -> None: ...
def unpack_directory(filename, extract_dir, progress_filter=...) -> None: ...
def unpack_zipfile(filename, extract_dir, progress_filter=...) -> None: ...
def unpack_tarfile(filename, extract_dir, progress_filter=...): ...

extraction_drivers: tuple[Callable[..., Incomplete], ...]
