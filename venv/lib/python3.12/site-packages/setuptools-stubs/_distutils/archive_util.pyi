from _typeshed import StrOrBytesPath, StrPath
from typing import Literal, overload

@overload
def make_archive(
    base_name: str,
    format: str,
    root_dir: StrOrBytesPath | None = None,
    base_dir: str | None = None,
    verbose: bool = False,
    dry_run: bool = False,
    owner: str | None = None,
    group: str | None = None,
) -> str: ...
@overload
def make_archive(
    base_name: StrPath,
    format: str,
    root_dir: StrOrBytesPath,
    base_dir: str | None = None,
    verbose: bool = False,
    dry_run: bool = False,
    owner: str | None = None,
    group: str | None = None,
) -> str: ...
def make_tarball(
    base_name: str,
    base_dir: StrPath,
    compress: Literal["gzip", "bzip2", "xz"] | None = "gzip",
    verbose: bool = False,
    dry_run: bool = False,
    owner: str | None = None,
    group: str | None = None,
) -> str: ...
def make_zipfile(base_name: str, base_dir: StrPath, verbose: bool = False, dry_run: bool = False) -> str: ...
