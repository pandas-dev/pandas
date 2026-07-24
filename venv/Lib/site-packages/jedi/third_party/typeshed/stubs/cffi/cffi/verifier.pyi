import io
import os
from _typeshed import Incomplete, StrPath
from typing import AnyStr
from typing_extensions import TypeAlias

NativeIO: TypeAlias = io.StringIO

class Verifier:
    ffi: Incomplete
    preamble: Incomplete
    flags: int | None
    kwds: dict[str, list[str] | tuple[str]]
    tmpdir: StrPath
    sourcefilename: str
    modulefilename: str
    ext_package: str | None
    def __init__(
        self,
        ffi,
        preamble,
        tmpdir: StrPath | None = None,
        modulename: str | None = None,
        ext_package: str | None = None,
        tag: str = "",
        force_generic_engine: bool = False,
        source_extension: str = ".c",
        flags: int | None = None,
        relative_to: os.PathLike[AnyStr] | None = None,
        **kwds: list[str] | tuple[str],
    ) -> None: ...
    def write_source(self, file=None) -> None: ...
    def compile_module(self) -> None: ...
    def load_library(self): ...
    def get_module_name(self) -> str: ...
    def get_extension(self): ...
    def generates_python_module(self) -> bool: ...
    def make_relative_to(
        self, kwds: dict[str, list[str] | tuple[str]], relative_to: os.PathLike[AnyStr] | None
    ) -> dict[str, list[str] | tuple[str]]: ...

def set_tmpdir(dirname: StrPath) -> None: ...
def cleanup_tmpdir(tmpdir: StrPath | None = None, keep_so: bool = False) -> None: ...
