from _typeshed import StrPath
from collections.abc import Iterable
from typing import ClassVar, Literal

from . import unix

class Compiler(unix.Compiler):
    src_extensions: ClassVar[list[str]]
    zos_compiler: Literal["ibm-openxl", "ibm-xlclang", "ibm-xlc"]
    def __init__(self, verbose: bool = False, force: bool = False) -> None: ...
    def runtime_library_dir_option(self, dir: str) -> str: ...
    def link(
        self,
        target_desc: str,
        objects: list[str] | tuple[str, ...],
        output_filename: str,
        output_dir: str | None = None,
        libraries: list[str] | tuple[str, ...] | None = None,
        library_dirs: list[str] | tuple[str, ...] | None = None,
        runtime_library_dirs: list[str] | tuple[str, ...] | None = None,
        export_symbols: Iterable[str] | None = None,
        debug: bool = False,
        extra_preargs: list[str] | None = None,
        extra_postargs: list[str] | None = None,
        build_temp: StrPath | None = None,
        target_lang: str | None = None,
    ) -> None: ...
