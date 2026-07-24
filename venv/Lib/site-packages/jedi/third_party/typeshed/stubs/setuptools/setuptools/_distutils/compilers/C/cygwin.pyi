from _typeshed import StrPath
from collections.abc import Callable, Iterable
from shlex import _ShlexInstream
from typing import ClassVar, Final, Literal, NoReturn
from typing_extensions import Never, deprecated

from ...version import LooseVersion
from . import unix

def get_msvcr() -> list[str]: ...

class Compiler(unix.Compiler):
    compiler_type: ClassVar[str]
    obj_extension: ClassVar[str]
    static_lib_extension: ClassVar[str]
    shared_lib_extension: ClassVar[str]
    dylib_lib_extension: ClassVar[str]
    static_lib_format: ClassVar[str]
    shared_lib_format: ClassVar[str]
    dylib_lib_format: ClassVar[str]
    exe_extension: ClassVar[str]
    cc: str
    cxx: str
    linker_dll: str
    linker_dll_cxx: str
    dll_libraries: list[str]
    def __init__(self, verbose: bool = False, force: bool = False) -> None: ...
    @property
    @deprecated(
        "gcc_version attribute of CygwinCCompiler is deprecated. "
        "Instead of returning actual gcc version a fixed value 11.2.0 is returned."
    )
    def gcc_version(self) -> LooseVersion: ...
    # `objects` and `libraries` uses list methods
    def link(
        self,
        target_desc: str,
        objects: list[str],  # type: ignore[override]
        output_filename: str,
        output_dir: str | None = None,
        libraries: list[str] | None = None,  # type: ignore[override]
        library_dirs: list[str] | tuple[str, ...] | None = None,
        runtime_library_dirs: list[str] | tuple[str, ...] | None = None,
        export_symbols: Iterable[str] | None = None,
        debug: bool = False,
        extra_preargs: list[str] | None = None,
        extra_postargs: list[str] | None = None,
        build_temp: StrPath | None = None,
        target_lang: str | None = None,
    ) -> None: ...
    # cygwin doesn't support rpath; prints a warning and returns an empty list
    def runtime_library_dir_option(self, dir: str) -> list[Never]: ...  # type: ignore[override]
    @property
    def out_extensions(self) -> dict[str, str]: ...

class MinGW32Compiler(Compiler):
    compiler_type: ClassVar[str]
    def __init__(self, verbose: bool = False, force: bool = False) -> None: ...
    def runtime_library_dir_option(self, dir: str) -> NoReturn: ...

CONFIG_H_OK: Final = "ok"
CONFIG_H_NOTOK: Final = "not ok"
CONFIG_H_UNCERTAIN: Final = "uncertain"

def check_config_h() -> tuple[Literal["ok", "not ok", "uncertain"], str]: ...
def is_cygwincc(cc: str | _ShlexInstream) -> bool: ...

get_versions: Callable[[], tuple[LooseVersion | None, ...]] | None
