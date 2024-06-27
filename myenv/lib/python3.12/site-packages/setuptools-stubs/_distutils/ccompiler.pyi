from _typeshed import BytesPath, StrPath
from collections.abc import Callable, Iterable
from typing import Any, ClassVar, Literal, TypeVar, overload
from typing_extensions import TypeAlias

_Macro: TypeAlias = tuple[str] | tuple[str, str | None]
_StrPathT = TypeVar("_StrPathT", bound=StrPath)
_BytesPathT = TypeVar("_BytesPathT", bound=BytesPath)

def gen_lib_options(
    compiler: CCompiler, library_dirs: list[str], runtime_library_dirs: list[str], libraries: list[str]
) -> list[str]: ...
def gen_preprocess_options(macros: list[_Macro], include_dirs: list[str]) -> list[str]: ...
def get_default_compiler(osname: str | None = ..., platform: str | None = ...) -> str: ...
def new_compiler(
    plat: str | None = ...,
    compiler: str | None = ...,
    verbose: bool | Literal[0, 1] = 0,
    dry_run: bool | Literal[0, 1] = 0,
    force: bool | Literal[0, 1] = 0,
) -> CCompiler: ...
def show_compilers() -> None: ...

class CCompiler:
    src_extensions: ClassVar[list[str] | None]
    obj_extension: ClassVar[str | None]
    static_lib_extension: ClassVar[str | None]
    shared_lib_extension: ClassVar[str | None]
    static_lib_format: ClassVar[str | None]
    shared_lib_format: ClassVar[str | None]
    exe_extension: ClassVar[str | None]
    language_map: ClassVar[dict[str, str]]
    language_order: ClassVar[list[str]]
    dry_run: bool
    force: bool
    verbose: bool
    output_dir: str | None
    macros: list[_Macro]
    include_dirs: list[str]
    libraries: list[str]
    library_dirs: list[str]
    runtime_library_dirs: list[str]
    objects: list[str]
    def __init__(
        self, verbose: bool | Literal[0, 1] = 0, dry_run: bool | Literal[0, 1] = 0, force: bool | Literal[0, 1] = 0
    ) -> None: ...
    def add_include_dir(self, dir: str) -> None: ...
    def set_include_dirs(self, dirs: list[str]) -> None: ...
    def add_library(self, libname: str) -> None: ...
    def set_libraries(self, libnames: list[str]) -> None: ...
    def add_library_dir(self, dir: str) -> None: ...
    def set_library_dirs(self, dirs: list[str]) -> None: ...
    def add_runtime_library_dir(self, dir: str) -> None: ...
    def set_runtime_library_dirs(self, dirs: list[str]) -> None: ...
    def define_macro(self, name: str, value: str | None = ...) -> None: ...
    def undefine_macro(self, name: str) -> None: ...
    def add_link_object(self, object: str) -> None: ...
    def set_link_objects(self, objects: list[str]) -> None: ...
    def detect_language(self, sources: str | list[str]) -> str | None: ...
    def find_library_file(self, dirs: list[str], lib: str, debug: bool | Literal[0, 1] = 0) -> str | None: ...
    def has_function(
        self,
        funcname: str,
        includes: list[str] | None = ...,
        include_dirs: list[str] | None = ...,
        libraries: list[str] | None = ...,
        library_dirs: list[str] | None = ...,
    ) -> bool: ...
    def library_dir_option(self, dir: str) -> str: ...
    def library_option(self, lib: str) -> str: ...
    def runtime_library_dir_option(self, dir: str) -> str: ...
    def set_executables(self, **args: str) -> None: ...
    def compile(
        self,
        sources: list[str],
        output_dir: str | None = ...,
        macros: list[_Macro] | None = ...,
        include_dirs: list[str] | None = ...,
        debug: bool | Literal[0, 1] = 0,
        extra_preargs: list[str] | None = ...,
        extra_postargs: list[str] | None = ...,
        depends: list[str] | None = ...,
    ) -> list[str]: ...
    def create_static_lib(
        self,
        objects: list[str],
        output_libname: str,
        output_dir: str | None = ...,
        debug: bool | Literal[0, 1] = 0,
        target_lang: str | None = ...,
    ) -> None: ...
    def link(
        self,
        target_desc: str,
        objects: list[str],
        output_filename: str,
        output_dir: str | None = ...,
        libraries: list[str] | None = ...,
        library_dirs: list[str] | None = ...,
        runtime_library_dirs: list[str] | None = ...,
        export_symbols: list[str] | None = ...,
        debug: bool | Literal[0, 1] = 0,
        extra_preargs: list[str] | None = ...,
        extra_postargs: list[str] | None = ...,
        build_temp: str | None = ...,
        target_lang: str | None = ...,
    ) -> None: ...
    def link_executable(
        self,
        objects: list[str],
        output_progname: str,
        output_dir: str | None = ...,
        libraries: list[str] | None = ...,
        library_dirs: list[str] | None = ...,
        runtime_library_dirs: list[str] | None = ...,
        debug: bool | Literal[0, 1] = 0,
        extra_preargs: list[str] | None = ...,
        extra_postargs: list[str] | None = ...,
        target_lang: str | None = ...,
    ) -> None: ...
    def link_shared_lib(
        self,
        objects: list[str],
        output_libname: str,
        output_dir: str | None = ...,
        libraries: list[str] | None = ...,
        library_dirs: list[str] | None = ...,
        runtime_library_dirs: list[str] | None = ...,
        export_symbols: list[str] | None = ...,
        debug: bool | Literal[0, 1] = 0,
        extra_preargs: list[str] | None = ...,
        extra_postargs: list[str] | None = ...,
        build_temp: str | None = ...,
        target_lang: str | None = ...,
    ) -> None: ...
    def link_shared_object(
        self,
        objects: list[str],
        output_filename: str,
        output_dir: str | None = ...,
        libraries: list[str] | None = ...,
        library_dirs: list[str] | None = ...,
        runtime_library_dirs: list[str] | None = ...,
        export_symbols: list[str] | None = ...,
        debug: bool | Literal[0, 1] = 0,
        extra_preargs: list[str] | None = ...,
        extra_postargs: list[str] | None = ...,
        build_temp: str | None = ...,
        target_lang: str | None = ...,
    ) -> None: ...
    def preprocess(
        self,
        source: str,
        output_file: str | None = ...,
        macros: list[_Macro] | None = ...,
        include_dirs: list[str] | None = ...,
        extra_preargs: list[str] | None = ...,
        extra_postargs: list[str] | None = ...,
    ) -> None: ...
    @overload
    def executable_filename(self, basename: str, strip_dir: Literal[0, False] = 0, output_dir: StrPath = ...) -> str: ...
    @overload
    def executable_filename(self, basename: StrPath, strip_dir: Literal[1, True], output_dir: StrPath = ...) -> str: ...
    def library_filename(
        self, libname: str, lib_type: str = "static", strip_dir: bool | Literal[0, 1] = 0, output_dir: StrPath = ""
    ) -> str: ...
    def object_filenames(
        self, source_filenames: Iterable[StrPath], strip_dir: bool | Literal[0, 1] = 0, output_dir: StrPath | None = ...
    ) -> list[str]: ...
    @overload
    def shared_object_filename(self, basename: str, strip_dir: Literal[0, False] = 0, output_dir: StrPath = ...) -> str: ...
    @overload
    def shared_object_filename(self, basename: StrPath, strip_dir: Literal[1, True], output_dir: StrPath = ...) -> str: ...
    def execute(self, func: Callable[..., object], args: tuple[Any, ...], msg: str | None = ..., level: int = ...) -> None: ...
    def spawn(self, cmd: list[str]) -> None: ...
    def mkpath(self, name: str, mode: int = ...) -> None: ...
    @overload
    def move_file(self, src: StrPath, dst: _StrPathT) -> _StrPathT | str: ...
    @overload
    def move_file(self, src: BytesPath, dst: _BytesPathT) -> _BytesPathT | bytes: ...
    def announce(self, msg: str, level: int = ...) -> None: ...
    def warn(self, msg: str) -> None: ...
    def debug_print(self, msg: str) -> None: ...
