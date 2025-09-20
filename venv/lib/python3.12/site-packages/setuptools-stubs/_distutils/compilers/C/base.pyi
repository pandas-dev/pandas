from _typeshed import BytesPath, Incomplete, StrPath, Unused
from collections.abc import Callable, Iterable, MutableSequence, Sequence
from typing import ClassVar, Final, Literal, TypeVar, overload
from typing_extensions import TypeAlias, TypeVarTuple, Unpack

_Macro: TypeAlias = tuple[str] | tuple[str, str | None]
_StrPathT = TypeVar("_StrPathT", bound=StrPath)
_BytesPathT = TypeVar("_BytesPathT", bound=BytesPath)
_Ts = TypeVarTuple("_Ts")

class Compiler:
    compiler_type: ClassVar[str]
    executables: ClassVar[dict[str, Incomplete]]

    # Subclasses that rely on the standard filename generation methods
    # implemented below should override these
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

    SHARED_OBJECT: Final = "shared_object"
    SHARED_LIBRARY: Final = "shared_library"
    EXECUTABLE: Final = "executable"
    def __init__(self, verbose: bool = False, dry_run: bool = False, force: bool = False) -> None: ...
    def add_include_dir(self, dir: str) -> None: ...
    def set_include_dirs(self, dirs: list[str]) -> None: ...
    def add_library(self, libname: str) -> None: ...
    def set_libraries(self, libnames: list[str]) -> None: ...
    def add_library_dir(self, dir: str) -> None: ...
    def set_library_dirs(self, dirs: list[str]) -> None: ...
    def add_runtime_library_dir(self, dir: str) -> None: ...
    def set_runtime_library_dirs(self, dirs: list[str]) -> None: ...
    def define_macro(self, name: str, value: str | None = None) -> None: ...
    def undefine_macro(self, name: str) -> None: ...
    def add_link_object(self, object: str) -> None: ...
    def set_link_objects(self, objects: list[str]) -> None: ...
    def detect_language(self, sources: str | list[str]) -> str | None: ...
    def find_library_file(self, dirs: Iterable[str], lib: str, debug: bool = False) -> str | None: ...
    def has_function(
        self,
        funcname: str,
        includes: Iterable[str] | None = None,
        include_dirs: list[str] | tuple[str, ...] | None = None,
        libraries: list[str] | None = None,
        library_dirs: list[str] | tuple[str, ...] | None = None,
    ) -> bool: ...
    def library_dir_option(self, dir: str) -> str: ...
    def library_option(self, lib: str) -> str: ...
    def runtime_library_dir_option(self, dir: str) -> str: ...
    def set_executables(self, **kwargs: str) -> None: ...
    def set_executable(self, key: str, value) -> None: ...
    def compile(
        self,
        sources: Sequence[StrPath],
        output_dir: str | None = None,
        macros: list[_Macro] | None = None,
        include_dirs: list[str] | tuple[str, ...] | None = None,
        debug: bool = False,
        extra_preargs: list[str] | None = None,
        extra_postargs: list[str] | None = None,
        depends: list[str] | tuple[str, ...] | None = None,
    ) -> list[str]: ...
    def create_static_lib(
        self,
        objects: list[str] | tuple[str, ...],
        output_libname: str,
        output_dir: str | None = None,
        debug: bool = False,
        target_lang: str | None = None,
    ) -> None: ...
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
    def link_executable(
        self,
        objects: list[str] | tuple[str, ...],
        output_progname: str,
        output_dir: str | None = None,
        libraries: list[str] | tuple[str, ...] | None = None,
        library_dirs: list[str] | tuple[str, ...] | None = None,
        runtime_library_dirs: list[str] | tuple[str, ...] | None = None,
        debug: bool = False,
        extra_preargs: list[str] | None = None,
        extra_postargs: list[str] | None = None,
        target_lang: str | None = None,
    ) -> None: ...
    def link_shared_lib(
        self,
        objects: list[str] | tuple[str, ...],
        output_libname: str,
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
    def link_shared_object(
        self,
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
    def preprocess(
        self,
        source: StrPath,
        output_file: StrPath | None = None,
        macros: list[_Macro] | None = None,
        include_dirs: list[str] | tuple[str, ...] | None = None,
        extra_preargs: list[str] | None = None,
        extra_postargs: Iterable[str] | None = None,
    ) -> None: ...
    @overload
    def executable_filename(self, basename: str, strip_dir: Literal[False] = False, output_dir: StrPath = "") -> str: ...
    @overload
    def executable_filename(self, basename: StrPath, strip_dir: Literal[True], output_dir: StrPath = "") -> str: ...
    def library_filename(
        self, libname: str, lib_type: str = "static", strip_dir: bool = False, output_dir: StrPath = ""
    ) -> str: ...
    @property
    def out_extensions(self) -> dict[str, str]: ...
    def object_filenames(
        self, source_filenames: Iterable[StrPath], strip_dir: bool = False, output_dir: StrPath | None = ""
    ) -> list[str]: ...
    @overload
    def shared_object_filename(self, basename: str, strip_dir: Literal[False] = False, output_dir: StrPath = "") -> str: ...
    @overload
    def shared_object_filename(self, basename: StrPath, strip_dir: Literal[True], output_dir: StrPath = "") -> str: ...
    def execute(
        self, func: Callable[[Unpack[_Ts]], Unused], args: tuple[Unpack[_Ts]], msg: str | None = None, level: int = 1
    ) -> None: ...
    def spawn(self, cmd: MutableSequence[bytes | StrPath]) -> None: ...
    def mkpath(self, name: str, mode: int = 0o777) -> None: ...
    @overload
    def move_file(self, src: StrPath, dst: _StrPathT) -> _StrPathT | str: ...
    @overload
    def move_file(self, src: BytesPath, dst: _BytesPathT) -> _BytesPathT | bytes: ...
    def announce(self, msg: str, level: int = 1) -> None: ...
    def warn(self, msg: str) -> None: ...
    def debug_print(self, msg: str) -> None: ...

def get_default_compiler(osname: str | None = None, platform: str | None = None) -> str: ...

compiler_class: dict[str, tuple[str, str, str]]

def show_compilers() -> None: ...
def new_compiler(
    plat: str | None = None, compiler: str | None = None, verbose: bool = False, dry_run: bool = False, force: bool = False
) -> Compiler: ...
def gen_preprocess_options(macros: Iterable[_Macro], include_dirs: Iterable[str]) -> list[str]: ...
def gen_lib_options(
    compiler: Compiler, library_dirs: Iterable[str], runtime_library_dirs: Iterable[str], libraries: Iterable[str]
) -> list[str]: ...
