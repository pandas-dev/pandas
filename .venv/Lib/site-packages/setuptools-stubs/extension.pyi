from _typeshed import StrPath
from collections.abc import Iterable

from ._distutils.extension import Extension as _Extension

def have_pyrex() -> bool: ...

class Extension(_Extension):
    py_limited_api: bool
    def __init__(
        self,
        name: str,
        sources: Iterable[StrPath],
        include_dirs: list[str] | None = None,
        define_macros: list[tuple[str, str | None]] | None = None,
        undef_macros: list[str] | None = None,
        library_dirs: list[str] | None = None,
        libraries: list[str] | None = None,
        runtime_library_dirs: list[str] | None = None,
        extra_objects: list[str] | None = None,
        extra_compile_args: list[str] | None = None,
        extra_link_args: list[str] | None = None,
        export_symbols: list[str] | None = None,
        swig_opts: list[str] | None = None,
        depends: list[str] | None = None,
        language: str | None = None,
        optional: bool | None = None,
        *,
        py_limited_api: bool = False,
    ) -> None: ...

class Library(Extension): ...
