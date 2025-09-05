from typing import ClassVar

from . import base

class Compiler(base.Compiler):
    src_extensions: ClassVar[list[str]]
    obj_extension: ClassVar[str]
    static_lib_extension: ClassVar[str]
    shared_lib_extension: ClassVar[str]
    dylib_lib_extension: ClassVar[str]
    xcode_stub_lib_extension: ClassVar[str]
    static_lib_format: ClassVar[str]
    shared_lib_format: ClassVar[str]
    dylib_lib_format: ClassVar[str]
    xcode_stub_lib_format: ClassVar[str]
    def runtime_library_dir_option(self, dir: str) -> str | list[str]: ...  # type: ignore[override]
