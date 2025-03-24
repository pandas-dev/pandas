from _typeshed import Incomplete
from typing import ClassVar, Final

from . import base

PLAT_SPEC_TO_RUNTIME: Final[dict[str, str]]

class Compiler(base.Compiler):
    compiler_type: ClassVar[str]
    executables: ClassVar[dict[str, Incomplete]]
    src_extensions: ClassVar[list[str]]
    res_extension: ClassVar[str]
    obj_extension: ClassVar[str]
    static_lib_extension: ClassVar[str]
    shared_lib_extension: ClassVar[str]
    # This was accidentally removed upstream and should be back pretty soon.
    # shared_lib_format: ClassVar[str]
    # static_lib_format = shared_lib_format
    static_lib_format: ClassVar[str]
    exe_extension: ClassVar[str]
    initialized: bool
    def initialize(self, plat_name: str | None = None) -> None: ...
    @property
    def out_extensions(self) -> dict[str, str]: ...
