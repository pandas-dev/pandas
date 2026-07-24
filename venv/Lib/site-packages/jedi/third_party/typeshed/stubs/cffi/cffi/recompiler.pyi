import io
from _typeshed import Incomplete, StrPath
from typing import Final
from typing_extensions import TypeAlias

from .cffi_opcode import *
from .error import VerificationError as VerificationError

VERSION_BASE: Final = 9729
VERSION_EMBEDDED: Final = 9985
VERSION_CHAR16CHAR32: Final = 10241
USE_LIMITED_API: Final = True

class GlobalExpr:
    name: Incomplete
    address: Incomplete
    type_op: Incomplete
    size: Incomplete
    check_value: Incomplete
    def __init__(self, name, address, type_op, size: int = 0, check_value: int = 0) -> None: ...
    def as_c_expr(self) -> str: ...
    def as_python_expr(self) -> str: ...

class FieldExpr:
    name: Incomplete
    field_offset: Incomplete
    field_size: Incomplete
    fbitsize: Incomplete
    field_type_op: Incomplete
    def __init__(self, name, field_offset, field_size, fbitsize, field_type_op) -> None: ...
    def as_c_expr(self) -> str: ...
    def as_python_expr(self) -> None: ...
    def as_field_python_expr(self) -> str: ...

class StructUnionExpr:
    name: Incomplete
    type_index: Incomplete
    flags: Incomplete
    size: Incomplete
    alignment: Incomplete
    comment: Incomplete
    first_field_index: Incomplete
    c_fields: Incomplete
    def __init__(self, name, type_index, flags, size, alignment, comment, first_field_index, c_fields) -> None: ...
    def as_c_expr(self) -> str: ...
    def as_python_expr(self) -> str: ...

class EnumExpr:
    name: Incomplete
    type_index: Incomplete
    size: Incomplete
    signed: Incomplete
    allenums: Incomplete
    def __init__(self, name, type_index, size, signed, allenums) -> None: ...
    def as_c_expr(self) -> str: ...
    def as_python_expr(self) -> str: ...

class TypenameExpr:
    name: Incomplete
    type_index: Incomplete
    def __init__(self, name, type_index) -> None: ...
    def as_c_expr(self) -> str: ...
    def as_python_expr(self) -> str: ...

class Recompiler:
    ffi: Incomplete
    module_name: str
    target_is_python: bool
    def __init__(self, ffi, module_name: str, target_is_python: bool = False) -> None: ...
    def needs_version(self, ver: int) -> None: ...
    cffi_types: list[Incomplete]
    def collect_type_table(self) -> None: ...
    ALL_STEPS: list[str]
    def collect_step_tables(self) -> None: ...
    def write_source_to_f(self, f, preamble: str) -> None: ...
    def write_c_source_to_f(self, f, preamble: str) -> None: ...
    def write_py_source_to_f(self, f) -> None: ...

NativeIO: TypeAlias = io.StringIO

def make_c_source(ffi, module_name: str, preamble: str, target_c_file, verbose: bool = False): ...
def make_py_source(ffi, module_name: str, target_py_file, verbose: bool = False): ...
def recompile(
    ffi,
    module_name: str | bytes,
    preamble: str | None,
    tmpdir: str = ".",
    call_c_compiler: bool = True,
    c_file=None,
    source_extension: str = ".c",
    extradir: StrPath | None = None,
    compiler_verbose: int = 1,
    target: str | None = None,
    debug: int | None = None,
    uses_ffiplatform: bool = True,
    **kwds,
): ...
